use core::num;
use std::{marker::PhantomData, iter::Map, error::Error, fmt::Display};

use ark_ec::{Group, pairing::Pairing};
use ark_ff::{PrimeField, UniformRand, Zero, One};
use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, Evaluations, univariate::DensePolynomial, Polynomial, DenseUVPolynomial};
use nalgebra::{SMatrix, DMatrix, DVector};
use rand::RngCore;

use crate::data_structures::{VectorCommitment, VCUniversalParams};

/// KZGKey represents the universal parameters, AKA reference string, for both
/// committing polynomials and verifying commitments
#[derive(Clone, Debug)]
struct KZGKey<G1: Group> {
    // The max degree this reference string supports
    degree: usize,

    // Reference string in G1 group
    ref_string_g1: Vec<G1>,
}

impl<G1, F> KZGKey<G1>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
{
    fn new_from_secret(secret: F, max_degree: usize) -> Self {
        let g1 = G1::generator();
        let mut params: Self = Self {
            degree: max_degree,
            ref_string_g1: vec![],
        };
        params.ref_string_g1.push(g1);

        let mut sec_cur: F = F::one();
        for i in 1..max_degree {
            sec_cur = sec_cur * secret; //α^i
            params.ref_string_g1.push(g1 * sec_cur);
        }

        params
    }

    /// Multiplies the array of coefficients with the reference string in G1
    /// and sums the group elements. 
    fn commit_g1(&self, coeffs: &[F]) -> (G1, usize) {
        let mut sum = G1::zero();
        let mut i = 0;
        for c in coeffs {
            sum += self.element_at_g1(i) * c;
            i += 1;
        }

        (sum, i)
    }

    fn element_at_g1(&self, index: usize) -> G1 {
        self.ref_string_g1[index]
    }
}

impl<G1, F> VCUniversalParams for KZGKey<G1>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
{
    fn max_size(&self) -> usize {
        self.degree
    }
}

/// KZGPreparedData will create a polynomial to commit to by interpolating the dataset.
/// Instead of the evaluations occuring at [0,1,...,n-1] they instead occur at [0,w,w^2,...w^n-1].
/// w is the n-th root of unity
struct KZGPreparedData<F: PrimeField> {
    /// The polynomial representing the set of data this originated from
    data: DensePolynomial<F>,

    /// The n-th root of unity
    w: F,

    /// The size of the data-set (not necessairly equal to n-th root of unity)
    size: usize,

    /// The Toeplitz matrix of the coefficients of `data`. This matrix is used, along with the reference
    /// string, to create a vector of proofs for *all* points in the evaluation domain.
    matrix: DMatrix<F>
}

impl<F: PrimeField> FromIterator<F> for KZGPreparedData<F> {
    fn from_iter<T: IntoIterator<Item = F>>(iter: T) -> Self {
        let points: Vec<F> = iter.into_iter().collect();
        let len = points.len();
        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::<F>::new(points.len()).unwrap();
        let poly = Evaluations::from_vec_and_domain(points, domain).interpolate();
        let num_coeffs = poly.degree()+1;
        Self {
            data: poly,
            //TODO: Need to make sure this is a *primitive* root of unity aka every i < d, w^i != 0
            w: domain.group_gen(),
            size: len,
            //matrix: DMatrix::<F>::from_element(0, 0, F::zero())
            //matrix: DMatrix::zeros(len, len)
            matrix: DMatrix::from_element(num_coeffs, num_coeffs, F::zero())
        }
    }
}

impl<F: PrimeField> KZGPreparedData<F> {
    /// Prepare the Toeplitz matrix based on the coefficients of the data polynomial.
    fn toeplitz_matrix(&mut self) {
        let coeffs = self.data.coeffs();
        let n = coeffs.len();
        let elements = DVector::from_row_slice(coeffs).transpose();

        println!("Elements {}x{}", elements.nrows(), elements.ncols());
        self.matrix.row_iter_mut().enumerate().for_each(|(i, mut row)| {
            let new_row = elements.view((0,i), (1, n-i));
            row.view_mut((0, i), (1, n-i)).copy_from(&new_row);
        })
        /*
        for i in 0..coeffs.len() {
            for j in 0..coeffs.len() {
                if i >= j {
                    self.matrix[(i, j)] = coeffs[i -j];
                }
            }
        }
        */
    }

    /// Evaluate the prepared data polynomial at an `index` in the range [0,...,n]. This will translate it to
    /// the corresponding root of unity raised to `index`
    fn evaluate_by_index(&self, index: usize) -> F {
        // TODO: Not fully-safe as we are ignoring error of F::from_str
        let eval_index = self.w.pow(&F::from_str(&index.to_string()).ok().unwrap().into_bigint());
        self.data.evaluate(&eval_index)
    }

    /// Return an array of the coefficients of this data's polynomial
    fn coeffs(&self) -> &[F] {
        self.data.coeffs()
    }
}

struct KZGCommitment<G: Group> {
    commit: G
}

impl<G: Group> KZGCommitment<G> {
    fn new(commit: G) -> Self {
        Self {
            commit
        }
    }
}

struct KZGProof<F: PrimeField, G: Group> {
    commit: G,
    data: F
}

struct KZGBatchProof<F: PrimeField, G:Group> {
    commit: G,
    data: Map<F,F>
}

#[derive(Clone, Debug)]
struct KZGError {}

impl KZGError {
    fn new() -> Self {
        Self {}
    }
}

impl Display for KZGError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", "")
    }
}

impl Error for KZGError {}

/// Implementation of the Feist-Khovratovich technique of "Fast Amortized KZG proofs".
/// 
struct KZGAmortized<E: Pairing> {
    _engine: PhantomData<E>
}

impl<E: Pairing> VectorCommitment<E::ScalarField> for KZGAmortized<E>
{
    type UniversalParams = KZGKey<E::G1>;
    type PreparedData = KZGPreparedData<E::ScalarField>;
    type Commitment = KZGCommitment<E::G1>;
    type Proof = KZGProof<E::ScalarField, E::G1>;
    type Error = KZGError;

    fn setup<R: rand::RngCore>(
        max_items: usize,
        rng: &mut R
    ) -> Result<Self::UniversalParams, Self::Error> {
        let secret = <E::ScalarField as UniformRand>::rand(rng);
        Ok( KZGKey::new_from_secret(secret, max_items) )
    }

    fn commit(
        key: &Self::UniversalParams,
        data: &Self::PreparedData,
    ) -> Result<Self::Commitment, Self::Error> {
        let commit = key.commit_g1(data.coeffs());
        Ok(KZGCommitment::new(commit.0))
    }

    fn prove(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        data: &Self::PreparedData,
    ) -> Result<Self::Proof, Self::Error> {
        todo!()    
    }

    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof
    ) -> Result<bool, Self::Error> {
        todo!()    
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::{Bn254};

    type F = <Bn254 as Pairing>::ScalarField;
    type G1 = <Bn254 as Pairing>::G1;
    type poly = DensePolynomial<<Bn254 as Pairing>::ScalarField>;

    fn gen_data(num: usize) -> Vec<F> {
        let mut data: Vec<F> = vec![];
        let mut rng = rand::thread_rng();
        for _i in 0..num {
            data.push(F::rand(&mut rng));
        }
        data
    }

    #[test]
    fn test_interpolation() {
        let data = gen_data(10);
        let mut prepared_data = KZGPreparedData::from_iter(data.to_vec());

        for (i, d) in data.iter().enumerate() {
            assert_eq!(prepared_data.evaluate_by_index(i), d.clone());
        }

        prepared_data.toeplitz_matrix();

        let matrix = prepared_data.matrix.clone_owned();
        let crs: KZGKey<<Bn254 as Pairing>::G1> = KZGKey::new_from_secret(F::from(10), 16);
        let crs_vec = DVector::from_row_slice(&crs.ref_string_g1).transpose();

        //let his = matrix * crs_vec;

        //let mut to_dft: Vec<F> = his.into();
    }

    #[test]
    fn test_blah() {
        // https://docs.rs/ark-poly/0.4.2/ark_poly/domain/trait.EvaluationDomain.html#method.fft

        let data = gen_data(2);
        let prepared_data = KZGPreparedData::from_iter(data.to_vec());
        let crs: KZGKey<<Bn254 as Pairing>::G1> = KZGKey::new_from_secret(F::from(2), 16);

        let g1_zeros = vec![G1::zero(); 2];
        let mut s_hat:Vec<G1> = crs.clone().ref_string_g1;
        s_hat.reverse();
        s_hat.extend(g1_zeros.to_owned());

        let f_zeros = vec![F::zero(); 2];
        let mut c_hat: Vec<F> = prepared_data.coeffs().to_vec();
        c_hat.reverse();
        c_hat.extend(f_zeros);

        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::<F>::new(4).unwrap();
        let y: Vec<G1> = domain.fft(&s_hat);
        let v: Vec<F> = domain.fft(&c_hat);

        let v_unity: F = domain.group_gen();
        let mut t = v_unity.clone();
        let mut vs: Vec<F> = vec![F::one()];
        for i in 1..4 {
            vs.push(t);
            t *= t;
        }

        //let y_vec = DVector::from_vec(y);
        //let v_vec = DVector::from_vec(v);
        //let vs_vec = DVector::from_vec(vs);

        //let u = y_vec.component_mul(&v_vec);
        let mut u = elementwise_mul(&y, &v);
        u = elementwise_mul(&u, &vs);

        let mut h_hat: Vec<G1> = domain.ifft(&u);
        println!("{:?}\n", h_hat);

        h_hat.extend(g1_zeros); //TODO: We only need n-m zeros (m is the degree of poly, n is 2^x proofs generated)
        let proofs = domain.fft(&h_hat);

        println!("proofs {:?}", proofs);
    }

    /*
    fn dot<G1: Group<ScalarField = F>, F: PrimeField>(v1: &Vec<G1>, v2: &Vec<F>) -> G1 {
        let mut result = G1::zero();
        for (i, &element) in v1.iter().enumerate() {
            result += element*v2[i];
        }

        result
    }
    */
    fn elementwise_mul<F: PrimeField, G1: Group<ScalarField = F>>(g_vec: &Vec<G1>, f_vec: &Vec<F>) -> Vec<G1> {
        let mut result: Vec<G1> = vec![];
        for (g, f) in g_vec.iter().zip(f_vec.iter()) {
            result.push(*g * *f);
        }

        result
    }
}