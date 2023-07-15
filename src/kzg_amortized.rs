use std::{marker::PhantomData, error::Error, fmt::Display, collections::HashMap};

use ark_ec::{Group, pairing::Pairing};
use ark_ff::{PrimeField, UniformRand, Zero, One, Field};
use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, Evaluations, univariate::DensePolynomial, Polynomial, DenseUVPolynomial};

use crate::data_structures::{VectorCommitment, VCUniversalParams};

/// KZGKey represents the universal parameters, AKA reference string, for both
/// committing polynomials and verifying commitments
#[derive(Clone, Debug)]
struct KZGKey<G1: Group, G2: Group> {
    // The max degree this reference string supports
    degree: usize,

    // Reference string in G1 group
    ref_string_g1: Vec<G1>,

    // Reference string in G2 group
    ref_string_g2: Vec<G2>
}

impl<F, G1, G2> KZGKey<G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>
{
    fn new_from_secret(secret: F, max_degree: usize) -> Self {
        let g1 = G1::generator();
        let g2 = G2::generator();
        let mut params: Self = Self {
            degree: max_degree,
            ref_string_g1: vec![],
            ref_string_g2: vec![]
        };
        params.ref_string_g1.push(g1);
        params.ref_string_g2.push(g2);

        let mut sec_cur: F = F::one();
        for i in 1..max_degree {
            sec_cur = sec_cur * secret; //α^i
            params.ref_string_g1.push(g1 * sec_cur);
            params.ref_string_g2.push(g2 * sec_cur);
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

    fn element_at_g2(&self, index: usize) -> G2 {
        self.ref_string_g2[index]
    }
}

impl<F, G1, G2> VCUniversalParams for KZGKey<G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>
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
        }
    }
}

impl<F: PrimeField> KZGPreparedData<F> {
    /// Evaluate the prepared data polynomial at an `index` in the range [0,...,n]. This will translate it to
    /// the corresponding root of unity raised to `index`
    fn evaluate_by_index(&self, index: usize) -> F {
        // TODO: Not fully-safe as we are ignoring error of F::from_str
        let eval_index = self.index_to_evaluation_point(index);
        self.data.evaluate(&eval_index)
    }

    /// Return an array of the coefficients of this data's polynomial
    fn coeffs(&self) -> &[F] {
        self.data.coeffs()
    }

    /// Return polynomial degree size
    fn degree(&self) -> usize {
        self.data.degree()
    }

    /// Returns w^i
    fn index_to_evaluation_point(&self, index: usize) -> F {
        // TODO: Look at Radix2EvaluationDomain.element()
        self.w.pow(&F::from_str(&index.to_string()).ok().unwrap().into_bigint())
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

    /// index is the position in the vector, NOT the evaluation point
    index: usize,

    /// Evaluation point
    /// TODO: pass root of unity with proof/batch proof?
    point: F,
    
    data: F
}

struct KZGBatchProof<F: PrimeField, G:Group> {
    //commit: G,

    // For now, there is no proof compression. This is a map of index -> (eval point, output, proof)
    proofs: HashMap<usize,(F, F, G)>
}

impl<F: PrimeField, G: Group> Default for KZGBatchProof<F, G> {
    fn default() -> Self {
        Self {
            proofs: HashMap::new()
        }
    }
}

#[derive(Clone, Debug)]
enum KZGError {
    DefaultError,
    DataExceedsMaxSize,
    InvalidDomain
}

impl KZGError {
    fn new() -> Self {
        Self::DefaultError
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
    type UniversalParams = KZGKey<E::G1, E::G2>;
    type PreparedData = KZGPreparedData<E::ScalarField>;
    type Commitment = KZGCommitment<E::G1>;
    type Proof = KZGProof<E::ScalarField, E::G1>;
    type BatchProof = KZGBatchProof<E::ScalarField, E::G1>;
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
        let commit = Self::single_proof(&key, &data, index)?;
        Ok(
            Self::Proof {
                commit,
                index,
                data: data.evaluate_by_index(index),
                point: data.index_to_evaluation_point(index)
            }
        )
    }

    fn prove_all(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        data: &Self::PreparedData
    ) -> Result<Self::BatchProof, Self::Error> {
        Self::get_all_proofs(&key, &data)
    }

    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof
    ) -> Result<bool, Self::Error> {
        // TODO: Need to be able to access domain's root of unity to be able to evaluate properly
        let pairing1 = E::pairing(proof.commit, key.element_at_g2(1) - (E::G2::generator() * proof.point));
        let pairing2 = E::pairing(commitment.commit - E::G1::generator()*proof.data, E::G2::generator());
        
        Ok ( pairing1 == pairing2 )
    }

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof
    ) -> Result<bool, Self::Error> {
        let s = key.element_at_g2(1);
        for (index, p) in proof.proofs.iter() {
            let pairing1 = E::pairing(p.2, s - (E::G2::generator() * p.0));
            let pairing2 = E::pairing(commitment.commit - E::G1::generator() * p.1, E::G2::generator());

            if pairing1 != pairing2 {
                return Ok ( false ) 
            }
        }

        Ok( true )
    }
}

impl<E: Pairing> KZGAmortized<E> {
    /// Get all the proofs for every coefficient of `data`'s polynomial representation
    fn get_all_proofs(
        key: &KZGKey<E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>,
    ) -> Result<KZGBatchProof<E::ScalarField, E::G1>, KZGError> {
        let all_proofs = Self::build_witness_matrix(key, data)?;
        let mut res = KZGBatchProof::default();

        let domain = GeneralEvaluationDomain::<E::ScalarField>::new(data.degree()).unwrap();
        let proofs = domain.fft(&all_proofs);

        for index in 0..data.size {
            let point = data.index_to_evaluation_point(index);
            res.proofs.insert(index, (point, data.evaluate_by_index(index), proofs[index]));
        }

        Ok(res)
    }

    /// Build the group elements from the FFT of the polynomial coefficients multiplied by the reference string
    fn build_witness_matrix(
        key: &KZGKey<E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>
    ) -> Result<Vec<E::G1>, KZGError> {
        let degree = data.degree();
        let coeffs = data.coeffs().to_vec();
        let domain = GeneralEvaluationDomain::<E::ScalarField>::new(degree*2).ok_or(KZGError::InvalidDomain)?;
        let domain_size = domain.size();

        let mut c_hat: Vec<E::ScalarField> = vec![coeffs[degree]];
        c_hat.extend(vec![E::ScalarField::zero(); degree+1]);
        c_hat.extend(&coeffs[0..degree]);

        let mut s_hat = key.ref_string_g1[0..degree].to_vec();
        s_hat.reverse();
        s_hat.extend(vec![E::G1::zero(); domain_size-degree]);

        let y = domain.fft(&c_hat);
        let v = domain.fft(&s_hat);
        let u = Self::elementwise_mul(&v, &y);
        let h_hat = domain.ifft(&u);

        Ok ( h_hat[0..degree].to_vec() )
    }

    fn single_proof(
        key: &KZGKey<E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>,
        index: usize
    ) -> Result<E::G1, KZGError> {
        let data_piece = data.evaluate_by_index(index);
        let point = data.index_to_evaluation_point(index);
        let mut poly = data.data.clone();
        poly -= &DensePolynomial::<E::ScalarField>::from_coefficients_slice(&[data_piece]);

        let divisor = DensePolynomial::<E::ScalarField>::from_coefficients_slice(&[E::ScalarField::zero() - point, E::ScalarField::one()]);
        let q = &poly / &divisor;
        let commit = key.commit_g1(q.coeffs());

        Ok( commit.0 )
    }

    fn elementwise_mul<F: PrimeField, G1: Group<ScalarField = F>>(g_vec: &Vec<G1>, f_vec: &Vec<F>) -> Vec<G1> {
        let mut result: Vec<G1> = vec![];
        for (g, f) in g_vec.iter().zip(f_vec.iter()) {
            result.push(*g * *f);
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate test;
    use test::Bencher;
    use ark_bn254::{Bn254};

    type F = <Bn254 as Pairing>::ScalarField;
    type G1 = <Bn254 as Pairing>::G1;
    type G2 = <Bn254 as Pairing>::G2;
    type poly = DensePolynomial<<Bn254 as Pairing>::ScalarField>;
    type KZG = KZGAmortized<Bn254>;

    const DATA_SIZE: usize = 20;
    const MAX_CRS: usize = 32;

    fn gen_data(num: usize) -> Vec<F> {
        let mut data: Vec<F> = vec![];
        let mut rng = rand::thread_rng();
        for _i in 0..num {
            data.push(F::rand(&mut rng));
        }
        data
    }

    fn setup(n: usize, max_degree: usize) -> (KZGPreparedData<F>, KZGKey<G1,G2>) {
        let data = gen_data(n);
        let prep = KZGPreparedData::from_iter(data);
        let crs = KZG::setup(max_degree, &mut rand::thread_rng()).unwrap();

        (prep, crs)
    }

    #[test]
    fn test_interpolation() {
        let data = gen_data(10);
        let mut prepared_data = KZGPreparedData::from_iter(data.to_vec());

        for (i, d) in data.iter().enumerate() {
            assert_eq!(prepared_data.evaluate_by_index(i), d.clone());
        }
    }

    #[test]
    fn test_single_proof() {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();

        for i in 0..DATA_SIZE {
            let proof = KZG::prove(&crs, &commit, i, &data).unwrap();
            assert!(KZG::verify(&crs, &commit, &proof).unwrap());
        }
    }

    #[test]
    fn test_batch_proof() {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();

        let proofs = KZG::prove_all(&crs, &commit, &data).unwrap();
        assert!(KZG::verify_batch(&crs, &commit, &proofs).unwrap())
    }


    #[bench]
    fn bench_single_proof(b: &mut Bencher) {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();
        //let mut rng = rand::thread_rng();
        b.iter(|| KZG::prove(&crs, &commit, 0, &data));
    }

    #[bench]
    fn bench_multi_proof(b: &mut Bencher) {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();

        b.iter(|| KZG::prove_all(&crs, &commit, &data));
    }

    fn vec_to_str<T: Display>(v: &Vec<T>) -> String {
        let mut s = "[\n".to_owned();
        for e in v {
            s.push_str(&format!("\t{}\n", e));
        }
        s.push_str("]");
        s
    }
}