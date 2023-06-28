use std::{marker::PhantomData, iter::Map, error::Error, fmt::Display};

use ark_ec::{Group, pairing::Pairing};
use ark_ff::{PrimeField, UniformRand};
use ark_poly::{GeneralEvaluationDomain, EvaluationDomain, Evaluations, univariate::DensePolynomial, Polynomial};
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
    size: usize
}

impl<F: PrimeField> FromIterator<F> for KZGPreparedData<F> {
    fn from_iter<T: IntoIterator<Item = F>>(iter: T) -> Self {
        let points: Vec<F> = iter.into_iter().collect();
        let len = points.len();
        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::<F>::new(points.len()).unwrap();
        let poly = Evaluations::from_vec_and_domain(points, domain).interpolate();
        
        Self {
            data: poly,
            w: domain.group_gen(),
            size: len
        }
    }
}

impl<F: PrimeField> KZGPreparedData<F> {
    
    /// Evaluate the prepared data polynomial at an `index` in [0,...,n]. This will translate it to
    /// the corresponding root of unity raised to `index`
    fn evaluate_by_index(&self, index: usize) -> F {
        // TODO: Not fully-safe as we are ignoring error of F::from_str
        let eval_index = self.w.pow(&F::from_str(&index.to_string()).ok().unwrap().into_bigint());
        self.data.evaluate(&eval_index)
    }
}

struct KZGCommitment {}

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
    type Commitment = KZGCommitment;
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
        todo!()
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
        let prepared_data = KZGPreparedData::from_iter(data.to_vec());
        println!("{:?}", prepared_data.data);

        for (i, d) in data.iter().enumerate() {
            assert_eq!(prepared_data.evaluate_by_index(i), d.clone());
        }
    }
}