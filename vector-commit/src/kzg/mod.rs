use std::marker::PhantomData;

use ark_ec::{pairing::Pairing, Group};
use ark_ff::{field_hashers::HashToField, FftField, Field, One, PrimeField, Zero};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain};
use thiserror::Error;

use crate::{
    precompute::PrecomputedLagrange,
    transcript::TranscriptHasher,
    utils::{inner_product, to_usize},
    PointGenerator, VCCommitment, VCUniversalParams, VectorCommitment,
};

use self::kzg_point_generator::KZGRandomPointGenerator;

pub mod kzg_point_generator;

pub type KZGCommitment<G: Group> = G;

/// KZGKey represents the universal parameters, AKA reference string, for both
/// committing polynomials and verifying commitments
#[derive(Clone, Debug)]
pub struct KZGKey<F: FftField, G1: Group, G2: Group> {
    /// The max number of elements this reference string supports
    size: usize,

    /// TEMP
    g1: Vec<G1>,

    /// The corresponding `PointGenerator` should commit directly to the lagrange polynomials
    /// as we work in evaluation form.
    lagrange_commitments: Vec<G1>,

    /// For G2, we only need α*g
    g2: G2,

    precompute: PrecomputedLagrange<F>,
}

impl<F, G1, G2> KZGKey<F, G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>,
{
    fn from_lagrange_vec(g1: Vec<G1>, lagrange_g1: Vec<G1>, g2: G2, unity: F) -> Self {
        let size = g1.len();
        Self {
            size,
            g1,
            lagrange_commitments: lagrange_g1,
            g2,
            precompute: PrecomputedLagrange::new(size),
        }
    }
}

impl<F, G1, G2> VCUniversalParams<F> for KZGKey<F, G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>,
{
    fn max_size(&self) -> usize {
        self.size
    }

    fn precompute(&self) -> &crate::precompute::PrecomputedLagrange<F> {
        &self.precompute
    }
}

pub struct KZGProof<F: Field, G: Group> {
    proof: KZGCommitment<G>,
    y: F,
}

#[derive(Error, Clone, Debug)]
pub enum KZGError {
    #[error("An unspecified error occurred")]
    DefaultError,
    //InvalidDomain,
    //OutOfDomainBounds,
}

/// Implementation of the Feist-Khovratovich technique of "Fast Amortized KZG proofs".
#[derive(PartialEq, Clone)]
pub struct KZG<E, H, D> {
    _engine: PhantomData<E>,
    _hasher: PhantomData<H>,
    _domain: PhantomData<D>,
}

impl<E: Pairing, D: EvaluationDomain<E::ScalarField>, H: HashToField<E::ScalarField>>
    VectorCommitment<E::ScalarField, D> for KZG<E, H, D>
{
    type UniversalParams = KZGKey<E::ScalarField, E::G1, E::G2>;
    type Commitment = KZGCommitment<E::G1>;
    type Proof = KZGProof<E::ScalarField, E::G1>;
    type BatchProof = Vec<E::G1>;
    type Error = KZGError;
    type PointGenerator = KZGRandomPointGenerator<E::G1>;
    type Transcript = TranscriptHasher<E::ScalarField, H>;

    fn setup(
        max_items: usize,
        gen: &Self::PointGenerator,
    ) -> Result<Self::UniversalParams, crate::PointGeneratorError> {
        let g1_points = gen.gen(max_items)?;
        let domain = D::new(max_items).unwrap();
        let points = domain.ifft(&g1_points);
        let g2 = E::G2::generator() * gen.secret().unwrap();
        Ok(KZGKey::from_lagrange_vec(
            g1_points,
            points,
            g2,
            domain.group_gen(),
        ))
    }

    fn commit(
        key: &Self::UniversalParams,
        data: &crate::lagrange_basis::LagrangeBasis<E::ScalarField, D>,
    ) -> Result<Self::Commitment, Self::Error> {
        Ok(inner_product(
            &key.lagrange_commitments,
            data.elements_ref(),
        ))
    }

    fn prove_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: E::ScalarField,
        data: &crate::lagrange_basis::LagrangeBasis<E::ScalarField, D>,
        transcript: Option<Self::Transcript>,
    ) -> Result<Self::Proof, Self::Error> {
        let evaluation = data.evaluate(key.precompute(), point);
        let q = if point <= E::ScalarField::from(key.max_size() as u64) {
            data.divide_by_vanishing(key.precompute(), to_usize(&point))
        } else {
            data.divive_by_vanishing_outside_domain(key.precompute(), point)
        };

        Ok(KZGProof {
            proof: inner_product(&key.lagrange_commitments, &q),
            y: evaluation,
        })
    }

    fn prove_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        indexes: Vec<usize>,
        data: &crate::lagrange_basis::LagrangeBasis<E::ScalarField, D>,
    ) -> Result<Self::BatchProof, Self::Error> {
        todo!()
    }

    fn verify_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: E::ScalarField,
        proof: &Self::Proof,
        transcript: Option<Self::Transcript>,
    ) -> Result<bool, Self::Error> {
        let p = if point < E::ScalarField::from(key.max_size() as u64) {
            //key.precompute().unity().pow(&[to_usize(point) as u64])
            key.precompute()
                .domain()
                .group_gen()
                .pow(&[to_usize(&point) as u64])
        } else {
            point
        };

        let pairing1 = E::pairing(proof.proof, key.g2 - (E::G2::generator() * p));
        let pairing2 = E::pairing(
            *commitment - (E::G1::generator() * proof.y),
            E::G2::generator(),
        );

        Ok(pairing1 == pairing2)
    }

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::lagrange_basis::LagrangeBasis;

    use super::*;
    use ark_bn254::Bn254;
    use ark_ff::{field_hashers::DefaultFieldHasher, PrimeField, UniformRand};
    use ark_poly::GeneralEvaluationDomain;
    use sha2::Sha256;

    type Hasher = DefaultFieldHasher<Sha256>;

    type F = <Bn254 as Pairing>::ScalarField;
    type G1 = <Bn254 as Pairing>::G1;
    type G2 = <Bn254 as Pairing>::G2;
    type D = GeneralEvaluationDomain<F>;

    type TKZG = KZG<Bn254, Hasher, GeneralEvaluationDomain<F>>;

    const DATA_SIZE: usize = 8;
    const MAX_CRS: usize = 16;

    fn gen_data(num: usize) -> Vec<F> {
        let mut data: Vec<F> = vec![];
        let mut rng = rand::thread_rng();
        for _i in 0..num {
            data.push(F::rand(&mut rng));
        }
        data
    }

    fn setup(n: usize, max_degree: usize) -> (LagrangeBasis<F, D>, KZGKey<F, G1, G2>) {
        let data = gen_data(n);
        let point_gen = KZGRandomPointGenerator::<G1>::default();

        let crs = TKZG::setup(max_degree, &point_gen).unwrap();
        let prep = LagrangeBasis::from_vec_and_domain(data, *crs.precompute().domain());

        (prep, crs)
    }

    #[test]
    fn test_single_proof() {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = TKZG::commit(&crs, &data).unwrap();

        for i in 0..DATA_SIZE {
            let proof = TKZG::prove(&crs, &commit, i, &data).unwrap();
            assert!(TKZG::verify(&crs, &commit, i, &proof).unwrap());
        }

        for i in DATA_SIZE..MAX_CRS {
            let proof = TKZG::prove(&crs, &commit, i, &data).unwrap();
            assert!(TKZG::verify(&crs, &commit, i, &proof).unwrap());
            assert!(proof.y == F::zero());
        }

        let outside_index = MAX_CRS + 1;
        let outside_proof = TKZG::prove(&crs, &commit, outside_index, &data).unwrap();
        assert!(TKZG::verify(&crs, &commit, outside_index, &outside_proof).unwrap());
    }
}
