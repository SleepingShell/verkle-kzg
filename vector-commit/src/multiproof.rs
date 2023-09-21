use std::{
    collections::HashMap,
    hash::Hash,
    iter::Sum,
    ops::{Mul, Sub},
};

use ark_ec::{pairing::Pairing, CurveGroup, Group};
use ark_ff::{field_hashers::HashToField, One, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;

use itertools::Itertools;
use rayon::prelude::*;

use crate::{
    ipa::IPA,
    kzg::KZG,
    lagrange_basis::LagrangeBasis,
    transcript::Transcript,
    utils::{invert_domain_at, powers_of},
    HasPrecompute, VCData, VCUniversalParams, VectorCommitment,
};

#[derive(Clone)]
pub struct MultiproofProverQuery<'a, C, D, F: Clone> {
    data: &'a D,
    commit: &'a C,
    z: usize,
    y: F,
}

impl<'a, C, D, F: Clone> MultiproofProverQuery<'a, C, D, F> {
    pub fn new(data: &'a D, commit: &'a C, z: usize, y: F) -> Self {
        Self { data, commit, z, y }
    }

    pub fn to_verifier_query(&self) -> MultiproofVerifierQuery<'a, C, F> {
        MultiproofVerifierQuery::<'a, C, F>::new(self.commit, self.z, self.y.clone())
    }
}

pub struct MultiproofVerifierQuery<'a, C, F> {
    commit: &'a C,
    z: usize,
    y: F,
}

impl<'a, C, F> MultiproofVerifierQuery<'a, C, F> {
    pub fn new(commit: &'a C, z: usize, y: F) -> Self {
        Self { commit, z, y }
    }
}

pub struct Multiproof<P, D> {
    proof: P,
    d: D,
}

pub trait VectorCommitmentMultiproof<G, D>:
    VectorCommitment<Data = LagrangeBasis<G::ScalarField, D>>
where
    G: Group,
    D: EvaluationDomain<G::ScalarField> + Sync + Send,
    <Self as VectorCommitment>::Commitment: CanonicalSerialize
        + Sub<Output = <Self as VectorCommitment>::Commitment>
        + Eq
        + Hash
        + Mul<G::ScalarField, Output = Self::Commitment>
        + Sum
        + Copy
        + Sync,
    <Self as VectorCommitment>::UniversalParams: HasPrecompute<G::ScalarField> + Sync,
{
    /// Create a multiproof that proves multiple datasets at (possibly) multiple different evaluation points
    fn prove_multiproof<'a>(
        key: &Self::UniversalParams,
        queries: &[MultiproofProverQuery<
            'a,
            Self::Commitment,
            LagrangeBasis<G::ScalarField, D>,
            G::ScalarField,
        >],
    ) -> Result<Multiproof<Self::Proof, Self::Commitment>, Self::Error> {
        let mut transcript = <Self as VectorCommitment>::Transcript::new("multiproof");
        for query in queries.iter() {
            transcript.append(query.commit, "C");
            transcript.append(&query.z, "z");
            transcript.append(&query.y, "y");
        }

        let r = transcript.digest("r", true);
        let r_pows = powers_of(r, queries.len());

        // Scale queries by their challenge
        let scaled_queries: Vec<(usize, LagrangeBasis<G::ScalarField, D>)> = queries
            .par_iter()
            .zip(r_pows.par_iter())
            .map(|(q, r)| (q.z, q.data * *r))
            .collect();

        // Group queries by their evaluation point
        let queries_by_point = scaled_queries.iter().into_group_map_by(|q| q.0);

        // Compute g(x)
        let mut g = LagrangeBasis::new_zero(key.max_size());
        let quotients: Vec<LagrangeBasis<G::ScalarField, D>> = queries_by_point
            .iter()
            .par_bridge()
            .map(|(point, queries)| {
                let mut total = LagrangeBasis::new_zero(key.max_size());
                queries.iter().for_each(|q| {
                    total += &q.1;
                });

                LagrangeBasis::from_vec_and_domain(
                    total.divide_by_vanishing(key.precompute(), *point),
                    D::new(key.max_size()).unwrap(),
                )
            })
            .collect();

        for quotient in quotients {
            g += &quotient;
        }

        // Commitment to g(x)
        let d = Self::commit(&key, &g).unwrap();
        transcript.append(&d, "D");

        // We will evaluate g(x) at the unknown-before-commit point t
        let t = transcript.digest("t", true);

        // Calculate all the t-z_i inversions at once
        let inversions = invert_domain_at::<G::ScalarField>(t, key.max_size());

        // Calculate h(x)
        let mut h = LagrangeBasis::new_zero(key.max_size());
        for (point, queries) in queries_by_point.iter() {
            for q in queries {
                h += &(&q.1 * inversions[*point]);
            }
        }

        let e = Self::commit(&key, &h).unwrap();
        transcript.append(&e, "E");

        let h_minus_g = h - g;

        let multiproof_commit = e - d.clone();
        let proof = Self::prove_point(key, &multiproof_commit, t, &h_minus_g, Some(transcript))?;
        Ok(Multiproof { proof, d })
    }

    fn verify_multiproof<'a>(
        key: &Self::UniversalParams,
        queries: &[MultiproofVerifierQuery<'a, Self::Commitment, G::ScalarField>],
        proof: &Multiproof<Self::Proof, Self::Commitment>,
    ) -> Result<bool, Self::Error> {
        let mut transcript = <Self as VectorCommitment>::Transcript::new("multiproof");
        for query in queries {
            transcript.append(query.commit, "C");
            transcript.append(&query.z, "z");
            transcript.append(&query.y, "y");
        }

        let r = transcript.digest("r", true);
        transcript.append(&proof.d, "D");
        let t = transcript.digest("t", true);

        let mut g2_of_t = G::ScalarField::zero();
        let mut r_pow = G::ScalarField::one();
        let mut e_coeffs = HashMap::<&Self::Commitment, G::ScalarField>::new();

        let inversions = invert_domain_at::<G::ScalarField>(t, key.max_size());

        for query in queries {
            let e_coeff = r_pow * inversions[query.z];
            e_coeffs
                .entry(query.commit)
                .and_modify(|c| *c += e_coeff)
                .or_insert(e_coeff);

            g2_of_t += e_coeff * query.y;
            r_pow *= r;
        }

        let e: Self::Commitment = e_coeffs.into_iter().map(|(c, coeff)| *c * coeff).sum();
        transcript.append(&e, "E");

        Self::verify_point(key, &(e - proof.d), t, &proof.proof, Some(transcript))
    }
}

impl<const N: usize, G, H, D> VectorCommitmentMultiproof<G, D> for IPA<N, G, H, D>
where
    G: CurveGroup,
    H: HashToField<G::ScalarField> + Sync,
    D: EvaluationDomain<G::ScalarField> + Sync + Send,
{
}

impl<E, H, D> VectorCommitmentMultiproof<E::G1, D> for KZG<E, H, D>
where
    E: Pairing,
    H: HashToField<E::ScalarField> + Sync,
    D: EvaluationDomain<E::ScalarField> + Sync + Send,
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        ipa::{IPACommitment, IPAPointGenerator, IPA},
        kzg::kzg_point_generator::KZGRandomPointGenerator,
    };

    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;
    use ark_ec::Group;
    use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};

    use ark_poly::GeneralEvaluationDomain;
    use rand::{thread_rng, Rng};
    use sha2::Sha256;

    type F = <Bn254 as Pairing>::ScalarField;
    type G = <Bn254 as Pairing>::G1;
    type Hasher = DefaultFieldHasher<Sha256>;

    const SIZE: usize = 32;
    type IPAT = IPA<SIZE, G, Hasher, GeneralEvaluationDomain<F>>;
    type KZGT = KZG<Bn254, Hasher, GeneralEvaluationDomain<F>>;

    #[test]
    fn test_ipa_multiproof() {
        let NUM_MULTIPROOF = 20;
        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();

        let all_data: Vec<(
            LagrangeBasis<F, GeneralEvaluationDomain<F>>,
            IPACommitment<G>,
        )> = (0..NUM_MULTIPROOF)
            .map(|_| {
                let r = F::rand(&mut thread_rng());
                let data = LagrangeBasis::<F, GeneralEvaluationDomain<F>>::from_vec(
                    (0..SIZE).map(|i| r + F::from(i as u64)).collect(),
                );
                let commit = IPAT::commit(&crs, &data).unwrap();

                (data, commit)
            })
            .collect();

        let queries: Vec<_> = all_data
            .iter()
            .map(|(data, commit)| {
                let z = thread_rng().gen_range(0..SIZE);
                MultiproofProverQuery {
                    commit: commit,
                    data: data,
                    z: z,
                    y: data[z],
                }
            })
            .collect();
        let mut verifier_queries: Vec<_> = queries.iter().map(|q| q.to_verifier_query()).collect();

        let mut proof = IPAT::prove_multiproof(&crs, &queries).unwrap();

        assert!(IPAT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        proof.d += G::generator();
        assert!(!IPAT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        proof.d -= G::generator();
        verifier_queries[0].y += F::one();
        assert!(!IPAT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        verifier_queries[0].y -= F::one();
        //proof.proof.l[0] += G::generator();
        //assert!(!IPAT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        //proof.proof.l[0] -= G::generator();
    }

    #[test]
    fn test_kzg_multiproof() {
        let NUM_MULTIPROOF = 20;
        let point_gen = KZGRandomPointGenerator::default();
        let crs = KZGT::setup(SIZE, &point_gen).unwrap();

        let all_data: Vec<(
            LagrangeBasis<F, GeneralEvaluationDomain<F>>,
            IPACommitment<G>,
        )> = (0..NUM_MULTIPROOF)
            .map(|_| {
                let r = F::rand(&mut thread_rng());
                let data = LagrangeBasis::<F, GeneralEvaluationDomain<F>>::from_vec(
                    (0..SIZE).map(|i| r + F::from(i as u64)).collect(),
                );
                let commit = KZGT::commit(&crs, &data).unwrap();

                (data, commit)
            })
            .collect();

        let queries: Vec<_> = all_data
            .iter()
            .map(|(data, commit)| {
                let z = thread_rng().gen_range(0..SIZE);
                MultiproofProverQuery {
                    commit: commit,
                    data: data,
                    z: z,
                    y: data[z],
                }
            })
            .collect();
        let mut verifier_queries: Vec<_> = queries.iter().map(|q| q.to_verifier_query()).collect();

        let mut proof = KZGT::prove_multiproof(&crs, &queries).unwrap();

        assert!(KZGT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        proof.d += G::generator();
        assert!(!KZGT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        proof.d -= G::generator();
        verifier_queries[0].y += F::one();
        assert!(!KZGT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        verifier_queries[0].y -= F::one();
        //proof.proof.l[0] += G::generator();
        //assert!(!IPAT::verify_multiproof(&crs, &verifier_queries, &proof).unwrap());
        //proof.proof.l[0] -= G::generator();
    }
}
