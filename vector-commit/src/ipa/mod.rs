use std::marker::PhantomData;

use ark_ec::{CurveGroup, Group};
use ark_ff::{field_hashers::HashToField, Field, One, PrimeField, Zero};

use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use thiserror::Error;

use crate::{
    lagrange_basis::LagrangeBasis,
    precompute::PrecomputedLagrange,
    transcript::{Transcript, TranscriptError, TranscriptHasher},
    utils::*,
    HasPrecompute, PointGenerator, VCCommitment, VCUniversalParams, VectorCommitment,
};

mod ipa_point_generator;
pub use ipa_point_generator::IPAPointGenerator;

use self::ipa_point_generator::EthereumHashToCurve;

pub struct IPAUniversalParams<const N: usize, G: Group, D: HashToField<G::ScalarField>> {
    g: [G; N], // Gens to commit the evaluations of the dataset
    q: G,      // Gen to commit to the inner product of the dataset with it's b vector
    precompute: PrecomputedLagrange<G::ScalarField>,

    digest: PhantomData<D>,
}

impl<const N: usize, G: Group, D: HashToField<G::ScalarField>> IPAUniversalParams<N, G, D> {
    fn new(g: [G; N], q: G) -> Self {
        Self {
            g,
            q,
            precompute: PrecomputedLagrange::new(N),
            digest: PhantomData,
        }
    }

    fn new_from_vec(all: Vec<G>) -> Self {
        let mut real_g = [G::zero(); N];
        for i in 0..N {
            real_g[i] = all[i];
        }
        Self {
            g: real_g,
            q: all[N],
            precompute: PrecomputedLagrange::new(N),
            digest: PhantomData,
        }
    }
}

impl<const N: usize, G: Group, D: HashToField<G::ScalarField>> VCUniversalParams
    for IPAUniversalParams<N, G, D>
{
    fn max_size(&self) -> usize {
        N
    }
}

impl<const N: usize, G: Group, D: HashToField<G::ScalarField>> HasPrecompute<G::ScalarField>
    for IPAUniversalParams<N, G, D>
{
    fn precompute(&self) -> &PrecomputedLagrange<G::ScalarField> {
        &self.precompute
    }
}

/// A commitment to the set of data
pub type IPACommitment<G> = G;

pub struct IPACommitProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
}

pub struct IPAProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
    y: G::ScalarField,
}

#[derive(Error, Clone, Debug)]
pub enum IPAError {
    #[error("Attempting to use an in-domain function outside of the domain")]
    OutOfDomain,

    #[error("Attempting to create an IPA commitment greater than the CRS size")]
    OutOfCRS,

    #[error("Transcript error")]
    TranscriptError(#[from] TranscriptError),
}

pub struct IPA<const N: usize, G, H, D> {
    _g: PhantomData<G>,
    _h: PhantomData<H>,
    _d: PhantomData<D>,
}

impl<const N: usize, G, H, D> VectorCommitment for IPA<N, G, H, D>
where
    G: CurveGroup,
    H: HashToField<G::ScalarField> + Sync,
    D: EvaluationDomain<G::ScalarField>,
{
    type UniversalParams = IPAUniversalParams<N, G, H>;
    //type PreparedData = IPAPreparedData<N, G::ScalarField>;
    //type PreparedData = LagrangeBasis<G::ScalarField, GeneralEvaluationDomain<G::ScalarField>>;
    type Commitment = IPACommitment<G>;
    type Data = LagrangeBasis<G::ScalarField, D>;
    type Proof = IPAProof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = IPAError;
    type PointGenerator = IPAPointGenerator<G, EthereumHashToCurve>;
    type Transcript = TranscriptHasher<G::ScalarField, H>;

    fn setup(
        max_items: usize,
        gen: &Self::PointGenerator,
    ) -> Result<Self::UniversalParams, crate::PointGeneratorError> {
        let gens = gen.gen(max_items + 1)?;
        // TODO: Perhaps the PointGenerator should also have a generic bound on its max size
        Ok(Self::UniversalParams::new_from_vec(gens))
    }

    fn commit(
        key: &Self::UniversalParams,
        data: &LagrangeBasis<G::ScalarField, D>,
    ) -> Result<Self::Commitment, Self::Error> {
        Ok(inner_product(&key.g, data.elements_ref()))
    }

    fn prove_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: G::ScalarField,
        data: &LagrangeBasis<G::ScalarField, D>,
        transcript: Option<Self::Transcript>,
    ) -> Result<Self::Proof, Self::Error> {
        let mut b = key.precompute.compute_barycentric_coefficients(point);
        low_level_ipa::<G, G::ScalarField, Self::Transcript>(
            &key.g,
            &key.q,
            data.elements_ref(),
            &b,
            commitment,
            point,
            transcript,
        )
    }

    fn prove_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        indexes: Vec<usize>,
        data: &LagrangeBasis<G::ScalarField, D>,
    ) -> Result<Self::BatchProof, Self::Error> {
        todo!()
    }

    fn verify_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: G::ScalarField,
        proof: &Self::Proof,
        transcript: Option<Self::Transcript>,
    ) -> Result<bool, Self::Error> {
        low_level_verify_ipa::<G, G::ScalarField, Self::Transcript>(
            &key.g,
            &key.q,
            &key.precompute.compute_barycentric_coefficients(point),
            commitment,
            point,
            proof,
            transcript,
        )
    }

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error> {
        todo!()
    }
}

impl<const N: usize, G, H, D> IPA<N, G, H, D>
where
    G: CurveGroup,
    H: HashToField<G::ScalarField> + Sync,
    D: EvaluationDomain<G::ScalarField>,
{
    /// Prove that we have made a valid commitment
    pub fn prove_commitment(
        key: &IPAUniversalParams<N, G, H>,
        commitment: &IPACommitment<G>,
        data: &LagrangeBasis<G::ScalarField, GeneralEvaluationDomain<G::ScalarField>>,
    ) -> IPACommitProof<G> {
        let max = data.max();
        let mut data = data.elements_ref()[0..max + 1].to_vec();
        let mut gens = key.g[0..max + 1].to_vec();
        let mut l = Vec::<G>::new();
        let mut r = Vec::<G>::new();

        let mut transcript = <Self as VectorCommitment>::Transcript::new("ipa");
        transcript.append(commitment, "C");
        let mut ra = transcript.digest("x", true);

        while data.len() > 1 {
            let (data_l, data_r) = split(&data);
            let (gens_l, gens_r) = split(&gens);

            let y_l = inner_product(&gens_r, &data_l);
            let y_r = inner_product(&gens_l, &data_r);
            l.push(y_l);
            r.push(y_r);

            transcript.append(&y_l, "L");
            transcript.append(&y_r, "R");
            ra = transcript.digest("x", true);

            data = vec_add_and_distribute(&data_l, &data_r, ra);
            gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
        }
        IPACommitProof {
            l: l,
            r: r,
            tip: data[0],
        }
    }

    /// Verify that a commitment is valid
    pub fn verify_commitment_proof(
        key: &IPAUniversalParams<N, G, H>,
        commitment: &IPACommitment<G>,
        proof: &IPACommitProof<G>,
    ) -> bool {
        let gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
        let mut c = commitment.clone();
        let mut points_coeffs = vec![G::ScalarField::one()];
        let mut transcript = <Self as VectorCommitment>::Transcript::new("ipa");
        transcript.append(commitment, "C");
        let mut ra = transcript.digest("x", true);

        for i in 0..proof.l.len() {
            transcript.append(&proof.l[i], "L");
            transcript.append(&proof.r[i], "R");
            ra = transcript.digest("x", true);

            c = proof.l[i] + c * ra + proof.r[i] * ra.square();
            points_coeffs = points_coeffs
                .into_iter()
                .map(|x| vec![x * ra, x])
                .flatten()
                .collect();
        }

        let combined_point = inner_product(&gens, &points_coeffs);
        c == combined_point * proof.tip
    }
}

fn low_level_ipa<G: Group<ScalarField = F>, F: PrimeField, T: Transcript<F>>(
    gens: &[G],
    q: &G,
    a: &[F],
    b: &[F],
    commitment: &IPACommitment<G>,
    input_point: F,
    prev_transcript: Option<T>,
) -> Result<IPAProof<G>, IPAError> {
    let eval = inner_product(&a, &b);

    let mut gens = gens[0..a.len()].to_vec();
    let mut data = a.clone().to_vec();
    let mut other = b.clone().to_vec();
    let mut transcript = match prev_transcript {
        Some(t) => t,
        None => T::new("ipa"),
    };
    transcript.append(commitment, "C")?;
    transcript.append(&input_point, "input point")?;
    transcript.append(&eval, "output point")?;

    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();
    let mut ra = transcript.digest("w", true);

    let q = *q * ra;
    while data.len() > 1 {
        let (data_l, data_r) = split(&data);
        let (gens_l, gens_r) = split(&gens);
        let (b_l, b_r) = split(&other);
        let y_l = inner_product(&gens_r, &data_l) + q * inner_product(&data_l, &b_r);
        let y_r = inner_product(&gens_l, &data_r) + q * inner_product(&data_r, &b_l);

        l.push(y_l);
        r.push(y_r);
        transcript.append(&y_l, "L")?;
        transcript.append(&y_r, "R")?;
        ra = transcript.digest("x", true);

        data = vec_add_and_distribute(&data_l, &data_r, ra);
        gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
        other = vec_add_and_distribute(&b_r, &b_l, ra);
    }

    Ok(IPAProof {
        l,
        r,
        tip: data[0],
        y: eval,
    })
}

fn low_level_verify_ipa<G: Group<ScalarField = F>, F: PrimeField, T: Transcript<F>>(
    gens: &[G],
    q: &G,
    b: &[F],
    commitment: &IPACommitment<G>,
    input_point: F,
    proof: &IPAProof<G>,
    prev_transcript: Option<T>,
) -> Result<bool, IPAError> {
    let mut c = commitment.clone();
    let mut transcript = match prev_transcript {
        Some(t) => t,
        None => T::new("ipa"),
    };
    transcript.append(commitment, "C")?;
    transcript.append(&input_point, "input point")?;
    transcript.append(&proof.y, "output point")?;
    let mut ra = transcript.digest("w", true);
    let mut points_coeffs = vec![G::ScalarField::one()];
    let q = *q * ra;
    c += q * proof.y;

    for i in 0..proof.l.len() {
        transcript.append(&proof.l[i], "L")?;
        transcript.append(&proof.r[i], "R")?;
        ra = transcript.digest("x", true);

        c = proof.l[i] + c * ra + proof.r[i] * ra.square();
        points_coeffs = points_coeffs
            .into_iter()
            .map(|x| vec![x * ra, x])
            .flatten()
            .collect();
    }

    let combined_point = inner_product(&gens, &points_coeffs);
    let combined_b = inner_product(&b, &points_coeffs);

    Ok(c == combined_point * proof.tip + q * (proof.tip * combined_b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VCData;

    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;
    use ark_ec::Group;
    use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};

    use rand::{thread_rng, Rng};
    use sha2::Sha256;

    type F = <Bn254 as Pairing>::ScalarField;
    type G = <Bn254 as Pairing>::G1;
    type Hasher = DefaultFieldHasher<Sha256>;

    const SIZE: usize = 32;
    type IPAT = IPA<SIZE, G, Hasher, GeneralEvaluationDomain<F>>;

    #[test]
    fn test_commit_evaluations() {
        let mut gens = [G::zero(); SIZE];
        (0..SIZE).for_each(|i| gens[i] = G::generator() * F::from(i as u64 + 1));
        let q = G::generator() * F::from(SIZE as u64 + 1);

        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();
        //let data = IPAPreparedData::<SIZE, F>::new_incremental(data_raw);
        let data = LagrangeBasis::<F, GeneralEvaluationDomain<F>>::from_vec(data_raw);

        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let proof = IPAT::prove_commitment(&crs, &commit, &data);
        assert!(IPAT::verify_commitment_proof(&crs, &commit, &proof));

        commit += G::generator();
        assert!(!IPAT::verify_commitment_proof(&crs, &commit, &proof));
    }

    #[test]
    fn test_eval_proof() {
        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();
        let data = LagrangeBasis::<F, GeneralEvaluationDomain<F>>::from_vec(data_raw);
        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let index = thread_rng().gen_range(0..SIZE) as usize;
        let proof = IPAT::prove(&crs, &commit, index, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, index, &proof).unwrap());

        let index_outside = SIZE * 2;
        let proof_outside = IPAT::prove(&crs, &commit, index_outside, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, index_outside, &proof_outside).unwrap());
        assert!(!IPAT::verify(&crs, &commit, index, &proof_outside).unwrap());
    }
}
