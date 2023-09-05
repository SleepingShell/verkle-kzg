use std::{
    collections::HashMap,
    error::Error,
    fmt::Display,
    iter::Sum,
    marker::PhantomData,
    ops::{Add, Mul},
};

use ark_ec::{CurveGroup, Group};
use ark_ff::{
    batch_inversion_and_mul, field_hashers::HashToField, BigInteger, FftField, Field, One,
    PrimeField, Zero,
};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_serialize::CanonicalSerialize;
use thiserror::Error;

use crate::{
    transcript::{Transcript, TranscriptError},
    PointGenerator, PointGeneratorError, VCPreparedData, VCUniversalParams, VectorCommitment,
};

mod ipa_point_generator;
pub use ipa_point_generator::IPAPointGenerator;

use self::ipa_point_generator::EthereumHashToCurve;

/// Precomputed barycentric elements in a domain of size `N`.
struct PrecomputedLagrange<const N: usize, F: Field> {
    // l(x) = ∏(X-x_i)
    vanishing_polynomial: DensePolynomial<F>,

    // l'(x) = ∑ l(x)/(X-x_i)
    vanishing_derivative: DensePolynomial<F>,

    // w_j = ∏_{m!=j}(x_j-x_m)^-1
    barycentric_weights: [F; N],
}

impl<const N: usize, F: PrimeField> PrecomputedLagrange<N, F> {
    fn new() -> Self {
        let domain = (0..N as u64).map(|x| F::from(x)).collect::<Vec<F>>();
        let vanishing = Self::compute_vanishing_polynomial(&domain);
        let derivative = Self::compute_vanishing_derivative(&domain, &vanishing);
        let weights = Self::compute_barycentric_weights(domain);

        Self {
            vanishing_polynomial: vanishing,
            vanishing_derivative: derivative,
            barycentric_weights: weights,
        }
    }

    fn compute_vanishing_polynomial(x_points: &Vec<F>) -> DensePolynomial<F> {
        let mut res = DensePolynomial::<F>::from_coefficients_slice(&[F::one()]);
        x_points.into_iter().for_each(|x| {
            res = &res * &DensePolynomial::<F>::from_coefficients_slice(&[F::zero() - x, F::one()]);
        });

        res
    }

    fn compute_vanishing_derivative(
        x_points: &Vec<F>,
        vanishing_roots: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let mut res = DensePolynomial::<F>::from_coefficients_slice(&[F::zero()]);
        x_points.into_iter().for_each(|x| {
            res += &(vanishing_roots
                / &DensePolynomial::<F>::from_coefficients_slice(&[F::zero() - x, F::one()]));
        });

        res
    }

    fn compute_barycentric_weights(points: Vec<F>) -> [F; N] {
        let mut weights = [F::zero(); N];
        points.iter().enumerate().for_each(|(i, x)| {
            let mut sum = F::one();
            points.iter().for_each(|x_inner| {
                if (x != x_inner) {
                    sum *= *x - *x_inner;
                }
            });
            weights[i] = sum.inverse().unwrap();
        });

        weights
    }

    /// Computes the b vector in IPA. When this vector is inner product'd by the evaluations in the domain,
    /// the result is the evaluation F(point).
    ///
    /// b_i = l(point) / l'(x_i)(z-x_i)
    fn compute_barycentric_coefficients(&self, point: F) -> Vec<F> {
        if point < F::from(N as u64) {
            let mut res = vec![F::zero(); N];
            //FIXME: THIS IS SO BAD OH MY GOD
            let mut point_usize = 0usize;
            point
                .into_bigint()
                .to_bytes_be()
                .into_iter()
                .for_each(|b| point_usize = (point_usize << 8) | b as usize);
            res[point_usize] = F::one();
            return res;
        }
        let vanishing_eval = self.vanishing_polynomial.evaluate(&point);
        let mut coeffs = Vec::new();
        (0..N as u64).into_iter().for_each(|i| {
            let f_i = F::from(i);
            coeffs.push(self.vanishing_derivative.evaluate(&f_i) * (point - f_i));
        });

        batch_inversion_and_mul(&mut coeffs, &vanishing_eval);

        coeffs
    }
}

pub struct IPAUniversalParams<const N: usize, G: Group, D: HashToField<G::ScalarField>> {
    g: [G; N], // Gens to commit the evaluations of the dataset
    q: G,      // Gen to commit to the inner product of the dataset with it's b vector
    precompute: PrecomputedLagrange<N, G::ScalarField>,

    digest: PhantomData<D>,
}

impl<const N: usize, G: Group, D: HashToField<G::ScalarField>> IPAUniversalParams<N, G, D> {
    fn new(g: [G; N], q: G) -> Self {
        Self {
            g,
            q,
            precompute: PrecomputedLagrange::new(),
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
            precompute: PrecomputedLagrange::new(),
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

/// A commitment to the set of data
pub type IPACommitment<G: Group> = G;

pub struct IPACommitProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
}

pub struct IPAProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
    x: usize,
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

pub struct IPAPreparedData<const N: usize, F: Field> {
    data: [F; N],

    // The highest x value that this dataset contains
    max: usize,
}

impl<const N: usize, F: FftField> IPAPreparedData<N, F> {
    /// An ORDERED list of evaluation points
    pub fn new(data_list: Vec<(usize, F)>) -> Self {
        let max = data_list[data_list.len() - 1].0;
        let mut data = [F::zero(); N];
        data_list.into_iter().for_each(|(i, x)| data[i] = x);
        Self { data, max }
    }

    /// The data domain is simply at 0..length-1
    pub fn new_incremental(data_list: Vec<F>) -> Self {
        let max = data_list.len() - 1;
        let mut data = [F::zero(); N];
        data_list
            .into_iter()
            .enumerate()
            .for_each(|(i, x)| data[i] = x);
        Self { data, max }
    }
}

impl<const N: usize, F: FftField> VCPreparedData for IPAPreparedData<N, F> {
    type Item = F;
    type Error = IPAError;

    fn from_vec(data: Vec<Self::Item>) -> Self {
        Self::new_incremental(data)
    }

    fn get(&self, index: usize) -> Option<&Self::Item> {
        if index < self.max {
            Some(&self.data[index])
        } else {
            None
        }
    }

    fn get_all(&self) -> Vec<(usize, Self::Item)> {
        self.data
            .into_iter()
            .enumerate()
            .filter(|(i, x)| *x != F::zero())
            .collect()
    }

    fn max_size(&self) -> usize {
        N
    }

    fn set_evaluation(&mut self, index: usize, value: Self::Item) -> Result<(), Self::Error> {
        if index > N {
            Err(IPAError::OutOfDomain)
        } else {
            self.data[index] = value;
            Ok(())
        }
    }
}

pub struct IPA<const N: usize, G, D> {
    _g: PhantomData<G>,
    _d: PhantomData<D>,
}

impl<const N: usize, G, D> VectorCommitment for IPA<N, G, D>
where
    G: CurveGroup,
    D: HashToField<G::ScalarField>,
{
    type UniversalParams = IPAUniversalParams<N, G, D>;
    type PreparedData = IPAPreparedData<N, G::ScalarField>;
    type Commitment = IPACommitment<G>;
    type Proof = IPAProof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = IPAError;
    type PointGenerator = IPAPointGenerator<G, EthereumHashToCurve>;

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
        data: &Self::PreparedData,
    ) -> Result<Self::Commitment, Self::Error> {
        Ok(inner_product(&key.g, &data.data))
    }

    fn convert_commitment_to_data(
        commit: &Self::Commitment,
    ) -> <Self::PreparedData as VCPreparedData>::Item {
        if commit.is_zero() {
            G::ScalarField::ZERO
        } else {
            let mut bytes = Vec::new();
            //TODO: Check
            commit.serialize_compressed(&mut bytes).unwrap();
            G::ScalarField::from_le_bytes_mod_order(&bytes)
        }
    }

    fn prove(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        data: &Self::PreparedData,
    ) -> Result<Self::Proof, Self::Error> {
        let mut b = key
            .precompute
            .compute_barycentric_coefficients(G::ScalarField::from(index as u64));
        let eval = if index <= data.max {
            data.data[index]
        } else {
            inner_product(&data.data, &b)
        };

        let mut gens = key.g[0..data.max + 1].to_vec();
        let mut data = data.data.to_vec();
        let mut transcript = Transcript::<G::ScalarField, D>::new("ipa");
        transcript.append(commitment, "C")?;
        transcript.append(&G::ScalarField::from(index as u64), "input point")?;
        transcript.append(&eval, "output point")?;

        //let hasher = D::new(&[]);
        let mut l: Vec<G> = Vec::new();
        let mut r: Vec<G> = Vec::new();
        let mut ra = transcript.hash("w", true);

        let q = key.q * ra;
        while data.len() > 1 {
            let (data_l, data_r) = split(data);
            let (gens_l, gens_r) = split(gens);
            let (b_l, b_r) = split(b);
            let y_l = inner_product(&gens_r, &data_l) + q * inner_product(&data_l, &b_r);
            let y_r = inner_product(&gens_l, &data_r) + q * inner_product(&data_r, &b_l);

            l.push(y_l);
            r.push(y_r);
            transcript.append(&y_l, "L")?;
            transcript.append(&y_r, "R")?;
            ra = transcript.hash("x", true);

            data = vec_add_and_distribute(&data_l, &data_r, ra);
            gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
            b = vec_add_and_distribute(&b_r, &b_l, ra);
        }

        Ok(IPAProof {
            l,
            r,
            tip: data[0],
            x: index,
            y: eval,
        })
    }

    fn prove_all(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        data: &Self::PreparedData,
    ) -> Result<Self::BatchProof, Self::Error> {
        todo!()
    }

    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
        let mut c = commitment.clone();
        let mut b = key
            .precompute
            .compute_barycentric_coefficients(G::ScalarField::from(proof.x as u64));
        let mut transcript = Transcript::<G::ScalarField, D>::new("ipa");
        transcript.append(commitment, "C")?;
        transcript.append(&G::ScalarField::from(proof.x as u64), "input point")?;
        transcript.append(&proof.y, "output point")?;
        let mut ra = transcript.hash("w", true);
        let mut points_coeffs = vec![G::ScalarField::one()];
        let q = key.q * ra;
        c += q * proof.y;

        for i in 0..proof.l.len() {
            transcript.append(&proof.l[i], "L")?;
            transcript.append(&proof.r[i], "R")?;
            ra = transcript.hash("x", true);

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

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error> {
        todo!()
    }
}
impl<const N: usize, G, D> IPA<N, G, D>
where
    G: Group,
    D: HashToField<G::ScalarField>,
{
    /// Prove that we have made a valid commitment
    pub fn prove_commitment(
        key: &IPAUniversalParams<N, G, D>,
        commitment: &IPACommitment<G>,
        data: &IPAPreparedData<N, G::ScalarField>,
    ) -> IPACommitProof<G> {
        let max = data.max;
        let mut data = data.data[0..max + 1].to_vec();
        let mut gens = key.g[0..max + 1].to_vec();
        let mut l = Vec::<G>::new();
        let mut r = Vec::<G>::new();

        let mut transcript = Transcript::<G::ScalarField, D>::new("ipa");
        transcript.append(commitment, "C");
        let mut ra = transcript.hash("x", true);

        while data.len() > 1 {
            let (data_l, data_r) = split(data);
            let (gens_l, gens_r) = split(gens);

            let y_l = inner_product(&gens_r, &data_l);
            let y_r = inner_product(&gens_l, &data_r);
            l.push(y_l);
            r.push(y_r);

            transcript.append(&y_l, "L");
            transcript.append(&y_r, "R");
            ra = transcript.hash("x", true);

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
        key: &IPAUniversalParams<N, G, D>,
        commitment: &IPACommitment<G>,
        proof: &IPACommitProof<G>,
    ) -> bool {
        let gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
        let mut c = commitment.clone();
        let mut points_coeffs = vec![G::ScalarField::one()];
        let mut transcript = Transcript::<G::ScalarField, D>::new("ipa");
        transcript.append(commitment, "C");
        let mut ra = transcript.hash("x", true);

        for i in 0..proof.l.len() {
            transcript.append(&proof.l[i], "L");
            transcript.append(&proof.r[i], "R");
            ra = transcript.hash("x", true);

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

fn serialize<T: CanonicalSerialize>(x: &T) -> Vec<u8> {
    let mut b = Vec::new();
    x.serialize_compressed(&mut b);

    b
}

fn inner_product<R: Copy, T: Mul<R, Output = T> + Sum<T> + Copy>(a: &[T], b: &[R]) -> T {
    //a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
    b.iter().enumerate().map(|(i, r)| a[i] * *r).sum()
}

//res_i = a_i + x*b_i
fn vec_add_and_distribute<R: Copy, T: Copy + Add<T, Output = T> + Mul<R, Output = T>>(
    a: &[T],
    b: &[T],
    x: R,
) -> Vec<T> {
    //TODO: Remove
    assert!(a.len() == b.len());
    a.iter().zip(b.iter()).map(|(a, b)| *a + (*b * x)).collect()
}

fn split<T: Clone>(a: Vec<T>) -> (Vec<T>, Vec<T>) {
    // TODO: Remove
    assert!(a.len() % 2 == 0);
    (a[0..a.len() / 2].to_vec(), a[a.len() / 2..].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;
    use ark_ec::Group;
    use ark_ff::field_hashers::DefaultFieldHasher;

    use rand::{thread_rng, Rng};
    use sha2::Sha256;

    type F = <Bn254 as Pairing>::ScalarField;
    type G = <Bn254 as Pairing>::G1;
    type Hasher = DefaultFieldHasher<Sha256>;

    const SIZE: usize = 32;
    type IPAT = IPA<SIZE, G, Hasher>;

    #[test]
    fn test_commit_evaluations() {
        let mut gens = [G::zero(); SIZE];
        (0..SIZE).for_each(|i| gens[i] = G::generator() * F::from(i as u64 + 1));
        let q = G::generator() * F::from(SIZE as u64 + 1);

        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();
        let data = IPAPreparedData::<SIZE, F>::new_incremental(data_raw);

        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let proof = IPAT::prove_commitment(&crs, &commit, &data);
        assert!(IPAT::verify_commitment_proof(&crs, &commit, &proof));

        commit += G::generator();
        assert!(!IPAT::verify_commitment_proof(&crs, &commit, &proof));
    }

    #[test]
    fn test_eval_proof() {
        let mut gens = [G::zero(); SIZE];
        (0..SIZE).for_each(|i| gens[i] = G::generator() * F::from(i as u64 + 1));
        let q = G::generator() * F::from(SIZE as u64 + 1);

        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();
        let data = IPAPreparedData::<SIZE, F>::new_incremental(data_raw);
        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let index = thread_rng().gen_range(0..SIZE) as usize;
        let proof = IPAT::prove(&crs, &commit, index, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, &proof).unwrap());

        let index_outside = SIZE * 2;
        let proof_outside = IPAT::prove(&crs, &commit, index_outside, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, &proof_outside).unwrap());
    }
}
