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
    transcript::{self, Transcript, TranscriptError},
    MultiProofQuery, PointGenerator, PointGeneratorError, VCPreparedData, VCUniversalParams,
    VectorCommitment,
};

mod ipa_point_generator;
pub use ipa_point_generator::IPAPointGenerator;

use self::ipa_point_generator::EthereumHashToCurve;

/// Precomputed barycentric elements in a domain of size `N`.
struct PrecomputedLagrange<const N: usize, F: Field> {
    // l(x) = ∏(X-x_i) AKA A(x)
    vanishing_polynomial: DensePolynomial<F>,

    // l'(x) = ∑ l(x)/(X-x_i) AKA A'(x)
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
    y: G::ScalarField,
}

pub struct IPAMultiProof<G: Group> {
    ipa: IPAProof<G>,
    d: G,
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

#[derive(Clone)]
pub struct IPAPreparedData<const N: usize, F: Field> {
    data: [F; N],

    // The highest x value that this dataset contains
    max: usize,
}

impl<const N: usize, F: PrimeField> IPAPreparedData<N, F> {
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

    /// Evaluate both inside and outside the domain
    fn evaluate(&self, precompute: &PrecomputedLagrange<N, F>, index: usize) -> F {
        match self.get(index) {
            Some(res) => *res,
            None => inner_product(
                &precompute.compute_barycentric_coefficients(F::from(index as u64)),
                &self.data,
            ),
        }
    }

    /// Compute the quotient polynomial q(x) = [f(X) - f(x_i)] / [X-x_i]
    fn divide_by_vanishing(&self, precompute: &PrecomputedLagrange<N, F>, index: usize) -> Vec<F> {
        let mut q = vec![F::zero(); N];
        let index_f = F::from(index as u64);
        let eval = self.data[index];

        for i in 0..N {
            if i == index {
                continue;
            }

            let denom = F::from(i as u64) - index_f;

            let sub = self.data[i] - eval;
            q[i] = sub * denom.inverse().unwrap();
            q[index] += sub
                * denom
                * precompute.vanishing_derivative.evaluate(&index_f)
                * precompute
                    .vanishing_derivative
                    .evaluate(&F::from(i as u64))
                    .inverse()
                    .unwrap();
        }

        q
    }
}

impl<const N: usize, F: PrimeField> VCPreparedData for IPAPreparedData<N, F> {
    type Item = F;
    type Error = IPAError;

    fn from_vec(data: Vec<Self::Item>) -> Self {
        Self::new_incremental(data)
    }

    fn get(&self, index: usize) -> Option<&Self::Item> {
        if index <= self.max {
            Some(&self.data[index])
        } else {
            None
        }
    }

    fn get_all(&self) -> Vec<(usize, Self::Item)> {
        self.data
            .into_iter()
            .enumerate()
            .filter(|(_i, x)| *x != F::zero())
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
    type MultiProof = IPAMultiProof<G>;

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
        low_level_ipa::<G, G::ScalarField, D>(
            &key.g,
            &key.q,
            &data.data,
            &b,
            commitment,
            G::ScalarField::from(index as u64),
            None,
        )
    }

    fn prove_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        indexes: Vec<usize>,
        data: &Self::PreparedData,
    ) -> Result<Self::BatchProof, Self::Error> {
        todo!()
    }

    fn prove_multiproof<'a>(
        key: &Self::UniversalParams,
        queries: &[MultiProofQuery<
            'a,
            Self::Commitment,
            Self::PreparedData,
            <Self::PreparedData as VCPreparedData>::Item,
        >],
    ) -> Result<Self::MultiProof, Self::Error> {
        let mut transcript = Transcript::<G::ScalarField, D>::new("multiproof");
        for query in queries {
            transcript.append(query.commit, "C");
            transcript.append(&query.z, "z");
            transcript.append(&query.y, "y");
        }

        let r = transcript.hash("r", true);
        let mut r_pow = G::ScalarField::one();

        let mut g = [G::ScalarField::zero(); N];
        for query in queries {
            let quotient = query.data.divide_by_vanishing(&key.precompute, query.z);
            for i in 0..N {
                g[i] += r_pow * quotient[i];
            }
            r_pow *= r;
        }

        // Commitment to g(x)
        let d = inner_product(&key.g, &g);
        transcript.append(&d, "D");

        let t = transcript.hash("t", true);
        r_pow = G::ScalarField::one();

        let mut h = [G::ScalarField::zero(); N];
        for query in queries {
            let denom = (t - G::ScalarField::from(query.z as u64))
                .inverse()
                .unwrap();
            for i in 0..N {
                h[i] += r_pow * query.data.data[i] * denom;
            }
            r_pow *= r;
        }

        let e = inner_product(&key.g, &h);
        transcript.append(&e, "E");

        let h_minus_g: Vec<G::ScalarField> = g
            .into_iter()
            .zip(h.into_iter())
            .map(|(g_i, h_i)| h_i - g_i)
            .collect();

        let multiproof_commit = e - d;
        let proof = low_level_ipa::<G, G::ScalarField, D>(
            &key.g,
            &key.q,
            &h_minus_g,
            &key.precompute.compute_barycentric_coefficients(t),
            &multiproof_commit,
            t,
            Some(transcript),
        )?;
        Ok(IPAMultiProof { ipa: proof, d: d })
    }

    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let input_point = G::ScalarField::from(index as u64);
        low_level_verify_ipa::<G, G::ScalarField, D>(
            &key.g,
            &key.q,
            &key.precompute.compute_barycentric_coefficients(input_point),
            commitment,
            input_point,
            proof,
            None,
        )
    }

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error> {
        todo!()
    }

    fn verify_multiproof<'a>(
        key: &Self::UniversalParams,
        queries: &[MultiProofQuery<
            'a,
            Self::Commitment,
            Self::PreparedData,
            <Self::PreparedData as VCPreparedData>::Item,
        >],
        proof: &Self::MultiProof,
    ) -> Result<bool, Self::Error> {
        let mut transcript = Transcript::<G::ScalarField, D>::new("multiproof");
        for query in queries {
            transcript.append(query.commit, "C");
            transcript.append(&query.z, "z");
            transcript.append(&query.y, "y");
        }

        let r = transcript.hash("r", true);
        transcript.append(&proof.d, "D");
        let t = transcript.hash("t", true);

        let mut g2_of_t = G::ScalarField::zero();
        let mut r_pow = G::ScalarField::one();
        let mut e_coeffs = HashMap::<&IPACommitment<G>, G::ScalarField>::new();

        for query in queries {
            let e_coeff = r_pow / (t - G::ScalarField::from(query.z as u64));
            e_coeffs
                .entry(query.commit)
                .and_modify(|c| *c += e_coeff)
                .or_insert(e_coeff);

            g2_of_t += e_coeff * query.y;
            r_pow *= r;
        }

        let e: G = e_coeffs.into_iter().map(|(c, coeff)| *c * coeff).sum();
        transcript.append(&e, "E");

        low_level_verify_ipa(
            &key.g,
            &key.q,
            &key.precompute.compute_barycentric_coefficients(t),
            &(e - proof.d),
            t,
            &proof.ipa,
            Some(transcript),
        )
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
            let (data_l, data_r) = split(&data);
            let (gens_l, gens_r) = split(&gens);

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

fn split<T: Clone>(a: &[T]) -> (Vec<T>, Vec<T>) {
    // TODO: Remove
    assert!(a.len() % 2 == 0);
    (a[0..a.len() / 2].to_vec(), a[a.len() / 2..].to_vec())
}

fn powers_of<T: Mul<T, Output = T> + One + Copy>(a: T, n: usize) -> Vec<T> {
    let mut res = Vec::with_capacity(n);
    let mut cur = T::one();
    res.push(cur);

    (1..n).for_each(|_| {
        cur = cur * a;
        res.push(cur);
    });

    res
}

fn low_level_ipa<G: Group<ScalarField = F>, F: PrimeField, D: HashToField<F>>(
    gens: &[G],
    q: &G,
    a: &[F],
    b: &[F],
    commitment: &IPACommitment<G>,
    input_point: F,
    prev_transcript: Option<Transcript<F, D>>,
) -> Result<IPAProof<G>, IPAError> {
    let eval = inner_product(&a, &b);

    let mut gens = gens[0..a.len()].to_vec();
    let mut data = a.clone().to_vec();
    let mut other = b.clone().to_vec();
    let mut transcript = match prev_transcript {
        Some(t) => t,
        None => Transcript::<G::ScalarField, D>::new("ipa"),
    };
    transcript.append(commitment, "C")?;
    transcript.append(&input_point, "input point")?;
    transcript.append(&eval, "output point")?;

    //let hasher = D::new(&[]);
    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();
    let mut ra = transcript.hash("w", true);

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
        ra = transcript.hash("x", true);

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

fn low_level_verify_ipa<G: Group<ScalarField = F>, F: PrimeField, D: HashToField<F>>(
    gens: &[G],
    q: &G,
    b: &[F],
    commitment: &IPACommitment<G>,
    input_point: F,
    proof: &IPAProof<G>,
    prev_transcript: Option<Transcript<F, D>>,
) -> Result<bool, IPAError> {
    let mut c = commitment.clone();
    let mut transcript = match prev_transcript {
        Some(t) => t,
        None => Transcript::<G::ScalarField, D>::new("ipa"),
    };
    transcript.append(commitment, "C")?;
    transcript.append(&input_point, "input point")?;
    transcript.append(&proof.y, "output point")?;
    let mut ra = transcript.hash("w", true);
    let mut points_coeffs = vec![G::ScalarField::one()];
    let q = *q * ra;
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

#[cfg(test)]
mod tests {
    use super::*;

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
        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();
        let data = IPAPreparedData::<SIZE, F>::new_incremental(data_raw);
        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let index = thread_rng().gen_range(0..SIZE) as usize;
        let proof = IPAT::prove(&crs, &commit, index, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, index, &proof).unwrap());

        let index_outside = SIZE * 2;
        let proof_outside = IPAT::prove(&crs, &commit, index_outside, &data).unwrap();
        assert!(IPAT::verify(&crs, &commit, index_outside, &proof_outside).unwrap());
    }

    #[test]
    fn test_multiproof() {
        let NUM_MULTIPROOF = 20;
        let point_gen = IPAPointGenerator::default();
        let crs = IPAT::setup(SIZE, &point_gen).unwrap();

        let all_data: Vec<(IPAPreparedData<SIZE, F>, IPACommitment<G>)> = (0..NUM_MULTIPROOF)
            .map(|_| {
                let data = IPAPreparedData::<SIZE, F>::new_incremental(
                    (0..SIZE).map(|_| F::rand(&mut thread_rng())).collect(),
                );
                let commit = IPAT::commit(&crs, &data).unwrap();

                (data, commit)
            })
            .collect();

        let queries: Vec<_> = all_data
            .iter()
            .map(|(data, commit)| {
                let z = thread_rng().gen_range(0..SIZE);
                MultiProofQuery {
                    commit: commit,
                    data: data,
                    z: z,
                    y: data.data[z],
                }
            })
            .collect();

        let proof = IPAT::prove_multiproof(&crs, &queries).unwrap();

        assert!(IPAT::verify_multiproof(&crs, &queries, &proof).unwrap());
    }
}
