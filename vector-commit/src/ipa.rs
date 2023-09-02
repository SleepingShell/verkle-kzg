use std::{
    collections::HashMap,
    error::Error,
    fmt::Display,
    iter::Sum,
    marker::PhantomData,
    ops::{Add, Mul},
};

use ark_ec::Group;
use ark_ff::{
    batch_inversion_and_mul, field_hashers::HashToField, FftField, Field, One, PrimeField, Zero,
};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_serialize::CanonicalSerialize;

use crate::{VCPreparedData, VCUniversalParams, VectorCommitment};

/// Precomputed barycentric elements in a domain of size `N`.
struct PrecomputedLagrange<const N: usize, F: Field> {
    // l(x) = ∏(X-x_i)
    vanishing_polynomial: DensePolynomial<F>,

    // ∑ l(x)/(X-x_i)
    vanishing_derivative: DensePolynomial<F>,

    // w_j = ∏_{m!=j}(x_j-x_m)^-1
    barycentric_weights: [F; N],
}

impl<const N: usize, F: FftField> PrecomputedLagrange<N, F> {
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
    /// the result is F(point)
    fn compute_barycentric_coefficients(&self, point: F) -> Vec<F> {
        if point < F::from(N as u64) {
            let mut res = vec![F::zero(); N];
            //res[point as usize] = F::one();
            return res;
        }
        let vanishing_eval = self.vanishing_polynomial.evaluate(&point);
        //let mut res_poly = DensePolynomial::<F>::from_coefficients_slice(&[F::zero()]);
        let mut coeffs = Vec::new();
        (0..N as u64).into_iter().for_each(|i| {
            let f_i = F::from(i);
            coeffs.push(self.vanishing_derivative.evaluate(&f_i) * (point - f_i));
        });

        batch_inversion_and_mul(&mut coeffs, &vanishing_eval);

        coeffs
    }
}

struct IPAUniversalParams<const N: usize, G: Group, D: HashToField<G::ScalarField>> {
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
}

impl<const N: usize, G: Group, D: HashToField<G::ScalarField>> VCUniversalParams
    for IPAUniversalParams<N, G, D>
{
    fn max_size(&self) -> usize {
        N
    }
}

/// A commitment to the set of data
type IPACommitment<G: Group> = G;

struct IPACommitProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
}

struct IPAProof<G: Group> {
    l: Vec<G>,
    r: Vec<G>,
    tip: G::ScalarField,
    x: usize,
    y: G::ScalarField,
}

#[derive(Clone, Debug)]
enum IPAError {
    OutOfDomain,
    OutOfCRS,
}

impl Display for IPAError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", "")
    }
}

impl Error for IPAError {}

struct IPAPreparedData<const N: usize, F: Field> {
    data: [F; N],

    // The highest x value that this dataset contains
    max: usize,
}

impl<const N: usize, F: FftField> IPAPreparedData<N, F> {
    /// An ORDERED list of evaluation points
    fn new(data_list: Vec<(usize, F)>) -> Self {
        let max = data_list[data_list.len() - 1].0;
        let mut data = [F::zero(); N];
        data_list.into_iter().for_each(|(i, x)| data[i] = x);
        Self { data, max }
    }

    /// The data domain is simply at 0..length-1
    fn new_incremental(data_list: Vec<F>) -> Self {
        let max = data_list.len() - 1;
        let mut data = [F::zero(); N];
        data_list
            .into_iter()
            .enumerate()
            .for_each(|(i, x)| data[i] = x);
        Self { data, max }
    }

    /*
    fn compute_b_vector(&self, index: usize, precompute: &PrecomputedLagrange<N, F>) -> Vec<F> {
        let index_f = F::from(index as u64);
        let universal_eval = precompute.vanishing_polynomial.evaluate(&index_f);
        let mut res = Vec::new();
        self.data.iter().for_each(|(x, _y)| {
            if index == *x {
                res.push(F::one());
            } else {
                let x_f = F::from(*x as u64);
                let rhs = index_f - x_f;
                let deriv_eval: F = precompute.vanishing_derivative.evaluate(&x_f);
                res.push(universal_eval / (deriv_eval * rhs));
            }
        });

        res
    }
    */
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

struct IPA<const N: usize, G, D> {
    _g: PhantomData<G>,
    _d: PhantomData<D>,
}

impl<const N: usize, G, D> VectorCommitment for IPA<N, G, D>
where
    G: Group,
    D: HashToField<G::ScalarField>,
{
    type UniversalParams = IPAUniversalParams<N, G, D>;
    type PreparedData = IPAPreparedData<N, G::ScalarField>;
    type Commitment = IPACommitment<G>;
    type Proof = IPAProof<G>;
    type BatchProof = Vec<Self::Proof>;
    type Error = IPAError;

    fn setup<R: rand::RngCore>(
        max_items: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        todo!()
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
        todo!()
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
        todo!()
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

        let hasher = D::new(&[]);
        let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];

        while data.len() > 1 {
            let (data_l, data_r) = split(data);
            let (gens_l, gens_r) = split(gens);

            let y_l = inner_product(&gens_r, &data_l);
            let y_r = inner_product(&gens_l, &data_r);
            l.push(y_l);
            r.push(y_r);

            ra = hasher.hash_to_field(
                &[serialize(&ra), serialize(&y_l), serialize(&y_r)].concat(),
                1,
            )[0];

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
        let hasher = D::new(&[]);
        let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];
        let mut points_coeffs = vec![G::ScalarField::one()];

        for i in 0..proof.l.len() {
            ra = hasher.hash_to_field(
                &[
                    serialize(&ra),
                    serialize(&proof.l[i]),
                    serialize(&proof.r[i]),
                ]
                .concat(),
                1,
            )[0];

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
/*
/// Commit to the dataset
fn commit_data<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    data: &IPAPreparedData<G::ScalarField>,
) -> IPACommitment<G> {
    let max = data.max;
    let gens = key.g[0..max + 1].to_vec();

    data.data.iter().map(|(i, v)| gens[*i] * v).sum()
}

/// Generate a proof that the dataset is committed to a valid polynomial
fn prove_data_commit<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    data: &IPAPreparedData<G::ScalarField>,
) -> IPACommitProof<G> {
    let max = data.max;
    let mut data = data.evaluations();
    let mut gens = key.g[0..max + 1].to_vec();
    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();

    let hasher = D::new(&[]);
    let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];

    while data.len() > 1 {
        let (data_l, data_r) = split(data);
        let (gens_l, gens_r) = split(gens);

        let y_l = inner_product(&gens_r, &data_l);
        let y_r = inner_product(&gens_l, &data_r);
        l.push(y_l);
        r.push(y_r);

        ra = hasher.hash_to_field(
            &[serialize(&ra), serialize(&y_l), serialize(&y_r)].concat(),
            1,
        )[0];

        data = vec_add_and_distribute(&data_l, &data_r, ra);
        gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
    }

    IPACommitProof {
        l: l,
        r: r,
        tip: data[0],
    }
}

fn verify_data_commit<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    proof: &IPACommitProof<G>,
) -> bool {
    let gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
    let mut c = commitment.clone();
    let hasher = D::new(&[]);
    let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];
    let mut points_coeffs = vec![G::ScalarField::one()];

    for i in 0..proof.l.len() {
        ra = hasher.hash_to_field(
            &[
                serialize(&ra),
                serialize(&proof.l[i]),
                serialize(&proof.r[i]),
            ]
            .concat(),
            1,
        )[0];

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

fn prove_eval<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    index: usize,
    data: &IPAPreparedData<G::ScalarField>,
) -> IPAProof<G> {
    let mut x_pows = data.compute_b_vector(index);
    let mut gens = key.g[0..data.max + 1].to_vec();
    let mut data = data.expanded_data();
    let y = data[index];

    let hasher = D::new(&[]);

    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();
    let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];

    let q = key.q * ra;
    while data.len() > 1 {
        let (data_l, data_r) = split(data);
        let (gens_l, gens_r) = split(gens);
        let (x_pows_l, x_pows_r) = split(x_pows);
        let y_l = inner_product(&gens_r, &data_l) + q * inner_product(&data_l, &x_pows_r);
        let y_r = inner_product(&gens_l, &data_r) + q * inner_product(&data_r, &x_pows_l);

        l.push(y_l);
        r.push(y_r);
        ra = hasher.hash_to_field(
            &[serialize(&ra), serialize(&y_l), serialize(&y_r)].concat(),
            1,
        )[0];

        data = vec_add_and_distribute(&data_l, &data_r, ra);
        gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
        x_pows = vec_add_and_distribute(&x_pows_r, &x_pows_l, ra);
    }

    IPAProof {
        l,
        r,
        tip: data[0],
        x: index,
        y,
    }
}

fn verify_eval<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    proof: &IPAProof<G>,
) -> bool {
    let gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
    let mut c = commitment.clone();
    let mut x_pows = vec![G::ScalarField::zero(); gens.len()];
    x_pows[proof.x] = G::ScalarField::one();
    let hasher = D::new(&[]);
    let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];
    let mut points_coeffs = vec![G::ScalarField::one()];
    let q = key.q * ra;
    c += q * proof.y;

    for i in 0..proof.l.len() {
        ra = hasher.hash_to_field(
            &[
                serialize(&ra),
                serialize(&proof.l[i]),
                serialize(&proof.r[i]),
            ]
            .concat(),
            1,
        )[0];

        c = proof.l[i] + c * ra + proof.r[i] * ra.square();
        points_coeffs = points_coeffs
            .into_iter()
            .map(|x| vec![x * ra, x])
            .flatten()
            .collect() //let coeffs: Vec<F> = (0..size).map(|i| F::from(i)).collect();ct();
    }

    let combined_point = inner_product(&gens, &points_coeffs);
    let combined_x_pow = inner_product(&x_pows, &points_coeffs);

    c == combined_point * proof.tip + q * (proof.tip * combined_x_pow)
}
*/

fn serialize<T: CanonicalSerialize>(x: &T) -> Vec<u8> {
    let mut b = Vec::new();
    x.serialize_compressed(&mut b);

    b
}

fn inner_product<R: Copy, T: Mul<R, Output = T> + Sum<T> + Copy>(a: &[T], b: &[R]) -> T {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
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
    fn test_commit_poly() {
        let mut gens = [G::zero(); SIZE];
        (0..SIZE).for_each(|i| gens[i] = G::generator() * F::from(i as u64 + 1));
        let q = G::generator() * F::from(SIZE as u64 + 1);

        let data_raw: Vec<F> = (0..SIZE as u64).map(|i| F::from(i)).collect();

        let crs = IPAUniversalParams::<SIZE, G, Hasher>::new(gens, q);
        let data = IPAPreparedData::<SIZE, F>::new_incremental(data_raw);

        let mut commit = IPAT::commit(&crs, &data).unwrap();

        let proof = IPAT::prove_commitment(&crs, &commit, &data);
        assert!(IPAT::verify_commitment_proof(&crs, &commit, &proof));

        commit += G::generator();
        assert!(!IPAT::verify_commitment_proof(&crs, &commit, &proof));
    }

    /*
    #[test]
    fn test_eval_proof() {
        let size: u64 = 32;
        let gens: Vec<G> = (0..size).map(|i| G::generator() * F::from(i + 1)).collect();
        let q = G::generator() * F::from(size + 1);

        let data_raw: Vec<F> = (0..size).map(|i| F::from(i)).collect();

        let crs = IPAUniversalParams::<G, Hasher>::new(gens, q);
        let data = IPAPreparedData::<F>::new_incremental(data_raw);
        let commit = commit_data(&crs, &data);

        let index = thread_rng().gen_range(0..size) as usize;
        let e_proof = prove_eval(&crs, &commit, index, &data);
        assert!(verify_eval(&crs, &commit, &e_proof));
    }
    */
}
