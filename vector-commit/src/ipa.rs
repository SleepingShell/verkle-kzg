use std::{
    collections::HashMap,
    iter::Sum,
    marker::PhantomData,
    ops::{Add, Mul, MulAssign},
};

use ark_ec::Group;
use ark_ff::{field_hashers::HashToField, FftField, Field, One};
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_serialize::CanonicalSerialize;

use crate::{VCPreparedData, VCUniversalParams, VectorCommitment};

struct IPAUniversalParams<G: Group, D: HashToField<G::ScalarField>> {
    g: Vec<G>, // Gens to commit the coefficients of the polynomial
    q: G,      // Gen to commit to the evaluation of the polynomial
    digest: PhantomData<D>,
}

impl<G: Group, D: HashToField<G::ScalarField>> IPAUniversalParams<G, D> {
    fn new(g: Vec<G>, q: G) -> Self {
        Self {
            g,
            q,
            digest: PhantomData,
        }
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
    x: G::ScalarField,
    y: G::ScalarField,
}

struct IPAPreparedData<F: Field> {
    // An ordered set of data, of (x, eval)
    data: Vec<(usize, F)>,

    // We keep the weights for barycentric interpolation w_j = ∏_{m != j}(x_j - x_m)^-1
    weights: Vec<F>,

    // TODO: Move outside of this so we only store once l(x) = ∏(X-x_i)
    universal_basis: DensePolynomial<F>,

    //
    universal_basis_derivative: DensePolynomial<F>,

    // The highest x value that this dataset contains
    max: usize,
}

impl<F: FftField> IPAPreparedData<F> {
    /// An ORDERED list of evaluation points
    fn new(data: Vec<(usize, F)>) -> Self {
        let max = data[data.len() - 1].0;
        let mut d_calc = Vec::<(F, F)>::new();
        let mut x_points = Vec::<F>::new();
        data.iter().for_each(|point| {
            x_points.push(F::from(point.0 as u64));
            d_calc.push((F::from(point.0 as u64), point.1));
        });
        let basis = Self::compute_universal_basis(x_points.clone());
        Self {
            data,
            weights: Self::compute_barycentric_weights(d_calc),
            universal_basis: basis.clone(),
            universal_basis_derivative: Self::compute_universal_basis_derivative(x_points, basis),
            max,
        }
    }

    /// The data domain is simply at 0..length-1
    fn new_incremental(data: Vec<F>) -> Self {
        let max = data.len();
        let mut data_stored = Vec::<(usize, F)>::new();
        data.into_iter().enumerate().for_each(|(i, v)| {
            data_stored.push((i, v));
        });

        Self::new(data_stored)
    }

    fn compute_universal_basis(x_points: Vec<F>) -> DensePolynomial<F> {
        let mut res = DensePolynomial::<F>::from_coefficients_slice(&[F::one()]);
        x_points.into_iter().for_each(|x| {
            res = &res * &DensePolynomial::<F>::from_coefficients_slice(&[F::zero() - x, F::one()]);
        });

        res
    }

    fn compute_universal_basis_derivative(
        x_points: Vec<F>,
        universal_basis: DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let mut res = DensePolynomial::<F>::from_coefficients_slice(&[F::zero()]);
        x_points.into_iter().for_each(|x| {
            res += &(&universal_basis
                / &DensePolynomial::<F>::from_coefficients_slice(&[F::zero() - x, F::one()]));
        });

        res
    }

    fn compute_barycentric_weights(points: Vec<(F, F)>) -> Vec<F> {
        let mut weights: Vec<F> = Vec::new();
        points.iter().for_each(|point| {
            let mut sum = F::one();
            points.iter().for_each(|p_inner| {
                if point.0 != p_inner.0 {
                    sum *= point.0 - p_inner.0;
                }
            });
            weights.push(sum.inverse().unwrap());
        });

        weights
    }

    fn compute_b_vector(&self, index: usize) -> Vec<F> {
        let index_f = F::from(index as u64);
        let universal_eval = self.universal_basis.evaluate(&index_f);
        let mut res = Vec::new();
        self.data.iter().for_each(|(x, _y)| {
            if index == *x {
                res.push(F::one());
            } else {
                let x_f = F::from(*x as u64);
                let rhs = index_f - x_f;
                let deriv_eval: F = self.universal_basis_derivative.evaluate(&x_f);
                res.push(universal_eval / (deriv_eval * rhs));
            }
        });

        res
    }

    /// Return the data of size `max` with zero entries included
    fn expanded_data(&self) -> Vec<F> {
        let mut res = vec![F::zero(); self.max + 1];
        for point in self.data.iter() {
            res[point.0] = point.1;
        }

        res
    }

    /// Simply return the stored evaluations without their index's (densely packed)
    fn evaluations(&self) -> Vec<F> {
        self.data.iter().map(|(_i, p)| *p).collect()
    }
}

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
        let (poly_l, poly_r) = split(data);
        let (gens_l, gens_r) = split(gens);

        let y_l = inner_product(&gens_r, &poly_l);
        let y_r = inner_product(&gens_l, &poly_r);
        l.push(y_l);
        r.push(y_r);

        ra = hasher.hash_to_field(
            &[serialize(&ra), serialize(&y_l), serialize(&y_r)].concat(),
            1,
        )[0];

        data = vec_add_and_distribute(&poly_l, &poly_r, ra);
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
    //let x = G::ScalarField::from(index as u64);
    //let mut x_pows: Vec<G::ScalarField> = (0..data.len()).map(|i| x.pow([i as u64])).collect();

    //TODO: Temporary evaluation at 'index'. Will be migrating to evluation form
    let y = inner_product(&data, &x_pows);
    let hasher = D::new(&[]);

    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();
    let mut ra = hasher.hash_to_field(&serialize(commitment), 1)[0];

    let q = key.q * ra;
    while data.len() > 1 {
        let (poly_l, poly_r) = split(data);
        let (gens_l, gens_r) = split(gens);
        let (x_pows_l, x_pows_r) = split(x_pows);

        let y_l = inner_product(&gens_r, &poly_l) + q * inner_product(&poly_l, &x_pows_r);
        let y_r = inner_product(&gens_l, &poly_r) + q * inner_product(&poly_r, &x_pows_l);

        l.push(y_l);
        r.push(y_r);

        ra = hasher.hash_to_field(
            &[serialize(&ra), serialize(&y_l), serialize(&y_r)].concat(),
            1,
        )[0];

        data = vec_add_and_distribute(&poly_l, &poly_r, ra);
        gens = vec_add_and_distribute(&gens_r, &gens_l, ra);
        x_pows = vec_add_and_distribute(&x_pows_r, &x_pows_l, ra);
    }

    IPAProof {
        l,
        r,
        tip: data[0],
        x: G::ScalarField::from(index as u64),
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
    let mut x_pows: Vec<G::ScalarField> =
        (0..gens.len()).map(|i| proof.x.pow([i as u64])).collect();
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
            .collect();
    }

    let combined_point = inner_product(&gens, &points_coeffs);
    let combined_x_pow = inner_product(&x_pows, &points_coeffs);

    c == combined_point * proof.tip + q * (proof.tip * combined_x_pow)
}

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

    use sha2::Sha256;

    type F = <Bn254 as Pairing>::ScalarField;
    type G = <Bn254 as Pairing>::G1;
    type Hasher = DefaultFieldHasher<Sha256>;

    #[test]
    fn test_commit_poly() {
        let size = 32;
        let gens: Vec<G> = (0..size).map(|i| G::generator() * F::from(i + 1)).collect();
        let q = G::generator() * F::from(size + 1);

        //let coeffs: Vec<F> = (0..size).map(|i| F::from(i)).collect();
        let data_raw: Vec<F> = (0..size).map(|i| F::from(i)).collect();

        let crs = IPAUniversalParams::<G, Hasher>::new(gens, q);
        let data = IPAPreparedData::<F>::new_incremental(data_raw);

        let mut commit = commit_data(&crs, &data);
        let poly_proof = prove_data_commit(&crs, &commit, &data);
        assert!(verify_data_commit(&crs, &commit, &poly_proof));

        commit += G::generator();
        assert!(!verify_data_commit(&crs, &commit, &poly_proof));
    }

    #[test]
    fn test_eval_proof() {
        let size = 32;
        let gens: Vec<G> = (0..size).map(|i| G::generator() * F::from(i + 1)).collect();
        let q = G::generator() * F::from(size + 1);

        let data_raw: Vec<F> = (0..size).map(|i| F::from(i)).collect();

        let crs = IPAUniversalParams::<G, Hasher>::new(gens, q);
        let data = IPAPreparedData::<F>::new_incremental(data_raw);
        let commit = commit_data(&crs, &data);

        let mut e_proof = prove_eval(&crs, &commit, 0, &data);
        assert!(verify_eval(&crs, &commit, &e_proof));
    }
}
