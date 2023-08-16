use std::{
    iter::Sum,
    marker::PhantomData,
    ops::{Add, Mul},
};

use ark_ec::Group;
use ark_ff::{field_hashers::HashToField, Field, One, Zero};
use ark_serialize::CanonicalSerialize;
use digest::Digest;

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
    data: Vec<F>,
}

impl<F: Field> IPAPreparedData<F> {
    fn new(data: Vec<F>) -> Self {
        Self { data }
    }
}

/// Commit to the dataset
fn commit_data<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    data: &IPAPreparedData<G::ScalarField>,
) -> IPACommitment<G> {
    let len = data.data.len();
    let gens = key.g[0..len].to_vec();

    inner_product(&gens, &data.data)
}

/// Generate a proof that the dataset is committed to a valid polynomial
fn prove_data_commit<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    data: &IPAPreparedData<G::ScalarField>,
) -> IPACommitProof<G> {
    let mut data = data.data.clone();
    let mut gens = key.g[0..data.len()].to_vec();
    let mut l: Vec<G> = Vec::new();
    let mut r: Vec<G> = Vec::new();

    let hasher = D::new(&[]);
    let mut com_bytes = Vec::new();
    commitment.serialize_compressed(&mut com_bytes);
    let mut ra = hasher.hash_to_field(&com_bytes, 1)[0];

    while data.len() > 1 {
        let (poly_l, poly_r) = split(data);
        let (gens_l, gens_r) = split(gens);

        let y_l = inner_product(&gens_r, &poly_l);
        let y_r = inner_product(&gens_l, &poly_r);
        l.push(y_l);
        r.push(y_r);

        let mut yl_bytes = Vec::new();
        let mut yr_bytes = Vec::new();
        y_l.serialize_compressed(&mut yl_bytes);
        y_r.serialize_compressed(&mut yr_bytes);

        let mut ra_bytes = Vec::new();
        ra.serialize_compressed(&mut ra_bytes);
        ra = hasher.hash_to_field(&[ra_bytes, yl_bytes, yr_bytes].concat(), 1)[0];

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
    let mut gens = key.g[0..(2usize).pow(proof.l.len() as u32)].to_vec();
    let mut c = commitment.clone();
    let hasher = D::new(&[]);
    let mut com_bytes = serialize(commitment);
    let mut ra = hasher.hash_to_field(&com_bytes, 1)[0];
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

fn prove<G: Group, D: HashToField<G::ScalarField>>(
    key: &IPAUniversalParams<G, D>,
    commitment: &IPACommitment<G>,
    index: usize,
    data: &IPAPreparedData<G::ScalarField>,
) {
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

    use std::ops::{Add, AddAssign};

    use ark_bn254::Bn254;
    use ark_ec::pairing::Pairing;
    use ark_ec::Group;
    use ark_ff::{field_hashers::DefaultFieldHasher, Field, One, UniformRand};

    use sha2::Sha256;

    type F = <Bn254 as Pairing>::ScalarField;
    type G = <Bn254 as Pairing>::G1;
    type hasher = DefaultFieldHasher<Sha256>;

    #[test]
    fn test_commit_poly() {
        let size = 8;
        let gens: Vec<G> = (0..size).map(|i| G::generator() * F::from(i + 1)).collect();
        let q = G::generator() * F::from(size + 1);

        let coeffs: Vec<F> = (0..size).map(|i| F::from(i)).collect();

        let crs = IPAUniversalParams::<G, hasher>::new(gens, q);
        let data = IPAPreparedData::<F>::new(coeffs);

        let commit = commit_data(&crs, &data);
        let mut poly_proof = prove_data_commit(&crs, &commit, &data);
        assert!(verify_data_commit(&crs, &commit, &poly_proof));
    }
}
