//! `vector-commit` is a collection of traits for use in a vector commitment (VC) scheme.
//! A vector of data (the dataset) is committed to with a binding property (you cannot change the dataset after commitment).
//! One can then generate proofs of inclusion for the data to the commitment. These proofs can then
//! verify that the data was indeed in the dataset committed to by the VC scheme.
//!
//! Most VC schemes aim to generate constant or logarithmic sized proofs with efficient verification.
//! Some VC scheme require a trusted setup in which parameters are generated for proving/verification.
//! The binding property of these schemes is reliant on no one knowing the secret used in the trusted setup.
use std::{collections::HashMap, error::Error, fmt::Debug};

use ark_ec::Group;
use ark_ff::{FftField, Field, PrimeField, Zero};
use ark_poly::EvaluationDomain;
use lagrange_basis::LagrangeBasis;
use precompute::PrecomputedLagrange;
use thiserror::Error;
use transcript::Transcript;

pub mod ipa;
pub mod kzg;
pub mod lagrange_basis;
pub mod multiproof;
pub mod precompute;
pub(crate) mod transcript;
pub(crate) mod utils;

/// The proving and verification parameters for the VC scheme
pub trait VCUniversalParams<F: FftField> {
    /// Outputs the maximum number of items that can be committed to with this key
    fn max_size(&self) -> usize;

    fn precompute(&self) -> &PrecomputedLagrange<F>;
}

/// The dataset that a VC scheme will work over. This should include any preprocessing that is required
/*
pub trait VCPreparedData {
    type Item: Clone + Zero;
    type Error: Debug;

    fn from_vec(data: Vec<Self::Item>) -> Self;

    fn set_evaluation(&mut self, index: usize, value: Self::Item) -> Result<(), Self::Error>;

    fn get(&self, index: usize) -> Option<&Self::Item>;

    fn get_all(&self) -> Vec<(usize, Self::Item)>;

    /// Return the max amount of data that can be stored in this data
    fn max_size(&self) -> usize;
}
*/

pub trait VCCommitment<F> {
    fn to_field(&self) -> F;
}

/// Default implementation when the proof is simply a group element
impl<G: Group> VCCommitment<G::ScalarField> for G {
    fn to_field(&self) -> G::ScalarField {
        if self.is_zero() {
            G::ScalarField::ZERO
        } else {
            let mut bytes: Vec<u8> = Vec::new();
            // TODO: Check
            self.serialize_compressed(&mut bytes).unwrap();
            G::ScalarField::from_le_bytes_mod_order(&bytes)
        }
    }
}

/// A vector commitment schemes allows committing to a vector of data and generating proofs of inclusion.
pub trait VectorCommitment<F: PrimeField, D: EvaluationDomain<F>> {
    /// The universal parameters for the vector commitment scheme.
    /// CURRENTLY this API does not support differing committing, proving and verifying keys
    type UniversalParams: VCUniversalParams<F>;

    /// The Commitment to a vector.
    type Commitment: VCCommitment<F> + PartialEq + Clone;

    /// The proof for a single member of a vector.
    type Proof;

    /// The proof for multiple members of a vector.
    type BatchProof;

    /// The error type for the scheme.
    type Error: Error + Debug;

    /// The type that will generate the CRS points of the scheme
    type PointGenerator;

    /// The challenge generator using the Fiat-Shamir technique
    type Transcript: Transcript<F>;

    /// Constructs the Universal parameters for the scheme, which allows committing
    /// and proving inclusion of vectors up to `max_items` items
    fn setup(
        max_items: usize,
        gen: &Self::PointGenerator,
    ) -> Result<Self::UniversalParams, PointGeneratorError>;

    /// Commit a prepared data vector (`data`) to the `key` UniversalParams.
    fn commit(
        key: &Self::UniversalParams,
        data: &LagrangeBasis<F, D>,
    ) -> Result<Self::Commitment, Self::Error>;

    /// Prove that a piece of data exists inside of `commitment`. The `index` represents the index
    /// of the data inside of `data`.
    fn prove(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        data: &LagrangeBasis<F, D>,
    ) -> Result<Self::Proof, Self::Error> {
        Self::prove_point(key, commitment, F::from(index as u64), data, None)
    }

    /// Perform the same operation as the `prove` method, but take in a `Self::Point` evaluation point
    fn prove_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: F,
        data: &LagrangeBasis<F, D>,
        transcript: Option<Self::Transcript>,
    ) -> Result<Self::Proof, Self::Error>;

    /// Generate a batch proof that proves all of the `indexes`.
    fn prove_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        indexes: Vec<usize>,
        data: &LagrangeBasis<F, D>,
    ) -> Result<Self::BatchProof, Self::Error>;

    /// Verify that the `proof` is valid with respect to the `key` and `commitment`
    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        Self::verify_point(key, commitment, F::from(index as u64), proof, None)
    }

    /// Perform the same operation as the `verify` method, but take in a `Self::Point` evaluation point
    fn verify_point(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        point: F,
        proof: &Self::Proof,
        transcript: Option<Self::Transcript>,
    ) -> Result<bool, Self::Error>;

    /// Verify the batch proof is valid
    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error>;
}

#[derive(Error, Debug)]
pub enum PointGeneratorError {
    #[error("Attempted to create generator outside of max allowed")]
    OutOfBounds,
    #[error("Attempt to serialize bytes into a non-exsistent point")]
    InvalidPoint,
}

pub trait PointGenerator {
    type Point;
    type Secret;

    fn gen(&self, num: usize) -> Result<Vec<Self::Point>, PointGeneratorError>;
    fn gen_at(&self, index: usize) -> Result<Self::Point, PointGeneratorError>;
    fn secret(&self) -> Option<Self::Secret>;
}
