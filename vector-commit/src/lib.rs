//! `vector-commit` is a collection of traits for use in a vector commitment (VC) scheme.
//! A vector of data (the dataset) is committed to with a binding property (you cannot change the dataset after commitment).
//! One can then generate proofs of inclusion for the data to the commitment. These proofs can then
//! verify that the data was indeed in the dataset committed to by the VC scheme.
//!
//! Most VC schemes aim to generate constant or logarithmic sized proofs with efficient verification.
//! Some VC scheme require a trusted setup in which parameters are generated for proving/verification.
//! The binding property of these schemes is reliant on no one knowing the secret used in the trusted setup.
use std::{error::Error, fmt::Debug};

use ark_ff::{PrimeField, Zero};
use rand::RngCore;
use thiserror::Error;

pub mod ipa;
pub mod kzg;
mod transcript;

/// The proving and verification parameters for the VC scheme
pub trait VCUniversalParams {
    /// Outputs the maximum number of items that can be committed to with this key
    fn max_size(&self) -> usize;
}

/// The dataset that a VC scheme will work over. This should include any preprocessing that is required
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

/// A vector commitment schemes allows committing to a vector of data and generating proofs of inclusion.
pub trait VectorCommitment {
    /// The universal parameters for the vector commitment scheme.
    /// CURRENTLY this API does not support differing committing, proving and verifying keys
    type UniversalParams: VCUniversalParams;

    /// The vector dataset that has gone through preparation to use with the Vector Commitment.
    type PreparedData: VCPreparedData;

    /// The Commitment to a vector.
    type Commitment: PartialEq + Clone;

    /// The proof for a single member of a vector.
    type Proof;

    /// The proof for multiple members of a vector.
    type BatchProof;

    /// The error type for the scheme.
    type Error: Error + Debug;

    /// The type that will generate the CRS points of the scheme
    type PointGenerator;

    /// Constructs the Universal parameters for the scheme, which allows committing
    /// and proving inclusion of vectors up to `max_items` items
    fn setup(
        max_items: usize,
        gen: &Self::PointGenerator,
    ) -> Result<Self::UniversalParams, PointGeneratorError>;

    /// Commit a prepared data vector (`data`) to the `key` UniversalParams.
    fn commit(
        key: &Self::UniversalParams,
        data: &Self::PreparedData,
    ) -> Result<Self::Commitment, Self::Error>;

    /// Prove that a piece of data exists inside of `commitment`. The `index` represents the index
    /// of the data inside of `data`.
    fn prove(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        data: &Self::PreparedData,
    ) -> Result<Self::Proof, Self::Error>;

    /// Generate all proofs of the dataset using the Feist-Khovratovich techique
    fn prove_all(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        data: &Self::PreparedData,
    ) -> Result<Self::BatchProof, Self::Error>;

    /// Verify that the `proof` is valid with respect to the `key` and `commitment`
    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error>;

    /// Verify multiple proofs are valid
    /// TODO: Keep this as boolean return value, or number of valid proofs? Once compression is implemeneted then will be boolean
    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error>;

    /// Converts a commitment to `PreparedData::Item`
    fn convert_commitment_to_data(
        commit: &Self::Commitment,
    ) -> <Self::PreparedData as VCPreparedData>::Item;
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
