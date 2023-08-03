use std::error::Error;

use ark_ff::{PrimeField};
use rand::RngCore;

pub trait VCUniversalParams {
    /// Outputs the maximum number of items that can be committed to with this key
    fn max_size(&self) -> usize;
}

pub trait VCPreparedData {
    type Item;
    type Error;

    fn from_vec(data: Vec<Self::Item>) -> Self;

    fn set_evaluation(&mut self, index: usize, value: Self::Item) -> Result<(), Self::Error>;

    /// Return the max amount of data that can be stored in this data
    fn max_size(&self) -> usize;
}

/// A vector commitment schemes allows committing to a vector of data over a Finite Field,
/// and generating proofs of inclusion.
/// TODO: This should be a separate module
pub trait VectorCommitment {
    /// The universal parameters for the vector commitment scheme.
    /// CURRENTLY this API does not support differing committing, proving and verifying keys
    /// TODO: Seprate keys^
    type UniversalParams: VCUniversalParams;

    /// The vector dataset that has gone through preparation to use with the Vector Commitment.
    type PreparedData: VCPreparedData;

    /// The Commitment to a vector.
    type Commitment;

    /// The proof for a single member of a vector.
    type Proof;

    /// The proof for multiple members of a vector.
    type BatchProof;

    /// The error type for the scheme.
    type Error: Error;

    /// Constructs the Universal parameters for the scheme, which allows committing
    /// and proving inclusion of vectors up to `max_items` items
    fn setup<R: RngCore>(
        max_items: usize,
        rng: &mut R
    ) -> Result<Self::UniversalParams, Self::Error>;

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
        data: &Self::PreparedData
    ) -> Result<Self::BatchProof, Self::Error>;

    /// Verify that the `proof` is valid with respect to the `key` and `commitment`
    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof
    ) -> Result<bool, Self::Error>;

    /// Verify multiple proofs are valid
    /// TODO: Keep this as boolean return value, or number of valid proofs? Once compression is implemeneted then will be boolean
    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof
    ) -> Result<bool, Self::Error>;

    /// Converts a commitment to `PreparedData::Item`
    fn convert_commitment_to_data(commit: &Self::Commitment) -> <Self::PreparedData as VCPreparedData>::Item;
    
}