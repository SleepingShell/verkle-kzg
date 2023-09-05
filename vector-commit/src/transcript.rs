use std::{error::Error, fmt::Display, marker::PhantomData};

use ark_ff::{field_hashers::HashToField, Field};
use ark_serialize::CanonicalSerialize;
use thiserror::Error;

#[derive(Error, Debug, Clone)]
pub enum TranscriptError {
    #[error("Invalid serialization")]
    InvalidSerialize,
}

pub(crate) struct Transcript<F: Field, D: HashToField<F>> {
    state: Vec<u8>,
    hasher: D,
    _f: PhantomData<F>,
}

impl<F: Field, D: HashToField<F>> Transcript<F, D> {
    pub(crate) fn new(label: &str) -> Self {
        Self {
            state: Vec::new(),
            hasher: D::new(label.as_bytes()),
            _f: PhantomData,
        }
    }

    pub(crate) fn append<T: CanonicalSerialize>(
        &mut self,
        value: &T,
        label: &str,
    ) -> Result<(), TranscriptError> {
        self.state.append(&mut label.as_bytes().to_vec());
        self.state.append(&mut Self::serialize(value)?);
        Ok(())
    }

    /// Hash the current transcript to a field element. `clear` will clear the state and append the hashed output
    pub(crate) fn hash(&mut self, label: &str, clear: bool) -> F {
        self.state.append(&mut label.as_bytes().to_vec());
        let res = self.hasher.hash_to_field(&self.state, 1)[0];
        if clear {
            self.state = Self::serialize(&res).unwrap();
            self.state.append(&mut label.as_bytes().to_vec());
        }
        res
    }

    fn serialize<T: CanonicalSerialize>(x: &T) -> Result<Vec<u8>, TranscriptError> {
        let mut b = Vec::new();
        if x.serialize_compressed(&mut b).is_err() {
            return Err(TranscriptError::InvalidSerialize);
        }

        Ok(b)
    }
}
