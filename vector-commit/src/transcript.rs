use std::marker::PhantomData;

use ark_ff::{field_hashers::HashToField, Field};
use ark_serialize::CanonicalSerialize;
use thiserror::Error;

#[derive(Error, Debug, Clone)]
pub enum TranscriptError {
    #[error("Invalid serialization")]
    InvalidSerialize,
}

pub trait Transcript<F> {
    /// Initialize a new hasher with a domain separator `label`
    fn new(label: &str) -> Self;

    /// Append an element and label to the state
    fn append<T: CanonicalSerialize>(
        &mut self,
        value: &T,
        label: &str,
    ) -> Result<(), TranscriptError>;

    /// Digest the current transcript to a field element. `clear` will clear the state and append the output to the state
    fn digest(&mut self, label: &str, clear: bool) -> F;
}

pub struct TranscriptHasher<F: Field, H: HashToField<F>> {
    state: Vec<u8>,
    hasher: H,
    _f: PhantomData<F>,
}

impl<F: Field, D: HashToField<F>> Transcript<F> for TranscriptHasher<F, D> {
    fn new(label: &str) -> Self {
        Self {
            state: Vec::new(),
            hasher: D::new(label.as_bytes()),
            _f: PhantomData,
        }
    }

    fn append<T: CanonicalSerialize>(
        &mut self,
        value: &T,
        label: &str,
    ) -> Result<(), TranscriptError> {
        self.state.append(&mut label.as_bytes().to_vec());
        self.state.append(&mut serialize(value)?);
        Ok(())
    }

    fn digest(&mut self, label: &str, clear: bool) -> F {
        self.state.append(&mut label.as_bytes().to_vec());
        let res = self.hasher.hash_to_field(&self.state, 1)[0];
        if clear {
            self.state = serialize(&res).unwrap();
            self.state.append(&mut label.as_bytes().to_vec());
        }
        res
    }
}

fn serialize<T: CanonicalSerialize>(x: &T) -> Result<Vec<u8>, TranscriptError> {
    let mut b = Vec::new();
    if x.serialize_compressed(&mut b).is_err() {
        return Err(TranscriptError::InvalidSerialize);
    }

    Ok(b)
}
