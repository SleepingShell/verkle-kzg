use std::marker::PhantomData;

use ark_ec::{hashing::HashToCurve, CurveGroup, Group};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use digest::Digest;
use num::traits::ToBytes;
use sha2::Sha256;

use crate::{PointGenerator, PointGeneratorError};

pub struct IPAPointGenerator<G, H> {
    max: usize,
    seed: Vec<u8>,
    _g: PhantomData<G>,
    _h: PhantomData<H>,
}

impl<G, H> IPAPointGenerator<G, H> {
    pub fn new(max: usize, seed: Vec<u8>) -> Self {
        Self {
            max,
            seed,
            _g: PhantomData,
            _h: PhantomData,
        }
    }
}

impl<G, H> Default for IPAPointGenerator<G, H> {
    fn default() -> Self {
        Self {
            max: 256,
            seed: "eth_verkle_oct_2021".to_owned().into_bytes(),
            _g: PhantomData,
            _h: PhantomData,
        }
    }
}

impl<G: CurveGroup, H: HashToCurve<G>> PointGenerator for IPAPointGenerator<G, H> {
    type Point = G;
    type Secret = Vec<u8>;

    fn gen(&self, num: usize) -> Result<Vec<Self::Point>, PointGeneratorError> {
        if num > self.max {
            return Err(PointGeneratorError::OutOfBounds);
        }
        let hasher = H::new(&self.seed).unwrap();
        let mut res: Vec<G> = Vec::with_capacity(num);
        for i in 0..num {
            res.push(
                hasher
                    .hash(&i.to_le_bytes())
                    .map_err(|_| PointGeneratorError::InvalidPoint)?
                    .into(),
            );
        }

        Ok(res)
    }

    fn gen_at(&self, index: usize) -> Result<Self::Point, PointGeneratorError> {
        if index > self.max {
            return Err(PointGeneratorError::OutOfBounds);
        }
        let hasher = H::new(&self.seed).unwrap();
        hasher
            .hash(&index.to_le_bytes())
            .map_err(|_| PointGeneratorError::InvalidPoint)
            .map(|p| p.into())
    }

    fn secret(&self) -> Option<Self::Secret> {
        Some(self.seed.clone())
    }
}
