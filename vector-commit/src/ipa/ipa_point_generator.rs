use std::marker::PhantomData;

use ark_ec::{
    hashing::{HashToCurve, HashToCurveError},
    AffineRepr, CurveGroup, Group,
};
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
        let mut i = 0usize;
        while res.len() < num {
            if let Ok(point) = hasher.hash(&i.to_le_bytes()) {
                res.push(point.into());
            }

            i += 1;
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

pub struct EthereumHashToCurve {
    domain: Vec<u8>,
}

impl<G: CurveGroup> HashToCurve<G> for EthereumHashToCurve {
    fn new(domain: &[u8]) -> Result<Self, ark_ec::hashing::HashToCurveError> {
        Ok(Self {
            domain: domain.to_vec(),
        })
    }

    fn hash(
        &self,
        message: &[u8],
    ) -> Result<<G as CurveGroup>::Affine, ark_ec::hashing::HashToCurveError> {
        let mut hasher = Sha256::new();
        hasher.update(&self.domain);
        hasher.update(message);
        let bytes = hasher.finalize();
        let point = <G::Affine as AffineRepr>::from_random_bytes(&bytes);
        point.ok_or(HashToCurveError::MapToCurveError(
            "Invalid point".to_owned(),
        ))
    }
}
