use ark_ec::{hashing::HashToCurve, CurveGroup, Group};
use ark_ff::{Field, One};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use digest::Digest;
use num::traits::ToBytes;
use sha2::Sha256;

use crate::{PointGenerator, PointGeneratorError};

pub struct KZGRandomPointGenerator<G: Group> {
    secret: G::ScalarField,
}

impl<G: Group> KZGRandomPointGenerator<G> {
    pub fn new(secret: G::ScalarField) -> Self {
        Self { secret }
    }
}

impl<G: Group> Default for KZGRandomPointGenerator<G> {
    fn default() -> Self {
        Self {
            secret: G::ScalarField::one(),
        }
    }
}

impl<G: Group> PointGenerator for KZGRandomPointGenerator<G> {
    type Point = G;
    type Secret = G::ScalarField;

    fn gen(&self, num: usize) -> Result<Vec<Self::Point>, PointGeneratorError> {
        let secret = self.secret;
        let gen = G::generator();
        let mut secret_cur = G::ScalarField::one();
        let mut res: Vec<G> = vec![gen; 1];
        for _i in 1..num {
            secret_cur *= secret;
            res.push(gen * secret_cur);
        }

        Ok(res)
    }

    fn gen_at(&self, index: usize) -> Result<Self::Point, PointGeneratorError> {
        Ok(G::generator() * self.secret.pow(&[index as u64]))
    }

    fn secret(&self) -> Option<Self::Secret> {
        Some(self.secret)
    }
}
