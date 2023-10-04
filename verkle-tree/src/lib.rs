//! Verkle Tree implementation over generic [vector commitment](https://docs.rs/vector-commit/latest/vector_commit/) schemes
//!
//! A Verkle Tree is a trie data structure in which children are committed to, and inclusion proven, using the alforementioned
//! vector commitment schemes.
//!
//! This crate is geared towards the Ethereum style of Verkle Trees, but aims to be generic as reasonably possible.
//! This style entails that there are two types of nodes: Internal and Extension.
//! - Internal nodes: A typical tree node that stores pointers to its' children nodes
//! - Extension nodes: The children of this node are the leaves, storing the actual data of the tree.
//!
//! Values in the tree are indexed by a `Key`. The implementation of the `Key` is an array of primitive integer types.
//! A single element of this array is referred to as a **Unit**. A stem refers to all units of a key except its terminating one.
//! This means that all leaves that share a stem will exist in the same Extension node.

use ark_ec::Group;
use ark_poly::EvaluationDomain;
use ark_serialize::CanonicalSerialize;
use bytemuck::{bytes_of, Pod};
use num::{One, Zero};
use std::{
    collections::{hash_map::Entry, HashMap},
    fmt::Debug,
    hash::Hash,
    marker::PhantomData,
};
use thiserror::Error;

use vector_commit::{
    multiproof::{VCCommitmentMultiProof, VectorCommitmentMultiproof},
    HasPrecompute, VCCommitment, VCData, VectorCommitment,
};

mod node;
use node::{Node, VerkleError};

/// KeyMethods defines methods that a key must implement
trait KeyMethods<const N: usize, UnitType> {
    /// Returns the index of where two keys differ. `cur_depth` is used as a hint for more efficient
    /// searching, and caller SHOULD assume that the keys are the same up to `cur_depth`
    fn next_diff_depth(&self, other: &Self, cur_depth: usize) -> usize;

    /// Splits a key into its stem and final unit
    fn split(self) -> ([UnitType; N], UnitType);

    /// Convert this key to a byte array
    fn to_bytes(&self) -> Vec<u8>;
}

/// A key represents how items are indexed inside of the verkle tree.
/// `N` represents the length of the key (and therefore a stem is of length N-1)
/// `UnitType` represents the data type for each unit of the key, which must be able to accomodate the arity of the tree.
type Key<const N: usize, UnitType> = [UnitType; N];

impl<const N: usize, UnitType: Zero + Copy + PartialEq + Pod> KeyMethods<N, UnitType>
    for Key<N, UnitType>
{
    fn next_diff_depth(&self, other: &Self, cur_depth: usize) -> usize {
        let mut d: usize = cur_depth + 1;
        while d < N {
            if self[d] != other[d] {
                break;
            }
            d += 1;
        }
        d
    }

    #[inline]
    fn split(self) -> ([UnitType; N], UnitType) {
        // TODO: Not too efficient
        let unit = self[N - 1];
        let stem = self;

        (stem, unit)
    }

    fn to_bytes(&self) -> Vec<u8> {
        self.map(|i| bytes_of(&i).to_owned()).concat()
    }
}

/// The Ethereum Verkle tree implementation requires that all values (which are 256-bits) must be split into
/// two 128-bit chunks. This is because many ECC orders (including the one Ethereum uses) are less than 256-bits.
pub trait SplittableValue {
    /// The type that is output from the split
    type Output;

    /// Split this type into a lower and upper half. Methodology is up to the implementor,
    /// but conventionally this will simply be split by the lower and upper bits of the type
    fn split(&self) -> (Self::Output, Self::Output);
}

/// A Verkle Tree implements convience functions that operate on a single root `Node`.
pub struct VerkleTree<const N: usize, K, VC, T, G, Domain>
where
    K: Eq + Hash,
    VC: VectorCommitment,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
    root: Node<N, K, VC, T>,
    _g: PhantomData<G>,
    _domain: PhantomData<Domain>,
}

/// Operations in this implementation block include all functionality outside of proving
impl<const N: usize, K, VC, T, G, Domain> VerkleTree<N, K, VC, T, G, Domain>
where
    K: Eq + Hash + Into<usize> + Zero + Copy + Pod + Into<usize>,
    VC: VectorCommitment,
    <VC::Data as VCData>::Item: Copy + One,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
    pub fn new() -> Self {
        Self {
            root: Node::new_internal(vec![]),
            _g: PhantomData,
            _domain: PhantomData,
        }
    }

    pub fn insert_single(&mut self, key: Key<N, K>, value: T) {
        let (stem, unit) = key.split();
        self.root.insert(stem, vec![(unit, value)], 0);
    }

    pub fn get_single(&self, key: &Key<N, K>) -> Option<&T> {
        let (stem, unit) = key.split();
        match self.root.get_stem(&stem, 0) {
            Some(stem) => stem.get_value(unit),
            None => None,
        }
    }

    pub fn commitment(&mut self, crs: &VC::UniversalParams) -> Result<VC::Commitment, VC::Error> {
        self.root.gen_commitment(crs).map(|c| c.clone())
    }

    fn path_to_stem<'a>(
        &'a self,
        stem: &Key<N, K>,
    ) -> Result<Vec<(Vec<K>, K, &'a Node<N, K, VC, T>)>, VerkleError> {
        let mut res = vec![];
        self.root.path_to_stem(stem, &mut res).map(move |_| res)
    }
}

/// This implementation block implements the proving functionality for the verkle tree
impl<const N: usize, K, VC, T, G, Domain> VerkleTree<N, K, VC, T, G, Domain>
where
    K: Eq + Hash + Into<usize> + Zero + Copy + Pod + Into<usize>,
    G: Group,
    Domain: EvaluationDomain<G::ScalarField> + Sync + Send,
    VC: VectorCommitmentMultiproof<G, Domain>,
    <VC::Data as VCData>::Item: Copy + One,
    VC::Commitment: VCCommitmentMultiProof<G::ScalarField>,
    VC::UniversalParams: HasPrecompute<G::ScalarField> + Sync,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
}

// TODO: Maybe publish a crate with a macro to allow derive(Default)
impl<const N: usize, K, VC, T, G, Domain> Default for VerkleTree<N, K, VC, T, G, Domain>
where
    K: Eq + Hash + Into<usize> + Zero + Copy + Pod + Into<usize>,
    VC: VectorCommitment,
    <VC::Data as VCData>::Item: Copy + One,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use ark_ec::pairing::Pairing;
    use ark_ff::{field_hashers::DefaultFieldHasher, PrimeField};
    use ark_poly::GeneralEvaluationDomain;
    use rand::{seq::SliceRandom, Fill, Rng};
    use sha2::Sha256;
    use std::ops::Add;

    use ark_bn254::Bn254;
    use vector_commit::kzg::{kzg_point_generator::KZGRandomPointGenerator, KZG};

    const KEY_LEN: usize = 3;
    type KEY_DATA_TYPE = u8;
    type KEYT = Key<KEY_LEN, KEY_DATA_TYPE>;

    type F = <Bn254 as Pairing>::ScalarField;
    type G1 = <Bn254 as Pairing>::G1;
    type Hasher = DefaultFieldHasher<Sha256>;

    #[derive(Debug, Clone, PartialEq)]
    struct U256([u8; 32]);
    type KZGT = KZG<Bn254, Hasher, GeneralEvaluationDomain<F>>;
    type TestTree = VerkleTree<KEY_LEN, KEY_DATA_TYPE, KZGT, U256>;

    impl SplittableValue for U256 {
        type Output = F;
        fn split(&self) -> (Self::Output, Self::Output) {
            (
                F::from_le_bytes_mod_order(&self.0[0..16]),
                F::from_le_bytes_mod_order(&self.0[16..32]),
            )
        }
    }

    impl Zero for U256 {
        fn zero() -> Self {
            U256([0; 32])
        }

        fn is_zero(&self) -> bool {
            self == &Self::zero()
        }

        fn set_zero(&mut self) {
            for i in 0..32 {
                self.0[i] = 0;
            }
        }
    }

    impl Add<Self> for U256 {
        type Output = Self;
        fn add(self, rhs: Self) -> Self::Output {
            let mut res = [0; 32];
            for i in 0..32 {
                res[i] = self.0[i] + rhs.0[i];
            }

            U256(res)
        }
    }

    impl Fill for U256 {
        fn try_fill<R: Rng + ?Sized>(&mut self, rng: &mut R) -> Result<(), rand::Error> {
            self.0.try_fill(rng)
        }
    }

    fn random_key(arity: KEY_DATA_TYPE, prefix: Option<&[KEY_DATA_TYPE]>) -> KEYT {
        let mut rng = rand::thread_rng();
        let mut res = [0; KEY_LEN];
        let mut p_size = 0;
        if let Some(p) = prefix {
            p_size = p.len();
            res[0..p_size].copy_from_slice(p);
        }

        for i in p_size..KEY_LEN {
            res[i] = rng.gen_range(0..arity);
        }

        res
    }

    fn random_u256() -> U256 {
        let mut res = U256::zero();
        res.0.try_fill(&mut rand::thread_rng());
        res
    }

    #[test]
    fn test_insert_get_leaves() {
        let num_leaves = 50;
        let mut rng = rand::thread_rng();

        let mut tree1 = TestTree::new();
        let mut tree2 = TestTree::new();

        // 1/4 of keys will share a stem
        let (stem, _) = random_key(255, None).split();
        let mut kvs: HashMap<KEYT, U256> = (0..num_leaves / 4)
            .map(|_| {
                let key = random_key(255, Some(&stem));

                (key.clone(), random_u256())
            })
            .collect();

        while kvs.len() < num_leaves {
            let key = random_key(255, None);
            kvs.insert(key.clone(), random_u256());
        }

        let keys: Vec<&KEYT> = kvs.keys().collect();
        let mut keys2 = keys.clone();
        keys2.shuffle(&mut rng);

        for (k1, k2) in keys.into_iter().zip(keys2.into_iter()) {
            let v1 = kvs.get(k1).unwrap();
            let v2 = kvs.get(k2).unwrap();
            tree1.insert_single(*k1, v1.clone());
            tree2.insert_single(*k2, v2.clone());
        }

        //assert!(tree1 == tree2);

        for k in kvs.keys().collect::<Vec<&KEYT>>() {
            let get1 = tree1.get_single(k);
            let get2 = tree2.get_single(k);

            assert!(get1 == get2);
            assert!(*get1.unwrap() == *kvs.get(k).unwrap());
        }
    }

    #[test]
    fn test_overwrite() {
        let mut tree = TestTree::new();

        let key = random_key(255, None);
        let val1 = random_u256();
        let val2 = random_u256();

        tree.insert_single(key, val1);
        tree.insert_single(key, val2.clone());

        assert!(tree.get_single(&key).unwrap() == &val2);
    }

    #[test]
    fn test_commitment() {
        let mut tree: TestTree = TestTree::new();
        let point_gen = KZGRandomPointGenerator::<G1>::default();
        let crs = KZGT::setup(256, &point_gen).unwrap();

        let key = random_key(255, None);
        let val = random_u256();
        tree.insert_single(key, val);

        let commit = tree.commitment(&crs);
        println!("{:?}", commit);
    }

    #[test]
    fn test_path_to_stem() {
        let mut tree: TestTree = TestTree::new();
        let point_gen = KZGRandomPointGenerator::<G1>::default();
        let crs = KZGT::setup(256, &point_gen).unwrap();

        let key = random_key(255, None);
        let val = random_u256();
        tree.insert_single(key, val);

        tree.insert_single(random_key(255, Some(&[key[0]])), random_u256());

        let path = tree.path_to_stem(&key);
        for (i, p) in path.unwrap().iter().enumerate() {
            assert!(p.0 == key[0..i + 1].to_vec());
            assert!(p.1 == key[i]);
        }
    }
}
