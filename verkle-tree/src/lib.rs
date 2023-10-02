use bytemuck::{bytes_of, Pod};
use num::{One, Zero};
use std::{
    collections::{hash_map::Entry, HashMap},
    fmt::Debug,
    hash::Hash,
};

use vector_commit::{VCCommitment, VCData, VectorCommitment};

trait KeyMethods<const N: usize, UnitType> {
    /// Returns the index of where two keys differ. `cur_depth` is used as a hint for more efficient
    /// searching, and caller SHOULD assume that the keys are the same up to `cur_depth`
    fn next_diff_depth(&self, other: &Self, cur_depth: usize) -> usize;

    /// Splits a key into its stem and final unit
    fn split(self) -> ([UnitType; N], UnitType);

    fn to_bytes(&self) -> Vec<u8>;
}

/// A key represents how items are indexed inside of the verkle tree.
/// `N` represents the length of the key (and therefore a stem is of length N-1)
/// `T` represents the data type for each unit of the key, which must be able to accomodate the arity of the tree.
type Key<const N: usize, T> = [T; N];

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

#[derive(Clone)]
enum Node<const N: usize, K, VC, T>
where
    VC: VectorCommitment,
{
    Internal {
        commit: Option<VC::Commitment>,
        children: HashMap<K, Node<N, K, VC, T>>,
    },
    Extension {
        stem: Key<N, K>, // TODO: Each stem has an extra unit because we cannot do N-1
        commit: Option<VC::Commitment>,
        leaves: HashMap<K, T>, // Sparse implementation...but at what cost :(
    },
}

impl<const N: usize, K, VC, T> Node<N, K, VC, T>
where
    K: Eq + Hash + Zero + Copy + Pod + Into<usize>,
    VC: VectorCommitment,
    <VC::Data as VCData>::Item: Copy + One,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq,
{
    fn new_extension(stem: Key<N, K>, values: Vec<(K, T)>) -> Self {
        Self::Extension {
            stem,
            commit: None,
            leaves: values.into_iter().collect(),
        }
    }

    fn new_internal(nodes: Vec<(K, Self)>) -> Self {
        Self::Internal {
            commit: None,
            children: nodes.into_iter().map(|v| (v.0, v.1)).collect(),
        }
    }

    /// Get the extension node equal to this stem. If it does not exist, return None
    fn get_stem(&self, stem: &Key<N, K>, mut cur_depth: usize) -> Option<&Self> {
        match self {
            Self::Extension {
                stem: self_stem, ..
            } => {
                if stem == self_stem {
                    Some(self)
                } else {
                    None
                }
            }
            Self::Internal { children, .. } => match children.get(&stem[cur_depth]) {
                Some(c) => {
                    cur_depth += 1;
                    c.get_stem(stem, cur_depth)
                }
                None => None,
            },
        }
    }

    /// Gets the value from an extension node
    fn get_value(&self, unit: K) -> Option<&T> {
        match self {
            Self::Extension { leaves, .. } => leaves.get(&unit),
            _ => panic!("Called get_value on non-extension node"),
        }
    }

    /// Recursively insert the `values` into an extension node that has its stem equal to `stem`
    /// This function will clear the stored commit of all touched nodes, as they are no longer valid
    fn insert(&mut self, stem: Key<N, K>, values: Vec<(K, T)>, cur_depth: usize) {
        match self {
            Self::Extension {
                stem: self_stem,
                commit,
                leaves,
            } => {
                // This function should only ever be called on an Extension node to insert values. I.e internal nodes
                // will never call into here. Therefore self_stem should ALWAYS be equal to stem
                if self_stem != &stem {
                    panic!("Traversed to extension node with differing stem");
                }
                *commit = None;
                values.into_iter().for_each(|v| {
                    leaves.insert(v.0, v.1);
                });
            }

            // If we are an internal node, then we will be inserting the value into a child node
            Self::Internal {
                commit, children, ..
            } => {
                // Set commit to None as the commit is no longer accurate
                *commit = None;
                let k = stem[cur_depth];
                let child = children.entry(k);
                match child {
                    // If we have a child that matches on the current unit of the key
                    Entry::Occupied(mut o) => {
                        let child_as_node = o.get_mut();
                        match child_as_node {
                            // If the child is an extension node, then we will either insert into it if the stem's match,
                            // or we must create a new branch resulting in a new inner node containing two extension nodes
                            Self::Extension {
                                stem: child_stem, ..
                            } => {
                                // The stems match, and therefore simply insert into the extension node
                                if &stem == child_stem || cur_depth == N - 2 {
                                    child_as_node.insert(stem, values, cur_depth + 1);

                                // The stems differ, so a new inner node is created with the children hashmap key
                                // equal to the differing unit
                                } else {
                                    let depth = child_stem.next_diff_depth(&stem, cur_depth);
                                    let nodes = vec![
                                        (stem[depth], Self::new_extension(stem, values)),
                                        (child_stem[depth], o.remove()),
                                    ];

                                    let new_internal = Self::new_internal(nodes);
                                    children.insert(k, new_internal);
                                }
                            }
                            Self::Internal { .. } => {
                                child_as_node.insert(stem, values, cur_depth + 1);
                            }
                        }
                    }
                    Entry::Vacant(v) => {
                        v.insert(Self::new_extension(stem, values));
                    }
                }
            }
        }
    }

    /// Generates a commitment recursively. If all of a Node's children have `Some` commitment,
    /// then there is no need to recurse further.
    fn gen_commitment(&mut self, crs: &VC::UniversalParams) -> Result<&VC::Commitment, VC::Error> {
        match self {
            Self::Extension {
                stem,
                commit,
                leaves,
            } => {
                if commit.is_some() {
                    return Ok(&commit.as_ref().unwrap());
                }

                let mut c1_values = [<VC::Data as VCData>::Item::zero(); N];
                let mut c2_values = [<VC::Data as VCData>::Item::zero(); N];

                for (&index, leaf) in leaves.iter() {
                    let (low, high) = leaf.split();
                    let index_low: usize = (2 * index.into()) % N;
                    let index_high: usize = (2 * index.into() + 1) % N;

                    if index.into() < (N / 2) {
                        c1_values[index_low] = low;
                        c1_values[index_high] = high;
                    } else {
                        c2_values[index_low] = low;
                        c2_values[index_high] = high;
                    }
                }

                let c1 = VC::commit(crs, &<VC::Data as VCData>::from_vec(c1_values.to_vec()))?;
                let c2 = VC::commit(crs, &<VC::Data as VCData>::from_vec(c2_values.to_vec()))?;

                let extension_data = vec![
                    <VC::Data as VCData>::Item::one(),
                    <VC::Data as VCData>::bytes_to_item(&stem.to_bytes()),
                    c1.to_data_item(),
                    c2.to_data_item(),
                ];

                let c = VC::commit(crs, &<VC::Data as VCData>::from_vec(extension_data))?;
                *commit = Some(c);

                Ok(commit.as_ref().unwrap())
            }
            Self::Internal { commit, children } => {
                if commit.is_some() {
                    return Ok(commit.as_ref().unwrap());
                }

                // FIXME FIXME THIS IS HARDCODED NEED TO FIX
                let mut vc_vec = vec![<VC::Data as VCData>::Item::zero(); 256];
                for (&k, child) in children.iter_mut() {
                    let cc = child.gen_commitment(crs)?;
                    vc_vec[k.into()] = cc.to_data_item();
                    //vc_vec.insert(*k.into(), cc.to_data_item());
                }

                let vc_data = <VC::Data as VCData>::from_vec(vc_vec);
                let c = VC::commit(crs, &vc_data)?;
                *commit = Some(c);

                Ok(commit.as_ref().unwrap())
            }
        }
    }
}

impl<const N: usize, K, VC, T> Debug for Node<N, K, VC, T>
where
    K: Debug,
    VC: VectorCommitment,
    VC::Commitment: Debug,
    T: Debug + Zero + PartialEq,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Extension {
                stem,
                commit,
                leaves,
            } => {
                f.write_fmt(format_args!("Extension ({:?}) {{\n", stem))?;
                f.write_fmt(format_args!("\tCommit: {:?}\n", commit))?;
                f.write_fmt(format_args!("\tChildren: {{\n"))?;
                for (i, v) in leaves.iter() {
                    if *v != T::zero() {
                        f.write_fmt(format_args!("({:?},{:?})\n", i, v))?;
                    }
                }
                f.write_str("\t}")
            }
            Self::Internal { commit, children } => {
                f.write_fmt(format_args!("Inner {{\n"))?;
                f.write_fmt(format_args!("\tCommit: {:?}\n", commit))?;
                f.write_fmt(format_args!("\tChildren: {{\n"))?;
                for c in children {
                    f.write_fmt(format_args!("\t\t{:?}\n", c.1))?;
                }
                f.write_str("\t}")
            }
        }
    }
}

/// A Verkle Tree implements convience functions that operate on a single root `Node`.
pub struct VerkleTree<const N: usize, K, VC, T>
where
    K: Eq + Hash,
    VC: VectorCommitment,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
    root: Node<N, K, VC, T>,
}

impl<const N: usize, K, VC, T> VerkleTree<N, K, VC, T>
where
    K: Eq + Hash + Into<usize> + Zero + Copy + Pod + Into<usize>,
    VC: VectorCommitment,
    <VC::Data as VCData>::Item: Copy + One,
    T: SplittableValue<Output = <VC::Data as VCData>::Item> + Zero + Clone + PartialEq + Debug,
{
    pub fn new() -> Self {
        Self {
            root: Node::new_internal(vec![]),
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

    fn random_key(arity: KEY_DATA_TYPE, prefix: Option<&KEYT>) -> Key<KEY_LEN, KEY_DATA_TYPE> {
        let mut rng = rand::thread_rng();
        let mut res = [0; KEY_LEN];
        let mut p_size = 0;
        if let Some(p) = prefix {
            p_size = p.len();
            res[0..KEY_LEN].copy_from_slice(p);
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
        let NUM_LEAVES = 50;
        let mut rng = rand::thread_rng();

        let mut tree1 = TestTree::new();
        let mut tree2 = TestTree::new();

        // 1/4 of keys will share a stem
        let (stem, _) = random_key(255, None).split();
        let mut kvs: HashMap<KEYT, U256> = (0..NUM_LEAVES / 4)
            .map(|_| {
                let key = random_key(255, Some(&stem));

                (key.clone(), random_u256())
            })
            .collect();

        while kvs.len() < NUM_LEAVES {
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
        let mut rng = rand::thread_rng();
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
        let mut rng = rand::thread_rng();
        let mut tree: TestTree = TestTree::new();
        let point_gen = KZGRandomPointGenerator::<G1>::default();
        let crs = KZGT::setup(256, &point_gen).unwrap();

        let key = random_key(255, None);
        let val = random_u256();
        tree.insert_single(key, val);

        let commit = tree.commitment(&crs);
        println!("{:?}", tree.root);
        println!("{:?}", commit);
    }
}
