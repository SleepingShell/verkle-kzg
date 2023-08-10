use std::{
    collections::{hash_map::Entry, HashMap},
    fmt::Debug,
    ops::Index,
};

use ark_ff::Zero;

use crate::{data_structures::VCPreparedData, VectorCommitment};

/// A Verkle Tree implementation, specifically using Ethereum's stem and extension model.
/// There are two types of nodes in this model, although I have attempted to leave open easy addition for more types.
/// - Inner nodes: Nodes whose children are either other inner nodes, or Extension nodes.
/// - Extension nodes: Nodes who contain the actual values stored in the tree, as leaves.

// TODO: Generic underlying unit of keys, will allow packing units into memory
trait Key {
    type FullKey;
    type Stem: Eq + Index<usize, Output = usize> + Debug;

    /// Return the full key length
    fn len() -> usize;

    fn diff_depth(cur_depth: usize, first: &Self::Stem, other: &Self::Stem) -> usize;

    /// Splits a full key into its stem and final unit
    fn split_full(full: &Self::FullKey) -> (Self::Stem, usize);
}

// N is the size of the key
// T is the leaf value type (a leaf does not need to know its own key)
enum Node<const N: usize, K, VC, T>
where
    K: Key,
    VC: VectorCommitment,
    VC::PreparedData: VCPreparedData<Item = T>,
{
    Internal {
        commit: Option<VC::Commitment>,
        vc_data: VC::PreparedData,
        children: HashMap<usize, Node<N, K, VC, T>>,
    },
    Extension {
        stem: K::Stem,
        commit: Option<VC::Commitment>,
        vc_data: VC::PreparedData,
    },
}

// The VerkleTree will receive every node that is now dirty, so it can later dispatch proof generation.
// If we want to parallelize proof generation, then need to keep track of the relation of nodes to one another, so different
// branches of equal sizing can be distributed.
// For non-parallel simpler implementation, we simply start from the largest keys (and therefore further down the tree).
//type ChangeCallbackFn<const N: usize, VC, T> = fn(Vec<usize>, &mut Node<N,VC,T>);

impl<const N: usize, K, VC, T> Node<N, K, VC, T>
where
    K: Key,
    VC: VectorCommitment,
    VC::PreparedData: VCPreparedData<Item = T>,
    T: Zero + Clone + PartialEq,
{
    fn new_extension(stem: K::Stem, values: Vec<(usize, T)>) -> Self {
        let mut vc_vec = vec![T::zero(); N + 1];
        values.into_iter().for_each(|v| vc_vec[v.0] = v.1);

        Self::Extension {
            stem,
            commit: None,
            vc_data: <VC::PreparedData as VCPreparedData>::from_vec(vc_vec),
        }
    }

    fn new_internal(nodes: Vec<(usize, Self)>) -> Self {
        Self::Internal {
            commit: None,
            vc_data: <VC::PreparedData as VCPreparedData>::from_vec(vec![]),
            children: nodes.into_iter().map(|v| (v.0, v.1)).collect(),
        }
    }

    //TODO: #[feature(generic_const_exprs)] would allow us to ensure the stem is N-1, but this is only available in nightly
    fn insert_into_stem(&mut self, cur_depth: usize, stem: K::Stem, values: Vec<(usize, T)>) {
        match self {
            Self::Extension {
                stem: self_stem,
                commit,
                vc_data,
                ..
            } => {
                // This function should only ever be called on an Extension node to insert values. I.e internal nodes
                // will never call into here. Therefore self_stem should ALWAYS be equal to stem
                if self_stem != &stem {
                    panic!("Traversed to extension node with differing stem");
                }
                *commit = None;
                values.into_iter().for_each(|v| {
                    vc_data.set_evaluation(v.0, v.1);
                });
            }
            Self::Internal {
                commit, children, ..
            } => {
                let k = stem[cur_depth];
                let child = children.entry(k);
                match child {
                    Entry::Occupied(mut o) => {
                        let child_as_node = o.get_mut();
                        match child_as_node {
                            Self::Extension {
                                stem: child_stem, ..
                            } => {
                                if &stem == child_stem || cur_depth == K::len() - 2 {
                                    child_as_node.insert_into_stem(cur_depth + 1, stem, values);
                                } else {
                                    let depth = K::diff_depth(cur_depth, child_stem, &stem);
                                    let nodes = vec![
                                        (stem[depth], Self::new_extension(stem, values)),
                                        (child_stem[depth], o.remove()),
                                    ];

                                    let new_internal = Self::new_internal(nodes);
                                    children.insert(k, new_internal);
                                }
                            }
                            Self::Internal { .. } => {
                                child_as_node.insert_into_stem(cur_depth + 1, stem, values);
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

    fn get_stem(&self, mut cur_depth: usize, stem: K::Stem) -> Option<&Self> {
        match self {
            Self::Extension {
                stem: self_stem, ..
            } => {
                if &stem == self_stem {
                    Some(self)
                } else {
                    None
                }
            }
            Self::Internal {
                commit, children, ..
            } => match children.get(&stem[cur_depth]) {
                Some(c) => {
                    cur_depth += 1;
                    c.get_stem(cur_depth, stem)
                }
                None => None,
            },
        }
    }

    fn get_value(&self, index: usize) -> Option<&T> {
        match self {
            Self::Extension { vc_data, .. } => vc_data.get(index),
            _ => panic!("Called get_value on non-extension node"),
        }
    }
}

struct VerkleTree<const N: usize, K, VC, T>
where
    K: Key,
    VC: VectorCommitment,
    VC::PreparedData: VCPreparedData<Item = T>,
    T: Zero + Clone + PartialEq,
{
    root: Node<N, K, VC, T>,
}

impl<const N: usize, K, VC, T> VerkleTree<N, K, VC, T>
where
    K: Key,
    VC: VectorCommitment,
    VC::PreparedData: VCPreparedData<Item = T>,
    T: Zero + Clone + PartialEq,
{
    fn new() -> Self {
        Self {
            root: Node::new_internal(vec![]),
        }
    }

    fn insert_single(&mut self, key: K::FullKey, value: T) {
        let (stem, unit) = K::split_full(&key);
        self.root.insert_into_stem(0, stem, vec![(unit, value)]);
    }

    fn get_single(&self, key: &K::FullKey) -> Option<&T> {
        let (stem, unit) = K::split_full(&key);
        match self.root.get_stem(0, stem) {
            Some(stem) => stem.get_value(unit),
            None => None,
        }
    }
}

#[derive(PartialEq)]
struct VanillaKey<const F: usize, const S: usize>;
impl<const F: usize, const S: usize> Key for VanillaKey<F, S> {
    type FullKey = [usize; F];
    type Stem = [usize; S];

    fn len() -> usize {
        F
    }

    fn diff_depth(cur_depth: usize, first: &Self::Stem, other: &Self::Stem) -> usize {
        let mut d: usize = cur_depth + 1;
        while d < S {
            if first[d] != other[d] {
                break;
            }
            d += 1;
        }
        d
    }

    fn split_full(full: &Self::FullKey) -> (Self::Stem, usize) {
        //full[0..S].try_into().unwrap()
        let (stem, unit) = full.split_at(S);
        (stem.try_into().unwrap(), unit[0])
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::pairing::Pairing;
    use rand::{seq::SliceRandom, Rng};

    use super::*;
    use crate::kzg_amortized::KZGAmortized;
    use ark_bn254::Bn254;

    const KEY_LEN: usize = 3;
    const ARITY: usize = 10;
    type TestKey = VanillaKey<KEY_LEN, { KEY_LEN - 1 }>;
    type TestFullKey = <TestKey as Key>::FullKey;
    type F = <Bn254 as Pairing>::ScalarField;
    type KZG = KZGAmortized<Bn254>;
    type TestTree = VerkleTree<ARITY, TestKey, KZG, F>;

    fn random_key(arity: usize, prefix: Option<&<TestKey as Key>::Stem>) -> TestFullKey {
        let mut rng = rand::thread_rng();
        let mut res = [0; KEY_LEN];
        let mut p_size = 0;
        if let Some(p) = prefix {
            p_size = p.len();
            res[0..KEY_LEN - 1].copy_from_slice(p);
        }
        //let real_len = len - res.len();

        for i in p_size..KEY_LEN {
            res[i] = rng.gen_range(0..arity);
        }

        res
    }

    #[test]
    fn test_insert_get_leaves() {
        let NUM_LEAVES = 50;
        let mut rng = rand::thread_rng();

        let mut tree1 = TestTree::new();
        let mut tree2 = TestTree::new();

        // 1/4 of keys will share a stem
        let (stem, _) = TestKey::split_full(&random_key(ARITY, None));
        let mut kvs: HashMap<TestFullKey, F> = (0..NUM_LEAVES / 4)
            .map(|_| {
                let key = random_key(ARITY, Some(&stem));

                (key.clone(), rng.gen())
            })
            .collect();

        while kvs.len() < NUM_LEAVES {
            let key = random_key(ARITY, None);
            kvs.insert(key.clone(), rng.gen());
        }

        let keys: Vec<&TestFullKey> = kvs.keys().collect();
        let mut keys2 = keys.clone();
        keys2.shuffle(&mut rng);

        for (k1, k2) in keys.into_iter().zip(keys2.into_iter()) {
            let v1 = kvs.get(k1).unwrap();
            let v2 = kvs.get(k2).unwrap();
            tree1.insert_single(*k1, v1.clone());
            tree2.insert_single(*k2, v2.clone());
        }

        //assert!(tree1 == tree2);

        for k in kvs.keys().collect::<Vec<&TestFullKey>>() {
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

        let key = random_key(3, None);
        let val1: F = rng.gen();
        let val2: F = rng.gen();

        tree.insert_single(key, val1);
        tree.insert_single(key, val2);

        assert!(tree.get_single(&key).unwrap() == &val2);
    }
}
