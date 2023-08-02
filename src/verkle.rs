use std::{
    collections::hash_map::Entry,
    collections::HashMap,
    fmt::Write,
    hash::Hash,
    marker::PhantomData,
    ops::{Deref, Index, IndexMut},
};

use crate::VectorCommitment;

trait Key: Eq + Sized {
    type Unit: Hash + Eq + Clone;
    type Action: ActionKey<Self, Unit = Self::Unit>;

    fn to_action(&self) -> Self::Action;
}

trait ActionKey<K> {
    type Unit;

    fn strip(&mut self);

    fn stripn(&mut self, n: usize);

    fn unstrip(&mut self);

    fn current_unit(&self) -> &Self::Unit;

    fn depth(&self) -> usize;

    /// The number of units left in the entire requested path
    fn units_left(&self) -> usize;

    fn key_up_to_here(&self) -> K;

    fn full_key(&self) -> K;
}

/// A Leaf node is the furthest API that can be reached when traversing the tree.
/// A vanilla implementation simply stores the key, value pair and return its own data.
/// With Ethereum, they have a concept of Extension nodes that enable 256-bit values to be committed to Bandersnatch,
/// and therefore extension nodes may be implemented as leaves
trait Leaf<K: Key, C>: Sized {
    type Value;

    fn new(key: K, value: Self::Value) -> Self;

    fn key(&self) -> &K;

    fn value(&self) -> &Self::Value;

    /// The value that the parent will be committing to for this node.
    /// For an Extension Node, this is the point_to_field of the commitment of its children
    /// For a Leaf, this is the hash of (key, value)
    fn commitment(&self) -> &C;

    ///// If we are a normal leaf, then self.key == key.
    //fn insert(&mut self, key: K, value: Self::Value) -> Option<Self::Value>;
}

#[derive(Clone, PartialEq)]
enum Node<K, C, L, VC>
where
    K: Key,
    VC: VectorCommitment + Clone,
{
    Internal {
        key: K,
        commit: Option<C>,
        vc_data: Option<VC::PreparedData>,
        is_dirty: bool,
        children: HashMap<K::Unit, Node<K, C, L, VC>>,
    },
    Leaf(L),
}

impl<K, C, L, VC> Node<K, C, L, VC>
where
    K: Key,
    L: Leaf<K, C> + Clone,
    VC: VectorCommitment + Clone,
{
    fn new_internal(key: K) -> Self {
        Self::Internal {
            key,
            commit: None,
            vc_data: None,
            is_dirty: true,
            children: HashMap::new(),
        }
    }

    fn new_leaf(l: L) -> Self {
        Self::Leaf(l)
    }

    fn key(&self) -> &K {
        match self {
            Self::Internal { key, .. } => key,
            Self::Leaf(leaf) => leaf.key(),
        }
    }

    fn insert(&mut self, mut q_key: K::Action, node: Self) -> Option<Self> {
        match self {
            Self::Internal {
                children, is_dirty, ..
            } => {
                *is_dirty = true;
                let child = children.entry(q_key.current_unit().clone());
                match child {
                    Entry::Occupied(mut o) => {
                        if q_key.units_left() == 1 {
                            Some(o.insert(node))
                        } else {
                            //o.get_mut().insert(q_key, node)
                            let c = o.get_mut();
                            match c {
                                Self::Internal { .. } => {
                                    q_key.strip();
                                    c.insert(q_key, node)
                                }
                                Self::Leaf(l) => {
                                    let mut l_key = l.key().to_action();

                                    if q_key.full_key() == l_key.full_key() {
                                        return Some(o.insert(node));
                                    }

                                    q_key.strip();
                                    l_key.stripn(q_key.depth());

                                    let mut int_node = Self::new_internal(q_key.key_up_to_here());
                                    int_node.insert(q_key, node);
                                    int_node.insert(l_key, Self::new_leaf(l.clone()));

                                    //children.insert(q_key.current_unit().clone(), int_node);
                                    o.insert(int_node);
                                    None
                                }
                            }
                        }
                    }
                    Entry::Vacant(v) => {
                        v.insert(node);
                        None
                    }
                }
            }

            Self::Leaf(l) => {
                // TODO: Make decision on if leaves will always be pure leaves, or if this is how extension
                //      nodes will be implemented
                panic!("Should not insert on leaf!")
            }
        }
    }

    fn get(&self, mut q_key: K::Action) -> Option<&Self> {
        match self {
            Self::Internal {
                key,
                commit,
                children,
                ..
            } => {
                let child = children.get(q_key.current_unit())?;

                if q_key.units_left() == 1 {
                    return Some(child);
                }

                match child {
                    Self::Internal { .. } => {
                        q_key.strip();
                        child.get(q_key)
                    }
                    Self::Leaf(l) => {
                        if *l.key() == q_key.full_key() {
                            Some(child)
                        } else {
                            None
                        }
                    }
                }
            }

            Self::Leaf(l) => {
                panic!("Should not reach here")
            }
        }
    }

    /// Generate the commitment of this node by recursively updating children that need to update their
    /// commitment
    /// Sets dirty to false
    ///
    fn gen_commitment(&mut self) -> C {
        match self {
            Self::Internal {
                key,
                commit,
                is_dirty,
                vc_data,
                children,
            } => {
                if commit.is_none() || *is_dirty {
                    let child_commits: Vec<(&K::Unit, C)> = children
                        .iter_mut()
                        .map(|(k, v)| (k, v.gen_commitment()))
                        .collect();

                    // If this is the first time initializing our VC data
                    if vc_data.is_none() {
                        for commit in child_commits {
                            
                        }
                    }
                } else {
                    commit.unwrap()
                }
            }
            Self::Leaf(l) => l.commitment().to_owned(),
        }
    }
}

impl<K, C, L, VC> core::fmt::Debug for Node<K, C, L, VC>
where
    C: core::fmt::Debug,
    K: Key + core::fmt::Debug,
    K::Unit: core::fmt::Debug,
    L: core::fmt::Debug,
    VC: VectorCommitment + Clone,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Internal {
                key,
                commit,
                children,
                ..
            } => {
                f.write_fmt(format_args!("Internal ({:?}) {{\n", key))?;
                f.write_fmt(format_args!("  Commit: {:?}\n", commit))?;
                if children.is_empty() {
                    return f.write_char('}');
                }

                f.write_fmt(format_args!("  Children: {{\n"))?;
                for c in children {
                    f.write_fmt(format_args!("    {:?}\n", c.1))?;
                }

                f.write_str("  }\n")?;
            }
            Self::Leaf(l) => {
                f.write_fmt(format_args!("{:?}", l))?;
            }
        }

        f.write_str("}")
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct VanillaLeaf<K, C, T> {
    key: K,
    value: T,
    commit: C,
}

impl<K, C, T> VanillaLeaf<K, C, T>
where
    C: Default,
{
    fn hash() -> C {
        C::default()
    }
}

impl<K, C, T> Leaf<K, C> for VanillaLeaf<K, C, T>
where
    K: Key,
    C: Default,
{
    type Value = T;

    fn new(key: K, value: Self::Value) -> Self {
        Self {
            key,
            value,
            commit: Self::hash(),
        }
    }

    fn key(&self) -> &K {
        &self.key
    }

    fn value(&self) -> &Self::Value {
        &self.value
    }

    fn commitment(&self) -> &C {
        &self.commit
    }
}

impl<U> Key for Vec<U>
where
    U: Hash + Clone + Eq,
{
    type Unit = U;
    type Action = AKey<Self>;

    fn to_action(&self) -> Self::Action {
        AKey::new(self.clone())
    }
}

/*
impl<const N: usize, U> Key for [U; N]
where
    U: Hash + Clone + Eq
{
    type Unit = U;
    type Action = AKey<Self>;

    fn to_action(&self) -> Self::Action {
        AKey::new(self.clone())
    }
}
*/

struct AKey<K> {
    full_key: K,
    cur_index: usize,
}

impl<K> AKey<K> {
    fn new(key: K) -> Self {
        Self {
            full_key: key,
            cur_index: 0,
        }
    }
}

impl<U> ActionKey<Vec<U>> for AKey<Vec<U>>
where
    U: Clone,
{
    type Unit = U;

    fn current_unit(&self) -> &Self::Unit {
        &self.full_key[self.cur_index]
    }

    fn depth(&self) -> usize {
        self.cur_index
    }

    fn strip(&mut self) {
        self.cur_index += 1;
    }

    fn stripn(&mut self, n: usize) {
        self.cur_index += n;
    }

    fn unstrip(&mut self) {
        self.cur_index -= 1;
    }

    fn units_left(&self) -> usize {
        self.full_key.len() - self.cur_index
    }

    fn key_up_to_here(&self) -> Vec<U> {
        self.full_key[0..self.cur_index].to_vec()
    }

    fn full_key(&self) -> Vec<U> {
        self.full_key.clone()
    }
}

#[cfg(test)]
mod tests {
    use rand::{seq::SliceRandom, Rng};

    use super::*;
    use crate::kzg_amortized::KZGAmortized;
    use ark_bn254::Bn254;

    #[derive(Clone, Debug, PartialEq)]
    struct Commit {}
    impl Default for Commit {
        fn default() -> Self {
            Self {}
        }
    }

    type K = Vec<usize>;
    type LeafT = VanillaLeaf<K, Commit, usize>;
    type N = Node<K, Commit, LeafT, KZGAmortized<Bn254>>;

    fn random_leaf(key: K) -> N {
        let val: usize = rand::random();
        N::new_leaf(LeafT::new(key, val))
    }

    fn random_key(len: usize, arity: usize, prefix: Option<&K>) -> K {
        let mut rng = rand::thread_rng();
        let mut res = Vec::new();
        if let Some(stem) = prefix {
            res.extend(stem);
        }
        let real_len = len - res.len();

        for _i in 0..real_len {
            res.push(rng.gen_range(0..arity));
        }

        res
    }

    #[test]
    fn test_insert_get_leaves() {
        let NUM_LEAVES = 50;
        let KEY_LEN = 3;
        let ARITY = 10;

        let mut root1 = N::new_internal(vec![]);
        let mut root2 = N::new_internal(vec![]);

        // 1/4 of keys will share a stem
        let stem = random_key(KEY_LEN - 1, ARITY, None);
        let mut kvs: HashMap<K, N> = (0..NUM_LEAVES / 4)
            .map(|_| {
                let key = random_key(KEY_LEN, ARITY, Some(&stem));

                (key.clone(), random_leaf(key))
            })
            .collect();

        while kvs.len() < NUM_LEAVES {
            let key = random_key(KEY_LEN, ARITY, None);
            kvs.insert(key.clone(), random_leaf(key));
        }

        let keys: Vec<&K> = kvs.keys().collect();
        let mut keys2 = keys.clone();
        keys2.shuffle(&mut rand::thread_rng());

        for (k1, k2) in keys.iter().zip(keys2.iter()) {
            let l1 = kvs.get(k1.clone()).unwrap();
            let l2 = kvs.get(k2.clone()).unwrap();
            root1.insert(k1.to_action(), l1.clone());
            root2.insert(k2.to_action(), l2.clone());
        }

        assert!(root1 == root2);

        for k in keys {
            let get1 = root1.get(k.to_action());
            let get2 = root2.get(k.to_action());

            assert!(get1 == get2);
            assert!(*get1.unwrap() == *kvs.get(k).unwrap());
        }
    }

    #[test]
    fn test_overwrite() {
        let mut root = N::new_internal(vec![]);

        let key = random_key(3, 10, None);
        let val1 = random_leaf(key.clone());
        let val2 = random_leaf(key.clone());

        root.insert(key.to_action(), val1);
        root.insert(key.to_action(), val2.clone());

        assert!(*root.get(key.to_action()).unwrap() == val2);
    }
}
