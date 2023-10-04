use std::{
    collections::{hash_map::Entry, HashMap},
    fmt::Debug,
    hash::Hash,
};

use bytemuck::Pod;
use num::{One, Zero};
use thiserror::Error;
use vector_commit::{VCCommitment, VCData, VectorCommitment};

use crate::{Key, KeyMethods, SplittableValue};

#[derive(Error, Debug)]
pub(crate) enum VerkleError {
    #[error("Invalid path requested")]
    InvalidPath,
}

/// The Node provides the recursive structure of the Verkle Tree.
///
/// Both node types store their optional cached commitments.
/// The below table provides some more information
/// on the node types:
/// | Type      | Description |
/// |-----------|-------------|
/// | Internal  | Stores children that are other Internal OR Extension nodes |
/// | Extension | Stores the actual data of the tree. Additionally stores the stem that all leaves of this node share |
///
/// Methods are implemented for:
/// - Storing and retrieving values
/// - Grabbing commitments
///
/// All of these methods are called recursively on children when required
#[derive(Clone)]
pub(crate) enum Node<const N: usize, K, VC, T>
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
    pub(crate) fn new_extension(stem: Key<N, K>, values: Vec<(K, T)>) -> Self {
        Self::Extension {
            stem,
            commit: None,
            leaves: values.into_iter().collect(),
        }
    }

    pub(crate) fn new_internal(nodes: Vec<(K, Self)>) -> Self {
        Self::Internal {
            commit: None,
            children: nodes.into_iter().map(|v| (v.0, v.1)).collect(),
        }
    }

    /// Get the extension node equal to this stem. If it does not exist, return None
    pub(crate) fn get_stem(&self, stem: &Key<N, K>, mut cur_depth: usize) -> Option<&Self> {
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

    /// Recursively finds the path that must be taken in the tree to get the `stem`
    ///
    /// Each node in the path will store:
    /// - The prefix of the key at the node
    /// - The index of itself in regard to it's parent
    /// - A reference to self
    pub(crate) fn path_to_stem<'a>(
        &'a self,
        stem: &Key<N, K>,
        path: &mut Vec<(Vec<K>, K, &'a Self)>,
    ) -> Result<(), VerkleError> {
        match self {
            Self::Extension { .. } => Ok(()),

            Self::Internal { children, .. } => {
                let depth = path.len();
                if let Some(child) = children.get(&stem[depth]) {
                    path.push((stem[0..depth + 1].to_vec(), stem[depth], self));
                    child.path_to_stem(stem, path)
                } else {
                    Err(VerkleError::InvalidPath)
                }
            }
        }
    }

    /// Gets the value from an extension node.
    ///
    /// ! Panics if called on an internal node
    pub(crate) fn get_value(&self, unit: K) -> Option<&T> {
        match self {
            Self::Extension { leaves, .. } => leaves.get(&unit),
            _ => panic!("Called get_value on non-extension node"),
        }
    }

    /// Recursively insert the `values` into an extension node that has its stem equal to `stem`
    /// This function will **clear the stored commitment** of all touched nodes, as they are no longer valid
    pub(crate) fn insert(&mut self, stem: Key<N, K>, values: Vec<(K, T)>, cur_depth: usize) {
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
    ///
    /// An Internal node will simply convert all children commitments to data items, and commit to that array
    ///
    /// The Extension node will generate commitments according to Ethereum's standard. This entails the following:
    /// 1. All leaves are split into their upper and lower halves (from `SplittableValue`)
    /// 2. The first half of leaves are committed to by `c1` with upper half preceding the lower half
    ///     - E.g a leaf at index 1 will set the `VC::Data` array with [..., leaf_1_low, leaf_1_upper, ...]
    /// 3. The upper half of leaves commit to `c2` using the same rules as above
    /// 4. `c1` and `c2` are encoded as `VC::Data::Item`s
    /// 5. The stem is encoded as a data item
    /// 6. Commit to the 4 data item array: `[1, stem, c1, c2]`
    pub(crate) fn gen_commitment(
        &mut self,
        crs: &VC::UniversalParams,
    ) -> Result<&VC::Commitment, VC::Error> {
        match self {
            Self::Extension {
                stem,
                commit,
                leaves,
            } => {
                if commit.is_some() {
                    return Ok(commit.as_ref().unwrap());
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

                // HACK FIXME THIS IS HARDCODED NEED TO FIX
                let mut vc_vec = vec![<VC::Data as VCData>::Item::zero(); 256];
                for (&k, child) in children.iter_mut() {
                    let cc = child.gen_commitment(crs)?;
                    vc_vec[k.into()] = cc.to_data_item();
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
                    f.write_fmt(format_args!("\t\t({:?}) {:?}\n", c.0, c.1))?;
                }
                f.write_str("\t}")
            }
        }
    }
}
