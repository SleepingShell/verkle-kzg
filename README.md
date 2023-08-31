This repository hosts a library that provides a VectorCommitment and several vector commitment scheme implementations, and a library for a verkle tree implementation that will become compatible with [Ethereum's spec](https://notes.ethereum.org/@vbuterin/verkle_tree_eip).

The verkle tree is similar to a merkle tree in that it is a tree structure that commits to its children, but instead of using hashes to provide binding, it uses a vector commitment scheme. The verkle tree itself is also a vector commitment scheme.

Below WIP documentation
# IPA
## Evaluation Form
- Used for efficient single proofs

## Polynomial Form
- Used when batch proofs will be required

# KZG