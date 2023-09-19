This repository hosts a library that provides a VectorCommitment and several vector commitment scheme implementations, and a library for a verkle tree implementation that will become compatible with [Ethereum's spec](https://notes.ethereum.org/@vbuterin/verkle_tree_eip).

The verkle tree is similar to a merkle tree in that it is a tree structure that commits to its children, but instead of using hashes to provide binding, it uses a vector commitment scheme. The verkle tree itself is also a vector commitment scheme.

# Cryptography Documentation
A Vector Commitment Scheme (VCS) enables one to commit to a vector of elements and then generate proofs of inclusion in regard to that commitment. The security of these schemes binds the vector of elements being proven to the initial commitment, so verifiers can be sure that proofs indeed prove that an element was contained in the commitment.

Many (and the ones in this repository) of VCS are built on Polynomial Commitment Schemes (PCS), utilizing the properties of polynomials to achieve the above goals.

Polynomials are usually presented in coefficient form $f(x) = a_0 + a_1x + a_2x^2...a_dx^{d-1}$. However, for the VCS it is more natural, and in fact more efficient, to represent these polynomials in *evaluation* form, which is to say simply the vector of data itself: $\overrightarrow{v} = <y_0, y_1,...,y_d>$ where we assume the domain (points of evaluation) are simply incrementing starting from 0, so $f(x_0) = y_0$. One can use Lagrange interpolation [[2]]([2]) to find the polynomial that intersects all points in $v$ with degree $d-1$. When the points of evaluation are at the $d$th roots of unity [[3]]([3]), the FFT [[4]]([4]) can be used to efficiently switch between these forms.

**A note on polynomial division:** A polynomial is divisible by $X-z$ if it has a zero at $z$, and the converse is true in that if a polynomial is divisible by $X-z$ with no remainder, than it is zero at $z$.

## Verkle Tree

## Vector Commitments (Underlying)
This repository contains two implementations of base vector commitments: IPA and KZG. Here is a nice comparison of their tradeoffs, taken from Dankrad [[1]]([1]).
| | Pedersen+IPA | KZG |
---|---|---|
Assumption | Discrete log | Bilinear group
Trusted Setup | No | Yes
Commitment size | 1 Group element | 1 Group Element
Proof Size | $2log_2(n)$ Group elements | 1 Group Element
Verification | $O(n)$ group operations | 1 Pairing  

# IPA
The Commom Reference String (CRS) of the IPA scheme is a set of ECC points in which the discrete log relation between them is unknown $\overrightarrow{g} = <g_0, g_1, ..., g_d>$ in addition to another random point $q$. You can think of $q$ as the "evaluation generator" in which the actual evaluation (piece of data) is committed to during proving. $$CRS = (\overrightarrow{g}, q)$$

Note that as the table above mentions, this **not** a trusted setup, as there must be no structure between these points (unlike KZG). For example, we could use a hash function that continously hashes its own output, and serialize each output to a ECC point. Given a sufficiently large curve, it should be probabilisticaly impossible to find $x$ such that $g_1 = xg_0$.

Let $\cdot$ refer to the inner product of two vectors, i.e $\overrightarrow{a} \cdot \overrightarrow{b} = \sum(a_i * b_i)$.

When working in evaluation form, we create a commitment to the dataset in the form of $$C = \overrightarrow{a} \cdot \overrightarrow{g}$$

Now, we want to prove that the element $y$ at index $z$ is indeed contained in commitment $C$. The bulk of this proving can be seen in the `low_level_ipa()` function. In order to prevent attackers from crafting their commitment to commit to false data, we utilize the Fiat-Shamir heuristic [[5]]([5]) to create a non-interactive protocol that produces a challenge $w$ based on the commitment, $z$ and $y$. This scales the evaluation generator $q$ to $q=wq$. 

Dankrad's article [[1]]([1]) does a much better job explaining the actual inner product argument than I can. Just note that because we are working in evaluation form, there is no commitment to $(1, z, z^2,...)$. As such, when we are generating our proof, we must generate the *barycentric coefficients* $\overrightarrow{b}$. These can take either two forms:
1. When we are proving a point inside of our domain (dataset), it is simply the vector of $0$s except for the index of $z$ set to one (as we only want the evaluation of that point). $\overrightarrow{a} \cdot \overrightarrow{b} = y$.
2. When we are proving outside of the domain, in the case of multiproofs, we utilize the barycentric coefficients from [[2]]([2]). These coefficients can be thought of as interpolating our domain of all $1$ values. When this is inner producted with our actual dataset, it produces what would be the output of the lagrange polynomial that interpolates our dataset. In otherwords, we can evaluate outside of our domain without actually interpolating a polynomial in coefficient form (which is computationally expensive)!

# KZG
Let $[s]_1$ refer to scalar multiplication of the generator ($G$) in G1 and $[s]_2$ to scalar multiplication of the generator ($H$) in G2. $[]_T$ refers to the pairing group

The CRS of the KZG scheme is of the form of a **structured** reference string (SRS), a set of ECC points in which they have the structure of a polynomial, i.e for some $\alpha$, $\overrightarrow{g} = g, [\alpha]_1, [\alpha^2]_1,...,[\alpha^{d-1}]_1$. The $\alpha$ is referred to as **toxic waste** and is the reason that this scheme requires a trusted setup. Trusted setups are usually ran as multi-party-computation (MPC) ceremonies, in which as long as a single participant destroys their secret the SRS can be considered secure (under its cryptographic assumptions). For example, Ethereum's KZG ceremony for EIP-4844 [[6]]([6]) saw over 140,000 participants! As long a single participant did not save their toxic waste, we can consider the SRS not compromised. We also need to commit to $\alpha$ in the $G2$ of our pairing function. So $$CRS = (\overrightarrow{g}, [\alpha]_2)$$.

With the SRS being structured in the form of a polynomial evaluation, what this enables us to do is to evaluate our data (which remember, **is** a polynomial, just in its evaluation form) at a random, unknown point. Because this point $\alpha$ is unknown to provers, they cannot use interpolation to find a polynomial that evaluates to a different set of data than their commitment $C$.

For efficiency reasons, instead of committing to $\alpha g$ (and their powers), we commit to the lagrange polynomials $l(\alpha)g$. This enables use the KZG scheme without having to switch to coefficient form! Similar to IPA, we can commit to our dataset $\overrightarrow{a}$ by $$C = \overrightarrow{a} \cdot \overrightarrow{g}$$

The proof of inclusion for KZG is simple: we just need to prove the following: $$q(X) = \frac{p(X) - y}{X-z}$$

Why? Remember the note on polynomial division above, a polynomial is only divisble by $(X-z)$ if it has a zero at $z$. And how do we convert our dataset to have a zero at $z$? We simply subtract it by $y$!. Similar to our commitment $C$, the proof ($\pi$) is the inner product of $q(X)$ with our SRS $$\pi = \overrightarrow{q} \cdot \overrightarrow{g}$$

Again, refer to Dankrad [[6]](6) for a better explanation, but the crux of verifying this proof is that the following must hold: $$e(\pi, [s-z]_2) = e(C - [y]_1, H)$$

Writing this out in terms of the pairing group, it becomes: $$[q(\alpha) * (\alpha-z)]_T = [p(\alpha) - y]_T$$

What this is proving is that the following equation holds, evaluated at the secret $\alpha$ point: $$q(X)(X-z) = p(X)-y$$


# References
<a id="1">[1]</a> https://dankradfeist.de/ethereum/2021/07/27/inner-product-arguments.html
<a id="2">[2]</a> https://en.wikipedia.org/wiki/Lagrange_polynomial
<a id="3">[3]</a> https://en.wikipedia.org/wiki/Root_of_unity
<a id="4">[4]</a> https://en.wikipedia.org/wiki/Fast_Fourier_transform
<a id="5">[5]</a> https://en.wikipedia.org/wiki/Fiat%E2%80%93Shamir_heuristic
<a id="6">[6]</a> https://ethereum.org/en/roadmap/danksharding/
<a id="7">[7]</a> https://dankradfeist.de/ethereum/2020/06/16/kate-polynomial-commitments.html
<a id="8">[8]</a> https://dankradfeist.de/ethereum/2021/06/18/pcs-multiproofs.html