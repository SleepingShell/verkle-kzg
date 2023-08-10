use std::{collections::HashMap, error::Error, fmt::Display, marker::PhantomData};

use ark_ec::{pairing::Pairing, Group};
use ark_ff::{Field, One, PrimeField, UniformRand, Zero};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Evaluations,
    GeneralEvaluationDomain, Polynomial,
};

use crate::{VCPreparedData, VCUniversalParams, VectorCommitment};

/// KZGKey represents the universal parameters, AKA reference string, for both
/// committing polynomials and verifying commitments
#[derive(Clone, Debug)]
pub struct KZGKey<F: PrimeField, G1: Group, G2: Group> {
    /// The max degree this reference string supports
    degree: usize,

    /// Reference string in G1 group
    ref_string_g1: Vec<G1>,

    /// For G2, we only need α*g
    g2_secret: G2,

    /// The lagrange polynomials are created and commited to.
    lagrange_commitments: Vec<G1>,

    /// domain.size / unity^i
    //unity_neg: Vec<F>,
    domain: GeneralEvaluationDomain<F>,
}

impl<F, G1, G2> KZGKey<F, G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>,
{
    fn new_from_secret(secret: F, max_degree: usize) -> Self {
        let g1 = G1::generator();
        let g2 = G2::generator();
        let domain = GeneralEvaluationDomain::<F>::new(max_degree).unwrap();

        let mut params: Self = Self {
            degree: max_degree,
            ref_string_g1: vec![],
            g2_secret: g2 * secret,
            lagrange_commitments: vec![], //TODO
            domain: domain,
        };
        params.ref_string_g1.push(g1);
        let mut sec_cur: F = F::one();
        for _i in 1..max_degree {
            sec_cur = sec_cur * secret; //α^i
            params.ref_string_g1.push(g1 * sec_cur);
        }
        //params.lagrange_commitments = domain.evaluate_all_lagrange_coefficients(secret).iter().map(|l| g1 * l).collect();
        params.lagrange_commitments = domain.ifft(&params.ref_string_g1);

        params
    }

    /// Multiplies the array of coefficients with the reference string in G1
    /// and sums the group elements.
    fn commit_g1(&self, coeffs: &[F]) -> G1 {
        let mut sum = G1::zero();
        for (i, c) in coeffs.iter().enumerate() {
            sum += self.element_at_g1(i) * c;
        }

        sum
    }

    fn commit_lagrange(&self, evaluations: &[F]) -> G1 {
        let mut sum = G1::zero();
        for (i, e) in evaluations.iter().enumerate() {
            sum += self.lagrange_commitments[i] * e;
        }

        sum
    }

    fn element_at_g1(&self, index: usize) -> G1 {
        self.ref_string_g1[index]
    }

    /// We only care about the α*g element in G2
    fn get_g2(&self) -> G2 {
        self.g2_secret
    }

    pub fn domain(&self) -> GeneralEvaluationDomain<F> {
        self.domain
    }

    fn unity(&self) -> F {
        self.domain.group_gen()
    }
}

impl<F, G1, G2> VCUniversalParams for KZGKey<F, G1, G2>
where
    F: PrimeField,
    G1: Group<ScalarField = F>,
    G2: Group<ScalarField = F>,
{
    fn max_size(&self) -> usize {
        self.degree
    }
}

/// KZGPreparedData will create a polynomial to commit to by interpolating the dataset.
/// Instead of the evaluations occuring at [0,1,...,n-1] they instead occur at [0,w,w^2,...w^n-1].
/// w is the n-th root of unity
#[derive(Clone, PartialEq)]
pub struct KZGPreparedData<F: PrimeField> {
    /// The evaluation domain of the dataset, aka the points that we will run polynomial evaluation at
    evaluations: Evaluations<F, GeneralEvaluationDomain<F>>,

    /// Because evaluations may be sparse, we cache which entries are filled
    filled_index: Vec<usize>,

    /// TODO: Remove
    poly: DensePolynomial<F>,

    /// The size of the data-set (not necessairly equal to n)
    size: usize,
}

impl<F: PrimeField> KZGPreparedData<F> {
    pub fn from_points_and_domain(mut points: Vec<F>, domain: GeneralEvaluationDomain<F>) -> Self {
        let len = points.len();
        if len < domain.size() {
            points.resize(domain.size(), F::zero());
        }
        let evals = Evaluations::from_vec_and_domain(points, domain);
        let poly = evals.interpolate_by_ref();

        let mut filled: Vec<usize> = Vec::new();
        for i in 0..len {
            filled.push(i);
        }

        Self {
            evaluations: evals,
            filled_index: filled,
            poly: poly,
            size: len,
        }
    }
    fn domain_size(&self) -> usize {
        self.evaluations.domain().size()
    }

    fn evaluate(&self, index: usize) -> F {
        self.evaluations[index]
    }

    /// Returns w^i
    fn index_to_point(&self, index: usize) -> F {
        self.evaluations.domain().element(index)
    }

    fn poly(&self) -> DensePolynomial<F> {
        self.evaluations.clone().interpolate()
    }
}

impl<F: PrimeField> VCPreparedData for KZGPreparedData<F> {
    type Item = F;
    type Error = KZGError;

    /// FIXIME: HACKY
    fn from_vec(data: Vec<Self::Item>) -> Self {
        let domain = GeneralEvaluationDomain::<F>::new(data.len()).unwrap();
        Self::from_points_and_domain(data, domain)
    }

    fn set_evaluation(&mut self, index: usize, value: Self::Item) -> Result<(), Self::Error> {
        // TODO: Domain expansion, although in Verkle case the domain size is the arity of the tree.
        if index > self.domain_size() {
            return Err(KZGError::OutOfDomainBounds);
        }
        self.evaluations.evals[index] = value;
        return Ok(());
    }

    fn get(&self, index: usize) -> Option<&Self::Item> {
        self.evaluations.evals.get(index)
    }

    fn get_all(&self) -> Vec<(usize, Self::Item)> {
        let mut res: Vec<(usize, Self::Item)> = Vec::new();
        for i in self.filled_index.iter() {
            res.push((*i, self.evaluations.evals[*i]));
        }
        res
    }

    fn max_size(&self) -> usize {
        self.domain_size()
    }
}

#[derive(PartialEq, Clone, Default, Debug)]
pub struct KZGCommitment<G: Group> {
    commit: G,
}

impl<G: Group> KZGCommitment<G> {
    fn new(commit: G) -> Self {
        Self { commit }
    }

    fn group_to_field(&self) -> G::ScalarField {
        if self.commit.is_zero() {
            <G::ScalarField as Field>::ZERO
        } else {
            let mut bytes: Vec<u8> = Vec::new();
            // TODO: Check
            self.commit.serialize_compressed(&mut bytes).unwrap();
            <G::ScalarField as PrimeField>::from_le_bytes_mod_order(&bytes)
        }
    }
}

pub struct KZGProof<F: PrimeField, G: Group> {
    commit: G,

    /// index is the position in the vector, NOT the evaluation point
    index: usize,

    /// Evaluation point
    /// TODO: pass root of unity with proof/batch proof?
    point: F,

    data: F,
}

pub struct KZGBatchProof<F: PrimeField, G: Group> {
    //commit: G,

    // For now, there is no proof compression. This is a map of index -> (eval point, output, proof)
    proofs: HashMap<usize, (F, F, G)>,
}

impl<F: PrimeField, G: Group> Default for KZGBatchProof<F, G> {
    fn default() -> Self {
        Self {
            proofs: HashMap::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub enum KZGError {
    DefaultError,
    DataExceedsMaxSize,
    InvalidDomain,
    OutOfDomainBounds,
}

impl KZGError {
    fn new() -> Self {
        Self::DefaultError
    }
}

impl Display for KZGError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", "")
    }
}

impl Error for KZGError {}

/// Implementation of the Feist-Khovratovich technique of "Fast Amortized KZG proofs".
#[derive(PartialEq, Clone)]
pub struct KZGAmortized<E: Pairing> {
    _engine: PhantomData<E>,
}

impl<E: Pairing> VectorCommitment for KZGAmortized<E> {
    type UniversalParams = KZGKey<E::ScalarField, E::G1, E::G2>;
    type PreparedData = KZGPreparedData<E::ScalarField>;
    type Commitment = KZGCommitment<E::G1>;
    type Proof = KZGProof<E::ScalarField, E::G1>;
    type BatchProof = KZGBatchProof<E::ScalarField, E::G1>;
    type Error = KZGError;

    fn setup<R: rand::RngCore>(
        max_items: usize,
        rng: &mut R,
    ) -> Result<Self::UniversalParams, Self::Error> {
        let secret = <E::ScalarField as UniformRand>::rand(rng);
        Ok(KZGKey::new_from_secret(secret, max_items))
    }

    fn commit(
        key: &Self::UniversalParams,
        data: &Self::PreparedData,
    ) -> Result<Self::Commitment, Self::Error> {
        Ok(KZGCommitment {
            commit: key.commit_lagrange(&data.evaluations.evals),
        })
    }

    fn prove(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        index: usize,
        data: &Self::PreparedData,
    ) -> Result<Self::Proof, Self::Error> {
        let commit = Self::single_proof(&key, &data, index)?;
        Ok(Self::Proof {
            commit,
            index,
            data: data.evaluate(index),
            point: data.index_to_point(index),
        })
    }

    fn prove_all(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        data: &Self::PreparedData,
    ) -> Result<Self::BatchProof, Self::Error> {
        Self::get_all_proofs(&key, &data)
    }

    fn verify(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        // TODO: Need to be able to access domain's root of unity to be able to evaluate properly
        let pairing1 = E::pairing(
            proof.commit,
            key.get_g2() - (E::G2::generator() * proof.point),
        );
        let pairing2 = E::pairing(
            commitment.commit - E::G1::generator() * proof.data,
            E::G2::generator(),
        );

        Ok(pairing1 == pairing2)
    }

    fn verify_batch(
        key: &Self::UniversalParams,
        commitment: &Self::Commitment,
        proof: &Self::BatchProof,
    ) -> Result<bool, Self::Error> {
        let s = key.get_g2();
        for (index, p) in proof.proofs.iter() {
            let pairing1 = E::pairing(p.2, s - (E::G2::generator() * p.0));
            let pairing2 = E::pairing(
                commitment.commit - E::G1::generator() * p.1,
                E::G2::generator(),
            );

            if pairing1 != pairing2 {
                return Ok(false);
            }
        }

        Ok(true)
    }

    fn convert_commitment_to_data(
        commit: &Self::Commitment,
    ) -> <Self::PreparedData as VCPreparedData>::Item {
        commit.group_to_field()
    }
}

impl<E: Pairing> KZGAmortized<E> {
    /// Get all the proofs for every coefficient of `data`'s polynomial representation
    fn get_all_proofs(
        key: &KZGKey<E::ScalarField, E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>,
    ) -> Result<KZGBatchProof<E::ScalarField, E::G1>, KZGError> {
        let all_proofs = Self::build_witness_matrix(key, data)?;
        let mut res = KZGBatchProof::default();

        let domain = data.evaluations.domain();
        let proofs = domain.fft(&all_proofs);

        for index in 0..data.size {
            let point = data.index_to_point(index);
            res.proofs
                .insert(index, (point, data.evaluate(index), proofs[index]));
        }

        Ok(res)
    }

    /// Build the group elements from the FFT of the polynomial coefficients multiplied by the reference string
    fn build_witness_matrix(
        key: &KZGKey<E::ScalarField, E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>,
    ) -> Result<Vec<E::G1>, KZGError> {
        let degree = data.poly.degree();
        let coeffs = data.poly.coeffs();
        let domain = GeneralEvaluationDomain::<E::ScalarField>::new(degree * 2)
            .ok_or(KZGError::InvalidDomain)?;
        let domain_size = domain.size();

        let mut c_hat: Vec<E::ScalarField> = vec![coeffs[degree]];
        c_hat.extend(vec![E::ScalarField::zero(); degree + 1]);
        c_hat.extend(&coeffs[0..degree]);

        let mut s_hat = key.ref_string_g1[0..degree].to_vec();
        s_hat.reverse();
        s_hat.extend(vec![E::G1::zero(); domain_size - degree]);

        let y = domain.fft(&c_hat);
        let v = domain.fft(&s_hat);
        let u = Self::elementwise_mul(&v, &y);
        let h_hat = domain.ifft(&u);

        Ok(h_hat[0..degree].to_vec())
    }

    fn single_proof(
        key: &KZGKey<E::ScalarField, E::G1, E::G2>,
        data: &KZGPreparedData<E::ScalarField>,
        index: usize,
    ) -> Result<E::G1, KZGError> {
        let data_piece = data.evaluate(index);
        let point = data.index_to_point(index);
        let mut poly = data.poly.clone();
        poly -= &DensePolynomial::<E::ScalarField>::from_coefficients_slice(&[data_piece]);

        let divisor = DensePolynomial::<E::ScalarField>::from_coefficients_slice(&[
            E::ScalarField::zero() - point,
            E::ScalarField::one(),
        ]);
        let q = &poly / &divisor;

        let commit = key.commit_g1(q.coeffs());

        Ok(commit)
    }

    fn elementwise_mul<F: PrimeField, G1: Group<ScalarField = F>>(
        g_vec: &Vec<G1>,
        f_vec: &Vec<F>,
    ) -> Vec<G1> {
        let mut result: Vec<G1> = vec![];
        for (g, f) in g_vec.iter().zip(f_vec.iter()) {
            result.push(*g * *f);
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bn254::Bn254;

    type F = <Bn254 as Pairing>::ScalarField;
    type G1 = <Bn254 as Pairing>::G1;
    type G2 = <Bn254 as Pairing>::G2;
    type KZG = KZGAmortized<Bn254>;

    const DATA_SIZE: usize = 4;
    const MAX_CRS: usize = 32;

    fn gen_data(num: usize) -> Vec<F> {
        let mut data: Vec<F> = vec![];
        let mut rng = rand::thread_rng();
        for _i in 0..num {
            data.push(F::rand(&mut rng));
        }
        data
    }

    fn setup(n: usize, max_degree: usize) -> (KZGPreparedData<F>, KZGKey<F, G1, G2>) {
        let data = gen_data(n);
        //let prep = KZGPreparedData::from_iter(data);
        let crs = KZG::setup(max_degree, &mut rand::thread_rng()).unwrap();
        let prep = KZGPreparedData::from_points_and_domain(data, crs.domain);

        (prep, crs)
    }

    #[test]
    fn test_single_proof() {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();

        for i in 0..DATA_SIZE {
            let proof = KZG::prove(&crs, &commit, i, &data).unwrap();
            assert!(KZG::verify(&crs, &commit, &proof).unwrap());
        }
    }

    #[test]
    fn test_batch_proof() {
        let (data, crs) = setup(DATA_SIZE, MAX_CRS);
        let commit = KZG::commit(&crs, &data).unwrap();

        let proofs = KZG::prove_all(&crs, &commit, &data).unwrap();
        assert!(KZG::verify_batch(&crs, &commit, &proofs).unwrap());
    }

    fn vec_to_str<T: Display>(v: &Vec<T>) -> String {
        let mut s = "[\n".to_owned();
        for e in v {
            s.push_str(&format!("\t{}\n", e));
        }
        s.push_str("]");
        s
    }
}
