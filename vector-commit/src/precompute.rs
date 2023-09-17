//! Working in a domain with d-th roots of unity enables a large computational efficiency increase
//! when working with polynomials in evaluation form.

use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::utils::to_usize;

/// Precomputes the evaluations (and inverses) of the derivative of the vanishing polynomial,
/// and the barycentric weights
#[derive(Clone, Debug)]
pub struct PrecomputedLagrange<F: FftField> {
    size: usize,

    unity: F,

    vanishing_evaluations: Vec<F>,

    vanishing_evaluations_inv: Vec<F>,
    //barycentric_weights: [F; N],
}

impl<F: PrimeField> PrecomputedLagrange<F> {
    pub(crate) fn new(size: usize) -> Self {
        let (evals, inv, unity) = Self::compute_vanishing_evaluations(size);
        Self {
            size,
            unity,
            vanishing_evaluations: evals,
            vanishing_evaluations_inv: inv,
        }
    }

    fn compute_vanishing_evaluations(size: usize) -> (Vec<F>, Vec<F>, F) {
        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::new(size).unwrap();
        let mut evals = vec![F::zero(); size];
        let mut inv = vec![F::zero(); size];

        let n_f = F::from(size as u64);
        let unity = domain.group_gen();
        for i in 0..size {
            evals[i] = n_f * unity.pow(&[i as u64]).inverse().unwrap();
            inv[i] = evals[i]; // Batch invert after loop
        }
        batch_inversion(&mut inv);

        (evals, inv, unity)
    }

    pub(crate) fn vanishing_at(&self, point: usize) -> F {
        self.vanishing_evaluations[point]
    }

    pub(crate) fn vanishing_inverse_at(&self, point: usize) -> F {
        self.vanishing_evaluations_inv[point]
    }

    /// Computes the b vector in IPA. When this vector is inner product'd by the evaluations in the domain,
    /// the result is the evaluation F(point).
    ///
    /// b_i = l(point) / l'(x_i)(z-x_i)
    pub(crate) fn compute_barycentric_coefficients(&self, point: F) -> Vec<F> {
        let mut res = vec![F::zero(); self.size];
        if point < F::from(self.size as u64) {
            let point_usize = to_usize(point);
            res[point_usize] = F::one();
            return res;
        }

        // t is the constant outside the summation in PCS multiproofs article
        let t = (point.pow(&[self.size as u64]) - F::one()) / F::from(self.size as u64);
        for i in 0..self.size {
            let pow = self.unity.pow(&[i as u64]);
            res[i] = (t * pow) / (point - pow);
        }

        res
    }
}
