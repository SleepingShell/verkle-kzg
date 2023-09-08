//! Working in a domain with d-th roots of unity enables a large computational efficiency increase
//! when working with polynomials in evaluation form.

use ark_ff::{BigInteger, FftField, Field, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

/// Precomputes the evaluations (and inverses) of the derivative of the vanishing polynomial,
/// and the barycentric weights
pub(crate) struct PrecomputedLagrange<const N: usize, F: FftField> {
    unity: F,

    vanishing_evaluations: [F; N],

    vanishing_evaluations_inv: [F; N],
    //barycentric_weights: [F; N],
}

impl<const N: usize, F: PrimeField> PrecomputedLagrange<N, F> {
    pub(crate) fn new() -> Self {
        let (evals, inv, unity) = Self::compute_vanishing_evaluations();
        Self {
            unity,
            vanishing_evaluations: evals,
            vanishing_evaluations_inv: inv,
        }
    }

    fn compute_vanishing_evaluations() -> ([F; N], [F; N], F) {
        let domain: GeneralEvaluationDomain<F> = GeneralEvaluationDomain::new(N).unwrap();
        let mut evals = [F::zero(); N];
        let mut inv = [F::zero(); N];

        let n_f = F::from(N as u64);
        let unity = domain.group_gen();
        for i in 0..N {
            evals[i] = n_f * unity.pow(&[i as u64]).inverse().unwrap();
            inv[i] = evals[i].inverse().unwrap();
        }

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
    pub(crate) fn compute_barycentric_coefficients(&self, point: F) -> [F; N] {
        let mut res = [F::zero(); N];
        if point < F::from(N as u64) {
            //FIXME: THIS IS SO BAD OH MY GOD
            let mut point_usize = 0usize;
            point
                .into_bigint()
                .to_bytes_be()
                .into_iter()
                .for_each(|b| point_usize = (point_usize << 8) | b as usize);
            res[point_usize] = F::one();
            return res;
        }

        let t = (point.pow(&[N as u64]) - F::one()) / F::from(N as u64);
        for i in 0..N {
            let pow = self.unity.pow(&[i as u64]);
            res[i] = (t * pow) / (point - pow);
        }

        res
    }
}
