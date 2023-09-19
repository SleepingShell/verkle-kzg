use std::ops::{AddAssign, Index, Mul, MulAssign, Sub};

use ark_ff::{batch_inversion, PrimeField};
use ark_poly::{univariate::DensePolynomial, EvaluationDomain, Evaluations};
use thiserror::Error;

use crate::{
    precompute::PrecomputedLagrange,
    utils::{inner_product, max, to_usize},
};

#[derive(Clone)]
pub struct LagrangeBasis<F: PrimeField, D: EvaluationDomain<F>> {
    /// The evaluations (data) stored in the vector
    evaluations: Evaluations<F, D>,

    /// The highest evaluation point
    max: usize,
}

impl<F: PrimeField, D: EvaluationDomain<F>> LagrangeBasis<F, D> {
    pub fn from_vec(data: Vec<F>) -> Self {
        let len = data.len();
        Self::from_vec_and_domain(data, D::new(len).unwrap())
    }

    pub fn from_vec_and_domain(data: Vec<F>, domain: D) -> Self {
        let max = data.len();
        Self {
            evaluations: Evaluations::from_vec_and_domain(data, domain),
            max,
        }
    }

    pub fn new_zero(size: usize) -> Self {
        Self {
            evaluations: Evaluations::from_vec_and_domain(
                vec![F::zero(); size],
                D::new(size).unwrap(),
            ),
            max: size,
        }
    }

    /// Returns the index of the highest evaluation point (can be smaller than the domain size)
    pub fn max(&self) -> usize {
        self.max - 1
    }

    pub fn elements(&self) -> impl Iterator<Item = &F> {
        self.evaluations.evals.iter()
    }

    pub fn elements_ref(&self) -> &[F] {
        &self.evaluations.evals
    }

    pub fn domain_size(&self) -> usize {
        self.evaluations.domain().size()
    }

    /// Evaluation has 3 paths:
    /// 1. We have stored the evaluation, return this
    /// 2. We have not stored the evaluation, but are within the domain. Return 0
    /// 3. We are outside of the domain. Evaluate using barycentric interpolation
    pub(crate) fn evaluate(&self, precompute: &PrecomputedLagrange<F>, point: F) -> F {
        if point <= F::from(self.max() as u64) {
            self[&point]
        } else if point > F::from(self.max() as u64) && point <= F::from(self.domain_size() as u64)
        {
            F::zero()
        } else {
            self.evaluate_outside_domain(precompute, point)
        }
    }

    pub(crate) fn evaluate_outside_domain(
        &self,
        precompute: &PrecomputedLagrange<F>,
        point: F,
    ) -> F {
        inner_product(
            &self.evaluations.evals,
            &precompute.compute_barycentric_coefficients(point),
        )
    }

    /// Returns w^i
    pub fn index_to_point(&self, index: usize) -> F {
        self.evaluations.domain().element(index)
    }

    /// Compute the quotient polynomial q(x) = [f(X) - f(x_i)] / [X-x_i]
    pub(crate) fn divide_by_vanishing(
        &self,
        precompute: &PrecomputedLagrange<F>,
        index: usize,
    ) -> Vec<F> {
        let mut q = vec![F::zero(); self.evaluations.domain().size()];
        let index_f = self.index_to_point(index);
        let eval = if index >= self.max {
            F::zero()
        } else {
            self[index]
        };
        let index_vanishing = precompute.vanishing_at(index);

        for i in 0..self.evaluations.domain().size() {
            if i == index {
                continue;
            }
            let i_f = self.index_to_point(i);
            let i_eval = if i >= self.max { F::zero() } else { self[i] };

            let sub = i_eval - eval;
            q[i] = sub / (i_f - index_f);
            q[index] +=
                sub * index_vanishing * precompute.vanishing_inverse_at(i) / (index_f - i_f);
        }

        q
    }

    pub(crate) fn divive_by_vanishing_outside_domain(
        &self,
        precompute: &PrecomputedLagrange<F>,
        point: F,
    ) -> Vec<F> {
        let mut q = vec![F::zero(); self.domain_size()];
        //let eval = self.evaluate_outside_domain(precompute, point);
        let eval = self.evaluate(precompute, point);

        let mut inversions = vec![F::zero(); self.domain_size()];
        for i in 0..self.domain_size() {
            inversions[i] = self.index_to_point(i) - point;
        }
        batch_inversion(&mut inversions);

        for i in 0..self.domain_size() {
            let i_eval = if i >= self.max { F::zero() } else { self[i] };
            q[i] = (i_eval - eval) * inversions[i];
        }

        q
    }

    /// Rarely would we want to go into coefficient form, but some computational speedups
    /// (i.e amortized KZG) only work in coefficient form
    pub(crate) fn interpolate(&self) -> DensePolynomial<F> {
        self.evaluations.interpolate_by_ref()
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> Index<usize> for LagrangeBasis<F, D> {
    type Output = F;

    fn index(&self, index: usize) -> &Self::Output {
        &self.evaluations[index]
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> Index<&F> for LagrangeBasis<F, D> {
    type Output = F;

    fn index(&self, index: &F) -> &Self::Output {
        &self.evaluations[to_usize(index)]
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> AddAssign<&Self> for LagrangeBasis<F, D> {
    fn add_assign(&mut self, rhs: &Self) {
        self.evaluations += &rhs.evaluations;
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> Sub for LagrangeBasis<F, D> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            evaluations: &self.evaluations - &rhs.evaluations,
            max: *max(&self.max, &rhs.max),
        }
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> Mul<F> for &LagrangeBasis<F, D> {
    type Output = LagrangeBasis<F, D>;
    fn mul(self, rhs: F) -> Self::Output {
        LagrangeBasis {
            evaluations: &self.evaluations * rhs,
            max: self.max,
        }
    }
}

impl<F: PrimeField, D: EvaluationDomain<F>> MulAssign<F> for LagrangeBasis<F, D> {
    fn mul_assign(&mut self, rhs: F) {
        self.evaluations = &self.evaluations * rhs;
    }
}

#[derive(Error, Clone, Debug)]
pub enum LagrangeError {}
