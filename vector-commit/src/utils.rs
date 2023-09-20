use std::{
    iter::Sum,
    ops::{Add, Mul},
};

use ark_ff::{batch_inversion, Field, One, PrimeField};
use ark_serialize::CanonicalSerialize;

pub(crate) fn serialize<T: CanonicalSerialize>(x: &T) -> Vec<u8> {
    let mut b = Vec::new();
    x.serialize_compressed(&mut b);

    b
}

pub(crate) fn inner_product<R: Copy, T: Mul<R, Output = T> + Sum<T> + Copy>(a: &[T], b: &[R]) -> T {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum()
    //b.iter().enumerate().map(|(i, r)| a[i] * *r).sum()
}

pub(crate) fn elementwise_mul<'a, R, T: Mul<&'a R, Output = T> + Copy>(
    a: &[T],
    b: &'a [R],
) -> Vec<T> {
    let mut res: Vec<T> = Vec::with_capacity(a.len());
    a.iter().zip(b.iter()).for_each(|(a, b)| res.push(*a * b));
    res
}

//res_i = a_i + x*b_i
pub(crate) fn vec_add_and_distribute<R: Copy, T: Copy + Add<T, Output = T> + Mul<R, Output = T>>(
    a: &[T],
    b: &[T],
    x: R,
) -> Vec<T> {
    assert!(a.len() == b.len());
    a.iter().zip(b.iter()).map(|(a, b)| *a + (*b * x)).collect()
}

pub(crate) fn split<T: Clone>(a: &[T]) -> (Vec<T>, Vec<T>) {
    (a[0..a.len() / 2].to_vec(), a[a.len() / 2..].to_vec())
}

pub(crate) fn powers_of<T: Mul<T, Output = T> + One + Copy>(a: T, n: usize) -> Vec<T> {
    let mut res = Vec::with_capacity(n);
    let mut cur = T::one();
    res.push(cur);

    (1..n).for_each(|_| {
        cur = cur * a;
        res.push(cur);
    });

    res
}

pub(crate) fn invert_domain_at<F: Field>(t: F, N: usize) -> Vec<F> {
    let mut res: Vec<F> = (0..N as u64).map(|i| t - F::from(i)).collect::<Vec<F>>();
    batch_inversion(&mut res);

    res
}

pub(crate) fn max<'a, T: Ord>(l: &'a T, r: &'a T) -> &'a T {
    if l < r {
        r
    } else {
        l
    }
}

pub(crate) fn to_usize<T: PrimeField>(x: &T) -> usize {
    x.into_bigint().as_ref()[0] as usize
}
