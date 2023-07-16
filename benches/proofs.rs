use rand::Rng;
use verkle_kzg::{kzg_amortized::*,VectorCommitment};

use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use criterion::{Criterion, criterion_group, criterion_main};
use ark_bn254::Bn254;

type F = <Bn254 as Pairing>::ScalarField;
type G1 = <Bn254 as Pairing>::G1;
type G2 = <Bn254 as Pairing>::G2;
type KZG = KZGAmortized<Bn254>;

const DATA_SIZE: usize = 20;
const MAX_CRS: usize = 32;

fn gen_data(num: usize) -> Vec<F> {
    let mut data: Vec<F> = vec![];
    let mut rng = rand::thread_rng();
    for _i in 0..num {
        data.push(F::rand(&mut rng));
    }
    data
}

fn setup(n: usize, max_degree: usize) -> (KZGPreparedData<F>, KZGKey<G1,G2>) {
    let data = gen_data(n);
    let prep = KZGPreparedData::from_iter(data);
    let crs = KZG::setup(max_degree, &mut rand::thread_rng()).unwrap();

    (prep, crs)
}

fn bench_single_proof(c: &mut Criterion) {
  let (data, crs) = setup(DATA_SIZE, MAX_CRS);
  let commit = KZG::commit(&crs, &data).unwrap();
  let mut rng = rand::thread_rng();
  c.bench_function("single proof", |b| b.iter(|| KZG::prove(&crs, &commit, rng.gen_range(0..DATA_SIZE), &data))); 
}

fn bench_multi_proof(c: &mut Criterion) {
  let (data, crs) = setup(DATA_SIZE, MAX_CRS);
  let commit = KZG::commit(&crs, &data).unwrap();
  c.bench_function("multi proof", |b| b.iter(|| KZG::prove_all(&crs, &commit, &data)));
}

criterion_group!(proofs, bench_single_proof, bench_multi_proof);
criterion_main!(proofs);