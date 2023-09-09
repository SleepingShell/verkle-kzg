use rand::Rng;
use vector_commit::{kzg::*, VectorCommitment};

use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::UniformRand;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

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

fn setup(n: usize, max_degree: usize) -> (KZGPreparedData<F>, KZGKey<F, G1, G2>) {
    let data = gen_data(n);
    let point_gen = KZGRandomPointGenerator::default();

    let crs = KZG::setup(max_degree, &point_gen).unwrap();
    let prep = KZGPreparedData::from_points_and_domain(data, crs.domain());

    (prep, crs)
}

fn bench_setup(c: &mut Criterion) {
    let base = 32;
    let point_gen = KZGRandomPointGenerator::default();

    let mut group = c.benchmark_group("kzg_crs_setup");
    group.sample_size(10);
    for size in [base, base * 64, base * 128, base * 512].iter() {
        group.throughput(criterion::Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| KZG::setup(size, &point_gen));
        });
    }

    group.finish();
}

fn bench_data_commitment(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);

    c.bench_function("kzg_commitment", |b| b.iter(|| KZG::commit(&crs, &data)));
}

fn bench_single_proof(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);
    let commit = KZG::commit(&crs, &data).unwrap();
    let mut rng = rand::thread_rng();
    c.bench_function("kzg_single_proof", |b| {
        b.iter(|| KZG::prove(&crs, &commit, rng.gen_range(0..DATA_SIZE), &data))
    });
}

fn bench_multi_proof(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);
    let commit = KZG::commit(&crs, &data).unwrap();
    c.bench_function("kzg_all_proof", |b| {
        b.iter(|| KZG::prove_batch(&crs, &commit, (0..DATA_SIZE).collect(), &data))
    });
}

criterion_group!(
    proofs,
    bench_single_proof,
    bench_multi_proof,
    bench_data_commitment,
    bench_setup
);
criterion_main!(proofs);
