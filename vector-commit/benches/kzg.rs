use ark_poly::GeneralEvaluationDomain;
use rand::Rng;
use sha2::Sha256;
use vector_commit::{
    kzg::{kzg_point_generator::KZGRandomPointGenerator, *},
    lagrange_basis::LagrangeBasis,
    VCUniversalParams, VectorCommitment,
};

use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

type F = <Bn254 as Pairing>::ScalarField;
type G1 = <Bn254 as Pairing>::G1;
type G2 = <Bn254 as Pairing>::G2;
type Hasher = DefaultFieldHasher<Sha256>;
type D = GeneralEvaluationDomain<F>;

type KZGT = KZG<Bn254, Hasher, D>;

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

fn setup(n: usize, max_degree: usize) -> (LagrangeBasis<F, D>, KZGKey<F, G1, G2>) {
    let data = gen_data(n);
    let point_gen = KZGRandomPointGenerator::<G1>::default();

    let crs = KZGT::setup(max_degree, &point_gen).unwrap();
    let prep = LagrangeBasis::from_vec_and_domain(data, *crs.precompute().domain());

    (prep, crs)
}

fn bench_setup(c: &mut Criterion) {
    let base = 32;
    let point_gen = KZGRandomPointGenerator::<G1>::default();

    let mut group = c.benchmark_group("kzg_crs_setup");
    group.sample_size(10);
    for size in [base, base * 64, base * 128, base * 512].iter() {
        group.throughput(criterion::Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| KZGT::setup(size, &point_gen));
        });
    }

    group.finish();
}

fn bench_data_commitment(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);

    c.bench_function("kzg_commitment", |b| b.iter(|| KZGT::commit(&crs, &data)));
}

fn bench_single_proof(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);
    let commit = KZGT::commit(&crs, &data).unwrap();
    let mut rng = rand::thread_rng();
    c.bench_function("kzg_single_proof", |b| {
        let i = rng.gen_range(0..DATA_SIZE);
        b.iter(|| KZGT::prove(&crs, &commit, i, &data))
    });
}

/*
fn bench_multi_proof(c: &mut Criterion) {
    let (data, crs) = setup(DATA_SIZE, MAX_CRS);
    let commit = KZGT::commit(&crs, &data).unwrap();
    c.bench_function("kzg_all_proof", |b| {
        b.iter(|| KZGT::prove_batch(&crs, &commit, (0..DATA_SIZE).collect(), &data))
    });
}
*/

criterion_group!(
    kzg_proofs,
    bench_single_proof,
    //bench_multi_proof,
    bench_data_commitment,
    bench_setup
);
criterion_main!(kzg_proofs);
