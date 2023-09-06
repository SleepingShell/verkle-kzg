use rand::{thread_rng, Rng};
use sha2::Sha256;
use vector_commit::{ipa::*, VectorCommitment};

use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

const SIZE: usize = 256;

type F = <Bn254 as Pairing>::ScalarField;
type G1 = <Bn254 as Pairing>::G1;
type Hasher = DefaultFieldHasher<Sha256>;

type IPAT = IPA<256, G1, Hasher>;

fn gen_data(num: usize) -> Vec<F> {
    let mut data: Vec<F> = vec![];
    let mut rng = thread_rng();
    for _i in 0..num {
        data.push(F::rand(&mut rng));
    }
    data
}

fn setup(
    n: usize,
    max_degree: usize,
) -> (
    IPAPreparedData<SIZE, F>,
    IPAUniversalParams<SIZE, G1, Hasher>,
) {
    let data = gen_data(n);
    let mut point_gen = IPAPointGenerator::default();
    point_gen.set_max(512);
    let crs = IPAT::setup(max_degree, &point_gen).unwrap();
    let prep = IPAPreparedData::<SIZE, F>::new_incremental(data);

    (prep, crs)
}

fn bench_commitment(c: &mut Criterion) {
    let (data, crs) = setup(SIZE, SIZE);

    c.bench_function("ipa_commitment", |b| b.iter(|| IPAT::commit(&crs, &data)));
}

fn bench_prove_single(c: &mut Criterion) {
    let (data, crs) = setup(SIZE, SIZE);
    let commit = IPAT::commit(&crs, &data).unwrap();

    c.bench_function("ipa_prove_single_in_domain", |b| {
        let index = thread_rng().gen_range(0..SIZE);
        b.iter(|| IPAT::prove(&crs, &commit, index, &data))
    });

    c.bench_function("ipa_prove_single_out_domain", |b| {
        let index = thread_rng().gen_range(SIZE..SIZE * 16);
        b.iter(|| IPAT::prove(&crs, &commit, index, &data))
    });
}

fn bench_verify_single(c: &mut Criterion) {
    let (data, crs) = setup(SIZE, SIZE);
    let commit = IPAT::commit(&crs, &data).unwrap();

    c.bench_function("ipa_verify_single_in_domain", |b| {
        let index = thread_rng().gen_range(0..SIZE);
        let proof = IPAT::prove(&crs, &commit, index, &data).unwrap();
        b.iter(|| IPAT::verify(&crs, &commit, index, &proof))
    });
}

criterion_group!(
    proofs,
    bench_commitment,
    bench_prove_single,
    bench_verify_single
);
criterion_main!(proofs);
