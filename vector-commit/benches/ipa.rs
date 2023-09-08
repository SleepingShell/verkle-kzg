use rand::{thread_rng, Rng};
use sha2::Sha256;
use vector_commit::{ipa::*, MultiProofQuery, VCPreparedData, VectorCommitment};

use ark_bn254::Bn254;
use ark_ec::pairing::Pairing;
use ark_ff::{field_hashers::DefaultFieldHasher, UniformRand};
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

const SIZE: usize = 256;

type F = <Bn254 as Pairing>::ScalarField;
type G1 = <Bn254 as Pairing>::G1;
type Hasher = DefaultFieldHasher<Sha256>;

type IPAT = IPA<256, G1, Hasher>;

fn gen_data(num: usize) -> IPAPreparedData<SIZE, F> {
    let mut data: Vec<F> = vec![];
    let mut rng = thread_rng();
    let r_f = F::rand(&mut rng);
    for i in 0..num {
        data.push(r_f + F::from(i as u64));
    }
    IPAPreparedData::<SIZE, F>::new_incremental(data)
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

    (data, crs)
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

fn bench_prove_multiproof(c: &mut Criterion) {
    let base = 32;
    let mut rng = rand::thread_rng();
    let (_, crs) = setup(SIZE, SIZE);

    let mut group = c.benchmark_group("ipa_multiproof_prove");
    group.sample_size(10);
    for size in [base, base * 64, base * 512].iter() {
        let all_data = (0..*size as usize)
            .map(|_| {
                let data = gen_data(SIZE);
                let commit = IPAT::commit(&crs, &data).unwrap();
                let challenge = rng.gen_range(0..SIZE);
                let eval = *data.get(challenge).unwrap();
                (data, commit, challenge, eval)
            })
            .collect::<Vec<_>>();
        let queries = all_data
            .iter()
            .map(|(d, c, z, y)| MultiProofQuery::new(d, c, *z, *y))
            .collect::<Vec<_>>();

        group.throughput(criterion::Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), &queries, |b, q| {
            b.iter(|| IPAT::prove_multiproof(&crs, q));
        });
    }
}

fn bench_verify_multiproof(c: &mut Criterion) {
    let base = 32;
    let mut rng = rand::thread_rng();
    let (_, crs) = setup(SIZE, SIZE);

    let mut group = c.benchmark_group("ipa_multiproof_verify");
    group.sample_size(10);
    for size in [base, base * 64, base * 512].iter() {
        let all_data = (0..*size as usize)
            .map(|_| {
                let data = gen_data(SIZE);
                let commit = IPAT::commit(&crs, &data).unwrap();
                let challenge = rng.gen_range(0..SIZE);
                let eval = *data.get(challenge).unwrap();
                (data, commit, challenge, eval)
            })
            .collect::<Vec<_>>();
        let queries = all_data
            .iter()
            .map(|(d, c, z, y)| MultiProofQuery::new(d, c, *z, *y))
            .collect::<Vec<_>>();
        let proof = IPAT::prove_multiproof(&crs, &queries).unwrap();

        group.throughput(criterion::Throughput::Elements(*size as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(size),
            &(&queries, &proof),
            |b, (q, p)| {
                b.iter(|| IPAT::verify_multiproof(&crs, *q, *p));
            },
        );
    }
}

criterion_group!(
    proofs,
    bench_commitment,
    bench_prove_single,
    bench_verify_single,
    bench_prove_multiproof,
    bench_verify_multiproof
);
criterion_main!(proofs);
