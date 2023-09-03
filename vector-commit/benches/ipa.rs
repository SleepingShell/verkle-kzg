use rand::Rng;
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
    let mut rng = rand::thread_rng();
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
    let point_gen = IPAPointGenerator::default();
    let crs = IPAT::setup(max_degree, &point_gen).unwrap();
    let prep = IPAPreparedData::<SIZE, F>::new_incremental(data);

    (prep, crs)
}

//fn bench_commitment(c: &mut Criterion) {
//    let (data, crs) = setup(SIZE, SIZE);
//}
