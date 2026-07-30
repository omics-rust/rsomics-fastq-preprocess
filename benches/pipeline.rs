use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_fastq_preprocess::transform::trim::trim_record;
use rsomics_fastq_preprocess::{FilterConfig, FixedTrim, PolyTailConfig};
use rsomics_seqio::OwnedRecord;

fn representative_record() -> OwnedRecord {
    let mut seq = b"ACGT".repeat(35);
    seq.extend_from_slice(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA");
    let qual = vec![b'I'; seq.len()];
    OwnedRecord {
        id: b"representative".to_vec(),
        seq,
        qual: Some(qual),
    }
}

fn bench_pipeline(c: &mut Criterion) {
    let source = representative_record();
    let filter = FilterConfig::default();
    c.bench_function("trim_filter_100k_records", |b| {
        b.iter(|| {
            for _ in 0..100_000 {
                let mut record = source.clone();
                trim_record(
                    &mut record,
                    FixedTrim::default(),
                    Some(PolyTailConfig::default()),
                    None,
                );
                let _ = filter
                    .check(
                        &record.seq,
                        record
                            .qual
                            .as_deref()
                            .expect("benchmark FASTQ record has quality"),
                    )
                    .expect("benchmark quality encoding is valid");
            }
        });
    });
}

criterion_group!(benches, bench_pipeline);
criterion_main!(benches);
