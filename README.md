# ternary-hash

Hashing, fingerprinting, and similarity estimation for ternary data `{-1, 0, +1}` — polynomial rolling hash, MinHash, Bloom filters, locality-sensitive hashing, and compact sketches.

## Why This Exists

When working with ternary data, you often need to quickly compare, deduplicate, or estimate similarity between sequences. General-purpose hashing doesn't exploit the three-valued alphabet. This crate provides hash structures and algorithms designed specifically for the ternary domain, from exact hashing to probabilistic similarity estimation.

## Core Concepts

- **TernaryHash** — Polynomial rolling hash with a Mersenne prime modulus, designed for ternary alphabets
- **TernaryMinHash** — MinHash sketch for estimating Jaccard similarity between ternary sets
- **TernaryBloomFilter** — Bloom filter with three-valued membership: `Yes`, `No`, `Maybe`
- **TernaryLSH** — Locality-sensitive hashing via banding for approximate nearest-neighbor search
- **TernarySketch** — Compact summary statistics (counts, fingerprint, mean) mergeable across streams

## Quick Start

```toml
# Cargo.toml
[dependencies]
ternary-hash = "0.1"
```

```rust
use ternary_hash::*;

// Exact hashing
let data: Vec<Trit> = vec![POS, ZERO, NEG, POS, ZERO];
let hash = TernaryHash::hash_slice(&data);
println!("Hash: {:016x}", hash);

// Incremental hashing
let mut h = TernaryHash::new();
h.update_slice(&data[..3]);
h.update_slice(&data[3..]);
assert_eq!(h.finish(), hash);

// MinHash similarity estimation
let mut a = TernaryMinHash::new(128);
let mut b = TernaryMinHash::new(128);
a.push_slice(&[POS, POS, ZERO, NEG]);
b.push_slice(&[POS, POS, ZERO, NEG, POS]);
println!("Estimated Jaccard: {:.3}", a.jaccard(&b));

// Bloom filter with ternary membership
let mut bf = TernaryBloomFilter::new(4096, 3);
bf.insert_slice(&[POS, NEG, ZERO]);
assert!(!bf.is_absent(POS));
assert!(bf.is_absent(NEG)); // single trits may collide; use larger filters
println!("Estimated FPR: {:.4}", bf.estimated_fpr());

// Locality-sensitive hashing
let lsh = TernaryLSH::new(8, 4);
let v1: Vec<Trit> = vec![1, 0, -1, 1, 0, -1, 1, 0];
let v2: Vec<Trit> = vec![1, 0, -1, 1, 0, -1, 1, 0];
let v3: Vec<Trit> = vec![-1, -1, -1, -1, -1, -1, -1, -1];
println!("v1~v2 candidates: {}", lsh.are_candidates(&v1, &v2));
println!("Matching bands: {}/{}", lsh.matching_bands(&v1, &v3), lsh.num_bands());

// Compact sketching
let mut s1 = TernarySketch::new();
let mut s2 = TernarySketch::new();
s1.add_slice(&[POS, POS, ZERO]);
s2.add_slice(&[POS, ZERO, NEG]);
println!("Similarity: {:.3}", s1.similarity(&s2));
s1.merge(&s2);
println!("Merged length: {}", s1.len());

// Fingerprinting utilities
let fp32 = fingerprint32(&data);
let fp64 = fingerprint64(&data);
println!("FNV-1a: {:08x} / {:016x}", fp32, fp64);
println!("Hamming similarity: {:.3}", hamming_similarity(&data, &data));
```

## API Overview

| Type / Function | Description |
|---|---|
| `TernaryHash` | Rolling polynomial hash (streaming or one-shot) |
| `TernaryMinHash` | Jaccard similarity estimation via MinHash sketches |
| `TernaryBloomFilter` | Probabilistic set membership with ternary results |
| `TernaryMembership` | `Yes` / `No` / `Maybe` enum |
| `TernaryLSH` | Band-based locality-sensitive hashing |
| `TernarySketch` | Mergeable compact summary with fingerprint |
| `fingerprint32` / `fingerprint64` | FNV-1a fingerprints |
| `hamming_similarity` | Position-match ratio between two sequences |

## How It Works

**TernaryHash** maps trits `{−1, 0, +1}` → `{0, 1, 2}` and applies polynomial rolling hash with base 31 and Mersenne prime modulus (2⁶¹ − 1).

**TernaryMinHash** maintains `k` independent hash functions and tracks the minimum value seen for each. Jaccard similarity is estimated by the fraction of matching minimum values between two sketches.

**TernaryBloomFilter** uses `k` hash positions per element. A ternary membership query returns `No` (definitely absent), `Yes` (small set, likely present), or `Maybe` (present but possibly a false positive).

**TernaryLSH** divides a trit vector's hash into bands. Two vectors are similarity candidates if any band matches — tunable sensitivity via band count and rows per band.

**TernarySketch** maintains per-value counts `[neg, zero, pos]`, a rolling fingerprint, length, and weighted sum. Sketches can be merged for distributed aggregation.

## Use Cases

1. **Duplicate detection in ternary datasets** — Hash sequences to quickly identify identical or near-identical entries
2. **Approximate nearest-neighbor search** — Use LSH to prune candidates before expensive exact comparisons
3. **Streaming aggregation** — Merge sketches from distributed ternary data sources without storing raw data
4. **Real-time deduplication** — Bloom filters for ternary event streams where exact membership isn't required

## Known Limitations

- **BloomFilter membership heuristic is fragile**: `TernaryBloomFilter::contains()` uses a hardcoded count threshold (≤3 → `Yes`, else → `Maybe`) to distinguish definite membership from possible false positives. This breaks down quickly — inserting more than 3 items means every match returns `Maybe`, even if the filter is sparse. The threshold should scale with filter capacity and item count.

- **MinHash operates on individual trits, not sets**: `TernaryMinHash::push()` hashes one trit at a time, so the sketch treats `[Pos, Pos, Pos]` as three separate items rather than a weighted set. This means Jaccard estimation reflects element frequency, not set membership — calling `push_slice(&[Pos, Pos])` twice gives a different result than `push(Pos)` four times.

- **LSH band hashing ignores vector length**: `TernaryLSH::hash()` computes band signatures from the first `rows_per_band` elements of each band's range. If the input vector is shorter than `num_bands × rows_per_band`, later bands hash fewer elements, making short vectors more likely to collide spuriously.

- **TernarySketch similarity is count-based, not order-aware**: `TernarySketch::similarity()` computes cosine similarity over per-value counts `[neg, zero, pos]`, discarding sequence order entirely. `[1, 0, -1]` and `[-1, 0, 1]` produce identical sketches (same counts), so similarity is 1.0 despite being opposite sequences.

- **Rolling hash has no removal operation**: `TernaryHash` supports incremental `update()` but cannot remove old elements, so it cannot maintain a sliding-window hash without recomputing from scratch.

## Ecosystem

Part of the **SuperInstance** ternary computing crate family:

- `ternary-compression-v2` — Multi-algorithm ternary compression
- `ternary-matrix` — Compact ternary matrix operations
- `ternary-pca` — Principal component analysis on ternary values
- `ternary-ga` — Genetic algorithms with ternary genomes
- `ternary-reservoir` — Echo state networks with ternary nodes
- `ternary-evolution-advanced` — Advanced evolutionary optimization
- `ternary-geometry` — Geometric algorithms in ternary space
- `ternary-causality` — Causal inference for ternary systems
- `ternary-consensus` — Distributed consensus for ternary agents

## License

MIT

## See Also
- **ternary-codes** — related
- **ternary-compression** — related
- **ternary-diff** — related
- **ternary-steganography** — related
- **ternary-database** — related

