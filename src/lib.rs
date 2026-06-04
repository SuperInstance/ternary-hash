#![forbid(unsafe_code)]

//! Hashing and fingerprinting for ternary data ({-1, 0, +1}).
//!
//! Provides TernaryHash, MinHash for ternary sets, BloomFilter with ternary
//! membership (yes/no/maybe), locality-sensitive hashing, and ternary sketching.

/// A trit value: -1, 0, or +1.
pub type Trit = i8;

/// Constants for ternary values.
pub const NEG: Trit = -1;
pub const ZERO: Trit = 0;
pub const POS: Trit = 1;

// ---------------------------------------------------------------------------
// TernaryHash
// ---------------------------------------------------------------------------

/// Hasher for ternary sequences using polynomial rolling hash with ternary bases.
#[derive(Clone, Debug)]
pub struct TernaryHash {
    state: u64,
    base: u64,
    modulus: u64,
}

impl TernaryHash {
    const DEFAULT_BASE: u64 = 31;
    const DEFAULT_MODULUS: u64 = (1u64 << 61) - 1; // Mersenne prime

    /// Create a new hasher with default parameters.
    pub fn new() -> Self {
        Self {
            state: 0,
            base: Self::DEFAULT_BASE,
            modulus: Self::DEFAULT_MODULUS,
        }
    }

    /// Create a hasher with custom base and modulus.
    pub fn with_params(base: u64, modulus: u64) -> Self {
        Self {
            state: 0,
            base,
            modulus,
        }
    }

    /// Update the hash with a single trit.
    pub fn update(&mut self, trit: Trit) {
        let mapped = (trit + 1) as u64; // map {-1,0,+1} -> {0,1,2}
        self.state = (self.state * self.base + mapped + 1) % self.modulus;
    }

    /// Update the hash with a slice of trits.
    pub fn update_slice(&mut self, trits: &[Trit]) {
        for &t in trits {
            self.update(t);
        }
    }

    /// Return the current hash value.
    pub fn finish(&self) -> u64 {
        self.state
    }

    /// Reset the hasher to initial state.
    pub fn reset(&mut self) {
        self.state = 0;
    }

    /// Convenience: hash a slice in one call.
    pub fn hash_slice(trits: &[Trit]) -> u64 {
        let mut h = Self::new();
        h.update_slice(trits);
        h.finish()
    }
}

impl Default for TernaryHash {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// MinHash for ternary sets
// ---------------------------------------------------------------------------

/// MinHash sketch for estimating Jaccard similarity of ternary sets.
#[derive(Clone, Debug)]
pub struct TernaryMinHash {
    num_hashes: usize,
    seeds: Vec<u64>,
    min_values: Vec<u64>,
}

impl TernaryMinHash {
    /// Create a new MinHash with `num_hashes` hash functions.
    pub fn new(num_hashes: usize) -> Self {
        let seeds: Vec<u64> = (0..num_hashes)
            .map(|i| {
                // Simple seed generation
                let s = (i as u64).wrapping_mul(0x9e3779b97f4a7c15);
                s ^ (s >> 30)
            })
            .collect();
        let min_values = vec![u64::MAX; num_hashes];
        Self {
            num_hashes,
            seeds,
            min_values,
        }
    }

    fn hash_trit_with_seed(&self, trit: Trit, seed: u64) -> u64 {
        let v = (trit + 1) as u64;
        let mut h = seed.wrapping_add(v.wrapping_mul(0x517cc1b727220a95));
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51afd7ed558ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
        h ^= h >> 33;
        h
    }

    /// Push a trit into the sketch.
    pub fn push(&mut self, trit: Trit) {
        for i in 0..self.num_hashes {
            let h = self.hash_trit_with_seed(trit, self.seeds[i]);
            if h < self.min_values[i] {
                self.min_values[i] = h;
            }
        }
    }

    /// Push a slice of trits.
    pub fn push_slice(&mut self, trits: &[Trit]) {
        for &t in trits {
            self.push(t);
        }
    }

    /// Estimate Jaccard similarity with another MinHash.
    pub fn jaccard(&self, other: &TernaryMinHash) -> f64 {
        if self.num_hashes != other.num_hashes {
            return 0.0;
        }
        let matches = self
            .min_values
            .iter()
            .zip(other.min_values.iter())
            .filter(|(a, b)| a == b)
            .count();
        matches as f64 / self.num_hashes as f64
    }

    /// Return the sketch values.
    pub fn sketch(&self) -> &[u64] {
        &self.min_values
    }
}

// ---------------------------------------------------------------------------
// BloomFilter with ternary membership
// ---------------------------------------------------------------------------

/// Ternary membership result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TernaryMembership {
    /// Definitely in the set.
    Yes,
    /// Definitely not in the set.
    No,
    /// Possibly in the set (false positive possible).
    Maybe,
}

/// Bloom filter supporting ternary membership queries.
#[derive(Clone, Debug)]
pub struct TernaryBloomFilter {
    bits: Vec<u64>,
    num_bits: usize,
    num_hashes: usize,
    count: usize,
}

impl TernaryBloomFilter {
    /// Create a new bloom filter with `num_bits` bits and `num_hashes` hash functions.
    pub fn new(num_bits: usize, num_hashes: usize) -> Self {
        let num_words = (num_bits + 63) / 64;
        Self {
            bits: vec![0u64; num_words],
            num_bits,
            num_hashes,
            count: 0,
        }
    }

    fn hash_trit(&self, trit: Trit, hash_idx: usize) -> usize {
        let v = (trit + 1) as u64;
        let seed = hash_idx as u64;
        let mut h = v.wrapping_add(seed.wrapping_mul(0x517cc1b727220a95));
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51afd7ed558ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
        (h as usize) % self.num_bits
    }

    fn set_bit(&mut self, idx: usize) {
        let word = idx / 64;
        let bit = idx % 64;
        self.bits[word] |= 1u64 << bit;
    }

    fn get_bit(&self, idx: usize) -> bool {
        let word = idx / 64;
        let bit = idx % 64;
        (self.bits[word] >> bit) & 1 == 1
    }

    /// Insert a trit into the filter.
    pub fn insert(&mut self, trit: Trit) {
        for i in 0..self.num_hashes {
            let idx = self.hash_trit(trit, i);
            self.set_bit(idx);
        }
        self.count += 1;
    }

    /// Insert a slice of trits.
    pub fn insert_slice(&mut self, trits: &[Trit]) {
        for &t in trits {
            self.insert(t);
        }
    }

    /// Check ternary membership of a trit.
    pub fn contains(&self, trit: Trit) -> TernaryMembership {
        for i in 0..self.num_hashes {
            if !self.get_bit(self.hash_trit(trit, i)) {
                return TernaryMembership::No;
            }
        }
        // Could be Yes or Maybe — use count heuristic
        if self.count == 0 {
            TernaryMembership::No
        } else if self.count <= 3 {
            TernaryMembership::Yes
        } else {
            TernaryMembership::Maybe
        }
    }

    /// Check if definitely not present.
    pub fn is_absent(&self, trit: Trit) -> bool {
        self.contains(trit) == TernaryMembership::No
    }

    /// Number of items inserted.
    pub fn count(&self) -> usize {
        self.count
    }

    /// Estimated false positive rate.
    pub fn estimated_fpr(&self) -> f64 {
        let k = self.num_hashes as f64;
        let n = self.count as f64;
        let m = self.num_bits as f64;
        (1.0 - (1.0 - 1.0 / m).powf(k * n)).powf(k)
    }

    /// Reset the filter.
    pub fn clear(&mut self) {
        for word in &mut self.bits {
            *word = 0;
        }
        self.count = 0;
    }
}

// ---------------------------------------------------------------------------
// Locality-Sensitive Hashing (LSH) for ternary similarity
// ---------------------------------------------------------------------------

/// LSH band for ternary vectors.
#[derive(Clone, Debug)]
pub struct TernaryLSH {
    num_bands: usize,
    rows_per_band: usize,
    hash_size: usize,
}

impl TernaryLSH {
    /// Create a new LSH with `num_bands` bands and `rows_per_band` rows per band.
    pub fn new(num_bands: usize, rows_per_band: usize) -> Self {
        Self {
            num_bands,
            rows_per_band,
            hash_size: num_bands * rows_per_band,
        }
    }

    /// Hash a ternary vector into band signatures.
    pub fn hash(&self, trits: &[Trit]) -> Vec<u64> {
        let mut bands = Vec::with_capacity(self.num_bands);
        for band in 0..self.num_bands {
            let mut h: u64 = (band as u64).wrapping_mul(0x9e3779b97f4a7c15);
            let start = band * self.rows_per_band;
            let end = start + self.rows_per_band;
            for i in start..end.min(trits.len()) {
                let v = (trits[i] + 1) as u64;
                h = h.wrapping_mul(31).wrapping_add(v + 1);
            }
            bands.push(h);
        }
        bands
    }

    /// Check if two ternary vectors are candidates for being similar.
    pub fn are_candidates(&self, a: &[Trit], b: &[Trit]) -> bool {
        let ha = self.hash(a);
        let hb = self.hash(b);
        ha.iter().zip(hb.iter()).any(|(x, y)| x == y)
    }

    /// Return the number of matching bands (approximate similarity).
    pub fn matching_bands(&self, a: &[Trit], b: &[Trit]) -> usize {
        let ha = self.hash(a);
        let hb = self.hash(b);
        ha.iter().zip(hb.iter()).filter(|(x, y)| x == y).count()
    }

    /// Return the hash size (total signature length).
    pub fn hash_size(&self) -> usize {
        self.hash_size
    }

    /// Number of bands.
    pub fn num_bands(&self) -> usize {
        self.num_bands
    }
}

// ---------------------------------------------------------------------------
// Ternary Sketching
// ---------------------------------------------------------------------------

/// A compact sketch (fingerprint) of a ternary sequence.
#[derive(Clone, Debug)]
pub struct TernarySketch {
    /// Feature counts: [count_neg, count_zero, count_pos]
    counts: [usize; 3],
    /// Running hash of the sequence.
    fingerprint: u64,
    /// Length of the sketched sequence.
    length: usize,
    /// Sum of mapped values (trit + 1).
    weighted_sum: i64,
}

impl TernarySketch {
    /// Create a new empty sketch.
    pub fn new() -> Self {
        Self {
            counts: [0; 3],
            fingerprint: 0,
            length: 0,
            weighted_sum: 0,
        }
    }

    /// Add a trit to the sketch.
    pub fn add(&mut self, trit: Trit) {
        let idx = (trit + 1) as usize; // 0, 1, or 2
        self.counts[idx] += 1;
        self.length += 1;
        self.weighted_sum += trit as i64;
        // Update fingerprint
        self.fingerprint = self
            .fingerprint
            .wrapping_mul(31)
            .wrapping_add((trit as u64).wrapping_add(2));
    }

    /// Add a slice of trits.
    pub fn add_slice(&mut self, trits: &[Trit]) {
        for &t in trits {
            self.add(t);
        }
    }

    /// Return the fingerprint hash.
    pub fn fingerprint(&self) -> u64 {
        self.fingerprint
    }

    /// Return the counts [neg, zero, pos].
    pub fn counts(&self) -> &[usize; 3] {
        &self.counts
    }

    /// Return the length of the sketched sequence.
    pub fn len(&self) -> usize {
        self.length
    }

    /// Check if the sketch is empty.
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Compute the mean trit value.
    pub fn mean(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.weighted_sum as f64 / self.length as f64
        }
    }

    /// Compute a simple similarity (normalized dot product) between two sketches.
    pub fn similarity(&self, other: &TernarySketch) -> f64 {
        if self.length == 0 || other.length == 0 {
            return 0.0;
        }
        // Use count overlap as proxy
        let dot: f64 = self
            .counts
            .iter()
            .zip(other.counts.iter())
            .map(|(a, b)| (*a as f64) * (*b as f64))
            .sum();
        let norm_a: f64 = self.counts.iter().map(|c| (*c as f64).powi(2)).sum::<f64>().sqrt();
        let norm_b: f64 = other.counts.iter().map(|c| (*c as f64).powi(2)).sum::<f64>().sqrt();
        if norm_a == 0.0 || norm_b == 0.0 {
            0.0
        } else {
            dot / (norm_a * norm_b)
        }
    }

    /// Merge another sketch into this one.
    pub fn merge(&mut self, other: &TernarySketch) {
        for i in 0..3 {
            self.counts[i] += other.counts[i];
        }
        self.length += other.length;
        self.weighted_sum += other.weighted_sum;
        // Combine fingerprints
        self.fingerprint = self
            .fingerprint
            .wrapping_mul(37)
            .wrapping_add(other.fingerprint);
    }

    /// Weighted sum of trits.
    pub fn weighted_sum(&self) -> i64 {
        self.weighted_sum
    }
}

impl Default for TernarySketch {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Fingerprint utilities
// ---------------------------------------------------------------------------

/// Compute a 32-bit fingerprint of a ternary sequence.
pub fn fingerprint32(trits: &[Trit]) -> u32 {
    let mut h: u32 = 0x811c9dc5; // FNV offset basis
    for &t in trits {
        let byte = (t + 1) as u8; // 0, 1, or 2
        h ^= byte as u32;
        h = h.wrapping_mul(0x01000193); // FNV prime
    }
    h
}

/// Compute a 64-bit fingerprint of a ternary sequence.
pub fn fingerprint64(trits: &[Trit]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325; // FNV-1a offset basis
    for &t in trits {
        let byte = (t + 1) as u8;
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3); // FNV-1a prime
    }
    h
}

/// Compute similarity between two ternary sequences using Hamming distance.
pub fn hamming_similarity(a: &[Trit], b: &[Trit]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return 0.0;
    }
    let matches = a.iter().zip(b.iter()).filter(|(x, y)| x == y).count();
    matches as f64 / a.len() as f64
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- TernaryHash tests ---

    #[test]
    fn test_hash_empty() {
        let h = TernaryHash::new();
        assert_eq!(h.finish(), 0);
    }

    #[test]
    fn test_hash_deterministic() {
        let data: Vec<Trit> = vec![1, 0, -1, 1, 1, 0, -1];
        let h1 = TernaryHash::hash_slice(&data);
        let h2 = TernaryHash::hash_slice(&data);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_hash_different_inputs() {
        let a: Vec<Trit> = vec![1, 1, 1];
        let b: Vec<Trit> = vec![-1, -1, -1];
        assert_ne!(TernaryHash::hash_slice(&a), TernaryHash::hash_slice(&b));
    }

    #[test]
    fn test_hash_order_matters() {
        let a: Vec<Trit> = vec![1, 0, -1];
        let b: Vec<Trit> = vec![-1, 0, 1];
        assert_ne!(TernaryHash::hash_slice(&a), TernaryHash::hash_slice(&b));
    }

    #[test]
    fn test_hash_reset() {
        let mut h = TernaryHash::new();
        h.update_slice(&[1, 0, -1]);
        assert_ne!(h.finish(), 0);
        h.reset();
        assert_eq!(h.finish(), 0);
    }

    // --- MinHash tests ---

    #[test]
    fn test_minhash_same_set() {
        let data: Vec<Trit> = vec![1, 0, -1, 1, 0];
        let mut a = TernaryMinHash::new(64);
        let mut b = TernaryMinHash::new(64);
        a.push_slice(&data);
        b.push_slice(&data);
        assert!((a.jaccard(&b) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_minhash_disjoint() {
        let mut a = TernaryMinHash::new(64);
        let mut b = TernaryMinHash::new(64);
        // Only +1 in a, only -1 in b
        for _ in 0..100 {
            a.push(1);
            b.push(-1);
        }
        let sim = a.jaccard(&b);
        assert!(sim < 0.5, "disjoint sets should have low similarity, got {}", sim);
    }

    #[test]
    fn test_minhash_partial_overlap() {
        let mut a = TernaryMinHash::new(128);
        let mut b = TernaryMinHash::new(128);
        // Both have some overlap (1 and 0), but different distributions
        for _ in 0..50 {
            a.push(1);
            a.push(0);
            b.push(0);
            b.push(-1);
        }
        let sim = a.jaccard(&b);
        assert!(sim > 0.0 && sim < 1.0);
    }

    #[test]
    fn test_minhash_sketch_length() {
        let mh = TernaryMinHash::new(32);
        assert_eq!(mh.sketch().len(), 32);
    }

    // --- BloomFilter tests ---

    #[test]
    fn test_bloom_insert_and_contains() {
        let mut bf = TernaryBloomFilter::new(1024, 3);
        bf.insert(1);
        assert_ne!(bf.contains(1), TernaryMembership::No);
        assert_eq!(bf.contains(-1), TernaryMembership::No);
    }

    #[test]
    fn test_bloom_no_false_negatives() {
        let mut bf = TernaryBloomFilter::new(512, 5);
        let trits: Vec<Trit> = vec![1, 0, -1, 1, -1, 0, 0, 1];
        for &t in &trits {
            bf.insert(t);
        }
        // All inserted values should not be No
        for &t in &trits {
            assert_ne!(bf.contains(t), TernaryMembership::No);
        }
    }

    #[test]
    fn test_bloom_count() {
        let mut bf = TernaryBloomFilter::new(256, 3);
        assert_eq!(bf.count(), 0);
        bf.insert(1);
        bf.insert(0);
        bf.insert(-1);
        assert_eq!(bf.count(), 3);
    }

    #[test]
    fn test_bloom_clear() {
        let mut bf = TernaryBloomFilter::new(256, 3);
        bf.insert(1);
        bf.insert(0);
        bf.clear();
        assert_eq!(bf.count(), 0);
        // After clear, all should be absent (no bits set)
        assert_eq!(bf.contains(1), TernaryMembership::No);
    }

    #[test]
    fn test_bloom_fpr_small() {
        let mut bf = TernaryBloomFilter::new(4096, 3);
        bf.insert(1);
        bf.insert(0);
        bf.insert(-1);
        assert!(bf.estimated_fpr() < 0.1);
    }

    // --- LSH tests ---

    #[test]
    fn test_lsh_identical_are_candidates() {
        let lsh = TernaryLSH::new(8, 4);
        let v: Vec<Trit> = vec![1, 0, -1, 1, 0, -1, 1, 0];
        assert!(lsh.are_candidates(&v, &v));
    }

    #[test]
    fn test_lhash_matching_bands_identical() {
        let lsh = TernaryLSH::new(4, 2);
        let v: Vec<Trit> = vec![1, 0, -1, 1, 0, -1, 1, 0];
        assert_eq!(lsh.matching_bands(&v, &v), 4);
    }

    #[test]
    fn test_lsh_different_vectors() {
        let lsh = TernaryLSH::new(8, 4);
        let a: Vec<Trit> = vec![1, 1, 1, 1, 1, 1, 1, 1];
        let b: Vec<Trit> = vec![-1, -1, -1, -1, -1, -1, -1, -1];
        // Very different vectors — may or may not be candidates but bands should differ
        let matches = lsh.matching_bands(&a, &b);
        assert!(matches < 8);
    }

    #[test]
    fn test_lsh_hash_size() {
        let lsh = TernaryLSH::new(5, 3);
        assert_eq!(lsh.hash_size(), 15);
    }

    // --- TernarySketch tests ---

    #[test]
    fn test_sketch_empty() {
        let s = TernarySketch::new();
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }

    #[test]
    fn test_sketch_counts() {
        let mut s = TernarySketch::new();
        s.add(1);
        s.add(1);
        s.add(0);
        s.add(-1);
        assert_eq!(s.counts(), &[1, 1, 2]); // [neg, zero, pos]
        assert_eq!(s.len(), 4);
    }

    #[test]
    fn test_sketch_mean() {
        let mut s = TernarySketch::new();
        s.add(1);
        s.add(-1);
        s.add(0);
        assert!((s.mean() - 0.0).abs() < 1e-9);
    }

    #[test]
    fn test_sketch_similarity_same() {
        let mut a = TernarySketch::new();
        let mut b = TernarySketch::new();
        for t in [1, 0, -1, 1] {
            a.add(t);
            b.add(t);
        }
        assert!((a.similarity(&b) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_sketch_merge() {
        let mut a = TernarySketch::new();
        let mut b = TernarySketch::new();
        a.add(1);
        b.add(-1);
        a.merge(&b);
        assert_eq!(a.len(), 2);
        assert_eq!(a.counts(), &[1, 0, 1]);
    }

    // --- Fingerprint tests ---

    #[test]
    fn test_fingerprint32_deterministic() {
        let d: Vec<Trit> = vec![1, 0, -1];
        assert_eq!(fingerprint32(&d), fingerprint32(&d));
    }

    #[test]
    fn test_fingerprint64_different() {
        let a: Vec<Trit> = vec![1, 1, 1];
        let b: Vec<Trit> = vec![-1, -1, -1];
        assert_ne!(fingerprint64(&a), fingerprint64(&b));
    }

    #[test]
    fn test_hamming_similarity_identical() {
        let v: Vec<Trit> = vec![1, 0, -1, 1];
        assert!((hamming_similarity(&v, &v) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_hamming_similarity_opposite() {
        let a: Vec<Trit> = vec![1, 0, -1];
        let b: Vec<Trit> = vec![-1, 0, 1];
        // Only middle element matches
        assert!((hamming_similarity(&a, &b) - (1.0 / 3.0)).abs() < 1e-9);
    }

    #[test]
    fn test_sketch_fingerprint_uniqueness() {
        let mut s1 = TernarySketch::new();
        let mut s2 = TernarySketch::new();
        s1.add_slice(&[1, 0, -1]);
        s2.add_slice(&[-1, 0, 1]);
        assert_ne!(s1.fingerprint(), s2.fingerprint());
    }

    #[test]
    fn test_bloom_is_absent() {
        let mut bf = TernaryBloomFilter::new(256, 3);
        bf.insert(1);
        assert!(!bf.is_absent(1));
        assert!(bf.is_absent(-1));
    }

    #[test]
    fn test_sketch_weighted_sum() {
        let mut s = TernarySketch::new();
        s.add(1);
        s.add(-1);
        s.add(1);
        assert_eq!(s.weighted_sum(), 1);
    }
}
