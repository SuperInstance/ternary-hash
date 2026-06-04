# Future Integration: ternary-hash

## Current State
Implements hashing for ternary data: `TernaryHash` (polynomial rolling hash with Mersenne prime modulus), `TernaryMinHash` for Jaccard similarity estimation, `TernaryBloomFilter` with three-state membership (yes/no/maybe), locality-sensitive hashing for approximate nearest neighbor, and ternary sketching.

## Integration Opportunities

### With ternary-cell / construct-core
Room identity verification. Each room's state vector is hashed via `TernaryHash::hash_slice()` to produce a fingerprint. If two rooms have identical fingerprints, they're in the same state. `TernaryBloomFilter` tracks which room states have been seen before — a room entering a novel state triggers a surprise signal. The "maybe" membership state (from the ternary bloom filter) handles uncertainty gracefully.

### With ternary-protocol
Message deduplication and integrity. Every protocol message carries a `TernaryHash` of its payload. Receivers verify the hash before processing. `TernaryMinHash` enables fast similarity checks between messages — if two messages have similar MinHash sketches, they carry similar information and can be merged or prioritized.

### With ternary-clustering
`TernaryMinHash` provides O(1) approximate similarity computation between room states. Instead of computing full Hamming distance between all room pairs (O(n²d)), compute MinHash sketches (O(nd)) and compare sketches (O(n²k) where k << d). This scales clustering to large fleets.

## Potential in Mature Systems
In PLATO, `TernaryHash` is the universal content-addressing scheme. Every construct state, every skill specification, every message payload has a hash. The `TernaryBloomFilter` runs at Layer 0 for fast novelty detection — a new state passes the bloom filter, and if it's not "yes" (seen before), it's escalated to Layer 1 for full analysis. Locality-sensitive hashing enables nearest-neighbor search across the fleet without centralized indexing.

## Cross-Pollination Ideas
**Music × Hashing:** Hash ternary chord progressions to identify similar passages across compositions. `TernaryMinHash` with locality-sensitive parameters groups chord sequences by voice-leading similarity. This enables music retrieval by structural similarity rather than exact match. Connects to `ternary-music`.

**Cryptography × Hashing:** The polynomial rolling hash with Mersenne prime is a simple cryptographic hash. Extend to a full ternary hash function (ternary-Merkle-Damgård?) for tamper-proof room state logging. Each room's state hash chains to the previous tick, creating an immutable audit trail.

## Dependencies for Next Steps
- Integration with `ternary-protocol` for message fingerprinting
- `TernaryBloomFilter` serialization for persistence across cell restarts
- Benchmark hash collision rates for typical room state vectors
