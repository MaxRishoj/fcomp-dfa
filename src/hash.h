#ifndef HASH_H_
#define HASH_H_

#include <climits>

#include "dfa.h"
#include "common.h"
#include "util.h"
#include "random.h"

const u32 MAX_BITS_FOR_ALPHABET = 8;

const u64 LARGE_MERSENNE_PRIME = ((u64)(1 << 31) - 1) * ((u64)(1 << 31) - 1);

inline void assert_dfa_size(const DFA *dfa) {
    assert(dfa->alphabet_size <= (1 << MAX_BITS_FOR_ALPHABET));
    assert(dfa->state_count   <= (1 << (31 - MAX_BITS_FOR_ALPHABET)));
}

inline u64 uniform_u64_for_hash() {
    u64 x;
    do { x = random_u64(); } while (x > LARGE_MERSENNE_PRIME);
    return x;
}

inline u32 pack_edge(state_t dest, u32 label) {
    return (label << (31 - MAX_BITS_FOR_ALPHABET)) | dest;
}

inline u32 hash_edge(state_t dest, u32 label, u64 a, u64 b) {
    const u32 v = pack_edge(dest, label);
    return (u64)(((u64)v*a + b) >> 32);
}

inline u32 calculate_single_min_hash(const DFA &dfa, state_t state, u64 a, u64 b) {
   u32 min_hash = UINT_MAX;
    for (i32 j = 0; j < dfa.alphabet_size; j++) {
        const state_t dest = at(dfa.txs, state, j);
        const u32 h = hash_edge(dest, j, a, b);
        if (h < min_hash) min_hash = h;
    }
    return min_hash;
}

inline void calculate_k_min_hashes(const DFA &dfa, state_t state, u64 a, u64 b, u32 min_hash_k, u32 *out_mins) {
    for (u32 k = 0; k < min_hash_k; k++) out_mins[k] = UINT_MAX;

    for (i32 j = 0; j < dfa.alphabet_size; j++) {
        const state_t dest = at(dfa.txs, state, j);
        const u32 h = hash_edge(dest, j, a, b);
        for (u32 k = 0; k < min_hash_k; k++) {
            if (h < out_mins[k]) {
                for (u32 l = min_hash_k - 1; l > k; l--) out_mins[l] = out_mins[l-1];
                out_mins[k] = h;
                break;
            }
        }
    }
    for (i32 k = 0; k < (i32)min_hash_k - 1; k++) assert(out_mins[k] <= out_mins[k+1]);
}

struct VectorHash {
    u32 size;
    u64 *seeds;
    u64 count = 0;
    u64 sum   = 0;
};

inline VectorHash vhash_init(u32 size) {
    VectorHash vhash;
    vhash.size  = size;
    vhash.seeds = alloc_array<u64>(size + 1); // b is final seed.
    for (u32 i = 0; i < vhash.size+1; i++) vhash.seeds[i] = random_u64();

    vhash.count = 0;
    vhash.sum   = 0;

    return vhash;
}

inline void vhash_clear(VectorHash &vhash) {
    vhash.count = 0;
    vhash.sum   = 0;
}

inline void vhash_free(VectorHash &vhash) {
    free_and_reset_array(vhash.seeds);
}

inline void vhash_add_pair(VectorHash &vhash, u32 x1, u32 x2) {
    const u64 a1 = vhash.seeds[vhash.count];
    const u64 a2 = vhash.seeds[vhash.count + 1];
    const u64 y1 = a1 + (u64)x2;
    const u64 y2 = a2 + (u64)x1;
    vhash.sum   += y1 * y2;
    vhash.count += 2;
}

inline void vhash_add_single(VectorHash &vhash, u32 x) {
    vhash_add_pair(vhash, x, x);
    vhash.count--;
}

inline u32 vhash_finalize(VectorHash &vhash) {
    assert(vhash.count == vhash.size);
    return (u64)(((u64)vhash.sum + vhash.seeds[vhash.size]) >> 32);
}

#endif // HASH_H_
