#ifndef RANDOM_H_
#define RANDOM_H_

#include "common.h"

// This code is based on a blog post by Pelle Evensen (public domain per comment).
// http://mostlymangling.blogspot.com/2018/07/on-mixing-functions-in-fast-splittable.html

extern u64  _RANDOM_STATE;

inline void seed_rng(u64 seed) {
    _RANDOM_STATE = seed;
}

inline u64 next_random_u64(u64 n) {
    const u64 z = 0x9fb21c651e98df25;

    n ^= ((n << 49) | (n >> 15)) ^ ((n << 24) | (n >> 40));
    n *= z;
    n ^= n >> 35;
    n *= z;
    n ^= n >> 28;

    return n;
}


inline u64 random_u64() {
    u64 out = next_random_u64(_RANDOM_STATE++);
    return out;
}

inline u32 random_u32() {
    return (u32)random_u64();
}

#endif // RANDOM_H_
