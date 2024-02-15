#ifndef COMMON_H_
#define COMMON_H_

#include "stdint.h"
#include <cassert>

typedef float    f32;
typedef double   f64;

typedef int64_t  i64;
typedef uint64_t u64;
typedef int32_t  i32;
typedef uint32_t u32;
typedef uint8_t  u8;

typedef u32 state_t;
typedef u8  label_t; // Type of symbols in the alphabet.

struct SmallEdge {
    state_t source;
    state_t dest;
};

struct Edge {
    state_t source;
    state_t dest;
    i32 heap_weight;
};

inline u64 to_edge_key(i32 state_a, i32 state_b) {
    assert(state_a != state_b);
    if (state_a < state_b) return ((u64)state_a << 32) | (u64)state_b;
    else                   return ((u64)state_b << 32) | (u64)state_a;
}

#endif // COMMON_H_
