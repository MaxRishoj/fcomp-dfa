#ifndef __HASHSET_H_
#define __HASHSET_H_

#include "common.h"

struct HashSet {
    u64 *data;
    u64 count;
    u64 capacity;
    u64 seeds[6]; // 6 = 3+3 for two 64->32 bit hash functions.
};

HashSet hset_init(u64 size);
void hset_free(HashSet &set);

void hset_insert(HashSet &set, u64 key);
void hset_insert(HashSet &set, u32 u, u32 v);

bool hset_contains(HashSet &set, u64 key);
bool hset_contains(HashSet &set, u32 u, u32 v);

#endif // __HASHSET_H_
