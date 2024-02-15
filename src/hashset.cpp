#include "hashset.h"

#include "common.h"
#include "hash.h"
#include "util.h"

const u64 KEY_EMPTY = 0;

const f32 max_load = 0.7; // Fill ratio before resizing.

inline u32 hash_key_u32(u64 key, u64 a1, u64 a2, u64 b) {
    return (u64)((a1 + key) * (a2 + (key >> 32)) + b) >> 32;
}
inline u64 hash_key(u64 key, u64* seeds) {
    const u64 h1 = (u64)hash_key_u32(key, seeds[0], seeds[1], seeds[2]);
    const u64 h2 = (u64)hash_key_u32(key, seeds[3], seeds[4], seeds[5]);
    return (h1 << 32) | h2;
}

void put_without_resize(u64 *data, u64 capacity, u64 *seeds, u64 key) {
    u64 i = hash_key(key, seeds) % capacity;
    while (data[i] != KEY_EMPTY) {
        i++;
        i -= capacity * (i >= capacity);
    }
    data[i] = key;
}

HashSet hset_init(u64 size) {
    HashSet set;
    set.count    = 0;
    set.capacity = size;
    set.data     = alloc_array<u64>(set.capacity);
    fill_array(set.data, set.capacity, KEY_EMPTY);
    for (u32 i = 0; i < 6; i++) set.seeds[i] = random_u64();

    return set;
}

void hset_free(HashSet &set) {
    free_and_reset_array(set.data);
}

void hset_insert(HashSet &set, u64 key) {
    assert(key != KEY_EMPTY);

    put_without_resize(set.data, set.capacity, set.seeds, key);
    set.count++;

    f64 load = (f64)set.count / (f64)set.capacity;
    if (load >= max_load) {
        u64 new_capacity = set.capacity * 2L;
        // printf("[HashSet] Doubling capacity: %ld -> %ld\n", set.capacity, new_capacity);

        u64 *new_data = alloc_array<u64>(new_capacity);
        fill_array(new_data, new_capacity, KEY_EMPTY);
        for (u32 i = 0; i < 6; i++) set.seeds[i] = random_u64();

        for (u64 i = 0; i < set.capacity; i++) {
            if (set.data[i] == KEY_EMPTY) continue;
            put_without_resize(new_data, new_capacity, set.seeds, set.data[i]);
        }

        free_and_reset_array(set.data);

        set.data = new_data;
        set.capacity = new_capacity;
    }
}
void hset_insert(HashSet &set, u32 u, u32 v) {
    return hset_insert(set, to_edge_key(u, v));
}

bool hset_contains(HashSet &set, u64 key) {
    assert(key != KEY_EMPTY);

    u64 i = hash_key(key, set.seeds) % set.capacity;
    while (set.data[i] != KEY_EMPTY) {
        if (set.data[i] == key) return true;
        i++;
        i -= set.capacity * (i >= set.capacity);
    }

    return false;
}
bool hset_contains(HashSet &set, u32 u, u32 v) {
    return hset_contains(set, to_edge_key(u, v));
}
