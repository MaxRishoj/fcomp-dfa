#include "bitarray.h"

#include <cassert>

#include "util.h"

// Allocates but does not initialize the bit array.
BitArray ba_make(i32 count) {
    const u32 words = (count / 32) + (count % 32 ? 1 : 0);
    BitArray ba;
    ba.count = count;
    ba.data  = alloc_array<u32>(words);
    return ba;
}

void ba_delete(BitArray &ba) {
    free_and_reset_array(ba.data);
    ba.count = 0;
}

bool ba_get(BitArray &ba, i32 index) {
    assert(0 <= index && index < (i32)ba.count);
    const u32 index_word  = index / 32;
    const u32 index_local = index - (index_word * 32); // We avoid modulo.
    return (ba.data[index_word] >> index_local) & 1;
}

void ba_set(BitArray &ba, i32 index, bool value) {
    assert(0 <= index && index < (i32)ba.count);
    const u32 index_word  = index / 32;
    const u32 index_local = index - (index_word * 32); // We avoid modulo.
    ba.data[index_word] &= ~(1 << index_local);  // Clear the bit.
    ba.data[index_word] |= value << index_local; // Set the bit.
}

void ba_zero(BitArray &ba) {
    const u32 words = (ba.count / 32) + (ba.count % 32 ? 1 : 0);
    memzero(ba.data, words);
}
