#ifndef __BITARRAY_H_
#define __BITARRAY_H_

#include "common.h"

struct BitArray {
    i32 count;
    u32 *data;
};

BitArray ba_make(i32 count);
void ba_delete(BitArray &ba);

bool ba_get(BitArray &ba, i32 index);
void ba_set(BitArray &ba, i32 index, bool value);

void ba_zero(BitArray &ba);

#endif // __BITARRAY_H_
