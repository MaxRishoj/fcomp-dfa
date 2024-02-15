#ifndef UTIL_H_
#define UTIL_H_

#include "common.h"
#include "cntalloc.h"

#include <cassert>
#include <string>
#include <cstring>

// Alias to avoid messing this call up.
template <typename T>
inline void memzero(T *array, u32 len) {
    memset(array, 0, sizeof(T) * len);
}

template <typename T>
inline void fill_array(T *array, u32 len, const T &value) {
    memset(array, value, sizeof(T) * len);
}

template <typename T>
void prefix_sum(T *array, u32 len) {
    T sum = 0;
    for (u32 i = 0; i < len; i++) {
        T pre = array[i];
        array[i] = sum;
        sum += pre;
    }
}

// Frees an array and sets the pointer to NULL.
template <typename T>
void free_and_reset_array(T *&array) { // Reference to pointer.
    if (array == NULL) return;
    CntAlloc<T>().deallocate(array);
    array = NULL;
}

// Allocates an array.
template <typename T>
T *alloc_array(u64 size) {
    T *ptr = CntAlloc<T>().allocate(size);
    assert(ptr != NULL);
    return ptr;
}

// Copy an array to another destination.
template <typename T>
void copy_array(T *dst, T* src, u32 count) {
    memcpy(dst, src, sizeof(T) * count);
}

// Swap two elements in an array.
template <typename T>
inline void swap_in_array(T *array, i64 i, i64 j) {
    T ti = array[i];
    array[i] = array[j];
    array[j] = ti;
}

#endif // UTIL_H_
