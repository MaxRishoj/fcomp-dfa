#ifndef CNTALLOC_H_
#define CNTALLOC_H_

#include <cstdlib>
#include <cstdio>
#include <new>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common.h"

/* Wrapper around malloc that tracks peak memory usage. */

// Update a memory count.
struct MemCount {
    i64 current = 0;
    i64 peak    = 0;
};

extern MemCount mem_count;

inline void reset_mem_count() {
    mem_count.current = 0;
    mem_count.peak    = 0;
}

extern std::unordered_map<void *, size_t> alloc_sizes;

template <typename T>
struct CntAlloc {
    using value_type = T;

  public:
    CntAlloc() = default;

    template <typename U>
    constexpr CntAlloc(const CntAlloc<U>&) noexcept {}

    T* allocate(size_t n) {
        T *p = (T *)malloc(n * sizeof(T));
        if (p == NULL) {
            fprintf(stderr, "[ALLOC] Failed to allocate %lu bytes!\n", n * sizeof(T));
            throw std::bad_alloc();
        }

        mem_count.current += sizeof(T) * n;
        if (mem_count.current > mem_count.peak) mem_count.peak = mem_count.current;

        alloc_sizes[(void *)p] = n;
        return p;
    }

    void deallocate(T *p, size_t n) noexcept {
        alloc_sizes.erase(p);
        free(p);
        mem_count.current -= sizeof(T) * n;
    }
    void deallocate(T *p) noexcept {
        size_t n = alloc_sizes[(void *)p];
        deallocate(p, n);
    }
};

template <typename T>
using Vector = std::vector<T, CntAlloc<T>>;

template <typename K, typename V>
using UnorderedMap = std::unordered_map<K, V, std::hash<K>, std::equal_to<K>, CntAlloc<std::pair<const K, V>>>;

template <typename K>
using UnorderedSet = std::unordered_set<K, std::hash<K>, std::equal_to<K>, CntAlloc<K>>;

#endif // CNTALLOC_H_
