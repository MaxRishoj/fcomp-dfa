#include "cntalloc.h"

MemCount mem_count;

std::unordered_map<void *, size_t> alloc_sizes;
