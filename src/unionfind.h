#ifndef __UNIONFIND_H_
#define __UNIONFIND_H_

#include "common.h"

struct UnionFind {
    i32 count;
    i32 *ids;
    i32 *szs;
};

UnionFind uf_init(i32 count);
void uf_delete(UnionFind &uf);

void uf_reset(UnionFind &uf);

i32 uf_find(UnionFind &uf, i32 a);
void uf_union(UnionFind &uf, i32 a, i32 b);

#endif // __UNIONFIND_H_
