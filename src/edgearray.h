#ifndef EDGEARRAY_H_
#define EDGEARRAY_H_

#include "common.h"

struct EdgeArray {
    Edge *data;
    i64 count;
    i64 size;
};

EdgeArray earray_init(i32 size);

void earray_free(EdgeArray &a);

void earray_clear(EdgeArray &a);

void earray_push(EdgeArray &a, const Edge &e);

#endif // EDGEARRAY_H_
