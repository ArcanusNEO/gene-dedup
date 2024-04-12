#include "cmacs.h"

#define Ki 1024l
#define Mi 1024l * 1024l
#define Gi 1024l * 1024l * 1024l
#define Ti 1024l * 1024l * 1024l * 1024l

#define A 0
#define T 3
#define C 1
#define G 2

#ifndef DEDUP_H
#define DEDUP_H 1

static const bool base_p[256]
    = { ['A'] = true, ['T'] = true, ['C'] = true, ['G'] = true,
        ['a'] = true, ['t'] = true, ['c'] = true, ['g'] = true };

static const byte bbmap[256] = { ['A'] = A, ['T'] = T, ['C'] = C, ['G'] = G,
                                 ['a'] = A, ['t'] = T, ['c'] = C, ['g'] = G };

#endif // DEDUP_H

#include "hashtab.h"
#include "basepack.h"
