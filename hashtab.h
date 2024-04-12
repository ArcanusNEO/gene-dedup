#include "cmacs.h"

#ifndef HASHTAB_H
#define HASHTAB_H 1

struct dense$ rnode
{
#include "rnode.def"
};

struct dense$ hnode
{
#include "rnode.def"
  struct hnode *next;
};

#endif // HASHTAB_H
