#include "cmacs.h"

#ifndef BASEPACK_H
#define BASEPACK_H 1

struct basepack
{
  byte _0 : 2;
  byte _1 : 2;
  byte _2 : 2;
  byte _3 : 2;
};

static inline$ byte
bpget (struct basepack bp[], word i)
{
  switch (i & 3)
    {
    case 0:
      return bp[i >> 2]._0;
    case 1:
      return bp[i >> 2]._1;
    case 2:
      return bp[i >> 2]._2;
    case 3:
      return bp[i >> 2]._3;
    }
}

static inline$ byte
bpset (struct basepack bp[], word i, byte base)
{
  switch (i & 3)
    {
    case 0:
      return bp[i >> 2]._0 = base;
    case 1:
      return bp[i >> 2]._1 = base;
    case 2:
      return bp[i >> 2]._2 = base;
    case 3:
      return bp[i >> 2]._3 = base;
    }
}

#endif /* BASEPACK_H */
