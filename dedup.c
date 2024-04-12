#include <unistd.h>
#include <omp.h>
#include "dedup.h"

#ifdef __USE_MISC
#define fread(...) fread_unlocked (__VA_ARGS__)
#define fwrite(...) fwrite_unlocked (__VA_ARGS__)
#endif

static struct basepack *bpool;
#define bpget(i) bpget (bpool, i)
#define bpset(i, base) bpset (bpool, i, base)

static struct rnode *rpool;

static struct hnode *hpool;
#define nullnr 0xffffffff
static omp_lock_t *hlock;

static long readlen = 256;
static long basenr = nullnr;
static long thread_nr = 2;
static long hashtab_len = 314606869;
static long recnr;

static long rolling_pow;
static long rolling_mi;

static inline$ byte
rpget (struct rnode *re, long i)
{
  long j = (i + re->offset) % readlen;
  if (re->dup & 2)
    j = readlen - j - 1;
  byte ret = bpget (re->index * readlen + j);
  if (re->dup & 1)
    ret = ~ret;
  return ret;
}

static void
fasta_parser (void)
{
  if (peekc () == '>')
    seekc ('\n');
  long i = 0;
  int c;
  while ((c = getc ()) != EOF)
    if (base_p[c])
      bpset (i++, bbmap[c]);
}

static u64
rhash (struct rnode *re)
{
  u64 ret = 0;
  const auto mod = hashtab_len;
  for (long j = 0; j < readlen; ++j)
    ret = (ret << 2 | rpget (re, readlen - j - 1)) % mod;
  return ret;
}

static bool
collision_p (u32 index, struct rnode *re)
{
  long u = index * readlen;
  for (long i = 0; i < readlen; ++i)
    if (bpget (u + i) != rpget (re, i))
      return false;
  return true;
}

static long
roll (struct rnode *re, u64 hash)
{
  const u64 mod = hashtab_len;
  const u64 rpow = rolling_pow;
  const u64 rmi = rolling_mi;

  auto h = hash;
  for (long i = 0; i < readlen; ++i)
    {
      /*    cycle shift    */
      /* 0 1 2 ... n - 1 n */
      /*         |         */
      /*         v         */
      /* n 0 1 2 ... n - 1 */
      auto cur = rpget (re, i);
      h = (h + readlen - cur + rpow * cur) % mod * rmi % mod;
      re->offset = (i + 1) % readlen;

      omp_set_lock (hlock + h);
      auto u = hpool + h;
      bool found = false;
      if (u->index == nullnr)
        u = null;
      while (u)
        {
          if (found = !collision_p (u->index, re))
            {
              re->index = u->index;
              break;
            }
          u = u->next;
        }
      omp_unset_lock (hlock + h);
      if (found)
        return true;
    }
  return false;
}

static void
dedup (void)
{
  const auto rem = basenr % readlen;
  const auto bound = basenr - rem;
#pragma omp parallel for
  for (long i = 0; i < bound; i += readlen)
    {
      u32 index = i / readlen;
      rpool[index].index = index;
      rpool[index].offset = rpool[index].dup = 0;
      u64 hsh[4];

      u16 j;
      for (j = 0; j < 4; ++j)
        {
          rpool[index].dup = j;
          hsh[j] = rhash (rpool + index);
          if (roll (rpool + index, hsh[j]))
            break;
        }
      if (j < 4)
        continue;

      rpool[index].index = index;
      rpool[index].offset = rpool[index].dup = 0;
      omp_set_lock (hlock + hsh[0]);
      struct hnode *hnd;
      if (hpool[hsh[0]].index == nullnr)
        hnd = hpool + hsh[0];
      else
        hnd = malloc (sizeof (hnd[0]));
      memcpy (hnd, rpool + index, sizeof (rpool[index]));
      hnd->next = hpool[hsh[0]].next;
      omp_unset_lock (hlock + hsh[0]);
    }
  if (rem)
    {
      u32 index = recnr - 1;
      rpool[index].index = index;
      rpool[index].offset = rpool[index].dup = 0;
      /* memset bpool */
    }
}

static void
cleanup (void)
{
  free (bpool);
  free (rpool);
  for (long i = 0; i < hashtab_len; ++i)
    {
      struct hnode *u = hpool[i].next;
      struct hnode *v;
      while (u)
        {
          v = u->next;
          free (u);
          u = v;
        }
      omp_destroy_lock (hlock + i);
    }
  free (hpool);
  free (hlock);
}

signed
main (signed argc, char *argv[])
{
  /* parse args */
  long opt;
  while ((opt = getopt (argc, argv, "j:l:n:s:")) != -1)
    {
      switch (opt)
        {
        case 'j':;
          thread_nr = atol (optarg);
          break;
        case 'l':;
          readlen = atol (optarg);
          break;
        case 'n':;
          basenr = atol (optarg);
          break;
        case 's':;
          hashtab_len = atol (optarg);
          break;
        default:;
          exit (1);
        }
    }
  /* check opt */
  if (thread_nr <= 0)
    exit (1);
  if (readlen <= 192 || readlen >= 16384 || readlen % 4 != 0)
    exit (1);
  if (basenr <= 0 || basenr == nullnr)
    exit (1);
  if (!isprime (hashtab_len))
    exit (1);

  /* setenv */
  recnr = basenr / readlen + !!(basenr % readlen);
  rolling_pow = qpow (4lu, readlen, hashtab_len);
  rolling_mi = qpow (4lu, hashtab_len - 2, hashtab_len);
  bpool = malloc (sizeof (bpool[0]) * (basenr / 4 + readlen / 4 + 4));
  rpool = malloc (sizeof (rpool[0]) * (recnr + 4));
  hpool = malloc (sizeof (hpool[0]) * hashtab_len + 4);
  hlock = malloc (sizeof (hlock[0]) * hashtab_len + 4);
  for (long i = 0; i < hashtab_len; ++i)
    {
      hpool[i].next = null;
      hpool[i].index = nullnr;
      omp_init_lock (hlock + i);
    }
  omp_set_num_threads (thread_nr);
  atexit (cleanup);

  /* input */
  if (optind == argc)
    fasta_parser ();
  for (long i = optind; i < argc; ++i)
    {
      freopen (argv[i], "r", stdin);
      fasta_parser ();
    }

  dedup ();

  /* output */
  u64 buf;
  /* write readlen */
  buf = readlen;
  fwrite (&buf, sizeof (buf), 1, stderr);
  /* write basenr */
  buf = basenr;
  fwrite (&buf, sizeof (buf), 1, stderr);
  /* write bitmap */
  long bitmapsz = recnr / 8 + !!(recnr % 8);
  smartptr byte *bitmap = calloc (1, bitmapsz);
  for (long i = 0; i < recnr; ++i)
    {
      if (rpool[i].index != i)
        bitput$ (bitmap, i);
    }
  fwrite (bitmap, bitmapsz, 1, stderr);
  /* write raw & index string */
  for (long i = 0; i < recnr; ++i)
    {
      if (rpool[i].index == i)
        fwrite (bpool + i * readlen / 4, readlen / 4, 1, stdout);
      else
        fwrite (rpool + i, sizeof (rpool[0]), 1, stderr);
    }
}
