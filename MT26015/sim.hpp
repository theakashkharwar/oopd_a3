// Core model shared by all parts (Phase 0). Header-only.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <algorithm>
enum : uint8_t { S = 0, I = 1, R = 2, V = 3, D = 4 };

struct Params {
  int N = 1000, days = 200, d = 7, I0 = 10, threads = 1, chunks = 1, snapK = 0, Q = 4, ioLatency = 0;
  uint64_t seed = 42;
  double beta = 0.08, mu = 0.01, v = 0, m = 1, lambda = 0.3, thOn = 0.02, thOff = 0.005;
  bool fixedDays = false, lockdown = true;   // lockdown=false: only masks apply                 // run all days even if epidemic ends
  std::string out = "out.csv", mode = "seq", split = "strips", snapDir = "snaps", snapMode = "async";
  // sweep options
  std::string betas = "0.04,0.06,0.08,0.10", vs = "0,0.2,0.4";
  int nseeds = 3, workers = 4; double timeLimit = 0;
};

inline Params parseArgs(int argc, char** argv) {
  Params p;
  for (int i = 1; i + 1 < argc; i += 2) {
    std::string k = argv[i], a = argv[i + 1];
    if (k.rfind("--", 0) != 0) { fprintf(stderr, "bad arg %s\n", k.c_str()); exit(1); }
    k = k.substr(2);
    if (k == "N") p.N = std::stoi(a);               else if (k == "days") p.days = std::stoi(a);
    else if (k == "d") p.d = std::stoi(a);          else if (k == "I0") p.I0 = std::stoi(a);
    else if (k == "threads") p.threads = std::stoi(a); else if (k == "chunks") p.chunks = std::stoi(a);
    else if (k == "snap-every") p.snapK = std::stoi(a); else if (k == "Q") p.Q = std::stoi(a);
    else if (k == "io-latency") p.ioLatency = std::stoi(a);
    else if (k == "seed") p.seed = std::stoull(a);  else if (k == "beta") p.beta = std::stod(a);
    else if (k == "mu") p.mu = std::stod(a);        else if (k == "v") p.v = std::stod(a);
    else if (k == "m") p.m = std::stod(a);          else if (k == "lambda") p.lambda = std::stod(a);
    else if (k == "theta-on") p.thOn = std::stod(a); else if (k == "theta-off") p.thOff = std::stod(a);
    else if (k == "lockdown") p.lockdown = std::stoi(a) != 0;
    else if (k == "fixed") p.fixedDays = std::stoi(a) != 0;
    else if (k == "out") p.out = a;                 else if (k == "mode") p.mode = a;
    else if (k == "snap-mode") p.snapMode = a;      else if (k == "betas") p.betas = a;
    else if (k == "vs") p.vs = a;                   else if (k == "seeds") p.nseeds = std::stoi(a);
    else if (k == "workers") p.workers = std::stoi(a); else if (k == "time-limit") p.timeLimit = std::stod(a);
    else if (k == "split") p.split = a;             else if (k == "snap-dir") p.snapDir = a;
    else { fprintf(stderr, "unknown option --%s\n", k.c_str()); exit(1); }
  }
  return p;
}

// ---- stateless RNG: depends only on (seed, day, cell, stream) ----
inline uint64_t mix64(uint64_t x) {  // splitmix64 finaliser
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
inline double rand01(uint64_t seed, uint64_t day, uint64_t cell, uint64_t stream) {
  uint64_t h = mix64(seed ^ mix64(day ^ mix64(cell ^ mix64(stream))));
  return (h >> 11) * (1.0 / 9007199254740992.0);  // [0,1)
}

struct Counts {
  uint64_t n[5] = {0, 0, 0, 0, 0};
  uint64_t chk = 0;  // order-independent grid hash (sum of per-cell hashes)
  Counts& operator+=(const Counts& o) { for (int i = 0; i < 5; i++) n[i] += o.n[i]; chk += o.chk; return *this; }
};

struct Grid {
  int N;
  std::vector<uint8_t> s, c;  // state, infected-day counter (0 unless I)
  explicit Grid(int n = 0) : N(n), s((size_t)n * n, S), c((size_t)n * n, 0) {}
};

inline void initGrid(Grid& g, const Params& p) {
  size_t T = (size_t)p.N * p.N;
  size_t nV = (size_t)std::llround(p.v * T), nI = std::min<size_t>(p.I0, T - nV);
  std::vector<uint32_t> perm(T);
  for (size_t i = 0; i < T; i++) perm[i] = (uint32_t)i;
  size_t need = nV + nI;  // partial Fisher-Yates: first `need` entries are a random sample
  for (size_t i = 0; i < need; i++) {
    size_t j = i + (size_t)(rand01(p.seed, 0xFFFFFFFFULL, i, 2) * (T - i));
    std::swap(perm[i], perm[j]);
  }
  for (size_t i = 0; i < nV; i++) g.s[perm[i]] = V;
  for (size_t i = nV; i < need; i++) { g.s[perm[i]] = I; g.c[perm[i]] = 0; }
}

inline uint64_t cellHash(size_t id, uint8_t s, uint8_t c) { return mix64(id * 1315423911ULL ^ ((uint64_t)s << 8) ^ c); }

inline Counts countGrid(const Grid& g) {  // used for the day-0 row
  Counts t;
  for (size_t i = 0; i < g.s.size(); i++) { t.n[g.s[i]]++; t.chk += cellHash(i, g.s[i], g.c[i]); }
  return t;
}

// Update rectangle rows [r0,r1) x cols [c0,c1): reads cur, writes only those cells of next.
inline Counts updateRegion(const Params& p, const Grid& cur, Grid& next, int r0, int r1, int c0, int c1, int rstep,
                           int day, const double* P) {
  const int N = cur.N;
  Counts t;
  for (int i = r0; i < r1; i += rstep)
    for (int j = c0; j < c1; j++) {
      size_t id = (size_t)i * N + j;
      uint8_t s = cur.s[id], c = cur.c[id], ns = s, nc = 0;
      if (s == S) {
        int k = 0;
        for (int a = std::max(i - 1, 0); a <= std::min(i + 1, N - 1); a++)
          for (int b = std::max(j - 1, 0); b <= std::min(j + 1, N - 1); b++)
            if ((a != i || b != j) && cur.s[(size_t)a * N + b] == I) k++;
        if (k > 0 && rand01(p.seed, day, id, 0) < P[k]) ns = I;
      } else if (s == I) {
        if (c + 1 == p.d) ns = (rand01(p.seed, day, id, 1) < p.mu) ? D : R;
        else { ns = I; nc = c + 1; }
      }
      next.s[id] = ns; next.c[id] = nc;
      t.n[ns]++; t.chk += cellHash(id, ns, nc);
    }
  return t;
}

// Intervention controller: called by the main thread between days.
struct Controller {
  const Params& p; bool locked = false;
  explicit Controller(const Params& q) : p(q) {}
  double beta() const { return p.beta * p.m * (locked ? p.lambda : 1.0); }
  void update(const Counts& t) {  // decision from day t applies from day t+1
    if (!p.lockdown) return;
    double f = (double)t.n[I] / ((double)p.N * p.N);
    if (!locked && f > p.thOn) locked = true;
    else if (locked && f < p.thOff) locked = false;
  }
};

inline void fillP(double* P, double beta) { for (int k = 0; k < 9; k++) P[k] = 1.0 - std::pow(1.0 - beta, k); }

inline void csvHeader(std::ostream& o) { o << "day,S,I,R,V,D,beta_eff,checksum\n"; }
inline void csvRow(std::ostream& o, int day, const Counts& t, double beta) {
  char b[256];
  snprintf(b, sizeof b, "%d,%llu,%llu,%llu,%llu,%llu,%.10g,%016llx\n", day, (unsigned long long)t.n[S],
           (unsigned long long)t.n[I], (unsigned long long)t.n[R], (unsigned long long)t.n[V],
           (unsigned long long)t.n[D], beta, (unsigned long long)t.chk);
  o << b;
}
inline double nowSec() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}