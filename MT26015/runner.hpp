// Shared by par.cpp and sweep.cpp: chunking, thread pool, one-run driver.
#pragma once
#include "sim.hpp"
#include <thread>
#include <future>
#include <functional>
#include <deque>
#include <mutex>
#include <condition_variable>
#include <memory>

struct Chunk { int r0, r1, c0, c1, rstep; };
inline std::vector<Chunk> makeChunks(const Params& p, int C) {
  int N = p.N; std::vector<Chunk> v;
  if (p.split == "strips") { C = std::min(C, N); for (int k = 0; k < C; k++) v.push_back({k * N / C, (k + 1) * N / C, 0, N, 1}); }
  else if (p.split == "interleaved") { C = std::min(C, N); for (int k = 0; k < C; k++) v.push_back({k, N, 0, N, C}); }
  else if (p.split == "blocks") {
    int gr = std::min(N, (int)std::ceil(std::sqrt((double)C))), gc = std::min(N, (C + gr - 1) / gr);
    for (int a = 0; a < gr; a++) for (int b = 0; b < gc; b++)
      v.push_back({a * N / gr, (a + 1) * N / gr, b * N / gc, (b + 1) * N / gc, 1});
  } else { fprintf(stderr, "unknown split\n"); exit(1); }
  return v;
}

// Thread pool created once. Uses mutex + condition variable for the task queue (justified in report).
class Pool {
  std::vector<std::thread> w; std::deque<std::function<void()>> q;
  std::mutex m; std::condition_variable cv; bool stop = false;
public:
  explicit Pool(int T) {
    for (int i = 0; i < T; i++) w.emplace_back([this] {
      for (;;) {
        std::function<void()> t;
        { std::unique_lock<std::mutex> lk(m); cv.wait(lk, [this] { return stop || !q.empty(); });
          if (stop && q.empty()) return; t = std::move(q.front()); q.pop_front(); }
        t();
      }});
  }
  template <class F> auto submit(F f) -> std::future<decltype(f())> {
    using R = decltype(f());
    auto t = std::make_shared<std::packaged_task<R()>>(std::move(f)); auto fu = t->get_future();
    { std::lock_guard<std::mutex> lk(m); q.push_back([t] { (*t)(); }); } cv.notify_one(); return fu;
  }
  size_t cancelPending() {  // drop tasks not yet started; their futures become broken_promise
    std::lock_guard<std::mutex> lk(m); size_t n = q.size(); q.clear(); return n;
  }
  ~Pool() { { std::lock_guard<std::mutex> lk(m); stop = true; } cv.notify_all(); for (auto& t : w) t.join(); }
};

struct Summary { double peakFrac = 0; int peakDay = 0; uint64_t totalInfected = 0, deaths = 0; int days = 0; };

// One complete run (used by the sweep). T>1 uses a private pool with T threads and T strips.
inline Summary runSim(const Params& p) {
  Grid cur(p.N), next(p.N); initGrid(cur, p);
  Controller ctl(p); Counts t = countGrid(cur); Summary s; double NN = (double)p.N * p.N;
  s.peakFrac = t.n[I] / NN;
  std::unique_ptr<Pool> pool; Params q = p; q.split = "strips";
  auto chunks = makeChunks(q, p.threads); if (p.threads > 1) pool.reset(new Pool(p.threads));
  double P[9]; int day = 1;
  for (; day <= p.days && t.n[I] > 0; day++) {
    fillP(P, ctl.beta()); t = Counts();
    if (pool) {
      std::vector<std::future<Counts>> fs;
      for (auto& c : chunks) fs.push_back(pool->submit([&, c, day] { return updateRegion(p, cur, next, c.r0, c.r1, c.c0, c.c1, 1, day, P); }));
      for (auto& f : fs) t += f.get();
    } else t = updateRegion(p, cur, next, 0, p.N, 0, p.N, 1, day, P);
    std::swap(cur, next); ctl.update(t);
    if (t.n[I] / NN > s.peakFrac) { s.peakFrac = t.n[I] / NN; s.peakDay = day; }
  }
  s.days = day - 1; s.totalInfected = t.n[I] + t.n[R] + t.n[D]; s.deaths = t.n[D];
  return s;
}