// Part B: threaded simulation. --mode seq|thread|pool  --split strips|blocks|interleaved
#include "sim.hpp"
#include <thread>
#include <future>
#include <functional>
#include <queue>
#include <mutex>
#include <condition_variable>

struct Chunk { int r0, r1, c0, c1, rstep; };

std::vector<Chunk> makeChunks(const Params& p, int C) {
  int N = p.N; std::vector<Chunk> v;
  if (p.split == "strips") {                       // C contiguous bands of rows
    C = std::min(C, N);
    for (int k = 0; k < C; k++) v.push_back({k * N / C, (k + 1) * N / C, 0, N, 1});
  } else if (p.split == "interleaved") {           // chunk k takes rows k, k+C, k+2C, ...
    C = std::min(C, N);
    for (int k = 0; k < C; k++) v.push_back({k, N, 0, N, C});
  } else if (p.split == "blocks") {                // gr x gc grid of rectangular blocks (>= C blocks)
    int gr = std::min(N, (int)std::ceil(std::sqrt((double)C)));
    int gc = std::min(N, (C + gr - 1) / gr);
    for (int a = 0; a < gr; a++)
      for (int b = 0; b < gc; b++) v.push_back({a * N / gr, (a + 1) * N / gr, b * N / gc, (b + 1) * N / gc, 1});
  } else { fprintf(stderr, "unknown split\n"); exit(1); }
  return v;
}

// Thread pool created once. NOTE: needs a mutex + condition variable for the task queue
// (justified in the report: idle workers must sleep until work arrives).
class Pool {
  std::vector<std::thread> w; std::queue<std::packaged_task<Counts()>> q;
  std::mutex m; std::condition_variable cv; bool stop = false;
public:
  explicit Pool(int T) {
    for (int i = 0; i < T; i++) w.emplace_back([this] {
      for (;;) {
        std::packaged_task<Counts()> t;
        { std::unique_lock<std::mutex> lk(m); cv.wait(lk, [this] { return stop || !q.empty(); });
          if (stop && q.empty()) return; t = std::move(q.front()); q.pop(); }
        t();
      }});
  }
  std::future<Counts> submit(std::function<Counts()> f) {
    std::packaged_task<Counts()> t(std::move(f)); auto fu = t.get_future();
    { std::lock_guard<std::mutex> lk(m); q.push(std::move(t)); } cv.notify_one(); return fu;
  }
  ~Pool() { { std::lock_guard<std::mutex> lk(m); stop = true; } cv.notify_all(); for (auto& t : w) t.join(); }
};

int main(int argc, char** argv) {
  Params p = parseArgs(argc, argv);
  int C = p.chunks <= 1 ? p.threads : p.chunks;
  std::vector<Chunk> chunks = makeChunks(p, p.mode == "seq" ? 1 : C);
  Grid cur(p.N), next(p.N);
  initGrid(cur, p);
  std::ofstream out(p.out); csvHeader(out);
  Controller ctl(p);
  Counts t = countGrid(cur); csvRow(out, 0, t, ctl.beta());
  std::unique_ptr<Pool> pool; if (p.mode == "pool") pool.reset(new Pool(p.threads));
  double P[9], t0 = nowSec(); int day = 1;
  for (; day <= p.days && (p.fixedDays || t.n[I] > 0); day++) {
    double b = ctl.beta(); fillP(P, b);
    auto work = [&, day](const Chunk& c) { return updateRegion(p, cur, next, c.r0, c.r1, c.c0, c.c1, c.rstep, day, P); };
    t = Counts();
    if (p.mode == "thread") {                       // one new thread per chunk, every day
      std::vector<Counts> res(chunks.size()); std::vector<std::thread> th;
      for (size_t k = 0; k < chunks.size(); k++) th.emplace_back([&, k] { res[k] = work(chunks[k]); });
      for (auto& x : th) x.join();
      for (auto& r : res) t += r;
    } else if (p.mode == "pool") {                  // C tasks on T persistent threads
      std::vector<std::future<Counts>> fs;
      for (auto& c : chunks) fs.push_back(pool->submit([&, c] { return work(c); }));
      for (auto& f : fs) t += f.get();
    } else for (auto& c : chunks) t += work(c);     // seq
    std::swap(cur, next); ctl.update(t); csvRow(out, day, t, b);
  }
  fprintf(stderr, "days=%d chunks=%zu sim_time=%.4f s\n", day - 1, chunks.size(), nowSec() - t0);
}