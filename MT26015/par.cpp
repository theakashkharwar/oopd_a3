// Part B: threaded simulation. --mode seq|thread|pool  --split strips|blocks|interleaved
#include "snap.hpp"

int main(int argc, char** argv) {
  Params p = parseArgs(argc, argv);
  int C = p.chunks <= 1 ? p.threads : p.chunks;
  std::vector<Chunk> chunks = makeChunks(p, p.mode == "seq" ? 1 : C);
  Grid cur(p.N), next(p.N);
  initGrid(cur, p);
  std::ofstream out(p.out); csvHeader(out);
  Controller ctl(p);
  Counts t = countGrid(cur); csvRow(out, 0, t, ctl.beta());
  std::unique_ptr<SnapWriter> snap; if (p.snapK > 0) snap.reset(new SnapWriter(p));
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
    if (snap && day % p.snapK == 0) snap->snapshot(cur, day);
  }
  if (snap) snap->finish();   // wait for all pending writes
  fprintf(stderr, "days=%d chunks=%zu sim_time=%.6f s\n", day - 1, chunks.size(), nowSec() - t0);
}