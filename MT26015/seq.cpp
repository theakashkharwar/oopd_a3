// Part A: single-threaded simulation.
#include "sim.hpp"
int main(int argc, char** argv) {
  Params p = parseArgs(argc, argv);
  Grid cur(p.N), next(p.N);
  initGrid(cur, p);
  std::ofstream out(p.out);
  csvHeader(out);
  Controller ctl(p);
  Counts t = countGrid(cur);
  csvRow(out, 0, t, ctl.beta());
  double P[9], t0 = nowSec();
  int day = 1;
  for (; day <= p.days && (p.fixedDays || t.n[I] > 0); day++) {
    double b = ctl.beta();
    fillP(P, b);
    t = updateRegion(p, cur, next, 0, p.N, 0, p.N, day, P);
    std::swap(cur, next);
    ctl.update(t);
    csvRow(out, day, t, b);
  }
  fprintf(stderr, "days=%d sim_time=%.4f s\n", day - 1, nowSec() - t0);
}