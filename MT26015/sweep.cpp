// C2: parameter sweep. Each run is an async task; results handled in the order they FINISH.
#include "runner.hpp"
#include <sstream>
struct Job { double beta, v; uint64_t seed; };
static std::vector<double> parseList(const std::string& s) { std::vector<double> r; std::stringstream ss(s); std::string x; while (std::getline(ss, x, ',')) r.push_back(std::stod(x)); return r; }
int main(int argc, char** argv) {
  Params base = parseArgs(argc, argv); if (base.out == "out.csv") base.out = "results.csv";
  std::vector<Job> jobs;
  for (double b : parseList(base.betas)) for (double v : parseList(base.vs)) for (int s = 1; s <= base.nseeds; s++) jobs.push_back({b, v, (uint64_t)s});
  size_t total = jobs.size(); double t0 = nowSec();
  Pool pool(std::max(1, base.workers));                 // R runs at once
  std::vector<std::future<Summary>> fut(total); std::vector<bool> handled(total, false);
  for (size_t i = 0; i < total; i++) {
    Params p = base; p.beta = jobs[i].beta; p.v = jobs[i].v; p.seed = jobs[i].seed;
    fut[i] = pool.submit([p] { return runSim(p); });
  }
  std::ofstream out(base.out); out << "beta,v,seed,peak_frac,peak_day,total_infected,deaths,days\n";
  size_t done = 0, handledN = 0; bool cancelled = false; std::vector<size_t> missing;
  while (handledN < total) {
    bool any = false;
    for (size_t i = 0; i < total; i++) {                // poll: take whichever future is ready (finish order)
      if (handled[i] || fut[i].wait_for(std::chrono::milliseconds(0)) != std::future_status::ready) continue;
      handled[i] = true; handledN++; any = true;
      try {
        Summary s = fut[i].get(); done++;
        char b[200]; snprintf(b, sizeof b, "%g,%g,%llu,%.6f,%d,%llu,%llu,%d\n", jobs[i].beta, jobs[i].v, (unsigned long long)jobs[i].seed,
                              s.peakFrac, s.peakDay, (unsigned long long)s.totalInfected, (unsigned long long)s.deaths, s.days);
        out << b; out.flush(); printf("%zu/%zu runs done (beta=%g v=%g seed=%llu)\n", done, total, jobs[i].beta, jobs[i].v, (unsigned long long)jobs[i].seed);
      } catch (const std::future_error&) { missing.push_back(i); }   // cancelled before it started
      catch (const std::exception& e) { missing.push_back(i); fprintf(stderr, "run %zu failed: %s\n", i, e.what()); }
    }
    if (base.timeLimit > 0 && !cancelled && nowSec() - t0 > base.timeLimit) {
      size_t n = pool.cancelPending(); cancelled = true; printf("time limit reached: cancelled %zu unstarted runs, waiting for running ones\n", n);
    }
    if (!any) std::this_thread::sleep_for(std::chrono::milliseconds(2));
  }
  printf("finished %zu/%zu runs in %.2f s\n", done, total, nowSec() - t0);
  if (!missing.empty()) { printf("missing runs:\n"); for (size_t i : missing) printf("  beta=%g v=%g seed=%llu\n", jobs[i].beta, jobs[i].v, (unsigned long long)jobs[i].seed); }
}