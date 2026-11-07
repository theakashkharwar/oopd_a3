// Part D benchmarks + SVG graphs, C++ only. Usage: ./bench d1|d2|d3|d4|d5|hw|all [--reps 5] [--quick] [--phys P] [--out results]
// Run from the folder that contains ./par and ./sweep (POSIX: Linux / macOS / WSL).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <thread>
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include "plot.hpp"
using namespace std;
#ifdef __APPLE__
const double DIV = 1048576.0;
#else
const double DIV = 1024.0;
#endif
int REPS = 5, PHYS, LOG; bool QK = false; string OUT = "results";
string sh(const string& c) {
  FILE* f = popen((c + " 2>/dev/null").c_str(), "r"); string s; char b[512]; size_t n;
  while (f && (n = fread(b, 1, sizeof b, f)) > 0) s.append(b, n);
  if (f) pclose(f); while (!s.empty() && (s.back() == '\n' || s.back() == ' ')) s.pop_back(); return s;
}
// ---- run a child process, capture stdout/stderr and its peak memory ----
struct Res { string out, err; double rss; };
Res run(const vector<string>& a) {
  char to[] = "/tmp/bo_XXXXXX", te[] = "/tmp/be_XXXXXX"; int fo = mkstemp(to), fe = mkstemp(te);
  pid_t pid = fork();
  if (pid == 0) { dup2(fo, 1); dup2(fe, 2); vector<char*> v; for (auto& s : a) v.push_back((char*)s.c_str()); v.push_back(nullptr); execvp(v[0], v.data()); _exit(127); }
  int st; rusage ru; wait4(pid, &st, 0, &ru);
  auto rd = [](int fd) { lseek(fd, 0, SEEK_SET); string s; char b[4096]; ssize_t n; while ((n = read(fd, b, sizeof b)) > 0) s.append(b, n); close(fd); return s; };
  Res r{rd(fo), rd(fe), ru.ru_maxrss / DIV}; unlink(to); unlink(te);
  if (!WIFEXITED(st) || WEXITSTATUS(st)) { fprintf(stderr, "FAILED: %s\n%s\n", a[0].c_str(), r.err.c_str()); exit(1); }
  return r;
}
struct M { vector<double> t; double mem = 0; double med() const { auto v = t; sort(v.begin(), v.end()); return v[v.size() / 2]; }
  double mn() const { return *min_element(t.begin(), t.end()); } double mx() const { return *max_element(t.begin(), t.end()); } };
M measure(const vector<string>& cmd, const string& key, int reps = 0) {
  M m; for (int i = 0; i < (reps ? reps : REPS); i++) {
    Res r = run(cmd); string all = r.err + r.out; size_t p = all.find(key);
    if (p == string::npos) { fprintf(stderr, "no '%s' in output\n", key.c_str()); exit(1); }
    m.t.push_back(atof(all.c_str() + p + key.size())); m.mem = max(m.mem, r.rss); }
  return m;
}
vector<string> parCmd(int N, int days, string mode, int T, int C = 0, string split = "strips", vector<string> extra = {}) {
  vector<string> c = {"./par", "--N", fmt(N), "--days", fmt(days), "--fixed", "1", "--out", "/dev/null", "--mode", mode, "--threads", fmt(T), "--chunks", fmt(C ? C : T), "--split", split};
  c.insert(c.end(), extra.begin(), extra.end()); return c;
}
double medv(vector<double> v) { sort(v.begin(), v.end()); return v.empty() ? 0 : v[v.size() / 2]; }
vector<int> pows(int top, vector<int> extra) { set<int> s(extra.begin(), extra.end()); for (int x = 1; x <= top; x *= 2) s.insert(x); return vector<int>(s.begin(), s.end()); }
struct Tab { vector<string> h; vector<vector<string>> r; };
void saveCsv(const string& name, const Tab& t) {
  ofstream f(OUT + "/" + name + ".csv"); for (size_t i = 0; i < t.h.size(); i++) f << (i ? "," : "") << t.h[i]; f << "\n";
  for (auto& r : t.r) { for (size_t i = 0; i < r.size(); i++) f << (i ? "," : "") << r[i]; f << "\n"; }
}
// ---------------- experiments ----------------
void hw() {
  string t = "Logical cores: " + fmt(LOG) + "\nPhysical cores: " + fmt(PHYS) + "\nOS: " + sh("uname -sr") + "\nCPU: " + sh("lscpu | grep 'Model name' | cut -d: -f2 | xargs") + sh("sysctl -n machdep.cpu.brand_string") +
    "\n" + sh("lscpu | grep -i cache") + "\nMemory: " + sh("free -h | grep Mem") + sh("sysctl -n hw.memsize") + "\nCompiler: " + sh("g++ --version | head -1") + "\n";
  ofstream(OUT + "/hardware.txt") << t; printf("%s", t.c_str());
}
void d1() {
  int N = QK ? 800 : 4000, days = QK ? 30 : 200; auto T = pows(2 * LOG, {PHYS, LOG}); vector<double> md, lo, hi, sp, ef, kf; Tab tb{{"T", "median", "min", "max", "speedup", "eff", "karp_flatt"}, {}};
  for (int t : T) { M m = measure(parCmd(N, days, "pool", t), "sim_time="); md.push_back(m.med()); lo.push_back(m.mn()); hi.push_back(m.mx()); printf("D1 T=%d %.3fs\n", t, m.med()); }
  vector<double> es;
  for (size_t i = 0; i < T.size(); i++) { double p = T[i], s = md[0] / md[i]; sp.push_back(s); ef.push_back(s / p); kf.push_back(p > 1 ? (1 / s - 1 / p) / (1 - 1 / p) : NAN); if (p > 1 && p <= PHYS) es.push_back(kf[i]);
    tb.r.push_back({fmt(p), fmt(md[i]), fmt(lo[i]), fmt(hi[i]), fmt(s), fmt(s / p), fmt(kf[i])}); }
  saveCsv("d1_strong", tb); double f = medv(es); vector<double> am, x = toD(T); for (double p : x) am.push_back(1 / (f + (1 - f) / p));
  Panel a{"Speedup", "threads", "speedup", true, true, {ser("measured", x, sp), ser("ideal", x, x, 'd'), ser("Amdahl f=" + num(f), x, am, 'd')}, {}, {(double)PHYS, (double)LOG}};
  Panel b{"Efficiency", "threads", "efficiency", true, false, {ser("eff", x, ef)}, {}, {(double)PHYS, (double)LOG}, 0};
  Panel c{"Karp-Flatt serial fraction", "threads", "e(p)", true, false, {ser("e(p)", vector<double>(x.begin() + 1, x.end()), vector<double>(kf.begin() + 1, kf.end()))}, {}, {(double)PHYS, (double)LOG}, 0};
  saveFig(OUT + "/d1_strong.svg", "D1 strong scaling N=" + fmt(N) + ", " + fmt(days) + " days (dotted: physical/logical cores)", {a, b, c}, 3);
}
void d2() {
  int base = QK ? 300 : 1000, days = QK ? 20 : 100; auto T = pows(LOG, {PHYS}); vector<double> md, ef, ss, N_; Tab tb{{"T", "N", "median", "eff", "scaled_speedup"}, {}};
  for (int t : T) { int N = lround(base * sqrt((double)t)); M m = measure(parCmd(N, days, "pool", t), "sim_time="); md.push_back(m.med()); N_.push_back(N); printf("D2 T=%d N=%d %.3fs\n", t, N, m.med()); }
  vector<double> fs; for (size_t i = 0; i < T.size(); i++) { ef.push_back(md[0] / md[i]); ss.push_back(T[i] * ef[i]); if (T[i] > 1) fs.push_back((T[i] - ss[i]) / (T[i] - 1)); tb.r.push_back({fmt(T[i]), fmt(N_[i]), fmt(md[i]), fmt(ef[i]), fmt(ss[i])}); }
  saveCsv("d2_weak", tb); double f = medv(fs); vector<double> x = toD(T), g; for (double p : x) g.push_back(p - f * (p - 1));
  Panel a{"Weak-scaling efficiency", "threads", "efficiency", true, false, {ser("eff", x, ef), ser("ideal", x, vector<double>(x.size(), 1), 'd')}, {}, {(double)PHYS}, 0};
  Panel b{"Scaled speedup", "threads", "scaled speedup", true, true, {ser("measured", x, ss), ser("ideal", x, x, 'd'), ser("Gustafson f=" + num(f), x, g, 'd')}, {}, {(double)PHYS}};
  saveFig(OUT + "/d2_weak.svg", "D2 weak scaling (" + fmt(base) + "x" + fmt(base) + " cells per thread)", {a, b}, 2);
}
void d3() {
  int T = PHYS, N = QK ? 600 : 2000, days = QK ? 15 : 50; Tab t1{{"C", "median", "min", "max"}, {}}, t2{{"split", "median"}, {}}, t3{{"N", "seq", "thread", "pool"}, {}};
  vector<double> cx, cm, cl, ch; double best = 1e300; int bc = T;
  for (int k : {1, 2, 4, 8, 16, 32, 64}) { M m = measure(parCmd(N, days, "pool", T, T * k), "sim_time="); cx.push_back(T * k); cm.push_back(m.med()); cl.push_back(m.mn()); ch.push_back(m.mx());
    t1.r.push_back({fmt(T * k), fmt(m.med()), fmt(m.mn()), fmt(m.mx())}); if (m.med() < best) { best = m.med(); bc = T * k; } }
  printf("D3.1 best C=%d (%.3fs)\n", bc, best);
  vector<string> names = {"strips", "blocks", "interleaved"}; Series sb; sb.mode = 'b';
  for (auto& s : names) { M m = measure(parCmd(N, days, "pool", T, bc, s), "sim_time="); sb.y.push_back(m.med()); t2.r.push_back({s, fmt(m.med())}); }
  vector<int> sizes = {50, 100, 200, 300, 500, 750, 1000}; int dd = QK ? 50 : 200; vector<double> sq, th, pl;
  for (int n : sizes) { sq.push_back(measure(parCmd(n, dd, "seq", 1), "sim_time=").med()); th.push_back(measure(parCmd(n, dd, "thread", T), "sim_time=").med()); pl.push_back(measure(parCmd(n, dd, "pool", T), "sim_time=").med());
    t3.r.push_back({fmt(n), fmt(sq.back()), fmt(th.back()), fmt(pl.back())}); }
  for (int m = 0; m < 2; m++) { auto& v = m ? pl : th; int c = -1; for (size_t i = 0; i < sizes.size(); i++) if (v[i] < sq[i]) { c = sizes[i]; break; } printf("D3.3 %s beats sequential from N=%d\n", m ? "pool" : "thread-per-day", c); }
  saveCsv("d3_chunks", t1); saveCsv("d3_split", t2); saveCsv("d3_sizes", t3);
  Series sc = ser("time", cx, cm); sc.lo = cl; sc.hi = ch;
  Panel a{"D3.1 chunk count (T=" + fmt(T) + ")", "chunks C", "time (s)", true, false, {sc}, {}, {}, 0}, b{"D3.2 grid division (C=" + fmt(bc) + ")", "", "time (s)", false, false, {sb}, names};
  Panel c{"D3.3 small grids", "N", "time (s)", true, true, {ser("seq", toD(sizes), sq), ser("thread/day", toD(sizes), th), ser("pool", toD(sizes), pl)}};
  saveFig(OUT + "/d3.svg", "D3 size of work units", {a, b, c}, 3);
}
void d4() {
  int P = PHYS, N = QK ? 80 : 300, days = QK ? 100 : 300; Tab tb{{"R", "T", "median", "runs_per_s", "peak_rss_mb"}, {}}; Series st, sr, sm; st.mode = sr.mode = sm.mode = 'b'; vector<string> L;
  for (int R = 1; R <= P; R++) if (P % R == 0) { int T = P / R;
    M m = measure({"./sweep", "--N", fmt(N), "--days", fmt(days), "--betas", "0.06,0.08,0.10,0.12", "--vs", "0,0.1", "--seeds", "4", "--workers", fmt(R), "--threads", fmt(T), "--out", "/dev/null"}, "runs in ", QK ? 2 : 0);
    L.push_back(fmt(R) + "x" + fmt(T)); st.y.push_back(m.med()); sr.y.push_back(32 / m.med()); sm.y.push_back(m.mem); tb.r.push_back({fmt(R), fmt(T), fmt(m.med()), fmt(32 / m.med()), fmt(m.mem)}); printf("D4 %dx%d %.2fs\n", R, T, m.med()); }
  saveCsv("d4_runs_vs_threads", tb);
  saveFig(OUT + "/d4.svg", "D4: 32 runs on P=" + fmt(P) + " cores (R runs x T threads)", {Panel{"Total time", "R x T", "seconds", false, false, {st}, L}, Panel{"Throughput", "R x T", "runs / second", false, false, {sr}, L}, Panel{"Peak memory", "R x T", "MB", false, false, {sm}, L}}, 3);
}
void d5() {
  int N = QK ? 300 : 1000, days = QK ? 30 : 100, T = PHYS, reps = QK ? 2 : REPS; string dir = OUT + "/_snap";
  auto go = [&](string mode, int K, int lat, int q) { M m = measure(parCmd(N, days, "pool", T, 0, "strips", {"--snap-every", fmt(K), "--snap-mode", mode, "--io-latency", fmt(lat), "--Q", fmt(q), "--snap-dir", dir}), "sim_time=", reps); sh("rm -rf " + dir); return m; };
  double base = measure(parCmd(N, days, "pool", T), "sim_time=", reps).med(); printf("D5 no snapshots %.3fs\n", base);
  Tab tb{{"K", "latency_ms", "none", "blocking", "async", "hidden_frac", "async_rss_mb"}, {}}; vector<int> lats = {0, 10, 50, 200}; vector<Panel> ps; vector<double> lx = toD(lats); Panel hid{"Write time hidden by async", "io-latency (ms)", "% hidden", false, false, {}, {}, {}, 0};
  for (int K : {5, 20}) { vector<double> bl, as; for (int lat : lats) { M b = go("blocking", K, lat, 4), a = go("async", K, lat, 4); double h = b.med() > base * 1.001 ? (b.med() - a.med()) / (b.med() - base) : NAN;
      bl.push_back(b.med()); as.push_back(a.med()); tb.r.push_back({fmt(K), fmt(lat), fmt(base), fmt(b.med()), fmt(a.med()), fmt(h), fmt(a.mem)}); printf("D5 K=%d lat=%d block=%.2f async=%.2f hidden=%.0f%%\n", K, lat, b.med(), a.med(), h * 100);
      if (hid.s.size() < 2 && lat == lats[0]) hid.s.push_back(ser("K=" + fmt(K), {}, {})); hid.s.back().x.push_back(lat); hid.s.back().y.push_back(h * 100); }
    ps.push_back(Panel{"Snapshot every K=" + fmt(K) + " days", "io-latency (ms)", "time (s)", false, false, {ser("none", lx, vector<double>(lats.size(), base)), ser("blocking", lx, bl), ser("async Q=4", lx, as)}, {}, {}, 0}); }
  saveCsv("d5_latency", tb); ps.push_back(hid); vector<double> qx, qt, qm; Tab tq{{"Q", "time", "rss_mb"}, {}};
  for (int q : {1, 2, 4, 8, 16}) { M m = go("async", 5, 50, q); qx.push_back(q); qt.push_back(m.med()); qm.push_back(m.mem); tq.r.push_back({fmt(q), fmt(m.med()), fmt(m.mem)}); printf("D5 Q=%d %.2fs %.0fMB\n", q, m.med(), m.mem); }
  saveCsv("d5_Q", tq); ps.push_back(Panel{"Effect of Q (K=5, 50 ms)", "Q", "time (s)", true, false, {ser("time", qx, qt)}, {}, {}, 0}); ps.push_back(Panel{"Effect of Q on memory", "Q", "peak RSS (MB)", true, false, {ser("rss", qx, qm)}, {}, {}, 0});
  saveFig(OUT + "/d5.svg", "D5 hiding I/O time", ps, 3);
}
int main(int argc, char** argv) {
  string exp = argc > 1 ? argv[1] : "all"; int ph = 0;
  for (int i = 2; i < argc; i++) { string a = argv[i]; if (a == "--quick") QK = true; else if (a == "--reps" && i + 1 < argc) REPS = atoi(argv[++i]); else if (a == "--phys" && i + 1 < argc) ph = atoi(argv[++i]); else if (a == "--out" && i + 1 < argc) OUT = argv[++i]; }
  LOG = max(1u, thread::hardware_concurrency()); string s = sh("lscpu -p=core,socket | grep -v '^#' | sort -u | wc -l"); if (s.empty() || atoi(s.c_str()) <= 0) s = sh("sysctl -n hw.physicalcpu");
  PHYS = ph ? ph : (atoi(s.c_str()) > 0 ? atoi(s.c_str()) : max(1, LOG / 2)); sh("mkdir -p " + OUT); hw();
  if (exp == "d1" || exp == "all") d1(); if (exp == "d2" || exp == "all") d2(); if (exp == "d3" || exp == "all") d3(); if (exp == "d4" || exp == "all") d4(); if (exp == "d5" || exp == "all") d5();
}