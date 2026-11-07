// Part E: epidemiology experiments (C++ only, SVG graphs).  ./epi e1|e2|e3|e4|e5|all [--N 500] [--seeds 10] [--maxdays 1000] [--workers W] [--out results_epi] [--quick]
#include "runner.hpp"
#include "plot.hpp"
struct Run { vector<Counts> v; };               // counts for day 0..end
Run simulate(Params p) {                        // sequential run until nobody is infected (or maxdays)
  Grid cur(p.N), next(p.N); initGrid(cur, p); Controller ctl(p); Counts t = countGrid(cur); Run r; r.v.push_back(t); double P[9];
  for (int day = 1; day <= p.days && t.n[I] > 0; day++) { fillP(P, ctl.beta()); t = updateRegion(p, cur, next, 0, p.N, 0, p.N, 1, day, P); swap(cur, next); ctl.update(t); r.v.push_back(t); }
  return r;
}
int WORKERS = 8, SEEDS = 10; string OUTD = "results_epi";
vector<vector<Run>> many(const vector<Params>& ps) {   // all runs of all settings in parallel, seeds 1..SEEDS
  Pool pool(WORKERS); vector<vector<future<Run>>> f(ps.size());
  for (size_t i = 0; i < ps.size(); i++) for (int s = 1; s <= SEEDS; s++) { Params p = ps[i]; p.seed = s; f[i].push_back(pool.submit([p] { return simulate(p); })); }
  vector<vector<Run>> out(ps.size()); for (size_t i = 0; i < ps.size(); i++) for (auto& x : f[i]) out[i].push_back(x.get()); return out;
}
struct Met { double peak = 0, peakDay = 0, attack = 0, deaths = 0, len = 0, major = 0; };
Met met(const Run& r, int N) {
  double NN = (double)N * N; Met m;
  for (size_t d = 0; d < r.v.size(); d++) { double f = r.v[d].n[I] / NN; if (f > m.peak) { m.peak = f; m.peakDay = d; } }
  auto& e = r.v.back(); m.attack = (e.n[I] + e.n[R] + e.n[D]) / NN; m.deaths = e.n[D] / NN; m.len = r.v.size() - 1; m.major = m.attack > 0.01; return m;
}
pair<double, double> ms(const vector<double>& x) { double m = 0, s = 0; for (double v : x) m += v; m /= x.size(); for (double v : x) s += (v - m) * (v - m); return {m, x.size() > 1 ? sqrt(s / (x.size() - 1)) : 0}; }
vector<double> col(const vector<Run>& rs, int N, double Met::*f) { vector<double> v; for (auto& r : rs) v.push_back(met(r, N).*f); return v; }
struct Curve { vector<double> x, mean, lo, hi; };
Curve curve(const vector<Run>& rs, int N, int state) {   // mean +- 1 sd over seeds, finished runs padded with their last day
  size_t T = 0; for (auto& r : rs) T = max(T, r.v.size()); Curve c; double NN = (double)N * N;
  for (size_t t = 0; t < T; t++) { vector<double> v; for (auto& r : rs) v.push_back(r.v[min(t, r.v.size() - 1)].n[state] / NN); auto q = ms(v); c.x.push_back(t); c.mean.push_back(q.first); c.lo.push_back(max(0.0, q.first - q.second)); c.hi.push_back(q.first + q.second); }
  return c;
}
Params base(int N, int days) { Params p; p.N = N; p.days = days; p.lockdown = false; return p; }
void csv(const string& n, const string& h, const vector<vector<double>>& rows) { ofstream f(OUTD + "/" + n + ".csv"); f << h << "\n"; for (auto& r : rows) { for (size_t i = 0; i < r.size(); i++) f << (i ? "," : "") << fmt(r[i]); f << "\n"; } }
Series withErr(string n, vector<double> x, vector<double> m, vector<double> sd) { Series s = ser(n, x, m); for (size_t i = 0; i < m.size(); i++) { s.lo.push_back(max(0.0, m[i] - sd[i])); s.hi.push_back(m[i] + sd[i]); } return s; }
int N = 500, MAXD = 1000;
vector<Run> baseRuns;  // default parameters, no lockdown (reused by E1 and E5)
void getBase() { if (baseRuns.empty()) baseRuns = many({base(N, MAXD)})[0]; }

void e1() {
  Params lk = base(N, MAXD); lk.lockdown = true; getBase(); auto lr = many({lk})[0]; vector<Panel> ps;
  const char* nm[] = {"S", "I", "R", "D"}; int st[] = {S, I, R, D}; vector<vector<double>> rows; int k = 0;
  for (auto* rs : {&baseRuns, &lr}) { Panel p{k ? "Default parameters + lockdown (lambda=0.3)" : "Default parameters, no lockdown", "day", "fraction of population", false, false, {}, {}, {}, 0};
    for (int j = 0; j < 4; j++) { Curve c = curve(*rs, N, st[j]); Series s = ser(nm[j], c.x, c.mean, 'c'); if (j == 1) { s.lo = c.lo; s.hi = c.hi; } p.s.push_back(s);
      if (!k && j == 0) for (size_t t = 0; t < c.x.size(); t += 5) { rows.push_back({c.x[t]}); } }
    ps.push_back(p); k++; }
  Curve s0 = curve(baseRuns, N, S), i0 = curve(baseRuns, N, I), r0 = curve(baseRuns, N, R), d0 = curve(baseRuns, N, D); rows.clear();
  for (size_t t = 0; t < s0.x.size(); t++) rows.push_back({s0.x[t], s0.mean[t], i0.mean[t], i0.lo[t], i0.hi[t], r0.mean[t], d0.mean[t]});
  csv("e1_default_curves", "day,S,I,I_lo,I_hi,R,D", rows); saveFig(OUTD + "/e1.svg", "E1 mean S/I/R/D over " + fmt(SEEDS) + " seeds (band: I +-1 sd), N=" + fmt(N), ps, 2);
  vector<double> pk; for (auto& r : baseRuns) pk.push_back(met(r, N).peak); auto q = ms(pk); auto pd = ms(col(baseRuns, N, &Met::peakDay)); auto at = ms(col(baseRuns, N, &Met::attack));
  printf("E1 baseline: peak I=%.3f+-%.3f at day %.0f+-%.0f, attack=%.3f+-%.3f\n", q.first, q.second, pd.first, pd.second, at.first, at.second);
}
void e2() {   // transmission threshold
  vector<double> bs = {0.01, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12}; vector<Params> ps; for (double b : bs) { Params p = base(N, MAXD); p.beta = b; ps.push_back(p); }
  auto R_ = many(ps); vector<double> am, asd, pm, len; vector<vector<double>> rows; double bc = -1;
  for (size_t i = 0; i < bs.size(); i++) { auto a = ms(col(R_[i], N, &Met::attack)); double pmj = ms(col(R_[i], N, &Met::major)).first, l = ms(col(R_[i], N, &Met::len)).first;
    am.push_back(a.first); asd.push_back(a.second); pm.push_back(pmj); len.push_back(l); rows.push_back({bs[i], a.first, a.second, pmj, l, 1 - pow(1 - bs[i], 7)}); if (bc < 0 && pmj >= 0.5) bc = bs[i];
    printf("E2 beta=%.3f attack=%.3f+-%.3f P(major)=%.1f mean_len=%.0f days\n", bs[i], a.first, a.second, pmj, l); }
  csv("e2_beta", "beta,attack_mean,attack_sd,p_major,mean_days,pair_prob_over_d", rows); printf("E2 estimated threshold: smallest beta with P(major)>=0.5 is %.3f (per-neighbour-pair infection prob over d=7 days: %.3f)\n", bc, bc > 0 ? 1 - pow(1 - bc, 7) : 0.0);
  saveFig(OUTD + "/e2.svg", "E2 effect of transmission probability (no lockdown, 10+ seeds)", {Panel{"Final fraction ever infected", "beta", "attack rate", false, false, {withErr("mean +-1sd", bs, am, asd)}, {}, {}, 0},
    Panel{"Fraction of runs with a major outbreak (>1% infected)", "beta", "P(major)", false, false, {ser("P", bs, pm)}, {}, {}, 0}, Panel{"Epidemic duration", "beta", "days", false, false, {ser("days", bs, len)}, {}, {}, 0}}, 3);
}
void e3() {   // vaccination / herd immunity
  vector<double> vs = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.9}; vector<Params> ps; for (double v : vs) { Params p = base(N, MAXD); p.v = v; ps.push_back(p); }
  auto R_ = many(ps); vector<double> am, asd, pm, pk; vector<vector<double>> rows; double vc = -1;
  for (size_t i = 0; i < vs.size(); i++) { auto a = ms(col(R_[i], N, &Met::attack)); double pmj = ms(col(R_[i], N, &Met::major)).first, p_ = ms(col(R_[i], N, &Met::peak)).first;
    am.push_back(a.first); asd.push_back(a.second); pm.push_back(pmj); pk.push_back(p_); rows.push_back({vs[i], a.first, a.second, pmj, p_}); if (vc < 0 && pmj <= 0.1) vc = vs[i];
    printf("E3 v=%.2f attack=%.3f+-%.3f P(major)=%.1f peak=%.3f\n", vs[i], a.first, a.second, pmj, p_); }
  double R0 = 8 * 0.08 * 7, mf = 1 - 1 / R0; csv("e3_vaccination", "v,attack_mean,attack_sd,p_major,peak_mean", rows);
  printf("E3 grid herd-immunity estimate (smallest v with P(major)<=0.1): %.2f   mean-field 1-1/R0 with R0=8*beta*d=%.2f: %.2f\n", vc, R0, mf);
  saveFig(OUTD + "/e3.svg", "E3 effect of vaccination (dotted: mean-field 1-1/R0)", {Panel{"Final fraction ever infected", "vaccinated fraction v", "attack rate", false, false, {withErr("mean +-1sd", vs, am, asd)}, {}, {mf}, 0},
    Panel{"Fraction of runs with a major outbreak", "vaccinated fraction v", "P(major)", false, false, {ser("P", vs, pm)}, {}, {mf}, 0}, Panel{"Peak infected fraction", "vaccinated fraction v", "peak I", false, false, {ser("peak", vs, pk)}, {}, {mf}, 0}}, 3);
}
void e4() {   // masks and lockdown
  struct Sc { string name; double m; bool lk; }; vector<Sc> sc = {{"none", 1, false}, {"mask .75", .75, false}, {"mask .5", .5, false}, {"lockdown", 1, true}, {"mask+lock", .5, true}};
  vector<Params> ps; for (auto& s : sc) { Params p = base(N, MAXD); p.m = s.m; p.lockdown = s.lk; ps.push_back(p); }
  auto R_ = many(ps); Panel pc{"Mean infected fraction over time", "day", "I / N^2", false, false, {}, {}, {}, 0}; Series hp, dp; hp.mode = dp.mode = 'b'; vector<string> L; vector<vector<double>> rows;
  for (size_t i = 0; i < sc.size(); i++) { Curve c = curve(R_[i], N, I); pc.s.push_back(ser(sc[i].name, c.x, c.mean, 'c')); auto pk = ms(col(R_[i], N, &Met::peak)), pd = ms(col(R_[i], N, &Met::peakDay)), at = ms(col(R_[i], N, &Met::attack));
    L.push_back(sc[i].name); hp.y.push_back(pk.first); hp.lo.push_back(max(0.0, pk.first - pk.second)); hp.hi.push_back(pk.first + pk.second); dp.y.push_back(pd.first); dp.lo.push_back(max(0.0, pd.first - pd.second)); dp.hi.push_back(pd.first + pd.second);
    rows.push_back({(double)i, pk.first, pk.second, pd.first, pd.second, at.first, at.second}); printf("E4 %-10s peak=%.4f+-%.4f day=%.0f+-%.0f attack=%.3f+-%.3f\n", sc[i].name.c_str(), pk.first, pk.second, pd.first, pd.second, at.first, at.second); }
  csv("e4_interventions", "scenario_index(none;mask.75;mask.5;lockdown;mask.5+lockdown),peak_mean,peak_sd,peak_day_mean,peak_day_sd,attack_mean,attack_sd", rows);
  saveFig(OUTD + "/e4.svg", "E4 masks and lockdown (bars: mean +-1sd)", {pc, Panel{"Peak height", "scenario", "peak infected fraction", false, false, {hp}, L}, Panel{"Peak timing", "scenario", "day of peak", false, false, {dp}, L}}, 3);
}
void e5() {   // grid vs standard SIR
  getBase(); double b = 8 * 0.08, g = 1.0 / 7, mu = 0.01, i0 = 10.0 / ((double)N * N), dt = 0.05; double y[4] = {1 - i0, i0, 0, 0};   // s,i,r,d
  auto f = [&](const double* u, double* d) { d[0] = -b * u[0] * u[1]; d[1] = b * u[0] * u[1] - g * u[1]; d[2] = (1 - mu) * g * u[1]; d[3] = mu * g * u[1]; };
  Curve ci = curve(baseRuns, N, I), cr = curve(baseRuns, N, R), cd = curve(baseRuns, N, D); size_t T = ci.x.size(); vector<double> sx, si, sr; double sp = 0, spd = 0;
  for (size_t t = 0; t < T; t++) { sx.push_back(t); si.push_back(y[1]); sr.push_back(y[2] + y[3]); if (y[1] > sp) { sp = y[1]; spd = t; }
    for (int k = 0; k < (int)(1 / dt + .5); k++) { double k1[4], k2[4], k3[4], k4[4], u[4]; f(y, k1); for (int j = 0; j < 4; j++) u[j] = y[j] + dt / 2 * k1[j]; f(u, k2); for (int j = 0; j < 4; j++) u[j] = y[j] + dt / 2 * k2[j]; f(u, k3);
      for (int j = 0; j < 4; j++) u[j] = y[j] + dt * k3[j]; f(u, k4); for (int j = 0; j < 4; j++) y[j] += dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]); } }
  vector<double> gc; for (size_t t = 0; t < T; t++) gc.push_back(cr.mean[t] + cd.mean[t]); double gp = 0, gpd = 0; for (size_t t = 0; t < T; t++) if (ci.mean[t] > gp) { gp = ci.mean[t]; gpd = t; }
  printf("E5 SIR: R0=%.2f peak=%.3f at day %.0f final size=%.3f | grid mean-curve: peak=%.3f at day %.0f final size=%.3f\n", b / g, sp, spd, sr.back(), gp, gpd, gc.back());
  Series gi = ser("grid mean +-1sd", ci.x, ci.mean, 'c'); gi.lo = ci.lo; gi.hi = ci.hi; vector<vector<double>> rows; for (size_t t = 0; t < T; t++) rows.push_back({(double)t, ci.mean[t], si[t], gc[t], sr[t]});
  csv("e5_sir_vs_grid", "day,grid_I,sir_I,grid_RplusD,sir_RplusD", rows);
  saveFig(OUTD + "/e5.svg", "E5 grid model vs standard SIR (beta_SIR=8*beta, gamma=1/d)", {Panel{"Infected fraction", "day", "I / N^2", false, false, {gi, ser("SIR", sx, si, 'd')}, {}, {}, 0}, Panel{"Recovered + dead (cumulative)", "day", "(R+D) / N^2", false, false, {ser("grid", ci.x, gc, 'c'), ser("SIR", sx, sr, 'd')}, {}, {}, 0}}, 2);
}
int main(int argc, char** argv) {
  string ex = argc > 1 ? argv[1] : "all"; WORKERS = max(1u, thread::hardware_concurrency());
  for (int i = 2; i < argc; i++) { string a = argv[i]; if (a == "--quick") { N = 100; SEEDS = 3; MAXD = 400; } else if (i + 1 < argc) { string v = argv[++i];
      if (a == "--N") N = stoi(v); else if (a == "--seeds") SEEDS = stoi(v); else if (a == "--maxdays") MAXD = stoi(v); else if (a == "--workers") WORKERS = stoi(v); else if (a == "--out") OUTD = v; } }
  if (SEEDS < 10 && ex != "all") fprintf(stderr, "note: assignment asks for >= 10 seeds\n"); system(("mkdir -p " + OUTD).c_str()); double t0 = nowSec();
  if (ex == "e1" || ex == "all") e1(); if (ex == "e2" || ex == "all") e2(); if (ex == "e3" || ex == "all") e3(); if (ex == "e4" || ex == "all") e4(); if (ex == "e5" || ex == "all") e5();
  printf("done in %.1f s\n", nowSec() - t0);
}