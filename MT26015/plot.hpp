// Minimal SVG plotting (no dependencies). Shared by bench.cpp and epi.cpp.
#pragma once
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
inline string fmt(double v) { char b[40]; snprintf(b, 40, "%.6g", v); return b; }
// ---- tiny SVG plotting ----
struct Series { string name; vector<double> x, y, lo, hi; char mode = 'l'; };  // l=line+markers(+error bars), d=dashed, c=curve(+shaded band lo..hi), b=bars(+error bars)
struct Panel { string title, xl, yl; bool logx = false, logy = false; vector<Series> s; vector<string> cats; vector<double> vl; double ymin = NAN; };
inline const char* COL[] = {"#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b"};
inline string num(double v) { char b[32]; snprintf(b, 32, "%.3g", v); return b; }
inline void saveFig(const string& file, const string& title, const vector<Panel>& ps, int cols) {
  const int W = 520, H = 360, ML = 70, MR = 20, MT = 40, MB = 55; int rows = ((int)ps.size() + cols - 1) / cols;
  ostringstream o; o << "<svg xmlns='http://www.w3.org/2000/svg' width='" << W * cols << "' height='" << H * rows + 30 << "' font-family='sans-serif' font-size='11'>\n<rect width='100%' height='100%' fill='white'/>\n"
    << "<text x='" << W * cols / 2 << "' y='18' text-anchor='middle' font-size='15' font-weight='bold'>" << title << "</text>\n";
  for (size_t k = 0; k < ps.size(); k++) {
    const Panel& p = ps[k]; int ox = (k % cols) * W, oy = (k / cols) * H + 30; bool bar = !p.cats.empty();
    double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
    for (auto& s : p.s) { for (double v : s.x) if (isfinite(v)) { x0 = min(x0, v); x1 = max(x1, v); }
      for (double v : s.y) if (isfinite(v)) { y0 = min(y0, v); y1 = max(y1, v); }
      for (double v : s.lo) if (isfinite(v)) y0 = min(y0, v); for (double v : s.hi) if (isfinite(v)) y1 = max(y1, v); }
    if (bar) { x0 = -0.5; x1 = p.cats.size() - 0.5; } for (double v : p.vl) { x0 = min(x0, v); x1 = max(x1, v); }
    if (x0 > x1 || y0 > y1) continue;
    if (bar) y0 = 0; else if (!isnan(p.ymin)) y0 = p.ymin;
    if (p.logy) { y0 = pow(10, floor(log10(y0))); y1 = pow(10, ceil(log10(y1))); if (y0 == y1) y1 *= 10; }
    else { if (y1 <= y0) y1 = y0 + 1; y1 += 0.06 * (y1 - y0); }
    if (p.logx) { x0 /= 1.2; x1 *= 1.2; } else if (!bar) { double d = x1 > x0 ? 0.03 * (x1 - x0) : 1; x0 -= d; x1 += d; }
    auto tx = [&](double v) { double a = p.logx ? log10(v) : v, l = p.logx ? log10(x0) : x0, h = p.logx ? log10(x1) : x1; return ox + ML + (a - l) / (h - l) * (W - ML - MR); };
    auto ty = [&](double v) { double a = p.logy ? log10(max(v, 1e-300)) : v, l = p.logy ? log10(y0) : y0, h = p.logy ? log10(y1) : y1; return oy + H - MB - (a - l) / (h - l) * (H - MT - MB); };
    o << "<rect x='" << ox + ML << "' y='" << oy + MT << "' width='" << W - ML - MR << "' height='" << H - MT - MB << "' fill='none' stroke='#444'/>\n";
    o << "<text x='" << ox + (ML + W - MR) / 2 << "' y='" << oy + 22 << "' text-anchor='middle' font-size='13'>" << p.title << "</text>\n"
      << "<text x='" << ox + (ML + W - MR) / 2 << "' y='" << oy + H - 10 << "' text-anchor='middle'>" << p.xl << "</text>\n"
      << "<text transform='translate(" << ox + 16 << "," << oy + (MT + H - MB) / 2 << ") rotate(-90)' text-anchor='middle'>" << p.yl << "</text>\n";
    vector<double> yt; if (p.logy) for (double e = log10(y0); e <= log10(y1) + 1e-9; e++) yt.push_back(pow(10, e)); else for (int i = 0; i < 5; i++) yt.push_back(y0 + (y1 - y0) * i / 4);
    for (double v : yt) o << "<line x1='" << ox + ML << "' x2='" << ox + W - MR << "' y1='" << ty(v) << "' y2='" << ty(v) << "' stroke='#ddd'/><text x='" << ox + ML - 5 << "' y='" << ty(v) + 4 << "' text-anchor='end'>" << num(v) << "</text>\n";
    if (bar) for (size_t i = 0; i < p.cats.size(); i++) o << "<text x='" << tx(i) << "' y='" << oy + H - MB + 15 << "' text-anchor='middle'>" << p.cats[i] << "</text>\n";
    else { vector<double> xt = p.s[0].x; sort(xt.begin(), xt.end()); xt.erase(unique(xt.begin(), xt.end()), xt.end());
      if (xt.size() > 10 && !p.logx) { xt.clear(); for (int i = 0; i < 5; i++) xt.push_back(x0 + (x1 - x0) * i / 4); }
      for (double v : xt) o << "<line x1='" << tx(v) << "' x2='" << tx(v) << "' y1='" << oy + H - MB << "' y2='" << oy + H - MB + 4 << "' stroke='#444'/><text x='" << tx(v) << "' y='" << oy + H - MB + 16 << "' text-anchor='middle'>" << num(v) << "</text>\n"; }
    for (double v : p.vl) o << "<line x1='" << tx(v) << "' x2='" << tx(v) << "' y1='" << oy + MT << "' y2='" << oy + H - MB << "' stroke='gray' stroke-dasharray='2,3'/>\n";
    for (size_t si = 0; si < p.s.size(); si++) { const Series& s = p.s[si]; const char* c = COL[si % 6];
      if (s.mode == 'b') { for (size_t i = 0; i < s.y.size(); i++) if (isfinite(s.y[i])) o << "<rect x='" << tx(i - 0.3) << "' y='" << ty(s.y[i]) << "' width='" << tx(i + 0.3) - tx(i - 0.3) << "' height='" << ty(y0) - ty(s.y[i]) << "' fill='" << c << "'/><text x='" << tx(i) << "' y='" << ty(s.y[i]) - 3 << "' text-anchor='middle'>" << num(s.y[i]) << "</text>\n"; for (size_t i = 0; i < s.lo.size() && i < s.y.size(); i++) o << "<line x1='" << tx(i) << "' x2='" << tx(i) << "' y1='" << ty(s.lo[i]) << "' y2='" << ty(s.hi[i]) << "' stroke='black' stroke-width='2'/>\n"; continue; }
      if (s.mode == 'c' && s.lo.size() == s.x.size()) { string bp; for (size_t i = 0; i < s.x.size(); i++) bp += fmt(tx(s.x[i])) + "," + fmt(ty(s.hi[i])) + " "; for (size_t i = s.x.size(); i-- > 0;) bp += fmt(tx(s.x[i])) + "," + fmt(ty(s.lo[i])) + " "; o << "<polygon points='" << bp << "' fill='" << c << "' fill-opacity='0.2'/>\n"; }
      string pts; for (size_t i = 0; i < s.x.size(); i++) if (isfinite(s.y[i])) pts += fmt(tx(s.x[i])) + "," + fmt(ty(s.y[i])) + " ";
      o << "<polyline points='" << pts << "' fill='none' stroke='" << c << "' stroke-width='2'" << (s.mode == 'd' ? " stroke-dasharray='6,4'" : "") << "/>\n";
      if (s.mode == 'l') for (size_t i = 0; i < s.x.size(); i++) if (isfinite(s.y[i])) {
        o << "<circle cx='" << tx(s.x[i]) << "' cy='" << ty(s.y[i]) << "' r='3' fill='" << c << "'/>\n";
        if (i < s.lo.size()) o << "<line x1='" << tx(s.x[i]) << "' x2='" << tx(s.x[i]) << "' y1='" << ty(s.lo[i]) << "' y2='" << ty(s.hi[i]) << "' stroke='" << c << "'/>\n"; } }
    if (p.s.size() > 1 && !bar) for (size_t si = 0; si < p.s.size(); si++)
      o << "<rect x='" << ox + ML + 8 << "' y='" << oy + MT + 8 + 15 * si << "' width='10' height='10' fill='" << COL[si % 6] << "'/><text x='" << ox + ML + 22 << "' y='" << oy + MT + 17 + 15 * si << "'>" << p.s[si].name << "</text>\n";
  }
  o << "</svg>\n"; ofstream(file) << o.str();
}
inline Series ser(string n, vector<double> x, vector<double> y, char m = 'l') { Series s; s.name = n; s.x = x; s.y = y; s.mode = m; return s; }
inline vector<double> toD(const vector<int>& v) { return vector<double>(v.begin(), v.end()); }