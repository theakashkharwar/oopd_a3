// C1: PPM snapshot writer, blocking or asynchronous with at most Q pending writes.
#pragma once
#include "runner.hpp"
#include <filesystem>
class SnapWriter {
  const Params& p; std::deque<std::future<std::string>> pend; int done = 0;
  std::vector<std::string> failures;
  static std::string write(std::vector<uint8_t> st, int N, std::string path, int latMs) {
    try {
      if (latMs > 0) std::this_thread::sleep_for(std::chrono::milliseconds(latMs));  // imitate slow storage
      static const uint8_t col[5][3] = {{60, 140, 220}, {220, 40, 40}, {60, 190, 90}, {240, 220, 60}, {30, 30, 30}};
      std::vector<uint8_t> img(st.size() * 3);
      for (size_t i = 0; i < st.size(); i++) for (int k = 0; k < 3; k++) img[3 * i + k] = col[st[i]][k];
      std::ofstream f(path, std::ios::binary); if (!f) return "cannot open " + path;
      f << "P6\n" << N << " " << N << "\n255\n"; f.write((const char*)img.data(), img.size());
      return f ? "" : "write failed " + path;
    } catch (std::exception& e) { return std::string(e.what()); }
  }
  void reap(std::future<std::string>& f) { std::string r = f.get(); done++; if (!r.empty()) failures.push_back(r); }
public:
  explicit SnapWriter(const Params& q) : p(q) { std::filesystem::create_directories(p.snapDir); }
  void snapshot(const Grid& g, int day) {
    char nm[64]; snprintf(nm, sizeof nm, "/day_%05d.ppm", day); std::string path = p.snapDir + nm;
    if (p.snapMode == "blocking") { std::string r = write(g.s, g.N, path, p.ioLatency); done++; if (!r.empty()) failures.push_back(r); return; }
    // Copy is required: the main loop reuses/overwrites this grid next day while the task is still reading it.
    std::vector<uint8_t> copy = g.s;
    while ((int)pend.size() >= std::max(1, p.Q)) { reap(pend.front()); pend.pop_front(); }  // wait for oldest
    pend.push_back(std::async(std::launch::async, write, std::move(copy), g.N, path, p.ioLatency));
  }
  void finish() {
    while (!pend.empty()) { reap(pend.front()); pend.pop_front(); }
    fprintf(stderr, "snapshots written=%d failed=%zu\n", done, failures.size());
    for (auto& f : failures) fprintf(stderr, "  FAILED: %s\n", f.c_str());
  }
};