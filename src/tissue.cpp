#include "tissue.hpp"

#include <algorithm>
#include <random>
#include <unordered_map>

namespace {

const double CELL_AREA = M_PI / 4.0; // aire visée d'une cellule (diamètre équivalent 1)

// Grille de recherche des germes (cases de côté h)
struct SeedGrid {
  double h;
  long nx, ny;
  std::vector<std::vector<int>> b;
  SeedGrid(double Lx, double Ly, double H) : h(H) {
    nx = std::max(1L, (long)std::ceil(Lx / h));
    ny = std::max(1L, (long)std::ceil(Ly / h));
    b.resize(nx * ny);
  }
  long ix(double x) const { return std::clamp((long)std::floor(x / h), 0L, nx - 1); }
  long iy(double y) const { return std::clamp((long)std::floor(y / h), 0L, ny - 1); }
  void add(int i, const V &p) { b[iy(p.y) * nx + ix(p.x)].push_back(i); }
  // Applique f(j) aux germes des cases de la couronne k autour de la case (cx, cy)
  template <class F> void ring(long cx, long cy, long k, F f) const {
    for (long jy = cy - k; jy <= cy + k; jy++) {
      for (long jx = cx - k; jx <= cx + k; jx++) {
        if (std::max(std::labs(jx - cx), std::labs(jy - cy)) != k || jx < 0 || jy < 0 || jx >= nx || jy >= ny) {
          continue;
        }
        for (int j : b[jy * nx + jx]) {
          f(j);
        }
      }
    }
  }
};

// Garde la partie de poly plus proche de p que de q
Polygon clipHalfPlane(const Polygon &poly, const V &p, const V &q) {
  V m = (p + q) * 0.5, d = q - p;
  Polygon out;
  size_t n = poly.size();
  for (size_t i = 0; i < n; i++) {
    const V &A = poly[i];
    const V &B = poly[(i + 1) % n];
    double fa = dot(A - m, d), fb = dot(B - m, d);
    if (fa <= 0.0) {
      out.push_back(A);
    }
    if ((fa < 0.0 && fb > 0.0) || (fa > 0.0 && fb < 0.0)) {
      out.push_back(A + (B - A) * (fa / (fa - fb)));
    }
  }
  return out;
}

// Voronoï par découpe de demi-plans : chaque cellule part du domaine (convexe) et est découpée par les germes
// de couronnes de cases croissantes, jusqu'à ce qu'aucun germe plus lointain ne puisse plus la couper.
std::vector<Polygon> voronoi(const std::vector<V> &seeds, const Domain &D, double h) {
  SeedGrid g(D.Lx, D.Ly, h);
  for (size_t i = 0; i < seeds.size(); i++) {
    g.add((int)i, seeds[i]);
  }
  const Polygon dom = D.polygon();
  std::vector<Polygon> cells(seeds.size());
  for (size_t i = 0; i < seeds.size(); i++) {
    const V &p = seeds[i];
    Polygon poly = dom;
    long cx = g.ix(p.x), cy = g.iy(p.y);
    for (long k = 0;; k++) {
      g.ring(cx, cy, k, [&](int j) {
        if (j != (int)i) {
          poly = clipHalfPlane(poly, p, seeds[j]);
        }
      });
      double R = 0.0;
      for (auto &q : poly) {
        R = std::max(R, norm(q - p));
      }
      if ((double)k * h >= 2.0 * R || (k > g.nx && k > g.ny)) {
        break;
      }
    }
    cells[i] = poly;
  }
  return cells;
}

// Indice du bundle le plus proche de q au sens de Mahalanobis (exclude : bundle à ignorer ; -1 si aucun)
int nearestBundle(const V &q, const std::vector<V> &seeds, const std::vector<int> &bundles, const Mahalanobis &M,
                  int exclude = -1) {
  int best = -1;
  double dBest = 1e300;
  for (int b : bundles) {
    if (b == exclude) {
      continue;
    }
    double d = M.d2(seeds[b] - q);
    if (d < dBest) {
      dBest = d;
      best = b;
    }
  }
  return best;
}

} // namespace

Mahalanobis::Mahalanobis(const Params &p) {
  double th = p.anisotropyAngle * M_PI / 180.0, c = std::cos(th), s = std::sin(th), K = p.anisotropy;
  double sxx = K * c * c + s * s, syy = K * s * s + c * c, sxy = (K - 1.0) * c * s;
  for (double *v : {&sxx, &sxy, &syy}) { // arrondis de cos(90°)
    if (std::fabs(*v) < 1e-12 * std::max(1.0, K)) {
      *v = 0.0;
    }
  }
  double det = sxx * syy - sxy * sxy;
  ixx = syy / det;
  ixy = -sxy / det;
  iyy = sxx / det;
}

double Tissue::spacing() { return std::sqrt(CELL_AREA); }

size_t Tissue::nBundles() const { return (size_t)std::count(isBundle.begin(), isBundle.end(), 1); }

void Tissue::init(const Params &p) {
  P = p;
  const Domain D(P);
  const double sp = spacing();
  const size_t N = (size_t)std::llround(D.area() / CELL_AREA);
  const double rmin = 0.7 * sp;

  std::mt19937_64 rng((unsigned)P.seed);
  std::uniform_real_distribution<double> ux(0.0, D.Lx), uy(0.0, D.Ly);
  seeds.clear();
  if (P.lattice == LATTICE_HEX || P.lattice == LATTICE_SQUARE) {
    // Réseau hexagonal ou carré de même densité, germes déplacés d'au plus disorder x le pas
    const bool hex = P.lattice == LATTICE_HEX;
    const double a = hex ? std::sqrt(2.0 * CELL_AREA / std::sqrt(3.0)) : std::sqrt(CELL_AREA);
    const double dy = hex ? 0.5 * std::sqrt(3.0) * a : a;
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    for (long k = -2; k <= (long)std::ceil(D.Ly / dy) + 2; k++) {
      double nk = ((double)k + 0.5) * dy;
      for (long i = -2; i <= (long)std::ceil(D.Lx / a) + 2; i++) {
        long parity = hex ? ((k % 2) + 2) % 2 : 0;
        double ak = ((double)i + 0.5) * a + (double)parity * 0.5 * a;
        double r = P.disorder * a * std::sqrt(u01(rng)), th = 2.0 * M_PI * u01(rng);
        V q(ak + r * std::cos(th), nk + r * std::sin(th));
        if (D.contains(q)) {
          seeds.push_back(q);
        }
      }
    }
  } else if (P.lattice == LATTICE_UNIFORM) {
    // Tirage uniforme (processus de Poisson, sans distance minimale)
    for (size_t i = 0; i < N; i++) {
      seeds.emplace_back(ux(rng), uy(rng));
    }
  } else {
    // Tirage avec distance minimale rmin (puis complément aléatoire si besoin)
    SeedGrid g(D.Lx, D.Ly, rmin);
    size_t attempts = 0;
    while (seeds.size() < N && attempts < 200 * N) {
      attempts++;
      V q(ux(rng), uy(rng));
      bool ok = true;
      long cx = g.ix(q.x), cy = g.iy(q.y);
      g.ring(cx, cy, 0, [&](int j) { ok = ok && norm(seeds[j] - q) >= rmin; });
      g.ring(cx, cy, 1, [&](int j) { ok = ok && norm(seeds[j] - q) >= rmin; });
      if (ok) {
        g.add((int)seeds.size(), q);
        seeds.push_back(q);
      }
    }
    while (seeds.size() < N) {
      seeds.emplace_back(ux(rng), uy(rng));
    }
  }

  // Bundles : une fraction bundleFraction des germes tirés au hasard, plus les germes les plus proches des points
  // imposés par 'bundle x y'
  isBundle.assign(seeds.size(), 0);
  if (P.bundleFraction > 0.0) {
    std::vector<int> idx(seeds.size());
    for (size_t i = 0; i < idx.size(); i++) {
      idx[i] = (int)i;
    }
    std::shuffle(idx.begin(), idx.end(), rng);
    size_t nb = std::min(idx.size(), (size_t)std::llround(P.bundleFraction * (double)seeds.size()));
    for (size_t i = 0; i < nb; i++) {
      isBundle[idx[i]] = 1;
    }
  }
  for (const V &b : P.bundleAt) {
    int i = nearestSeed(b);
    if (i >= 0) {
      isBundle[i] = 1;
    }
  }
  iteration = 0;
  lastMove = 0.0;
  updateVoronoi();
}

void Tissue::restore(const Params &p, std::vector<V> s, std::vector<char> bundle, int it) {
  P = p;
  seeds = std::move(s);
  isBundle = std::move(bundle);
  isBundle.resize(seeds.size(), 0);
  iteration = it;
  lastMove = 0.0;
  updateVoronoi();
}

void Tissue::updateVoronoi() { cells = voronoi(seeds, Domain(P), spacing()); }

// Itération de point fixe de l'article : chaque germe ordinaire va au centroïde c de sa cellule, déplacé vers le
// bundle b le plus proche (Mahalanobis) d'une fraction alpha / (1 + alpha) du segment [c, b]. Les bundles suivent
// bundleMotion : lloyd (centroïde), fixed (immobiles) ou attract (attirés comme les autres germes, vers le bundle
// le plus proche autre qu'eux-mêmes). Tous les germes sont mis à jour simultanément, avec les positions des
// bundles du début de l'itération.
void Tissue::step() {
  const Mahalanobis maha(P);
  std::vector<int> bundles;
  for (size_t i = 0; i < seeds.size(); i++) {
    if (isBundle[i]) {
      bundles.push_back((int)i);
    }
  }
  std::vector<V> next = seeds;
  double move = 0.0;
  for (size_t i = 0; i < seeds.size(); i++) {
    if (isBundle[i] && P.bundleMotion == MOTION_FIXED) {
      continue;
    }
    V c = cells[i].size() >= 3 ? centroid(cells[i]) : seeds[i];
    int b = -1;
    if (P.alpha > 0.0 && (!isBundle[i] || P.bundleMotion == MOTION_ATTRACT)) {
      b = nearestBundle(c, seeds, bundles, maha, isBundle[i] ? (int)i : -1);
    }
    next[i] = (b >= 0) ? (c + seeds[b] * P.alpha) * (1.0 / (1.0 + P.alpha)) : c;
    move += norm(next[i] - seeds[i]);
  }
  lastMove = seeds.empty() ? 0.0 : move / (double)seeds.size() / spacing();
  seeds = next;
  iteration++;
  updateVoronoi();
}

void Tissue::run() {
  while (iteration < P.lloydIterations) {
    step();
  }
}

int Tissue::nearestSeed(const V &q) const {
  int best = -1;
  double dBest = 1e300;
  for (size_t i = 0; i < seeds.size(); i++) {
    double d = norm(seeds[i] - q);
    if (d < dBest) {
      dBest = d;
      best = (int)i;
    }
  }
  return best;
}

int Tissue::toggleBundle(const V &q) {
  int i = nearestSeed(q);
  if (i >= 0) {
    isBundle[i] = !isBundle[i];
  }
  return i;
}

// Sommets partagés : les sommets calculés indépendamment par chaque cellule sont fusionnés à tolérance près
Mesh Tissue::mesh() const {
  const double sp = spacing();
  const double tolV = 1e-8 * sp, hb = 1e-6 * sp;
  const Domain D(P);
  Mesh m;
  std::unordered_map<long long, std::vector<int>> hash;
  auto key = [](long long i, long long j) { return i * 1000003LL + j; };
  for (const auto &poly : cells) {
    std::vector<int> ids;
    for (const auto &q : poly) {
      long long i = (long long)std::floor(q.x / hb), j = (long long)std::floor(q.y / hb);
      int found = -1;
      for (long long di = -1; di <= 1 && found < 0; di++) {
        for (long long dj = -1; dj <= 1 && found < 0; dj++) {
          auto it = hash.find(key(i + di, j + dj));
          if (it == hash.end()) {
            continue;
          }
          for (int v : it->second) {
            if (norm(m.verts[v] - q) < tolV) {
              found = v;
              break;
            }
          }
        }
      }
      if (found < 0) {
        found = (int)m.verts.size();
        m.verts.push_back(q);
        hash[key(i, j)].push_back(found);
      }
      if (ids.empty() || ids.back() != found) {
        ids.push_back(found);
      }
    }
    while (ids.size() > 1 && ids.front() == ids.back()) {
      ids.pop_back();
    }
    bool bnd = false;
    for (int v : ids) {
      bnd = bnd || D.onBoundary(m.verts[v], 1e-6 * sp);
    }
    m.cells.push_back(ids);
    m.boundary.push_back(bnd);
  }
  return m;
}

Metrics computeMetrics(const Tissue &t, const Mesh &m) {
  struct Acc {
    size_t n{0};
    double SI{0.0}, A{0.0}, A2{0.0}, dNN{0.0};
    double mean(double v) const { return n ? v / (double)n : 0.0; }
  };
  Acc all, reg, bun;
  const double sp = Tissue::spacing();
  const Mahalanobis maha(t.P);
  const size_t nc = m.cells.size();

  std::vector<Polygon> poly(nc);
  std::vector<V> cen(nc);
  std::vector<char> ok(nc, 0);
  SeedGrid g(t.P.Lx, t.P.Ly, sp);
  for (size_t c = 0; c < nc; c++) {
    if (m.cells[c].size() < 3) {
      continue;
    }
    for (int r : m.cells[c]) {
      poly[c].push_back(m.verts[r]);
    }
    cen[c] = centroid(poly[c]);
    ok[c] = 1;
    g.add((int)c, cen[c]);
  }
  std::vector<int> bundleCells;
  for (size_t c = 0; c < nc; c++) {
    if (ok[c] && t.isBundle[c]) {
      bundleCells.push_back((int)c);
    }
  }

  Metrics r;
  double dM = 0.0, sides = 0.0;
  size_t nM = 0, nInner = 0;
  for (size_t c = 0; c < nc; c++) {
    if (!ok[c]) {
      continue;
    }
    double A = std::fabs(signedArea(poly[c])), per = perimeter(poly[c]);
    // plus proche voisin : recherche par couronnes croissantes de la grille
    double d = 1e300;
    long cx = g.ix(cen[c].x), cy = g.iy(cen[c].y);
    for (long k = 0; k <= std::max(g.nx, g.ny); k++) {
      if (d < 1e300 && (double)(k - 1) * sp > d) {
        break;
      }
      g.ring(cx, cy, k, [&](int j) {
        if (j != (int)c) {
          d = std::min(d, norm(cen[j] - cen[c]));
        }
      });
    }
    for (Acc *a : {&all, t.isBundle[c] ? &bun : &reg}) {
      a->n++;
      a->SI += 4.0 * M_PI * A / (per * per);
      a->A += A;
      a->A2 += A * A;
      a->dNN += (d < 1e300) ? d : 0.0;
    }
    if (!t.isBundle[c] && !bundleCells.empty()) {
      double dBest = 1e300;
      for (int j : bundleCells) {
        dBest = std::min(dBest, maha.d2(cen[j] - cen[c]));
      }
      dM += std::sqrt(dBest);
      nM++;
    }
    if (!m.boundary[c]) {
      size_t ns = m.cells[c].size();
      if (r.sidesInner.size() <= ns) {
        r.sidesInner.resize(ns + 1, 0);
      }
      r.sidesInner[ns]++;
      sides += (double)ns;
      nInner++;
    }
  }
  r.nCells = nc;
  r.nBundles = t.nBundles();
  r.meanSI = all.mean(all.SI);
  double mA = all.mean(all.A);
  r.cvA = (all.n > 1 && mA > 0.0) ? std::sqrt(std::max(0.0, all.mean(all.A2) - mA * mA)) / mA : 0.0;
  r.dNN = all.mean(all.dNN) / sp;
  r.dM = nM ? dM / (double)nM / sp : 0.0;
  r.areaRegular = reg.mean(reg.A);
  r.areaBundle = bun.mean(bun.A);
  r.SIRegular = reg.mean(reg.SI);
  r.SIBundle = bun.mean(bun.SI);
  r.meanSides = nInner ? sides / (double)nInner : 0.0;
  r.lastMove = t.lastMove;
  return r;
}

std::vector<std::pair<std::string, double>> Metrics::items() const {
  return {
      {"nCells", (double)nCells},     {"nBundles", (double)nBundles}, {"meanSI", meanSI},
      {"cvA", cvA},                   {"dNN", dNN},                   {"dM", dM},
      {"areaRegular", areaRegular},   {"areaBundle", areaBundle},     {"SIRegular", SIRegular},
      {"SIBundle", SIBundle},         {"meanSides", meanSides},       {"lastMove", lastMove},
  };
}
