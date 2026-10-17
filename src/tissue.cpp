#include "tissue.hpp"
#include "expr.hpp"

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
std::vector<Polygon> voronoiCells(const std::vector<V> &seeds, const Domain &D, double h) {
  SeedGrid g(D.Lx, D.Ly, h);
  for (size_t i = 0; i < seeds.size(); i++) {
    g.add((int)i, seeds[i]);
  }
  // germes de chaque case triés par position : le résultat (arrondis compris) ne dépend pas de leur numérotation,
  // ce qui permet de reprendre exactement un calcul enregistré (germes extérieurs placés à la fin du fichier)
  for (auto &bucket : g.b) {
    std::sort(bucket.begin(), bucket.end(), [&](int i, int j) {
      return seeds[i].x < seeds[j].x || (seeds[i].x == seeds[j].x && seeds[i].y < seeds[j].y);
    });
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

// ================================================================================================
// Domaine : rectangle des germes et forme { f <= 0 }
// ================================================================================================

Domain::Domain(const Params &p) : Lx(p.Lx), Ly(p.Ly) {
  if (!p.shape.empty()) {
    shape = std::make_shared<const Expr>(p.shape, shapeVariables());
  }
}

double Domain::f(const V &q) const {
  if (!shape) {
    return -1.0;
  }
  const double cx = 0.5 * Lx, cy = 0.5 * Ly, dx = q.x - cx, dy = q.y - cy;
  const double v[8] = {q.x, q.y, std::sqrt(dx * dx + dy * dy), std::atan2(dy, dx), cx, cy, Lx, Ly};
  return shape->eval(v);
}

// Le contour de la cellule est échantillonné au pas h ; les points de passage du bord de la forme sont trouvés par
// dichotomie. Entre une sortie et l'entrée suivante, le bord de la forme est suivi : des points de la corde sont
// projetés sur f = 0 (méthode de Newton). Les points de passage situés sur une paroi commune à deux cellules sont
// calculés de la même façon des deux côtés : ils deviennent des sommets partagés.
Polygon Domain::clip(const Polygon &poly, double h, bool &cut) const {
  cut = false;
  if (!shape || poly.size() < 3) {
    return poly;
  }
  struct Sample {
    V p;
    double f;
    bool vertex;
  };
  std::vector<Sample> s;
  for (size_t i = 0; i < poly.size(); i++) {
    const V &A = poly[i], &B = poly[(i + 1) % poly.size()];
    int n = std::max(1, (int)std::ceil(norm(B - A) / h));
    for (int k = 0; k < n; k++) {
      V p = A + (B - A) * ((double)k / n);
      s.push_back({p, f(p), k == 0});
    }
  }
  const size_t n = s.size();
  int start = -1;
  bool anyOut = false;
  for (size_t i = 0; i < n; i++) {
    if (s[i].f <= 0.0) {
      start = start < 0 ? (int)i : start;
    } else {
      anyOut = true;
    }
  }
  if (!anyOut) {
    return poly;
  }
  cut = true;
  if (start < 0) {
    return {};
  }
  // point du bord entre a (intérieur) et b (extérieur)
  auto root = [&](V a, V b) {
    for (int it = 0; it < 60; it++) {
      V m = (a + b) * 0.5;
      (f(m) <= 0.0 ? a : b) = m;
    }
    return (a + b) * 0.5;
  };
  auto insideCell = [&](const V &q) {
    for (size_t i = 0; i < poly.size(); i++) {
      if (cross(poly[(i + 1) % poly.size()] - poly[i], q - poly[i]) < -1e-12 * h * h) {
        return false;
      }
    }
    return true;
  };
  Polygon out;
  V exitPt;
  const double eps = 1e-6 * h;
  for (size_t k = 0; k < n; k++) {
    const Sample &a = s[(start + k) % n], &b = s[(start + k + 1) % n];
    bool inA = a.f <= 0.0, inB = b.f <= 0.0;
    if (inA && a.vertex) {
      out.push_back(a.p);
    }
    if (inA && !inB) {
      exitPt = root(a.p, b.p);
      out.push_back(exitPt);
    } else if (!inA && inB) {
      V entry = root(b.p, a.p);
      // bord de la forme entre la sortie et cette entrée
      int m = (int)std::ceil(norm(entry - exitPt) / (0.5 * h)) - 1;
      for (int j = 1; j <= m; j++) {
        V p = exitPt + (entry - exitPt) * ((double)j / (m + 1));
        for (int it = 0; it < 10; it++) {
          double fp = f(p);
          V g((f(p + V(eps, 0)) - f(p - V(eps, 0))) / (2 * eps), (f(p + V(0, eps)) - f(p - V(0, eps))) / (2 * eps));
          double gg = dot(g, g);
          if (gg <= 0.0 || !std::isfinite(gg)) {
            break;
          }
          p = p - g * (fp / gg);
          if (std::fabs(fp) < 1e-12 * h) {
            break;
          }
        }
        if (insideCell(p)) {
          out.push_back(p);
        }
      }
      out.push_back(entry);
    }
  }
  return out.size() >= 3 ? out : Polygon();
}

// ================================================================================================
// Tissu
// ================================================================================================

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

size_t Tissue::nBundles() const {
  size_t n = 0;
  for (size_t i = 0; i < seeds.size(); i++) {
    n += isBundle[i] && inside(i);
  }
  return n;
}

void Tissue::init(const Params &p) {
  P = p;
  domain = Domain(P);
  const Domain &D = domain;
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

  // Bundles : une fraction bundleFraction des germes du tissu tirés au hasard, plus les germes les plus proches
  // des points imposés par 'bundle x y'
  isBundle.assign(seeds.size(), 0);
  updateVoronoi();
  if (P.bundleFraction > 0.0) {
    std::vector<int> idx;
    for (size_t i = 0; i < seeds.size(); i++) {
      if (inside(i)) {
        idx.push_back((int)i);
      }
    }
    std::shuffle(idx.begin(), idx.end(), rng);
    size_t nb = std::min(idx.size(), (size_t)std::llround(P.bundleFraction * (double)idx.size()));
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
}

void Tissue::restore(const Params &p, std::vector<V> s, std::vector<char> bundle, int it) {
  P = p;
  domain = Domain(P);
  seeds = std::move(s);
  isBundle = std::move(bundle);
  isBundle.resize(seeds.size(), 0);
  iteration = it;
  lastMove = 0.0;
  updateVoronoi();
}

void Tissue::updateVoronoi() {
  voronoi = voronoiCells(seeds, domain, spacing());
  cells.resize(seeds.size());
  areaOut.assign(seeds.size(), 0.0);
  for (size_t i = 0; i < seeds.size(); i++) {
    bool cut = false;
    cells[i] = domain.clip(voronoi[i], 0.25 * spacing(), cut);
    if (cut) {
      areaOut[i] = std::max(0.0, std::fabs(signedArea(voronoi[i])) - std::fabs(signedArea(cells[i])));
    }
  }
}

// Itération de point fixe de l'article : chaque germe ordinaire va au centroïde c de sa cellule, déplacé vers le
// bundle b le plus proche (Mahalanobis) d'une fraction alpha / (1 + alpha) du segment [c, b]. Les bundles suivent
// bundleMotion : lloyd (centroïde), fixed (immobiles) ou attract (attirés comme les autres germes, vers le bundle
// le plus proche autre qu'eux-mêmes). Tous les germes sont mis à jour simultanément, avec les positions des
// bundles du début de l'itération.
// Avec une forme, c est le centroïde de la partie intérieure de la cellule et une cellule coupée par le bord va
// en p' = (c + alpha b + w c) / (1 + alpha + w), w = beta * aire extérieure : avec beta grand, elle suit Lloyd sans
// être attirée. Les germes extérieurs (cellule vide) restent immobiles.
void Tissue::step() {
  const Mahalanobis maha(P);
  std::vector<int> bundles;
  for (size_t i = 0; i < seeds.size(); i++) {
    if (isBundle[i] && inside(i)) {
      bundles.push_back((int)i);
    }
  }
  std::vector<V> next = seeds;
  double move = 0.0;
  for (size_t i = 0; i < seeds.size(); i++) {
    if (!inside(i) || (isBundle[i] && P.bundleMotion == MOTION_FIXED)) {
      continue;
    }
    V c = centroid(cells[i]);
    int b = -1;
    if (P.alpha > 0.0 && (!isBundle[i] || P.bundleMotion == MOTION_ATTRACT)) {
      b = nearestBundle(c, seeds, bundles, maha, isBundle[i] ? (int)i : -1);
    }
    const double w = P.beta * areaOut[i];
    V num = w > 0.0 ? c * (1.0 + w) : c;
    double den = 1.0 + w;
    if (b >= 0) {
      num += seeds[b] * P.alpha;
      den += P.alpha;
    }
    next[i] = den != 1.0 ? num * (1.0 / den) : num;
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
    if (!inside(i)) {
      continue;
    }
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
CellMesh Tissue::mesh() const {
  const double sp = spacing();
  const double tolV = 1e-8 * sp, hb = 1e-6 * sp;
  CellMesh m;
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
    if (ids.size() < 3) {
      ids.clear();
    }
    m.cells.push_back(ids);
  }
  // Bord du tissu : parois qui n'appartiennent qu'à une cellule
  std::unordered_map<long long, int> uses;
  const long long nv = (long long)m.verts.size();
  auto edgeKey = [&](int a, int b) { return (long long)std::min(a, b) * nv + std::max(a, b); };
  for (const auto &ids : m.cells) {
    for (size_t k = 0; k < ids.size(); k++) {
      uses[edgeKey(ids[k], ids[(k + 1) % ids.size()])]++;
    }
  }
  m.boundary.assign(m.cells.size(), 0);
  for (size_t c = 0; c < m.cells.size(); c++) {
    const auto &ids = m.cells[c];
    for (size_t k = 0; k < ids.size(); k++) {
      int a = ids[k], b = ids[(k + 1) % ids.size()];
      if (uses[edgeKey(a, b)] == 1) {
        m.boundary[c] = 1;
        m.boundaryEdges.emplace_back(a, b);
      }
    }
  }
  return m;
}

Metrics computeMetrics(const Tissue &t, const CellMesh &m) {
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
  r.nCells = all.n;
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
