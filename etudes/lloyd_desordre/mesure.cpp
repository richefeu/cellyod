// mesure : indicateurs de désordre d'un pavage cellyod (fichiers de cellules).
//
// Usage : mesure cells1.txt [cells2.txt ...]
// Écrit une ligne par fichier (en-tête en première ligne). Seules les cellules « profondes » sont comptées : ni
// cellule de bord (paroi sur le bord du tissu), ni voisine d'une cellule de bord, pour ne pas mêler le désordre et
// l'effet du bord (les cellules de la deuxième couronne sont déformées par celles du bord).
//
//   iteration  nombre d'itérations de Lloyd du fichier
//   n          nombre de cellules profondes
//   cvA        coefficient de variation des aires
//   cvL        coefficient de variation des longueurs de paroi (parois entre deux cellules profondes)
//   sigmaN     écart type du nombre de côtés
//   p6         fraction d'hexagones
//   SI         indice de forme moyen 4 pi A / P^2 (hexagone régulier : 0.9069)
//   G          moment d'inertie normalisé : somme des int |x - p|^2 sur les cellules / (2 n <A>^2), p germe
//              (énergie de Lloyd adimensionnée ; pavage hexagonal régulier : 5 / (36 sqrt 3) = 0.080188)
//   psi6       ordre orientationnel local < |1/n_i sum_j exp(6 i theta_ij)| > sur les voisins de Voronoï
//              (1 : réseau hexagonal parfait)
//   cvNN       coefficient de variation de la distance d'un germe à son plus proche voisin
//   residu     distance moyenne germe - centroïde, rapportée à l'espacement visé (0 : point fixe de Lloyd)
#include "io.hpp"

#include <cmath>
#include <complex>
#include <cstdio>
#include <unordered_map>

namespace {

struct Stat {
  double s{0.0}, s2{0.0};
  size_t n{0};
  void add(double v) {
    s += v;
    s2 += v * v;
    n++;
  }
  double mean() const { return n ? s / (double)n : 0.0; }
  double sd() const { return n > 1 ? std::sqrt(std::max(0.0, s2 / (double)n - mean() * mean())) : 0.0; }
  double cv() const { return mean() != 0.0 ? sd() / mean() : 0.0; }
};

// Moment d'inertie polaire du polygone par rapport au point p : int |x - p|^2 dx
double polarMoment(const Polygon &poly, const V &p) {
  double I = 0.0;
  for (size_t k = 0; k < poly.size(); k++) {
    V a = poly[k] - p, b = poly[(k + 1) % poly.size()] - p;
    I += cross(a, b) * (a.x * a.x + a.x * b.x + b.x * b.x + a.y * a.y + a.y * b.y + b.y * b.y);
  }
  return std::fabs(I) / 12.0;
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage : mesure cells1.txt [cells2.txt ...]\n");
    return 1;
  }
  std::printf("iteration n cvA cvL sigmaN p6 SI G psi6 cvNN residu\n");
  for (int f = 1; f < argc; f++) {
    try {
      Tissue t = loadCells(argv[f]);
      CellMesh m = t.mesh();
      const size_t nc = m.cells.size();

      // voisins de Voronoï : cellules partageant une paroi
      std::unordered_map<long long, std::vector<int>> owners;
      const long long nv = (long long)m.verts.size();
      auto key = [&](int a, int b) { return (long long)std::min(a, b) * nv + std::max(a, b); };
      for (size_t c = 0; c < nc; c++) {
        const auto &ids = m.cells[c];
        for (size_t k = 0; k < ids.size(); k++) {
          owners[key(ids[k], ids[(k + 1) % ids.size()])].push_back((int)c);
        }
      }
      std::vector<std::vector<int>> neigh(nc);
      for (const auto &e : owners) {
        if (e.second.size() == 2) {
          neigh[e.second[0]].push_back(e.second[1]);
          neigh[e.second[1]].push_back(e.second[0]);
        }
      }
      // cellules profondes : ni de bord, ni voisines d'une cellule de bord
      std::vector<char> deep(nc, 0);
      for (size_t c = 0; c < nc; c++) {
        deep[c] = !m.cells[c].empty() && !m.boundary[c];
        for (int j : neigh[c]) {
          deep[c] = deep[c] && !m.boundary[j];
        }
      }
      Stat len;
      for (const auto &e : owners) {
        if (e.second.size() == 2 && deep[e.second[0]] && deep[e.second[1]]) {
          long long u = e.first / nv, w = e.first % nv;
          len.add(norm(m.verts[u] - m.verts[w]));
        }
      }

      Stat area, sides, si, nn, res;
      double inertia = 0.0, psi = 0.0;
      size_t hex = 0;
      for (size_t c = 0; c < nc; c++) {
        if (!deep[c]) {
          continue;
        }
        Polygon p;
        for (int v : m.cells[c]) {
          p.push_back(m.verts[v]);
        }
        double A = std::fabs(signedArea(p)), P = perimeter(p);
        area.add(A);
        sides.add((double)p.size());
        hex += p.size() == 6;
        si.add(4.0 * M_PI * A / (P * P));
        inertia += polarMoment(p, t.seeds[c]);
        res.add(norm(centroid(p) - t.seeds[c]) / Tissue::spacing());
        std::complex<double> z(0.0, 0.0);
        double dmin = 1e300;
        for (int j : neigh[c]) {
          V d = t.seeds[j] - t.seeds[c];
          z += std::polar(1.0, 6.0 * std::atan2(d.y, d.x));
          dmin = std::min(dmin, norm(d));
        }
        if (!neigh[c].empty()) {
          psi += std::abs(z) / (double)neigh[c].size();
          nn.add(dmin);
        }
      }
      const double n = (double)area.n;
      std::printf("%d %zu %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n", t.iteration, area.n, area.cv(), len.cv(),
                  sides.sd(), hex / n, si.mean(), inertia / (2.0 * n * area.mean() * area.mean()), psi / n, nn.cv(),
                  res.mean());
    } catch (const std::exception &e) {
      std::fprintf(stderr, "mesure : %s\n", e.what());
      return 1;
    }
  }
  return 0;
}
