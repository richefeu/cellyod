// mesure : indicateurs de désordre d'un pavage cellyod (fichiers de cellules).
//
// Usage : mesure cells1.txt [cells2.txt ...]     indicateurs, une ligne par fichier (en-tête en première ligne)
//         mesure -hist cells.txt                  histogrammes des angles aux sommets et des orientations de paroi
//
// Les indicateurs portent sur les cellules « profondes » : ni cellule de bord (paroi sur le bord du tissu), ni
// voisine d'une cellule de bord, pour ne pas mêler le désordre et l'effet du bord (les cellules de la deuxième
// couronne sont déformées par celles du bord). Pour le taux de transitions T1, les fichiers d'un même calcul doivent
// être donnés dans l'ordre des itérations.
//
//   iteration  nombre d'itérations de Lloyd du fichier
//   n          nombre de cellules profondes
// Tailles et formes :
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
// Orientation des parois (parois entre deux cellules profondes, pondérées par leur longueur ; theta : angle de la
// paroi avec l'axe x, dans [0, pi[) :
//   S6         ordre orientationnel global | sum l exp(6 i theta) | / sum l : 1 si toutes les parois suivent les
//              trois directions d'un même réseau hexagonal, proche de 0 si les orientations sont réparties au
//              hasard ou en domaines d'orientations différentes
//   S2         anisotropie | sum l exp(2 i theta) | / sum l : 0 sans direction privilégiée
//   sigmaV     écart type (degrés) des angles entre parois consécutives aux sommets profonds (120° partout pour
//              un pavage hexagonal régulier : 0)
// Connectivité :
//   deg3       fraction des sommets profonds de degré 3 (pavage de Voronoï générique : 1)
//   aboav      paramètre a de la loi d'Aboav-Weaire n m(n) = (6 - a) n + 6 a + mu2, m(n) nombre moyen de côtés des
//              voisins d'une cellule à n côtés (régression sur les cellules profondes ; 0 si tous les n sont égaux)
//   f57        fraction des parois séparant un pentagone et un heptagone
//   p5, p7     fractions de pentagones et d'heptagones (sans corrélation entre voisins, on attendrait
//              f57 = 2 q5 q7, avec q_n = n p_n / 6 la fraction des extrémités de paroi sur une cellule à n côtés)
//   paired     fraction des cellules non hexagonales ayant un voisin de défaut opposé (n < 6 contre n > 6)
//   t1         transitions T1 par cellule et par itération depuis le fichier précédent (voisinages perdus et
//              gagnés / 4 : une transition T1 en supprime un et en crée un ; 0 pour le premier fichier)
#include "io.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <unordered_map>
#include <unordered_set>

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

// Angle (radians, dans [0, pi[) de la paroi [a, b]
double edgeAngle(const V &a, const V &b) {
  double th = std::atan2(b.y - a.y, b.x - a.x);
  if (th < 0.0) {
    th += M_PI;
  }
  return th >= M_PI ? th - M_PI : th;
}

// Topologie du pavage : voisinages des cellules et des sommets, cellules et sommets profonds
struct Analysis {
  CellMesh m;
  std::vector<std::vector<int>> neigh;        // voisins de Voronoï (cellules partageant une paroi)
  std::vector<std::pair<int, int>> edges;     // parois entre deux cellules (indices de sommets)
  std::vector<std::pair<int, int>> edgeCells; // les deux cellules de chaque paroi
  std::vector<char> deep;                     // cellules profondes
  std::vector<std::vector<int>> vertNeigh;    // sommets reliés à chaque sommet par une paroi
  std::vector<char> deepVert;                 // sommets dont toutes les cellules sont profondes

  explicit Analysis(const Tissue &t) : m(t.mesh()) {
    const size_t nc = m.cells.size();
    std::unordered_map<long long, std::vector<int>> owners;
    const long long nv = (long long)m.verts.size();
    auto key = [&](int a, int b) { return (long long)std::min(a, b) * nv + std::max(a, b); };
    for (size_t c = 0; c < nc; c++) {
      const auto &ids = m.cells[c];
      for (size_t k = 0; k < ids.size(); k++) {
        owners[key(ids[k], ids[(k + 1) % ids.size()])].push_back((int)c);
      }
    }
    neigh.resize(nc);
    vertNeigh.resize(m.verts.size());
    for (const auto &e : owners) {
      int u = (int)(e.first / nv), w = (int)(e.first % nv);
      vertNeigh[u].push_back(w);
      vertNeigh[w].push_back(u);
      if (e.second.size() == 2) {
        neigh[e.second[0]].push_back(e.second[1]);
        neigh[e.second[1]].push_back(e.second[0]);
        edges.emplace_back(u, w);
        edgeCells.emplace_back(e.second[0], e.second[1]);
      }
    }
    deep.assign(nc, 0);
    for (size_t c = 0; c < nc; c++) {
      deep[c] = !m.cells[c].empty() && !m.boundary[c];
      for (int j : neigh[c]) {
        deep[c] = deep[c] && !m.boundary[j];
      }
    }
    deepVert.assign(m.verts.size(), 1);
    std::vector<char> used(m.verts.size(), 0);
    for (size_t c = 0; c < nc; c++) {
      for (int v : m.cells[c]) {
        used[v] = 1;
        deepVert[v] = deepVert[v] && deep[c];
      }
    }
    for (size_t v = 0; v < m.verts.size(); v++) {
      deepVert[v] = deepVert[v] && used[v];
    }
  }

  bool deepEdge(size_t e) const { return deep[edgeCells[e].first] && deep[edgeCells[e].second]; }

  // Angles (degrés) entre parois consécutives autour du sommet v
  std::vector<double> vertexAngles(int v) const {
    std::vector<double> th;
    for (int w : vertNeigh[v]) {
      th.push_back(std::atan2(m.verts[w].y - m.verts[v].y, m.verts[w].x - m.verts[v].x));
    }
    std::sort(th.begin(), th.end());
    std::vector<double> out;
    for (size_t k = 0; k < th.size(); k++) {
      double d = (k + 1 < th.size() ? th[k + 1] : th[0] + 2.0 * M_PI) - th[k];
      out.push_back(d * 180.0 / M_PI);
    }
    return out;
  }
};

// Histogrammes (densités par degré, classes de 2°) des angles aux sommets profonds et des orientations des parois
// profondes (pondérées par leur longueur)
int histograms(const char *file) {
  Tissue t = loadCells(file);
  Analysis an(t);
  const int nb = 90;
  std::vector<double> hv(nb, 0.0), he(nb, 0.0);
  double sv = 0.0, se = 0.0;
  for (size_t v = 0; v < an.m.verts.size(); v++) {
    if (an.deepVert[v]) {
      for (double a : an.vertexAngles((int)v)) {
        hv[std::min(nb - 1, (int)(a / 2.0))] += 1.0;
        sv += 1.0;
      }
    }
  }
  for (size_t e = 0; e < an.edges.size(); e++) {
    if (an.deepEdge(e)) {
      const V &a = an.m.verts[an.edges[e].first], &b = an.m.verts[an.edges[e].second];
      double l = norm(b - a);
      he[std::min(nb - 1, (int)(edgeAngle(a, b) * 180.0 / M_PI / 2.0))] += l;
      se += l;
    }
  }
  std::printf("angle sommets parois\n");
  for (int k = 0; k < nb; k++) {
    std::printf("%g %.6g %.6g\n", 2.0 * k + 1.0, hv[k] / (2.0 * sv), he[k] / (2.0 * se));
  }
  return 0;
}

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage : mesure cells1.txt [cells2.txt ...]\n        mesure -hist cells.txt\n");
    return 1;
  }
  try {
    if (!std::strcmp(argv[1], "-hist") && argc == 3) {
      return histograms(argv[2]);
    }
    std::printf("iteration n cvA cvL sigmaN p6 SI G psi6 cvNN residu S6 S2 sigmaV deg3 aboav f57 paired t1 p5 p7\n");
    std::unordered_set<long long> prevPairs;
    size_t prevNc = 0;
    int prevIt = -1;
    for (int f = 1; f < argc; f++) {
      Tissue t = loadCells(argv[f]);
      Analysis an(t);
      const CellMesh &m = an.m;
      const size_t nc = m.cells.size();

      // longueurs et orientations des parois profondes
      Stat len;
      std::complex<double> z6(0.0, 0.0), z2(0.0, 0.0);
      double sumL = 0.0;
      size_t nEdges = 0, n57 = 0;
      for (size_t e = 0; e < an.edges.size(); e++) {
        if (!an.deepEdge(e)) {
          continue;
        }
        const V &a = m.verts[an.edges[e].first], &b = m.verts[an.edges[e].second];
        double l = norm(b - a), th = edgeAngle(a, b);
        len.add(l);
        z6 += l * std::polar(1.0, 6.0 * th);
        z2 += l * std::polar(1.0, 2.0 * th);
        sumL += l;
        size_t na = m.cells[an.edgeCells[e].first].size(), nb = m.cells[an.edgeCells[e].second].size();
        nEdges++;
        n57 += (na == 5 && nb == 7) || (na == 7 && nb == 5);
      }

      // angles et degrés des sommets profonds
      Stat vang;
      size_t nVert = 0, nDeg3 = 0;
      for (size_t v = 0; v < m.verts.size(); v++) {
        if (an.deepVert[v]) {
          nVert++;
          nDeg3 += an.vertNeigh[v].size() == 3;
          for (double a : an.vertexAngles((int)v)) {
            vang.add(a);
          }
        }
      }

      // cellules profondes
      Stat area, sides, si, nn, res;
      double inertia = 0.0, psi = 0.0;
      size_t hex = 0, defects = 0, paired = 0, n5 = 0, n7 = 0;
      double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0; // Aboav-Weaire : y = n m(n) en fonction de x = n
      std::unordered_set<long long> pairs;              // voisinages entre cellules profondes
      for (size_t c = 0; c < nc; c++) {
        if (!an.deep[c]) {
          continue;
        }
        Polygon p;
        for (int v : m.cells[c]) {
          p.push_back(m.verts[v]);
        }
        const double ns = (double)p.size();
        double A = std::fabs(signedArea(p)), P = perimeter(p);
        area.add(A);
        sides.add(ns);
        hex += p.size() == 6;
        n5 += p.size() == 5;
        n7 += p.size() == 7;
        si.add(4.0 * M_PI * A / (P * P));
        inertia += polarMoment(p, t.seeds[c]);
        res.add(norm(centroid(p) - t.seeds[c]) / Tissue::spacing());
        std::complex<double> z(0.0, 0.0);
        double dmin = 1e300, sumSides = 0.0;
        bool opposite = false;
        for (int j : an.neigh[c]) {
          V d = t.seeds[j] - t.seeds[c];
          z += std::polar(1.0, 6.0 * std::atan2(d.y, d.x));
          dmin = std::min(dmin, norm(d));
          double nj = (double)m.cells[j].size();
          sumSides += nj;
          opposite = opposite || (ns - 6.0) * (nj - 6.0) < 0.0;
          if (an.deep[j]) {
            pairs.insert((long long)std::min<size_t>(c, j) * (long long)nc + (long long)std::max<size_t>(c, j));
          }
        }
        if (!an.neigh[c].empty()) {
          psi += std::abs(z) / (double)an.neigh[c].size();
          nn.add(dmin);
          sx += ns;
          sy += sumSides; // n m(n) : somme des côtés des voisins
          sxx += ns * ns;
          sxy += ns * sumSides;
        }
        if (p.size() != 6) {
          defects++;
          paired += opposite;
        }
      }
      const double n = (double)area.n;
      const double varx = sxx / n - (sx / n) * (sx / n);
      const double aboav = varx > 1e-12 ? 6.0 - (sxy / n - (sx / n) * (sy / n)) / varx : 0.0;

      // transitions T1 depuis le fichier précédent (même tissu, germes de même numéro)
      double t1 = 0.0;
      if (prevIt >= 0 && prevNc == nc && t.iteration > prevIt) {
        size_t changed = 0;
        for (long long q : prevPairs) {
          changed += !pairs.count(q);
        }
        for (long long q : pairs) {
          changed += !prevPairs.count(q);
        }
        t1 = (double)changed / 4.0 / n / (double)(t.iteration - prevIt);
      }
      prevPairs = pairs;
      prevNc = nc;
      prevIt = t.iteration;

      std::printf("%d %zu %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n",
                  t.iteration, area.n, area.cv(), len.cv(), sides.sd(), hex / n, si.mean(),
                  inertia / (2.0 * n * area.mean() * area.mean()), psi / n, nn.cv(), res.mean(),
                  sumL > 0 ? std::abs(z6) / sumL : 0.0, sumL > 0 ? std::abs(z2) / sumL : 0.0, vang.sd(),
                  nVert ? (double)nDeg3 / (double)nVert : 0.0, aboav, nEdges ? (double)n57 / (double)nEdges : 0.0,
                  defects ? (double)paired / (double)defects : 0.0, t1, n5 / n, n7 / n);
    }
  } catch (const std::exception &e) {
    std::fprintf(stderr, "mesure : %s\n", e.what());
    return 1;
  }
  return 0;
}
