// Géométrie 2D élémentaire : vecteurs et polygones.
#pragma once

#include <cmath>
#include <vector>

struct V {
  double x{0.0}, y{0.0};
  V() = default;
  V(double X, double Y) : x(X), y(Y) {}
  V operator+(const V &o) const { return V(x + o.x, y + o.y); }
  V operator-(const V &o) const { return V(x - o.x, y - o.y); }
  V operator*(double k) const { return V(x * k, y * k); }
  V &operator+=(const V &o) {
    x += o.x;
    y += o.y;
    return *this;
  }
};

inline double dot(const V &a, const V &b) { return a.x * b.x + a.y * b.y; }
inline double cross(const V &a, const V &b) { return a.x * b.y - a.y * b.x; }
inline double norm(const V &a) { return std::sqrt(dot(a, a)); }

using Polygon = std::vector<V>;

inline double signedArea(const Polygon &p) {
  double a = 0.0;
  for (size_t i = 0; i < p.size(); i++) {
    a += cross(p[i], p[(i + 1) % p.size()]);
  }
  return 0.5 * a;
}

inline double perimeter(const Polygon &p) {
  double l = 0.0;
  for (size_t i = 0; i < p.size(); i++) {
    l += norm(p[(i + 1) % p.size()] - p[i]);
  }
  return l;
}

inline V centroid(const Polygon &p) {
  double a = 0.0;
  V c;
  for (size_t i = 0; i < p.size(); i++) {
    const V &u = p[i];
    const V &w = p[(i + 1) % p.size()];
    double k = cross(u, w);
    a += k;
    c += (u + w) * k;
  }
  if (std::fabs(a) < 1e-300) {
    V m;
    for (auto &q : p) {
      m += q;
    }
    return m * (1.0 / (double)p.size());
  }
  return c * (1.0 / (3.0 * a));
}

inline bool isConvex(const Polygon &p) {
  const size_t n = p.size();
  for (size_t i = 0; i < n; i++) {
    V e1 = p[(i + 1) % n] - p[i], e2 = p[(i + 2) % n] - p[(i + 1) % n];
    if (cross(e1, e2) < -1e-12 * norm(e1) * norm(e2)) {
      return false;
    }
  }
  return true;
}

// Polygone (sens trigonométrique) décalé de d vers l'intérieur : chaque côté recule de d parallèlement à
// lui-même. Polygone convexe : intersection des demi-plans décalés (exacte ; les côtés trop courts disparaissent
// naturellement). Polygone non convexe : sommets décalés en onglet, le long de la bissectrice. Renvoie un
// polygone vide si la cellule disparaît ou se retourne (d trop grand pour elle).
inline Polygon insetPolygon(const Polygon &p, double d) {
  const size_t n = p.size();
  if (n < 3 || d == 0.0) {
    return p;
  }
  if (isConvex(p)) {
    Polygon out = p;
    for (size_t i = 0; i < n && out.size() >= 3; i++) {
      const V &A = p[i], &B = p[(i + 1) % n];
      double L = norm(B - A);
      if (L <= 0.0) {
        continue;
      }
      auto g = [&](const V &q) { return cross(B - A, q - A) / L - d; }; // >= 0 : côté conservé
      Polygon next;
      for (size_t k = 0; k < out.size(); k++) {
        const V &P = out[k], &Q = out[(k + 1) % out.size()];
        double gp = g(P), gq = g(Q);
        if (gp >= 0.0) {
          next.push_back(P);
        }
        if ((gp >= 0.0) != (gq >= 0.0)) {
          next.push_back(P + (Q - P) * (gp / (gp - gq)));
        }
      }
      out = next;
    }
    return out.size() >= 3 && signedArea(out) > 0.0 ? out : Polygon();
  }
  // Non convexe : chaque côté devient une droite décalée de d ; les sommets sont les intersections des droites
  // consécutives. Un côté qui se retourne (trop court pour le décalage) est supprimé et les sommets recalculés,
  // comme il disparaîtrait dans l'intersection des demi-plans d'un polygone convexe.
  struct Line {
    V p, u; // point et direction unitaire de la droite décalée
    V start; // sommet de départ du côté, décalé selon la normale (intersection avec une droite parallèle)
  };
  std::vector<Line> lines;
  for (size_t i = 0; i < n; i++) {
    V e = p[(i + 1) % n] - p[i];
    double L = norm(e);
    if (L > 0.0) {
      V u = e * (1.0 / L), nrm(-u.y, u.x); // normale intérieure
      lines.push_back({p[i] + nrm * d, u, p[i] + nrm * d});
    }
  }
  auto meet = [](const Line &a, const Line &b) {
    double c = cross(a.u, b.u);
    if (std::fabs(c) < 1e-9) { // côtés alignés
      return b.start;
    }
    return a.p + a.u * (cross(b.p - a.p, b.u) / c);
  };
  Polygon out;
  for (size_t pass = 0; pass <= n; pass++) {
    const size_t m = lines.size();
    if (m < 3) {
      return Polygon();
    }
    out.resize(m);
    for (size_t j = 0; j < m; j++) {
      out[j] = meet(lines[(j + m - 1) % m], lines[j]);
    }
    std::vector<Line> kept;
    for (size_t j = 0; j < m; j++) {
      if (dot(out[(j + 1) % m] - out[j], lines[j].u) > 0.0) {
        kept.push_back(lines[j]);
      }
    }
    if (kept.size() == m) {
      break;
    }
    lines = kept;
  }
  return out.size() >= 3 && signedArea(out) > 0.0 ? out : Polygon();
}
