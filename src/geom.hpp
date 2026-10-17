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
