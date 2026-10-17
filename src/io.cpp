#include "io.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

void saveCells(const std::string &file, const Tissue &t, const CellMesh &m) {
  std::ofstream os(file);
  if (!os) {
    throw std::runtime_error("impossible d'écrire " + file);
  }
  os.precision(17);
  os << "cellyod 1\n";
  writeParams(os, t.P, "param ");
  os << "domain 0 0 " << t.P.Lx << ' ' << t.P.Ly << '\n';
  os << "iteration " << t.iteration << '\n';
  os << "vertices " << m.verts.size() << '\n';
  for (const V &q : m.verts) {
    os << q.x << ' ' << q.y << '\n';
  }
  size_t nIn = 0;
  for (const auto &ids : m.cells) {
    nIn += !ids.empty();
  }
  os << "cells " << nIn << '\n';
  for (size_t c = 0; c < m.cells.size(); c++) {
    if (m.cells[c].empty()) {
      continue;
    }
    os << (int)t.isBundle[c] << ' ' << t.seeds[c].x << ' ' << t.seeds[c].y << ' ' << m.cells[c].size();
    for (int v : m.cells[c]) {
      os << ' ' << v;
    }
    os << '\n';
  }
  if (nIn < m.cells.size()) {
    os << "frozen " << m.cells.size() - nIn << '\n';
    for (size_t c = 0; c < m.cells.size(); c++) {
      if (m.cells[c].empty()) {
        os << t.seeds[c].x << ' ' << t.seeds[c].y << '\n';
      }
    }
  }
}

Tissue loadCells(const std::string &file) {
  std::ifstream in(file);
  if (!in) {
    throw std::runtime_error("impossible d'ouvrir " + file);
  }
  Params p;
  std::vector<V> seeds;
  std::vector<char> bundle;
  int it = 0;
  bool header = false;
  std::string line;
  // ligne suivante sans commentaire, ni ligne vide
  auto next = [&](std::istringstream &is) {
    while (std::getline(in, line)) {
      auto hash = line.find('#');
      if (hash != std::string::npos) {
        line.resize(hash);
      }
      if (line.find_first_not_of(" \t\r") != std::string::npos) {
        is.clear();
        is.str(line);
        return true;
      }
    }
    return false;
  };
  std::istringstream is;
  while (next(is)) {
    std::string key;
    is >> key;
    if (key == "cellyod") {
      header = true;
    } else if (key == "param") {
      std::string k, rest;
      is >> k;
      std::getline(is, rest);
      setParam(p, k, rest);
    } else if (key == "iteration") {
      is >> it;
    } else if (key == "domain") {
      // redondant avec Lx, Ly
    } else if (key == "vertices") {
      size_t n = 0;
      is >> n;
      for (size_t i = 0; i < n; i++) {
        if (!next(is)) {
          throw std::runtime_error(file + " : sommets manquants");
        }
      }
    } else if (key == "cells") {
      size_t n = 0;
      is >> n;
      for (size_t i = 0; i < n; i++) {
        int type;
        V s;
        if (!next(is) || !(is >> type >> s.x >> s.y)) {
          throw std::runtime_error(file + " : cellule " + std::to_string(i) + " illisible");
        }
        seeds.push_back(s);
        bundle.push_back(type == 1);
      }
    } else if (key == "frozen") {
      size_t n = 0;
      is >> n;
      for (size_t i = 0; i < n; i++) {
        V s;
        if (!next(is) || !(is >> s.x >> s.y)) {
          throw std::runtime_error(file + " : germe extérieur " + std::to_string(i) + " illisible");
        }
        seeds.push_back(s);
        bundle.push_back(0);
      }
    } else {
      throw std::runtime_error(file + " : mot-clé inconnu " + key);
    }
  }
  if (!header) {
    throw std::runtime_error(file + " n'est pas un fichier cellyod");
  }
  Tissue t;
  t.restore(p, seeds, bundle, it);
  return t;
}

// Couleur 0xRRGGBB au format SVG
static std::string svgColor(uint32_t c) {
  char buf[8];
  std::snprintf(buf, sizeof(buf), "#%06x", (unsigned)c);
  return buf;
}

// Rendu selon le thème (theme.hpp) : parois et bord du tissu, germes cerclés, bundles en rouge (optionnels)
void writeSVG(const std::string &file, const Tissue &t, const CellMesh &m, const SvgOptions &opt) {
  std::ofstream os(file);
  if (!os) {
    throw std::runtime_error("impossible d'écrire " + file);
  }
  const Theme &th = renderTheme(opt.dark);
  const double s = opt.widthPx / t.P.Lx, H = t.P.Ly * s;
  auto X = [&](double x) { return x * s; };
  auto Y = [&](double y) { return H - y * s; }; // axe y vers le haut
  char buf[64];
  auto f = [&](double v) {
    std::snprintf(buf, sizeof(buf), "%.2f", v);
    return std::string(buf);
  };
  os << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << f(opt.widthPx) << "\" height=\"" << f(H)
     << "\" viewBox=\"0 0 " << f(opt.widthPx) << ' ' << f(H) << "\">\n";
  os << "<rect width=\"100%\" height=\"100%\" fill=\"" << svgColor(th.background) << "\"/>\n";
  const double sw = std::max(0.5, 0.04 * s);
  os << "<g fill=\"" << svgColor(th.cell) << "\" stroke=\"" << svgColor(th.edge) << "\" stroke-width=\"" << f(sw) << "\" stroke-linejoin=\"round\">\n";
  for (const auto &ids : m.cells) {
    if (ids.empty()) {
      continue;
    }
    os << "<polygon points=\"";
    for (int v : ids) {
      os << f(X(m.verts[v].x)) << ',' << f(Y(m.verts[v].y)) << ' ';
    }
    os << "\"/>\n";
  }
  os << "</g>\n";
  // bord du tissu
  os << "<path fill=\"none\" stroke=\"" << svgColor(th.boundary) << "\" stroke-width=\"" << f(2.0 * sw) << "\" stroke-linecap=\"round\" d=\"";
  for (const auto &e : m.boundaryEdges) {
    os << 'M' << f(X(m.verts[e.first].x)) << ',' << f(Y(m.verts[e.first].y)) << 'L' << f(X(m.verts[e.second].x))
       << ',' << f(Y(m.verts[e.second].y));
  }
  os << "\"/>\n";
  const double rs = std::max(1.0, 0.07 * s), rb = std::max(2.0, 0.22 * s);
  if (opt.seeds) {
    os << "<g fill=\"" << svgColor(th.seed) << "\" stroke=\"" << svgColor(th.seedRing) << "\" stroke-width=\"" << f(0.4 * rs) << "\">\n";
    for (size_t i = 0; i < t.seeds.size(); i++) {
      if (t.inside(i) && !(opt.bundles && t.isBundle[i])) {
        os << "<circle cx=\"" << f(X(t.seeds[i].x)) << "\" cy=\"" << f(Y(t.seeds[i].y)) << "\" r=\"" << f(rs)
           << "\"/>\n";
      }
    }
    os << "</g>\n";
  }
  if (opt.bundles) {
    os << "<g fill=\"" << svgColor(th.bundle) << "\">\n";
    for (size_t i = 0; i < t.seeds.size(); i++) {
      if (t.inside(i) && t.isBundle[i]) {
        os << "<circle cx=\"" << f(X(t.seeds[i].x)) << "\" cy=\"" << f(Y(t.seeds[i].y)) << "\" r=\"" << f(rb)
           << "\"/>\n";
      }
    }
    os << "</g>\n";
  }
  os << "</svg>\n";
}

void appendMorphology(const std::string &file, const Tissue &t, const Metrics &mt) {
  bool fresh = !std::filesystem::exists(file);
  std::ofstream os(file, std::ios::app);
  if (!os) {
    throw std::runtime_error("impossible d'écrire " + file);
  }
  const auto items = mt.items();
  if (fresh) {
    os << "#";
    for (const auto &s : paramSpecs()) {
      if (s.kind != ParamKind::Text) { // la forme (texte avec espaces) n'est pas une colonne
        os << ' ' << s.name;
      }
    }
    os << " iteration";
    for (const auto &it : items) {
      os << ' ' << it.first;
    }
    os << '\n';
  }
  os.precision(8);
  bool first = true;
  for (const auto &s : paramSpecs()) {
    if (s.kind != ParamKind::Text) {
      os << (first ? "" : " ") << s.text(t.P);
      first = false;
    }
  }
  os << ' ' << t.iteration;
  for (const auto &it : items) {
    os << ' ' << it.second;
  }
  os << '\n';
}
