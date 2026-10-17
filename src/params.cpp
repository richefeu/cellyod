#include "params.hpp"
#include "expr.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

const std::vector<std::string> LATTICES = {"hex", "square", "random", "uniform"};
const std::vector<std::string> BUNDLE_MOTIONS = {"lloyd", "fixed", "attract"};

namespace {

ParamSpec real(const char *name, const char *group, double Params::*m, double lo, double hi, bool reseed,
               const char *help, bool logScale = false) {
  ParamSpec s{name, group, help, ParamKind::Real};
  s.real = m;
  s.lo = lo;
  s.hi = hi;
  s.reseed = reseed;
  s.logScale = logScale;
  return s;
}

ParamSpec integer(const char *name, const char *group, int Params::*m, double lo, double hi, bool reseed,
                  const char *help) {
  ParamSpec s{name, group, help, ParamKind::Int};
  s.integer = m;
  s.lo = lo;
  s.hi = hi;
  s.reseed = reseed;
  return s;
}

ParamSpec choice(const char *name, const char *group, int Params::*m, const std::vector<std::string> &c,
                 bool reseed, const char *help) {
  ParamSpec s{name, group, help, ParamKind::Choice};
  s.integer = m;
  s.choices = &c;
  s.lo = 0.0;
  s.hi = (double)c.size() - 1.0;
  s.reseed = reseed;
  return s;
}

ParamSpec text(const char *name, const char *group, std::string Params::*m, bool reseed, const char *help) {
  ParamSpec s{name, group, help, ParamKind::Text};
  s.str = m;
  s.reseed = reseed;
  return s;
}

std::string fmt(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.15g", v);
  return buf;
}

} // namespace

const std::vector<ParamSpec> &paramSpecs() {
  static const std::vector<ParamSpec> specs = {
      real("Lx", "Domaine", &Params::Lx, 2.0, 200.0, true, "largeur du domaine (en diamètres de cellule)"),
      real("Ly", "Domaine", &Params::Ly, 2.0, 200.0, true, "hauteur du domaine (en diamètres de cellule)"),
      text("shape", "Domaine", &Params::shape, true,
           "forme du tissu : expression de f(x, y), le tissu est la partie f <= 0 du rectangle (vide : tout le "
           "rectangle). Variables x y, r et t (rayon et angle depuis le centre), cx cy Lx Ly. Ex. : r - 18"),
      real("beta", "Domaine", &Params::beta, 1e-2, 1e8, false,
           "cellules coupées par le bord de la forme : p' = (c + alpha b + beta A_out c) / (1 + alpha + beta "
           "A_out), c centroïde de la partie intérieure, A_out aire extérieure (grand : pas d'attraction)",
           true),
      choice("lattice", "Germes", &Params::lattice, LATTICES, true,
             "état initial : hex (réseau hexagonal perturbé), square (grille carrée perturbée), random (tirage "
             "avec distance minimale, alias poisson), uniform (tirage uniforme)"),
      real("disorder", "Germes", &Params::disorder, 0.0, 1.0, true,
           "perturbation des germes des réseaux hex et square, en fraction du pas (0 : réseau parfait)"),
      integer("seed", "Germes", &Params::seed, 1, 1000, true, "graine du générateur aléatoire"),
      integer("lloydIterations", "Lloyd", &Params::lloydIterations, 0, 500, false,
              "nombre d'itérations de Lloyd (attraction vers les bundles comprise)"),
      real("bundleFraction", "Bundles", &Params::bundleFraction, 0.0, 0.3, true,
           "fraction des germes tirés au hasard comme bundles"),
      real("alpha", "Bundles", &Params::alpha, 0.0, 0.5, false,
           "intensité de l'attraction vers le bundle le plus proche : p' = (c + alpha b) / (1 + alpha)"),
      real("anisotropy", "Bundles", &Params::anisotropy, 1.0, 1e4, false,
           "K, anisotropie de la distance de Mahalanobis (1 : euclidienne ; >> 1 : colonnes). Forme "
           "« anisotropy K angle » acceptée",
           true),
      real("anisotropyAngle", "Bundles", &Params::anisotropyAngle, 0.0, 180.0, false,
           "orientation des colonnes, en degrés depuis l'axe x (90 : Sigma = diag(1, K))"),
      choice("bundleMotion", "Bundles", &Params::bundleMotion, BUNDLE_MOTIONS, false,
             "déplacement des bundles : lloyd (centroïde), fixed (immobiles), attract (attirés comme les autres)"),
  };
  return specs;
}

const ParamSpec *findParam(const std::string &name) {
  for (const auto &s : paramSpecs()) {
    if (name == s.name) {
      return &s;
    }
  }
  return nullptr;
}

double ParamSpec::get(const Params &p) const {
  return kind == ParamKind::Real ? p.*real : kind == ParamKind::Text ? 0.0 : (double)(p.*integer);
}

void ParamSpec::set(Params &p, double v) const {
  if (kind == ParamKind::Real) {
    p.*real = v;
  } else if (kind != ParamKind::Text) {
    p.*integer = (int)std::lround(v);
  }
}

std::string ParamSpec::text(const Params &p) const {
  switch (kind) {
  case ParamKind::Real:
    return fmt(p.*real);
  case ParamKind::Int:
    return std::to_string(p.*integer);
  case ParamKind::Choice:
    return (*choices)[p.*integer];
  case ParamKind::Text:
    return p.*str;
  }
  return "";
}

// Variables utilisables dans l'expression de la forme (cf. Domain, tissue.cpp)
const std::vector<std::string> &shapeVariables() {
  static const std::vector<std::string> v = {"x", "y", "r", "t", "cx", "cy", "Lx", "Ly"};
  return v;
}

void ParamSpec::parse(Params &p, const std::string &value) const {
  if (kind == ParamKind::Text) {
    size_t b = value.find_first_not_of(" \t\r"), e = value.find_last_not_of(" \t\r");
    std::string v = b == std::string::npos ? "" : value.substr(b, e - b + 1);
    if (str == &Params::shape && !v.empty()) {
      Expr check(v, shapeVariables()); // lève une exception si l'expression est invalide
    }
    p.*str = v;
    return;
  }
  std::istringstream is(value);
  std::string extra;
  if (kind == ParamKind::Choice) {
    std::string v;
    is >> v;
    if (v == "poisson" && choices == &LATTICES) {
      v = "random";
    }
    for (size_t i = 0; i < choices->size(); i++) {
      if ((*choices)[i] == v) {
        p.*integer = (int)i;
        return;
      }
    }
    std::string list;
    for (const auto &c : *choices) {
      list += (list.empty() ? "" : ", ") + c;
    }
    throw std::runtime_error(std::string(name) + " : valeur « " + v + " » inconnue (" + list + ")");
  }
  double v;
  if (!(is >> v) || (is >> extra)) {
    throw std::runtime_error(std::string(name) + " : valeur invalide « " + value + " »");
  }
  set(p, v);
}

void setParam(Params &p, const std::string &key, const std::string &value) {
  std::istringstream is(value);
  if (key == "anisotropy") { // anisotropy K [angle]
    double K, deg;
    if (!(is >> K)) {
      throw std::runtime_error("anisotropy : valeur invalide « " + value + " »");
    }
    p.anisotropy = K;
    if (is >> deg) {
      p.anisotropyAngle = deg;
    }
  } else if (key == "bundle") {
    V b;
    if (!(is >> b.x >> b.y)) {
      throw std::runtime_error("bundle : deux coordonnées attendues");
    }
    p.bundleAt.push_back(b);
  } else if (key == "svgSeeds" || key == "svgLegend") {
    int v;
    if (!(is >> v)) {
      throw std::runtime_error(key + " : 0 ou 1 attendu");
    }
    (key == "svgSeeds" ? p.svgSeeds : p.svgLegend) = v != 0;
  } else if (key == "output" || key == "svg" || key == "morphology") {
    std::string v;
    is >> v;
    (key == "output" ? p.output : key == "svg" ? p.svg : p.morphology) = v;
  } else if (const ParamSpec *s = findParam(key)) {
    s->parse(p, value);
  } else {
    throw std::runtime_error("mot-clé inconnu : " + key);
  }
  if (p.Lx <= 0.0 || p.Ly <= 0.0) {
    throw std::runtime_error("Lx et Ly doivent être positifs");
  }
  if (p.anisotropy <= 0.0) {
    throw std::runtime_error("anisotropy doit être positive");
  }
  if (p.beta < 0.0) {
    throw std::runtime_error("beta doit être positif");
  }
  if (p.alpha < 0.0 || p.bundleFraction < 0.0 || p.bundleFraction > 1.0) {
    throw std::runtime_error("alpha doit être positif et bundleFraction compris entre 0 et 1");
  }
}

Params readParams(const std::string &file) {
  std::ifstream in(file);
  if (!in) {
    throw std::runtime_error("impossible d'ouvrir " + file);
  }
  Params p;
  std::string line;
  int n = 0;
  while (std::getline(in, line)) {
    n++;
    auto hash = line.find('#');
    if (hash != std::string::npos) {
      line.resize(hash);
    }
    std::istringstream is(line);
    std::string key;
    if (!(is >> key)) {
      continue;
    }
    std::string rest;
    std::getline(is, rest);
    try {
      setParam(p, key, rest);
    } catch (const std::exception &e) {
      throw std::runtime_error(file + ":" + std::to_string(n) + " : " + e.what());
    }
  }
  return p;
}

void writeParams(std::ostream &os, const Params &p, const char *prefix) {
  for (const auto &s : paramSpecs()) {
    os << prefix << s.name << ' ' << s.text(p) << '\n';
  }
  for (const V &b : p.bundleAt) {
    os << prefix << "bundle " << fmt(b.x) << ' ' << fmt(b.y) << '\n';
  }
}

void printParamHelp(std::ostream &os) {
  Params def;
  const char *group = "";
  for (const auto &s : paramSpecs()) {
    if (std::string(group) != s.group) {
      group = s.group;
      os << "\n  [" << group << "]\n";
    }
    os << "  " << s.name << " (défaut " << s.text(def) << ")\n      " << s.help << "\n";
  }
  os << "\n  [Autres]\n"
        "  bundle x y\n      bundle imposé : germe le plus proche du point (répétable)\n"
        "  output fichier (défaut cells.txt)\n      fichier de cellules produit\n"
        "  svg fichier (défaut cells.svg)\n      vue SVG produite (none : pas de SVG)\n"
        "  svgSeeds 0|1, svgLegend 0|1 (défaut 1)\n      SVG : germes, légende des bundles\n"
        "  morphology fichier\n      ajoute une ligne de métriques morphologiques au fichier\n";
}
