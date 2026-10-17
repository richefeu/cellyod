// Paramètres de cellyod et leur description (registre).
//
// Chaque paramètre scalaire est décrit une seule fois dans paramSpecs() (params.cpp) : nom du mot-clé, groupe,
// aide, bornes du curseur et effet d'une modification. La lecture des fichiers de paramètres, l'écriture dans
// les fichiers de cellules, l'aide de la ligne de commande et le panneau de l'interface sont construits à partir
// de ce registre : ajouter un paramètre = un champ dans Params + une ligne dans paramSpecs().
//
// Toutes les grandeurs sont sans dimension : l'unité de longueur est le diamètre équivalent visé des cellules
// (aire visée d'une cellule : pi/4).
#pragma once

#include "geom.hpp"

#include <iosfwd>
#include <string>
#include <vector>

struct Params {
  // Domaine rectangulaire [0, Lx] x [0, Ly]
  double Lx{40.0}, Ly{40.0};
  // Germes initiaux
  int lattice{0};       // indice dans LATTICES
  double disorder{0.3}; // réseaux hex et square : perturbation des germes, en fraction du pas
  int seed{1};
  // Itérations de Lloyd (attraction vers les bundles comprise)
  int lloydIterations{10};
  // Tissu végétal (Trivino et al., morphogenèse de Voronoï variationnelle)
  double bundleFraction{0.0};    // fraction des germes tirés au hasard comme bundles
  double alpha{0.0};             // intensité de l'attraction : p' = (c + alpha b) / (1 + alpha)
  double anisotropy{1.0};        // K : rapport des valeurs propres de Sigma
  double anisotropyAngle{90.0};  // orientation des colonnes (degrés depuis l'axe x)
  int bundleMotion{0};           // indice dans BUNDLE_MOTIONS
  std::vector<V> bundleAt;       // bundles imposés : germe le plus proche de chaque point
  // Sorties de la ligne de commande
  std::string output{"cells.txt"};
  std::string svg{"cells.svg"};
  std::string morphology; // fichier auquel une ligne de métriques est ajoutée (balayages)
};

enum Lattice { LATTICE_HEX = 0, LATTICE_SQUARE, LATTICE_RANDOM, LATTICE_UNIFORM };
enum BundleMotion { MOTION_LLOYD = 0, MOTION_FIXED, MOTION_ATTRACT };
extern const std::vector<std::string> LATTICES;
extern const std::vector<std::string> BUNDLE_MOTIONS;

enum class ParamKind { Real, Int, Choice };

struct ParamSpec {
  const char *name;
  const char *group;
  const char *help;
  ParamKind kind;
  double Params::*real{nullptr};
  int Params::*integer{nullptr}; // Int et Choice (indice dans choices)
  const std::vector<std::string> *choices{nullptr};
  double lo{0.0}, hi{1.0}; // bornes conseillées (curseurs de l'interface)
  bool logScale{false};
  bool reseed{false}; // une modification impose de régénérer les germes (sinon : effet dès l'itération suivante)

  double get(const Params &p) const;
  void set(Params &p, double v) const;
  std::string text(const Params &p) const;
  void parse(Params &p, const std::string &value) const; // lève std::runtime_error
};

const std::vector<ParamSpec> &paramSpecs();
const ParamSpec *findParam(const std::string &name);

// Affecte un mot-clé (registre, ou mots-clés spéciaux : anisotropy K [angle], bundle x y, output, svg,
// morphology). Lève std::runtime_error avec un message en clair.
void setParam(Params &p, const std::string &key, const std::string &value);
// Lit un fichier de paramètres « mot-clé valeur(s) », # pour les commentaires
Params readParams(const std::string &file);
// Écrit tous les paramètres du tissu (pas les sorties), une ligne par mot-clé précédé de prefix
void writeParams(std::ostream &os, const Params &p, const char *prefix = "");
// Aide : liste des mots-clés avec défaut et description
void printParamHelp(std::ostream &os);
