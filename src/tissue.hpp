// Tissu cellulaire : germes, pavage de Voronoï et itérations de Lloyd modifiées.
//
// Modèle de Trivino et al. (morphogenèse de Voronoï variationnelle) : une fraction des germes sont des bundles
// (faisceaux vasculaires) ; à chaque itération, les autres germes vont en p' = (c + alpha b) / (1 + alpha), c
// centroïde de leur cellule et b bundle le plus proche au sens de la distance de Mahalanobis de matrice Sigma.
// Avec alpha = 0, c'est la régularisation de Lloyd. C'est une descente de E = E_CVL - alpha E_att.
#pragma once

#include "params.hpp"

#include <string>
#include <utility>
#include <vector>

// Domaine de pavage. Seul le rectangle [0, Lx] x [0, Ly] existe pour l'instant ; tout ce qui dépend de la forme
// du domaine passe par cette structure.
struct Domain {
  double Lx{1.0}, Ly{1.0};
  explicit Domain(const Params &p) : Lx(p.Lx), Ly(p.Ly) {}
  double area() const { return Lx * Ly; }
  bool contains(const V &q) const { return q.x > 0.0 && q.x < Lx && q.y > 0.0 && q.y < Ly; }
  bool onBoundary(const V &q, double tol) const {
    return std::fabs(q.x) < tol || std::fabs(q.x - Lx) < tol || std::fabs(q.y) < tol || std::fabs(q.y - Ly) < tol;
  }
  Polygon polygon() const { return {V(0.0, 0.0), V(Lx, 0.0), V(Lx, Ly), V(0.0, Ly)}; }
};

// Carré de la distance de Mahalanobis d² = u^T Sigma^-1 u, avec Sigma = K a a^T + b b^T (a : direction des
// colonnes à anisotropyAngle de l'axe x, b perpendiculaire). Toute matrice définie positive s'écrit ainsi à un
// facteur près, et seule la forme de Sigma compte.
struct Mahalanobis {
  double ixx{1.0}, ixy{0.0}, iyy{1.0}; // Sigma^-1
  explicit Mahalanobis(const Params &p);
  double d2(const V &u) const { return ixx * u.x * u.x + 2.0 * ixy * u.x * u.y + iyy * u.y * u.y; }
};

// Pavage à sommets partagés
struct Mesh {
  std::vector<V> verts;
  std::vector<std::vector<int>> cells; // sommets de chaque cellule, sens trigonométrique (indice = germe)
  std::vector<char> boundary;          // cellule touchant le bord du domaine
};

class Tissue {
public:
  Params P;
  std::vector<V> seeds;
  std::vector<char> isBundle;
  std::vector<Polygon> cells; // Voronoï des germes courants, découpé au domaine (indice = germe)
  int iteration{0};
  double lastMove{0.0}; // déplacement moyen des germes à la dernière itération, en fraction de spacing()

  // Germes et bundles initiaux à partir de p, puis pavage ; iteration = 0
  void init(const Params &p);
  // Reprend un état enregistré (fichier de cellules)
  void restore(const Params &p, std::vector<V> s, std::vector<char> bundle, int it);
  // Une itération de Lloyd modifiée (avec les paramètres courants de P)
  void step();
  // Itère jusqu'à P.lloydIterations
  void run();
  // Le germe le plus proche de q change de statut (bundle <-> ordinaire) ; renvoie son indice
  int toggleBundle(const V &q);
  int nearestSeed(const V &q) const;

  // Espacement visé entre germes : côté du carré de même aire qu'une cellule, sqrt(pi/4)
  static double spacing();
  size_t nBundles() const;
  Mesh mesh() const;

private:
  void updateVoronoi();
};

// Métriques morphologiques de l'article, calculées sur le pavage (toutes les cellules, bords compris) :
//   <SI>        indice de forme moyen, SI = 4 pi A / P² (compacité ; hexagone régulier : 0.907)
//   sigma_A/<A> coefficient de variation des aires (hétérogénéité)
//   <d_NN>      distance moyenne au plus proche voisin (centroïdes), rapportée à l'espacement visé
//   <d_M>       distance de Mahalanobis moyenne d'une cellule ordinaire au bundle le plus proche, rapportée à
//               l'espacement visé
struct Metrics {
  size_t nCells{0}, nBundles{0};
  double meanSI{0.0}, cvA{0.0}, dNN{0.0}, dM{0.0};
  double areaRegular{0.0}, areaBundle{0.0}, SIRegular{0.0}, SIBundle{0.0};
  double meanSides{0.0};      // cellules intérieures
  std::vector<int> sidesInner; // histogramme du nombre de côtés des cellules intérieures (indice = côtés)
  double lastMove{0.0};

  // Métriques scalaires nommées (affichage et courbes de l'interface, fichier morphology)
  std::vector<std::pair<std::string, double>> items() const;
};

Metrics computeMetrics(const Tissue &t, const Mesh &m);
