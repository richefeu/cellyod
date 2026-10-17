// Tissu cellulaire : germes, pavage de Voronoï et itérations de Lloyd modifiées.
//
// Modèle de Trivino et al. (morphogenèse de Voronoï variationnelle) : une fraction des germes sont des bundles
// (faisceaux vasculaires) ; à chaque itération, les autres germes vont en p' = (c + alpha b) / (1 + alpha), c
// centroïde de leur cellule et b bundle le plus proche au sens de la distance de Mahalanobis de matrice Sigma.
// Avec alpha = 0, c'est la régularisation de Lloyd. C'est une descente de E = E_CVL - alpha E_att.
#pragma once

#include "params.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

class Expr;

// Domaine. Les germes remplissent le rectangle [0, Lx] x [0, Ly] et le pavage de Voronoï y est découpé ; le tissu
// est la partie { f <= 0 } de ce rectangle, f étant l'expression shape (sans forme : tout le rectangle). Une
// cellule est découpée par le bord de la forme ; une cellule entièrement extérieure est vide et son germe reste
// immobile : les germes extérieurs forment la paroi qui contient le tissu.
struct Domain {
  double Lx{1.0}, Ly{1.0};
  std::shared_ptr<const Expr> shape;

  Domain() = default;
  explicit Domain(const Params &p); // lève std::runtime_error si l'expression est invalide
  double area() const { return Lx * Ly; } // aire du rectangle des germes
  bool contains(const V &q) const { return q.x > 0.0 && q.x < Lx && q.y > 0.0 && q.y < Ly; }
  Polygon polygon() const { return {V(0.0, 0.0), V(Lx, 0.0), V(Lx, Ly), V(0.0, Ly)}; }
  bool hasShape() const { return (bool)shape; }
  double f(const V &q) const; // < 0 dans la forme
  bool inShape(const V &q) const { return !shape || f(q) <= 0.0; }
  // Partie de la cellule convexe poly (sens trigonométrique) dans la forme, bord courbe échantillonné au pas h ;
  // cut : la cellule est coupée par le bord de la forme
  Polygon clip(const Polygon &poly, double h, bool &cut) const;
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
struct CellMesh {
  std::vector<V> verts;
  std::vector<std::vector<int>> cells;          // sommets de chaque cellule, sens trigonométrique (indice = germe ;
                                                // vide pour un germe extérieur à la forme)
  std::vector<char> boundary;                   // cellule touchant le bord du tissu
  std::vector<std::pair<int, int>> boundaryEdges; // parois du bord du tissu (n'appartenant qu'à une cellule)
};

class Tissue {
public:
  Params P;
  std::vector<V> seeds;
  std::vector<char> isBundle;
  Domain domain;
  std::vector<Polygon> voronoi; // Voronoï des germes courants, découpé au rectangle (indice = germe)
  std::vector<Polygon> cells;   // partie de chaque cellule dans la forme (vide : germe extérieur)
  std::vector<double> areaOut;  // aire de la partie extérieure à la forme
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
  // Germe le plus proche de q parmi ceux dont la cellule n'est pas vide (-1 si aucun)
  int nearestSeed(const V &q) const;
  bool inside(size_t i) const { return cells[i].size() >= 3; }

  // Espacement visé entre germes : côté du carré de même aire qu'une cellule, sqrt(pi/4)
  static double spacing();
  size_t nBundles() const; // bundles dans la forme
  CellMesh mesh() const;

private:
  void updateVoronoi();
};

// Métriques morphologiques de l'article, calculées sur le tissu (toutes les cellules, bords compris) :
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

Metrics computeMetrics(const Tissue &t, const CellMesh &m);
