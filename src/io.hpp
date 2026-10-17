// Entrées/sorties : fichier de cellules (format cellyod), SVG, fichier de métriques.
//
// Format de cellules (texte, # pour les commentaires) :
//
//   cellyod 1
//   param <mot-clé> <valeur>        tous les paramètres du tissu (cf. params.hpp)
//   domain xmin ymin xmax ymax
//   iteration <n>                   nombre d'itérations de Lloyd effectuées
//   vertices <N>                    puis N lignes : x y
//   cells <M>                       puis M lignes : type sx sy nv v0 v1 ... v(nv-1)
//
// Les sommets sont partagés entre cellules ; ceux d'une cellule sont donnés dans le sens trigonométrique.
// L'identifiant d'une cellule est son rang (à partir de 0). type : 0 cellule ordinaire, 1 bundle. (sx, sy) est
// le germe de la cellule : avec les paramètres, il permet de reprendre les itérations là où elles s'étaient
// arrêtées (le pavage est entièrement déterminé par les germes).
#pragma once

#include "tissue.hpp"

#include <string>

void saveCells(const std::string &file, const Tissue &t, const Mesh &m);
// Relit un fichier de cellules et reconstruit le tissu (lève std::runtime_error)
Tissue loadCells(const std::string &file);

struct SvgOptions {
  double widthPx{1000.0};
  bool seeds{true};
};
void writeSVG(const std::string &file, const Tissue &t, const Mesh &m, const SvgOptions &opt = {});

// Ajoute une ligne de paramètres et de métriques (en-tête écrit si le fichier est nouveau)
void appendMorphology(const std::string &file, const Tissue &t, const Metrics &mt);
