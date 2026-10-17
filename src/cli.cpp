// cellyod : génération de tissus cellulaires (pavage de Voronoï, Lloyd, attraction vers des bundles).
//
// Usage : cellyod [params.txt | -load cells.txt] [mot-clé=valeur ...]
//         cellyod -h
#include "io.hpp"

#include <cstring>
#include <iomanip>
#include <iostream>

static void usage() {
  std::cout << "Usage : cellyod [params.txt | -load cells.txt] [mot-clé=valeur ...]\n"
               "\n"
               "  params.txt          fichier de paramètres « mot-clé valeur(s) » (# : commentaire)\n"
               "  -load cells.txt     reprend un fichier de cellules et poursuit les itérations jusqu'à\n"
               "                      lloydIterations (les paramètres sont ceux du fichier)\n"
               "  mot-clé=valeur      remplace un paramètre (ex. alpha=0.05 seed=3 \"anisotropy=1000 90\")\n"
               "\nMots-clés (longueurs en diamètres de cellule) :\n";
  printParamHelp(std::cout);
}

int main(int argc, char *argv[]) {
  try {
    Params P;
    Tissue t;
    bool loaded = false;
    std::vector<std::string> overrides;
    for (int i = 1; i < argc; i++) {
      if (!std::strcmp(argv[i], "-h") || !std::strcmp(argv[i], "--help")) {
        usage();
        return 0;
      } else if (!std::strcmp(argv[i], "-load") && i + 1 < argc) {
        t = loadCells(argv[++i]);
        P = t.P;
        loaded = true;
      } else if (std::strchr(argv[i], '=')) {
        overrides.push_back(argv[i]);
      } else {
        P = readParams(argv[i]);
      }
    }
    for (const auto &o : overrides) {
      auto eq = o.find('=');
      setParam(P, o.substr(0, eq), o.substr(eq + 1));
    }
    if (loaded) {
      t.P = P; // paramètres « vivants » (alpha, anisotropie...) éventuellement modifiés ; les germes sont gardés
    } else {
      t.init(P);
    }
    t.run();

    CellMesh m = t.mesh();
    Metrics mt = computeMetrics(t, m);
    saveCells(P.output, t, m);
    if (!P.svg.empty() && P.svg != "none") {
      SvgOptions so;
      so.seeds = P.svgSeeds;
      so.legend = P.svgLegend;
      so.dark = P.svgDark;
      writeSVG(P.svg, t, m, so);
    }
    if (!P.morphology.empty()) {
      appendMorphology(P.morphology, t, mt);
    }

    std::cout << std::setprecision(4);
    std::cout << "cellyod : " << mt.nCells << " cellules (" << mt.nBundles << " bundles), " << m.verts.size()
              << " sommets, " << t.iteration << " itérations, domaine " << P.Lx << " x " << P.Ly << "\n";
    std::cout << "  morphologie : <SI> = " << mt.meanSI << ", sigma_A/<A> = " << mt.cvA
              << ", <d_NN>/espacement = " << mt.dNN;
    if (mt.nBundles) {
      std::cout << ", <d_M>/espacement = " << mt.dM;
    }
    std::cout << "\n  côtés (cellules intérieures, moyenne " << mt.meanSides << ") :";
    size_t nInner = 0;
    for (int n : mt.sidesInner) {
      nInner += (size_t)n;
    }
    for (size_t k = 0; k < mt.sidesInner.size(); k++) {
      if (mt.sidesInner[k]) {
        std::cout << "  " << k << " : " << std::setprecision(3) << 100.0 * mt.sidesInner[k] / (double)nInner
                  << " %";
      }
    }
    std::cout << "\n  -> " << P.output << (P.svg.empty() || P.svg == "none" ? "" : ", " + P.svg)
              << (P.morphology.empty() ? "" : ", " + P.morphology) << "\n";
  } catch (const std::exception &e) {
    std::cerr << "cellyod : " << e.what() << "\n";
    return 1;
  }
  return 0;
}
