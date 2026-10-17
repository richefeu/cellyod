// Thèmes de rendu, communs au SVG (io.cpp) et à l'interface : noir sur blanc (défaut) ou blanc sur noir.
#pragma once

#include <cstdint>

struct Theme {
  // couleurs 0xRRGGBB
  uint32_t background; // fond (hors du tissu)
  uint32_t cell;       // remplissage des cellules
  uint32_t edge;       // parois
  uint32_t boundary;   // bord du tissu (trait épais)
  uint32_t seed;       // germes
  uint32_t seedRing;   // contour des germes
  uint32_t bundle;     // bundles
  uint32_t text;       // légende
};

inline const Theme &renderTheme(bool dark) {
  static const Theme light{0xffffff, 0xffffff, 0x1a1a1a, 0x000000, 0x9a9aa3, 0x55565f, 0xe0524d, 0x1a1a1a};
  static const Theme darkTheme{0x000000, 0x0f111a, 0xf2f2f5, 0xffffff, 0x8e8f9c, 0xc9cad3, 0xf07470, 0xf2f2f5};
  return dark ? darkTheme : light;
}
