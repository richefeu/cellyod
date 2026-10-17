// cellyod-gui : interface raylib/raygui de cellyod.
//
// Le tissu évolue en direct (itérations de Lloyd animées) ; tous les paramètres du registre paramSpecs() sont
// réglables à la souris dans le panneau de droite. Un paramètre marqué * régénère les germes ; les autres
// s'appliquent dès l'itération suivante.
//
// Souris (vue)  : clic gauche = bundle on/off sur la cellule ; glisser (gauche, droit ou milieu) = déplacer ;
//                 molette = zoom
// Clavier       : Espace lancer/pause, S ou -> une itération, R régénérer, F recadrer, Ctrl+S enregistrer
//
// Usage : cellyod-gui [params.txt | cells.txt]
//         cellyod-gui [fichier] --screenshot image.png    (calcule, enregistre une capture et quitte)
#include "raylib.h"
#include "rlgl.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

#define RINI_IMPLEMENTATION
#include "rini.h"

// Fonte d'interface embarquée (Open Sans, OFL), générée par 'xxd -i guifont.ttf'
#include "guifont_ttf.h"

#include "../src/io.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

// ================================================================================================
// Thème et configuration
// ================================================================================================

static const Color COL_BG = {244, 248, 239, 255};       // fond de la vue
static const Color COL_CELL = {227, 239, 216, 255};     // cellule ordinaire
static const Color COL_BUNDLE = {127, 176, 105, 255};   // bundle
static const Color COL_EDGE = {90, 154, 79, 255};       // parois
static const Color COL_SEED = {90, 154, 79, 255};       // germes
static const Color COL_SEED_B = {47, 93, 39, 255};      // germes des bundles
static const Color COL_PANEL = {255, 255, 255, 255};
static const Color COL_PANEL_EDGE = {221, 225, 230, 255};
static const Color COL_TEXT = {43, 47, 54, 255};
static const Color COL_MUTED = {130, 137, 145, 255};
static const Color COL_ACCENT = {45, 125, 210, 255};
static const Color COL_ERROR = {200, 60, 50, 255};

static const char *CONFIG_FILE = "cellyod-gui.ini";
static int cfgWinW = 1280, cfgWinH = 820;
static int cfgPanelW = 330;
static int cfgFontSize = 16;

static Font uiFont;

static void uiText(const char *s, float x, float y, float size, Color c) {
  DrawTextEx(uiFont, s, {x, y}, size, 0.0f, c);
}

static void applyTheme() {
  GuiSetFont(uiFont);
  GuiSetStyle(DEFAULT, TEXT_SIZE, cfgFontSize);
  GuiSetStyle(DEFAULT, TEXT_SPACING, 0);
  GuiSetStyle(DEFAULT, BORDER_WIDTH, 1);
  GuiSetStyle(DEFAULT, BACKGROUND_COLOR, 0xffffffff);
  GuiSetStyle(DEFAULT, LINE_COLOR, 0xdde1e6ff);
  GuiSetStyle(DEFAULT, BORDER_COLOR_NORMAL, 0xc9ced4ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_NORMAL, 0xffffffff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_NORMAL, 0x2b2f36ff);
  GuiSetStyle(DEFAULT, BORDER_COLOR_FOCUSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_FOCUSED, 0xeaf2fbff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_FOCUSED, 0x1b4a7aff);
  GuiSetStyle(DEFAULT, BORDER_COLOR_PRESSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, BASE_COLOR_PRESSED, 0x2d7dd2ff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_PRESSED, 0xffffffff);
  GuiSetStyle(DEFAULT, BORDER_COLOR_DISABLED, 0xe6e9edff);
  GuiSetStyle(DEFAULT, BASE_COLOR_DISABLED, 0xf2f3f5ff);
  GuiSetStyle(DEFAULT, TEXT_COLOR_DISABLED, 0x9aa0a6ff);
  GuiSetStyle(TOGGLE, GROUP_PADDING, 2);
  GuiSetStyle(DEFAULT, TEXT_ALIGNMENT, TEXT_ALIGN_CENTER);
  GuiSetStyle(LABEL, TEXT_ALIGNMENT, TEXT_ALIGN_LEFT);
  GuiSetStyle(TEXTBOX, TEXT_ALIGNMENT, TEXT_ALIGN_RIGHT);
  GuiSetStyle(SLIDER, SLIDER_PADDING, 1);
}

// ================================================================================================
// Champs de coloration des cellules : pour en ajouter un, une ligne dans COLOR_FIELDS
// ================================================================================================

struct CellData {
  std::vector<double> area, SI;
  std::vector<int> sides;
};

struct ColorField {
  const char *name;
  bool discrete;                                    // nombre de côtés : palette discrète
  double (*value)(const CellData &, size_t c);      // nullptr : aplat (cellules / bundles)
};

static const ColorField COLOR_FIELDS[] = {
    {"Uni", false, nullptr},
    {"Côtés", true, [](const CellData &d, size_t c) { return (double)d.sides[c]; }},
    {"Aire", false, [](const CellData &d, size_t c) { return d.area[c]; }},
    {"SI", false, [](const CellData &d, size_t c) { return d.SI[c]; }},
};
static const int N_COLOR_FIELDS = (int)(sizeof(COLOR_FIELDS) / sizeof(COLOR_FIELDS[0]));

// Palette continue (approximation de viridis, éclaircie pour garder les parois lisibles)
static Color rampColor(double t) {
  static const Color stops[] = {{68, 1, 84, 255}, {59, 82, 139, 255}, {33, 145, 140, 255}, {94, 201, 98, 255},
                                {253, 231, 37, 255}};
  t = std::clamp(t, 0.0, 1.0) * 4.0;
  int i = std::min(3, (int)t);
  double f = t - i;
  Color a = stops[i], b = stops[i + 1];
  return {(unsigned char)(a.r + (b.r - a.r) * f), (unsigned char)(a.g + (b.g - a.g) * f),
          (unsigned char)(a.b + (b.b - a.b) * f), 255};
}

static Color sidesColor(int n) {
  static const Color pal[] = {{142, 68, 173, 255},  {214, 69, 65, 255},  {240, 150, 60, 255},
                              {205, 228, 190, 255}, {95, 170, 225, 255}, {45, 110, 180, 255},
                              {25, 60, 110, 255}};
  return pal[std::clamp(n, 3, 9) - 3];
}

// ================================================================================================
// État de l'application
// ================================================================================================

struct App {
  Params P; // paramètres édités
  Tissue T;
  CellMesh M;
  Metrics MT;
  CellData D;
  double vmin{0.0}, vmax{1.0}; // plage du champ de coloration continu
  bool dirty{true};

  // simulation
  bool playing{false};
  bool stopAtTarget{true};
  int speed{1}; // itérations par image

  // historique des métriques (courbe)
  std::vector<int> histIt;
  std::vector<std::vector<double>> hist;
  int plotMetric{3};

  // vue
  Rectangle view{};
  double zoom{10.0};
  V center;
  bool showSeeds{true}, showEdges{true};
  int colorField{0};
  bool dragging{false};
  Vector2 pressPos{};

  // panneau
  Vector2 scroll{};
  float contentH{1000.0f};
  int editing{-1}; // paramètre en cours de saisie (-1 : aucun ; -2 : nom de fichier)
  char editBuf[128]{};
  char file[256]{"cells.txt"};
  std::string message;
  bool messageError{false};
  std::string hoverHelp;
};

static App app;

static void say(const std::string &m, bool error = false) {
  app.message = m;
  app.messageError = error;
}

// CellMesh, métriques et données de coloration après une modification du tissu
static void refresh() {
  App &a = app;
  a.M = a.T.mesh();
  a.MT = computeMetrics(a.T, a.M);
  size_t nc = a.M.cells.size();
  a.D.area.assign(nc, 0.0);
  a.D.SI.assign(nc, 0.0);
  a.D.sides.assign(nc, 0);
  for (size_t c = 0; c < nc; c++) {
    Polygon p;
    for (int v : a.M.cells[c]) {
      p.push_back(a.M.verts[v]);
    }
    a.D.sides[c] = (int)p.size();
    if (p.size() >= 3) {
      double A = std::fabs(signedArea(p)), per = perimeter(p);
      a.D.area[c] = A;
      a.D.SI[c] = 4.0 * M_PI * A / (per * per);
    }
  }
  const ColorField &f = COLOR_FIELDS[a.colorField];
  if (f.value && !f.discrete) {
    a.vmin = 1e300;
    a.vmax = -1e300;
    for (size_t c = 0; c < nc; c++) {
      double v = f.value(a.D, c);
      a.vmin = std::min(a.vmin, v);
      a.vmax = std::max(a.vmax, v);
    }
    if (a.vmax <= a.vmin) {
      a.vmax = a.vmin + 1e-12;
    }
  }
  std::vector<double> vals;
  for (const auto &it : a.MT.items()) {
    vals.push_back(it.second);
  }
  if (!a.histIt.empty() && a.T.iteration < a.histIt.back()) {
    a.histIt.clear();
    a.hist.clear();
  }
  if (!a.histIt.empty() && a.histIt.back() == a.T.iteration) {
    a.hist.back() = vals;
  } else {
    a.histIt.push_back(a.T.iteration);
    a.hist.push_back(vals);
  }
  a.dirty = false;
}

static void fitView() {
  App &a = app;
  a.zoom = 0.95 * std::min(a.view.width / a.T.P.Lx, a.view.height / a.T.P.Ly);
  a.center = V(0.5 * a.T.P.Lx, 0.5 * a.T.P.Ly);
}

// Itère jusqu'à lloydIterations sans animation, en gardant les métriques de chaque itération (courbe)
static void runToTarget() {
  if (app.dirty) {
    refresh();
  }
  while (app.T.iteration < app.T.P.lloydIterations) {
    app.T.step();
    refresh();
  }
  app.playing = false;
}

static void regenerate() {
  App &a = app;
  bool resize = a.T.seeds.empty() || a.T.P.Lx != a.P.Lx || a.T.P.Ly != a.P.Ly;
  a.T.init(a.P);
  a.histIt.clear();
  a.hist.clear();
  a.dirty = true;
  if (resize) {
    fitView();
  }
}

// Paramètres modifiés dans le panneau : régénère les germes si un paramètre « reseed » a changé, sinon les
// nouvelles valeurs s'appliquent à l'itération suivante
static void applyParams() {
  App &a = app;
  bool reseed = false;
  for (const auto &s : paramSpecs()) {
    reseed = reseed || (s.reseed && s.get(a.P) != s.get(a.T.P));
  }
  if (reseed) {
    regenerate();
    a.playing = a.P.lloydIterations > 0;
  } else {
    a.T.P = a.P;
    if (a.stopAtTarget && a.T.iteration < a.P.lloydIterations) {
      a.playing = true;
    }
  }
}

static void loadFile(const std::string &name) {
  App &a = app;
  try {
    try {
      a.T = loadCells(name);
    } catch (const std::exception &) {
      // pas un fichier de cellules : fichier de paramètres ?
      a.P = readParams(name);
      regenerate();
      a.playing = a.P.lloydIterations > 0;
      fitView();
      say("Paramètres lus : " + name);
      return;
    }
    a.P = a.T.P;
    a.histIt.clear();
    a.hist.clear();
    a.dirty = true;
    a.playing = false;
    fitView();
    std::snprintf(a.file, sizeof(a.file), "%s", name.c_str());
    say("Cellules lues : " + name + " (itération " + std::to_string(a.T.iteration) + ")");
  } catch (const std::exception &e) {
    say(e.what(), true);
  }
}

static std::string stem(const std::string &f) {
  std::filesystem::path p(f);
  return (p.parent_path() / p.stem()).string();
}

static void saveAll() {
  App &a = app;
  try {
    std::string f = a.file;
    if (f.empty()) {
      f = "cells.txt";
    }
    if (a.dirty) {
      refresh();
    }
    saveCells(f, a.T, a.M);
    writeSVG(stem(f) + ".svg", a.T, a.M);
    say("Enregistré : " + f + " et " + stem(f) + ".svg");
  } catch (const std::exception &e) {
    say(e.what(), true);
  }
}

static void saveParams() {
  App &a = app;
  std::string f = stem(a.file) + "_params.txt";
  FILE *fp = std::fopen(f.c_str(), "w");
  if (!fp) {
    say("impossible d'écrire " + f, true);
    return;
  }
  std::ostringstream os;
  os << "# cellyod : paramètres enregistrés par cellyod-gui\n";
  writeParams(os, a.T.P);
  std::fputs(os.str().c_str(), fp);
  std::fclose(fp);
  say("Paramètres enregistrés : " + f);
}

// ================================================================================================
// Vue
// ================================================================================================

static Vector2 toScreen(const V &q) {
  const App &a = app;
  return {(float)(a.view.x + 0.5 * a.view.width + (q.x - a.center.x) * a.zoom),
          (float)(a.view.y + 0.5 * a.view.height - (q.y - a.center.y) * a.zoom)};
}

static V toWorld(Vector2 s) {
  const App &a = app;
  return V(a.center.x + (s.x - a.view.x - 0.5 * a.view.width) / a.zoom,
           a.center.y - (s.y - a.view.y - 0.5 * a.view.height) / a.zoom);
}

static Color cellColor(size_t c) {
  const App &a = app;
  const ColorField &f = COLOR_FIELDS[a.colorField];
  if (!f.value) {
    return a.T.isBundle[c] ? COL_BUNDLE : COL_CELL;
  }
  if (f.discrete) {
    Color col = sidesColor((int)f.value(a.D, c));
    if (a.M.boundary[c]) {
      col = ColorAlpha(col, 0.35f); // cellules de bord : nombre de côtés non significatif
    }
    return col;
  }
  return rampColor((f.value(a.D, c) - a.vmin) / (a.vmax - a.vmin));
}

static void drawTissue() {
  App &a = app;
  BeginScissorMode((int)a.view.x, (int)a.view.y, (int)a.view.width, (int)a.view.height);
  DrawRectangleRec(a.view, COL_BG);
  rlDisableBackfaceCulling();
  std::vector<Vector2> pts;
  for (size_t c = 0; c < a.M.cells.size(); c++) {
    const auto &cell = a.M.cells[c];
    if (cell.size() < 3) {
      continue;
    }
    pts.clear();
    for (int v : cell) {
      pts.push_back(toScreen(a.M.verts[v]));
    }
    DrawTriangleFan(pts.data(), (int)pts.size(), cellColor(c));
  }
  rlEnableBackfaceCulling();
  if (a.showEdges) {
    float th = (float)std::clamp(0.03 * a.zoom, 1.0, 3.0);
    Color edge = a.colorField == 0 ? COL_EDGE : Color{40, 44, 52, 150};
    for (const auto &cell : a.M.cells) {
      for (size_t k = 0; k < cell.size(); k++) {
        int u = cell[k], w = cell[(k + 1) % cell.size()];
        if (u < w || cell.size() < 3) { // chaque paroi intérieure est partagée : tracée une fois
          DrawLineEx(toScreen(a.M.verts[u]), toScreen(a.M.verts[w]), th, edge);
        }
      }
    }
    // les parois de bord n'apparaissent qu'une fois (sens quelconque) : contour du domaine
    Vector2 p0 = toScreen(V(0, 0)), p1 = toScreen(V(a.T.P.Lx, a.T.P.Ly));
    DrawRectangleLinesEx({p0.x, p1.y, p1.x - p0.x, p0.y - p1.y}, th, edge);
  }
  if (a.showSeeds) {
    float r = (float)std::clamp(0.06 * a.zoom, 1.0, 5.0);
    for (size_t i = 0; i < a.T.seeds.size(); i++) {
      bool b = a.T.isBundle[i];
      DrawCircleV(toScreen(a.T.seeds[i]), b ? 1.6f * r : r, b ? COL_SEED_B : COL_SEED);
    }
  } else if (a.colorField != 0) {
    // en mode coloré, les bundles restent repérés
    float r = (float)std::clamp(0.1 * a.zoom, 2.0, 6.0);
    for (size_t i = 0; i < a.T.seeds.size(); i++) {
      if (a.T.isBundle[i]) {
        DrawCircleV(toScreen(a.T.seeds[i]), r, COL_SEED_B);
      }
    }
  }

  // légende
  const ColorField &f = COLOR_FIELDS[a.colorField];
  float lx = a.view.x + 12, ly = a.view.y + a.view.height - 34;
  if (f.value && f.discrete) {
    DrawRectangle((int)lx - 6, (int)ly - 6, 7 * 34 + 12, 32, ColorAlpha(WHITE, 0.85f));
    for (int n = 3; n <= 9; n++) {
      float x = lx + (n - 3) * 34;
      DrawRectangle((int)x, (int)ly, 30, 20, sidesColor(n));
      uiText(TextFormat(n == 9 ? "%d+" : "%d", n), x + 8, ly + 1, 16, n >= 7 || n <= 4 ? WHITE : COL_TEXT);
    }
  } else if (f.value) {
    DrawRectangle((int)lx - 6, (int)ly - 6, 252, 32, ColorAlpha(WHITE, 0.85f));
    for (int k = 0; k < 120; k++) {
      DrawRectangle((int)lx + 60 + k, (int)ly, 1, 20, rampColor(k / 119.0));
    }
    uiText(TextFormat("%.3g", a.vmin), lx, ly + 1, 16, COL_TEXT);
    uiText(TextFormat("%.3g", a.vmax), lx + 186, ly + 1, 16, COL_TEXT);
  }
  EndScissorMode();
}

// Survol de la vue : informations sur la cellule sous la souris
static std::string hoverInfo() {
  App &a = app;
  Vector2 m = GetMousePosition();
  if (!CheckCollisionPointRec(m, a.view) || a.T.seeds.empty()) {
    return "";
  }
  V w = toWorld(m);
  if (w.x < 0 || w.y < 0 || w.x > a.T.P.Lx || w.y > a.T.P.Ly) {
    return TextFormat("x = %.3f   y = %.3f", w.x, w.y);
  }
  int c = a.T.nearestSeed(w); // cellule de Voronoï contenant w
  return TextFormat("x = %.3f   y = %.3f   |   cellule %d%s : %d côtés, aire %.4g, SI %.3f", w.x, w.y, c,
                    a.T.isBundle[c] ? " (bundle)" : "", a.D.sides[c], a.D.area[c], a.D.SI[c]);
}

static void handleViewInput() {
  App &a = app;
  Vector2 m = GetMousePosition();
  bool inView = CheckCollisionPointRec(m, a.view);
  float wheel = GetMouseWheelMove();
  if (inView && wheel != 0.0f) {
    V before = toWorld(m);
    a.zoom *= std::pow(1.15, wheel);
    a.zoom = std::clamp(a.zoom, 0.2, 5000.0);
    V after = toWorld(m);
    a.center = a.center + (before - after);
  }
  bool anyPress = IsMouseButtonPressed(MOUSE_BUTTON_LEFT) || IsMouseButtonPressed(MOUSE_BUTTON_RIGHT) ||
                  IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE);
  if (inView && anyPress) {
    a.dragging = true;
    a.pressPos = m;
  }
  if (a.dragging) {
    Vector2 d = GetMouseDelta();
    a.center = a.center - V(d.x / a.zoom, -d.y / a.zoom);
    if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
      // clic sans déplacement : bundle on/off
      if (std::fabs(m.x - a.pressPos.x) + std::fabs(m.y - a.pressPos.y) < 4.0f) {
        V w = toWorld(m);
        if (w.x > 0 && w.y > 0 && w.x < a.T.P.Lx && w.y < a.T.P.Ly) {
          int i = a.T.toggleBundle(w);
          // bundles imposés (bundle x y) : gardés à la régénération et enregistrés avec les paramètres
          auto &B = a.P.bundleAt;
          B.erase(std::remove_if(B.begin(), B.end(), [&](const V &b) { return a.T.nearestSeed(b) == i; }),
                  B.end());
          if (a.T.isBundle[i]) {
            B.push_back(a.T.seeds[i]);
          }
          a.T.P.bundleAt = B;
          a.dirty = true;
        }
      }
    }
    if (!IsMouseButtonDown(MOUSE_BUTTON_LEFT) && !IsMouseButtonDown(MOUSE_BUTTON_RIGHT) &&
        !IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)) {
      a.dragging = false;
    }
  }
}

// ================================================================================================
// Panneau
// ================================================================================================

struct Cursor {
  float x, y, w;
  Rectangle row(float h, float gap = 4.0f) {
    Rectangle r{x, y, w, h};
    y += h + gap;
    return r;
  }
};

static void section(Cursor &c, const char *title) {
  c.y += 8;
  uiText(title, c.x, c.y, (float)cfgFontSize, COL_ACCENT);
  c.y += cfgFontSize + 2;
  DrawLineEx({c.x, c.y}, {c.x + c.w, c.y}, 1.0f, COL_PANEL_EDGE);
  c.y += 6;
}

static void muted(Cursor &c, const char *text) {
  uiText(text, c.x, c.y, (float)cfgFontSize - 2, COL_MUTED);
  c.y += cfgFontSize + 2;
}

// Arrondit une valeur issue d'un curseur (fichiers de paramètres lisibles)
static double niceValue(const ParamSpec &s, double v) {
  if (s.kind != ParamKind::Real) {
    return std::round(v);
  }
  if (s.logScale) {
    if (v <= 0.0) {
      return v;
    }
    double e = std::pow(10.0, std::floor(std::log10(v)) - 2);
    return std::round(v / e) * e;
  }
  double step = (s.hi - s.lo) / 1000.0;
  double e = std::pow(10.0, std::floor(std::log10(step)));
  return std::round(v / e) * e;
}

static void hoverHelp(Rectangle r, const char *help) {
  if (CheckCollisionPointRec(GetMousePosition(), r)) {
    app.hoverHelp = help;
  }
}

// Ligne d'un paramètre : nom, valeur éditable au clavier, curseur ou choix
static bool paramWidget(Cursor &c, int idx, const ParamSpec &s) {
  App &a = app;
  bool changed = false;
  float h = (float)cfgFontSize + 6;
  Rectangle line = c.row(h, 2);
  hoverHelp({line.x, line.y, line.width, line.height + h}, s.help);
  GuiLabel({line.x, line.y, line.width - 100, h}, TextFormat("%s%s", s.name, s.reseed ? " *" : ""));
  if (s.kind != ParamKind::Choice) {
    Rectangle box{line.x + line.width - 96, line.y, 96, h};
    bool edit = a.editing == idx;
    if (!edit) {
      double v = s.get(a.P);
      std::snprintf(a.editBuf, sizeof(a.editBuf), "%s",
                    s.kind == ParamKind::Real ? TextFormat("%.5g", v) : s.text(a.P).c_str());
      static char tmp[128];
      std::snprintf(tmp, sizeof(tmp), "%s", a.editBuf);
      if (GuiTextBox(box, tmp, sizeof(tmp), false)) {
        a.editing = idx;
        std::snprintf(a.editBuf, sizeof(a.editBuf), "%s", s.text(a.P).c_str());
      }
    } else if (GuiTextBox(box, a.editBuf, sizeof(a.editBuf), true)) {
      try {
        s.parse(a.P, a.editBuf);
        changed = true;
      } catch (const std::exception &e) {
        say(e.what(), true);
      }
      a.editing = -1;
    }
    Rectangle bar = c.row(12, 8);
    double v = s.get(a.P);
    float t, t0;
    if (s.logScale) {
      t = (float)(std::log(std::max(v, s.lo) / s.lo) / std::log(s.hi / s.lo));
    } else {
      t = (float)((v - s.lo) / (s.hi - s.lo));
    }
    t = std::clamp(t, 0.0f, 1.0f);
    t0 = t;
    GuiSliderBar(bar, nullptr, nullptr, &t, 0.0f, 1.0f);
    if (t != t0) {
      double nv = s.logScale ? s.lo * std::pow(s.hi / s.lo, (double)t) : s.lo + (s.hi - s.lo) * (double)t;
      s.set(a.P, niceValue(s, nv));
      changed = s.get(a.P) != v;
    }
  } else {
    std::string items;
    for (const auto &ch : *s.choices) {
      items += (items.empty() ? "" : ";") + ch;
    }
    Rectangle row = c.row(h + 2, 8);
    int n = (int)s.choices->size();
    float w = (row.width - (n - 1) * 2) / n;
    int active = a.P.*(s.integer), before = active;
    GuiToggleGroup({row.x, row.y, w, row.height}, items.c_str(), &active);
    if (active != before) {
      a.P.*(s.integer) = active;
      changed = true;
    }
  }
  return changed;
}

static void drawPlot(Rectangle r) {
  App &a = app;
  DrawRectangleRec(r, {250, 251, 252, 255});
  DrawRectangleLinesEx(r, 1.0f, COL_PANEL_EDGE);
  if (a.hist.size() < 2) {
    uiText("(courbe : au moins deux itérations)", r.x + 6, r.y + r.height / 2 - 8, 14, COL_MUTED);
    return;
  }
  int k = a.plotMetric;
  double lo = 1e300, hi = -1e300;
  for (const auto &h : a.hist) {
    lo = std::min(lo, h[k]);
    hi = std::max(hi, h[k]);
  }
  if (hi - lo < 1e-12 * std::max(1.0, std::fabs(hi))) {
    hi = lo + 1.0;
    lo -= 1.0;
  }
  double i0 = a.histIt.front(), i1 = std::max((double)a.histIt.back(), i0 + 1);
  auto P = [&](size_t j) {
    return Vector2{(float)(r.x + 4 + (r.width - 8) * (a.histIt[j] - i0) / (i1 - i0)),
                   (float)(r.y + r.height - 4 - (r.height - 8) * (a.hist[j][k] - lo) / (hi - lo))};
  };
  for (size_t j = 1; j < a.hist.size(); j++) {
    DrawLineEx(P(j - 1), P(j), 1.5f, COL_ACCENT);
  }
  uiText(TextFormat("%.4g", hi), r.x + 4, r.y + 2, 13, COL_MUTED);
  uiText(TextFormat("%.4g", lo), r.x + 4, r.y + r.height - 16, 13, COL_MUTED);
  const char *it = TextFormat("it %d", a.histIt.back());
  uiText(it, r.x + r.width - 4 - MeasureTextEx(uiFont, it, 13, 0).x, r.y + r.height - 16, 13, COL_MUTED);
}

static void drawSidesHistogram(Rectangle r) {
  App &a = app;
  const auto &h = a.MT.sidesInner;
  int total = 0, top = 1;
  for (int n : h) {
    total += n;
    top = std::max(top, n);
  }
  if (!total) {
    return;
  }
  float bw = r.width / 8.0f;
  for (int n = 3; n <= 10; n++) {
    int v = n < (int)h.size() ? h[n] : 0;
    float bh = (r.height - 34) * v / (float)top;
    float x = r.x + (n - 3) * bw;
    DrawRectangle((int)x + 2, (int)(r.y + r.height - 18 - bh), (int)bw - 4, (int)bh, sidesColor(n));
    uiText(TextFormat("%d", n), x + bw / 2 - 4, r.y + r.height - 16, 13, COL_MUTED);
    if (100.0 * v / total >= 0.5) {
      uiText(TextFormat("%.0f%%", 100.0 * v / total), x + 4, r.y + r.height - 18 - bh - 14, 12, COL_TEXT);
    }
  }
}

static void drawPanel(float W, float H) {
  App &a = app;
  float pw = (float)cfgPanelW;
  Rectangle bounds{W - pw, 0, pw, H};
  DrawRectangleRec(bounds, COL_PANEL);
  DrawLineEx({bounds.x, 0}, {bounds.x, H}, 1.0f, COL_PANEL_EDGE);
  Rectangle content{bounds.x, 0, pw - 14, a.contentH};
  Rectangle sv{};
  GuiScrollPanel(bounds, nullptr, content, &a.scroll, &sv);
  BeginScissorMode((int)sv.x, (int)sv.y, (int)sv.width, (int)sv.height);
  Cursor c{sv.x + 12, sv.y + a.scroll.y + 4, sv.width - 24};
  float y0 = c.y;
  float bh = (float)cfgFontSize + 8;

  // --- Simulation
  section(c, "SIMULATION");
  {
    Rectangle r = c.row(bh);
    float w = (r.width - 8) / 3;
    if (GuiButton({r.x, r.y, w, bh}, a.playing ? "Pause" : "Lancer")) {
      a.playing = !a.playing;
    }
    if (GuiButton({r.x + w + 4, r.y, w, bh}, "Pas")) {
      a.playing = false;
      a.T.step();
      a.dirty = true;
    }
    if (GuiButton({r.x + 2 * (w + 4), r.y, w, bh}, "Régénérer")) {
      regenerate();
      a.playing = a.P.lloydIterations > 0;
    }
    Rectangle r2 = c.row(bh);
    if (GuiButton(r2, "Calcul direct jusqu'à lloydIterations")) {
      regenerate();
      runToTarget();
    }
    uiText(TextFormat("Itération %d / %d     déplacement %.2e", a.T.iteration, a.P.lloydIterations,
                      a.T.lastMove),
           c.x, c.y, (float)cfgFontSize, COL_TEXT);
    c.y += cfgFontSize + 6;
    Rectangle cb = c.row(16, 6);
    GuiCheckBox({cb.x, cb.y, 16, 16}, "Arrêt à lloydIterations", &a.stopAtTarget);
    Rectangle sl = c.row(cfgFontSize + 4, 2);
    GuiLabel({sl.x, sl.y, sl.width - 60, sl.height}, "Itérations par image");
    uiText(TextFormat("%d", a.speed), sl.x + sl.width - 30, sl.y + 2, (float)cfgFontSize, COL_TEXT);
    float sp = (float)a.speed;
    GuiSliderBar(c.row(12, 6), nullptr, nullptr, &sp, 1.0f, 50.0f);
    a.speed = std::max(1, (int)std::lround(sp));
  }

  // --- Paramètres (registre)
  section(c, "PARAMÈTRES");
  muted(c, "* : régénère les germes ; sinon effet immédiat");
  {
    const char *group = "";
    bool changed = false;
    const auto &specs = paramSpecs();
    for (size_t i = 0; i < specs.size(); i++) {
      if (std::strcmp(group, specs[i].group) != 0) {
        group = specs[i].group;
        c.y += 2;
        uiText(group, c.x, c.y, (float)cfgFontSize - 1, COL_MUTED);
        c.y += cfgFontSize + 2;
      }
      changed = paramWidget(c, (int)i, specs[i]) || changed;
    }
    muted(c, TextFormat("Clic sur une cellule : bundle on/off (%d imposés)", (int)a.P.bundleAt.size()));
    if (!a.P.bundleAt.empty() && GuiButton(c.row(bh), "Oublier les bundles imposés")) {
      a.P.bundleAt.clear();
      a.T.P.bundleAt.clear();
      say("Bundles imposés oubliés (effet à la prochaine régénération)");
    }
    if (changed) {
      applyParams();
    }
  }

  // --- Affichage
  section(c, "AFFICHAGE");
  {
    std::string items;
    for (int i = 0; i < N_COLOR_FIELDS; i++) {
      items += std::string(i ? ";" : "") + COLOR_FIELDS[i].name;
    }
    Rectangle r = c.row(bh, 8);
    float w = (r.width - (N_COLOR_FIELDS - 1) * 2) / N_COLOR_FIELDS;
    int before = a.colorField;
    GuiToggleGroup({r.x, r.y, w, bh}, items.c_str(), &a.colorField);
    if (before != a.colorField) {
      a.dirty = true;
    }
    Rectangle cb = c.row(16, 8);
    GuiCheckBox({cb.x, cb.y, 16, 16}, "Germes", &a.showSeeds);
    GuiCheckBox({cb.x + cb.width / 2, cb.y, 16, 16}, "Parois", &a.showEdges);
    if (GuiButton(c.row(bh), "Recadrer (F)")) {
      fitView();
    }
  }

  // --- Métriques
  section(c, "MÉTRIQUES");
  {
    const auto items = a.MT.items();
    std::string names;
    for (size_t i = 0; i < items.size(); i++) {
      const char *v = TextFormat("%.4g", items[i].second);
      bool sel = (int)i == a.plotMetric;
      Rectangle r = c.row((float)cfgFontSize + 2, 0);
      if (CheckCollisionPointRec(GetMousePosition(), r) && IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
          CheckCollisionPointRec(GetMousePosition(), sv)) {
        a.plotMetric = (int)i;
      }
      if (sel) {
        DrawRectangleRec(r, {234, 242, 251, 255});
      }
      uiText(items[i].first.c_str(), r.x + 4, r.y, (float)cfgFontSize, sel ? COL_ACCENT : COL_TEXT);
      uiText(v, r.x + r.width - 4 - MeasureTextEx(uiFont, v, (float)cfgFontSize, 0).x, r.y, (float)cfgFontSize,
             COL_TEXT);
    }
    c.y += 4;
    muted(c, "Clic sur une métrique : courbe en fonction de l'itération");
    drawPlot(c.row(120, 8));
    muted(c, "Côtés des cellules intérieures");
    drawSidesHistogram(c.row(90, 8));
  }

  // --- Fichiers
  section(c, "FICHIERS");
  {
    Rectangle r = c.row(bh);
    GuiSetStyle(TEXTBOX, TEXT_ALIGNMENT, TEXT_ALIGN_LEFT);
    if (GuiTextBox(r, a.file, sizeof(a.file), a.editing == -2)) {
      a.editing = a.editing == -2 ? -1 : -2;
    }
    GuiSetStyle(TEXTBOX, TEXT_ALIGNMENT, TEXT_ALIGN_RIGHT);
    Rectangle b = c.row(bh);
    float w = (b.width - 8) / 3;
    if (GuiButton({b.x, b.y, w, bh}, "Enregistrer")) {
      saveAll();
    }
    if (GuiButton({b.x + w + 4, b.y, w, bh}, "Charger")) {
      loadFile(a.file);
    }
    if (GuiButton({b.x + 2 * (w + 4), b.y, w, bh}, "Paramètres")) {
      saveParams();
    }
    muted(c, "Enregistrer : cellules + SVG ; Paramètres : _params.txt");
    muted(c, "Glisser-déposer : fichier de cellules ou de paramètres");
  }

  // --- Aide
  section(c, "CLAVIER");
  muted(c, "Espace lancer/pause   S ou -> un pas   R régénérer");
  muted(c, "F recadrer   Ctrl+S enregistrer");
  c.y += 8;

  a.contentH = c.y - y0;
  EndScissorMode();
}

// ================================================================================================
// Configuration persistante
// ================================================================================================

static void loadConfig() {
  if (!FileExists(CONFIG_FILE)) {
    return;
  }
  rini_data ini = rini_load(CONFIG_FILE);
  cfgWinW = rini_get_value_fallback(ini, "window_width", cfgWinW);
  cfgWinH = rini_get_value_fallback(ini, "window_height", cfgWinH);
  cfgPanelW = rini_get_value_fallback(ini, "panel_width", cfgPanelW);
  cfgFontSize = rini_get_value_fallback(ini, "font_size", cfgFontSize);
  app.showSeeds = rini_get_value_fallback(ini, "show_seeds", app.showSeeds) != 0;
  app.showEdges = rini_get_value_fallback(ini, "show_edges", app.showEdges) != 0;
  app.colorField = std::clamp(rini_get_value_fallback(ini, "color_field", app.colorField), 0, N_COLOR_FIELDS - 1);
  app.speed = std::max(1, rini_get_value_fallback(ini, "speed", app.speed));
  app.stopAtTarget = rini_get_value_fallback(ini, "stop_at_target", app.stopAtTarget) != 0;
  rini_unload(&ini);
}

static void saveConfig() {
  rini_data ini = rini_load(nullptr);
  rini_set_value(&ini, "window_width", GetScreenWidth(), "Largeur de la fenêtre (px)");
  rini_set_value(&ini, "window_height", GetScreenHeight(), "Hauteur de la fenêtre (px)");
  rini_set_value(&ini, "panel_width", cfgPanelW, "Largeur du panneau (px)");
  rini_set_value(&ini, "font_size", cfgFontSize, "Taille de police");
  rini_set_value(&ini, "show_seeds", app.showSeeds, "Afficher les germes (0/1)");
  rini_set_value(&ini, "show_edges", app.showEdges, "Afficher les parois (0/1)");
  rini_set_value(&ini, "color_field", app.colorField, "Champ de coloration (0 uni, 1 côtés, 2 aire, 3 SI)");
  rini_set_value(&ini, "speed", app.speed, "Itérations par image");
  rini_set_value(&ini, "stop_at_target", app.stopAtTarget, "Arrêt à lloydIterations (0/1)");
  rini_save(ini, CONFIG_FILE);
  rini_unload(&ini);
}

// ================================================================================================
// Boucle principale
// ================================================================================================

int main(int argc, char *argv[]) {
  std::string input, screenshot;
  for (int i = 1; i < argc; i++) {
    if (!std::strcmp(argv[i], "--screenshot") && i + 1 < argc) {
      screenshot = argv[++i];
    } else {
      input = argv[i];
    }
  }
  loadConfig();

  SetConfigFlags(FLAG_WINDOW_RESIZABLE | FLAG_MSAA_4X_HINT);
  SetTraceLogLevel(LOG_WARNING);
  InitWindow(cfgWinW, cfgWinH, "cellyod");
  SetExitKey(KEY_NULL);
  SetTargetFPS(60);

  // fonte : latin-1 + grec (alpha, sigma)
  std::vector<int> cps;
  for (int c = 32; c < 256; c++) {
    if (c < 127 || c >= 160) {
      cps.push_back(c);
    }
  }
  for (int c : {0x3B1, 0x3C3, 0x3A3}) {
    cps.push_back(c);
  }
  uiFont = LoadFontFromMemory(".ttf", guifont_ttf, (int)guifont_ttf_len, 48, cps.data(), (int)cps.size());
  SetTextureFilter(uiFont.texture, TEXTURE_FILTER_BILINEAR);
  applyTheme();

  auto layoutView = [&]() {
    float W = (float)GetScreenWidth(), H = (float)GetScreenHeight();
    app.view = {0, 0, W - cfgPanelW, H - 28};
  };
  layoutView();
  if (!input.empty()) {
    loadFile(input);
  } else {
    regenerate();
    app.playing = app.P.lloydIterations > 0;
  }
  fitView();
  if (!screenshot.empty()) {
    runToTarget();
  }

  int frame = 0;
  while (!WindowShouldClose()) {
    layoutView();
    app.hoverHelp.clear();

    // fichiers déposés
    if (IsFileDropped()) {
      FilePathList dropped = LoadDroppedFiles();
      if (dropped.count > 0) {
        loadFile(dropped.paths[0]);
      }
      UnloadDroppedFiles(dropped);
    }

    // clavier (hors saisie)
    if (app.editing == -1) {
      bool ctrl = IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL) || IsKeyDown(KEY_LEFT_SUPER) ||
                  IsKeyDown(KEY_RIGHT_SUPER);
      if (IsKeyPressed(KEY_SPACE)) {
        app.playing = !app.playing;
      }
      if ((IsKeyPressed(KEY_S) && !ctrl) || IsKeyPressed(KEY_RIGHT)) {
        app.playing = false;
        app.T.step();
        app.dirty = true;
      }
      if (IsKeyPressed(KEY_S) && ctrl) {
        saveAll();
      }
      if (IsKeyPressed(KEY_R)) {
        regenerate();
        app.playing = app.P.lloydIterations > 0;
      }
      if (IsKeyPressed(KEY_F)) {
        fitView();
      }
    }

    handleViewInput();

    // simulation : au plus 'speed' itérations par image, dans un budget de temps
    if (app.playing) {
      double t0 = GetTime();
      for (int k = 0; k < app.speed && GetTime() - t0 < 0.03; k++) {
        if (app.stopAtTarget && app.T.iteration >= app.P.lloydIterations) {
          app.playing = false;
          break;
        }
        app.T.step();
        app.dirty = true;
      }
    }
    if (app.dirty) {
      refresh();
    }

    BeginDrawing();
    ClearBackground(COL_PANEL);
    float W = (float)GetScreenWidth(), H = (float)GetScreenHeight();
    drawTissue();
    drawPanel(W, H);

    // barre d'état : aide du paramètre survolé, informations de la cellule survolée, ou dernier message
    Rectangle bar{0, H - 28, W - cfgPanelW, 28};
    DrawRectangleRec(bar, {248, 249, 250, 255});
    DrawLineEx({0, bar.y}, {bar.width, bar.y}, 1.0f, COL_PANEL_EDGE);
    std::string info = app.hoverHelp.empty() ? hoverInfo() : app.hoverHelp;
    bool err = false;
    if (info.empty()) {
      info = app.message;
      err = app.messageError;
    }
    if (info.empty()) {
      info = TextFormat("%zu cellules, %zu bundles", app.MT.nCells, app.MT.nBundles);
    }
    uiText(info.c_str(), 10, bar.y + 5, (float)cfgFontSize, err ? COL_ERROR : COL_TEXT);
    EndDrawing();

    if (!screenshot.empty() && ++frame == 3) {
      TakeScreenshot(screenshot.c_str());
      break;
    }
  }

  if (screenshot.empty()) {
    saveConfig();
  }
  UnloadFont(uiFont);
  CloseWindow();
  return 0;
}
