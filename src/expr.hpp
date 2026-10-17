// Expressions numériques compilées (sans dépendance), pour définir une forme par une fonction f(x, y).
//
// Opérateurs : + - * / ^ (puissance, associative à droite), parenthèses, moins unaire.
// Fonctions  : sin cos tan asin acos atan exp log sqrt abs floor ceil, atan2(y, x) min(a, b) max(a, b) pow(a, b)
// Constantes : pi, e. Les variables sont nommées à la compilation.
#pragma once

#include <string>
#include <vector>

class Expr {
public:
  // Compile text avec les variables vars (lève std::runtime_error avec un message en clair)
  Expr(const std::string &text, const std::vector<std::string> &vars);
  // Évalue avec les valeurs des variables, dans l'ordre donné à la compilation
  double eval(const double *values) const;

  enum Code { NUM, VAR, ADD, SUB, MUL, DIV, POW, NEG, FUN1, FUN2 };
  struct Op { // opération de la notation polonaise inverse
    Code code;
    double num{0.0};
    int index{0}; // variable ou fonction
  };
private:
  friend struct ExprParser;
  std::vector<Op> ops;
};
