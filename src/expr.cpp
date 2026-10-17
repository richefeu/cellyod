#include "expr.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

namespace {

const char *FUN1_NAMES[] = {"sin", "cos", "tan", "asin", "acos", "atan", "exp", "log", "sqrt", "abs", "floor", "ceil"};
double (*const FUN1_PTR[])(double) = {::sin, ::cos, ::tan, ::asin, ::acos, ::atan, ::exp, ::log, ::sqrt, ::fabs, ::floor, ::ceil};
const char *FUN2_NAMES[] = {"atan2", "min", "max", "pow"};
double fmin2(double a, double b) { return a < b ? a : b; }
double fmax2(double a, double b) { return a > b ? a : b; }
double (*const FUN2_PTR[])(double, double) = {::atan2, fmin2, fmax2, ::pow};
const int MAX_STACK = 64;

} // namespace

// Analyseur récursif descendant produisant une suite d'opérations en notation polonaise inverse
struct ExprParser {
  const std::string &s;
  const std::vector<std::string> &vars;
  size_t pos{0};
  std::vector<Expr::Op> &out;
  int depth{0}, maxDepth{0};

  [[noreturn]] void fail(const std::string &msg) const {
    throw std::runtime_error("expression « " + s + " » : " + msg + " (position " + std::to_string(pos + 1) + ")");
  }
  void push(Expr::Op op, int delta) {
    out.push_back(op);
    depth += delta;
    maxDepth = depth > maxDepth ? depth : maxDepth;
  }
  void skip() {
    while (pos < s.size() && std::isspace((unsigned char)s[pos])) {
      pos++;
    }
  }
  bool accept(char c) {
    skip();
    if (pos < s.size() && s[pos] == c) {
      pos++;
      return true;
    }
    return false;
  }
  void expect(char c) {
    if (!accept(c)) {
      fail(std::string("« ") + c + " » attendu");
    }
  }

  void expr() {
    term();
    for (;;) {
      if (accept('+')) {
        term();
        push({Expr::ADD}, -1);
      } else if (accept('-')) {
        term();
        push({Expr::SUB}, -1);
      } else {
        return;
      }
    }
  }
  void term() {
    unary();
    for (;;) {
      if (accept('*')) {
        unary();
        push({Expr::MUL}, -1);
      } else if (accept('/')) {
        unary();
        push({Expr::DIV}, -1);
      } else {
        return;
      }
    }
  }
  void unary() {
    if (accept('-')) {
      unary();
      push({Expr::NEG}, 0);
    } else if (accept('+')) {
      unary();
    } else {
      power();
    }
  }
  void power() {
    primary();
    if (accept('^')) {
      unary();
      push({Expr::POW}, -1);
    }
  }
  void primary() {
    skip();
    if (pos >= s.size()) {
      fail("expression incomplète");
    }
    if (accept('(')) {
      expr();
      expect(')');
      return;
    }
    char c = s[pos];
    if (std::isdigit((unsigned char)c) || c == '.') {
      const char *b = s.c_str() + pos;
      char *e = nullptr;
      double v = std::strtod(b, &e);
      pos += (size_t)(e - b);
      push({Expr::NUM, v}, 1);
      return;
    }
    if (!std::isalpha((unsigned char)c) && c != '_') {
      fail(std::string("caractère inattendu « ") + c + " »");
    }
    size_t b = pos;
    while (pos < s.size() && (std::isalnum((unsigned char)s[pos]) || s[pos] == '_')) {
      pos++;
    }
    std::string id = s.substr(b, pos - b);
    if (accept('(')) {
      for (int i = 0; i < (int)(sizeof(FUN1_NAMES) / sizeof(*FUN1_NAMES)); i++) {
        if (id == FUN1_NAMES[i]) {
          expr();
          expect(')');
          push({Expr::FUN1, 0.0, i}, 0);
          return;
        }
      }
      for (int i = 0; i < (int)(sizeof(FUN2_NAMES) / sizeof(*FUN2_NAMES)); i++) {
        if (id == FUN2_NAMES[i]) {
          expr();
          expect(',');
          expr();
          expect(')');
          push({Expr::FUN2, 0.0, i}, -1);
          return;
        }
      }
      fail("fonction inconnue « " + id + " »");
    }
    if (id == "pi") {
      push({Expr::NUM, M_PI}, 1);
      return;
    }
    if (id == "e") {
      push({Expr::NUM, M_E}, 1);
      return;
    }
    for (size_t i = 0; i < vars.size(); i++) {
      if (id == vars[i]) {
        push({Expr::VAR, 0.0, (int)i}, 1);
        return;
      }
    }
    std::string list;
    for (const auto &v : vars) {
      list += (list.empty() ? "" : ", ") + v;
    }
    fail("variable inconnue « " + id + " » (variables : " + list + ")");
  }
};

Expr::Expr(const std::string &text, const std::vector<std::string> &vars) {
  ExprParser p{text, vars, 0, ops};
  p.expr();
  p.skip();
  if (p.pos != text.size()) {
    p.fail("caractère inattendu « " + text.substr(p.pos, 1) + " »");
  }
  if (p.maxDepth > MAX_STACK) {
    p.fail("expression trop imbriquée");
  }
}

double Expr::eval(const double *v) const {
  double st[MAX_STACK];
  int n = 0;
  for (const Op &o : ops) {
    switch (o.code) {
    case NUM: st[n++] = o.num; break;
    case VAR: st[n++] = v[o.index]; break;
    case ADD: n--; st[n - 1] += st[n]; break;
    case SUB: n--; st[n - 1] -= st[n]; break;
    case MUL: n--; st[n - 1] *= st[n]; break;
    case DIV: n--; st[n - 1] /= st[n]; break;
    case POW: n--; st[n - 1] = std::pow(st[n - 1], st[n]); break;
    case NEG: st[n - 1] = -st[n - 1]; break;
    case FUN1: st[n - 1] = FUN1_PTR[o.index](st[n - 1]); break;
    case FUN2: n--; st[n - 1] = FUN2_PTR[o.index](st[n - 1], st[n]); break;
    }
  }
  return st[0];
}
