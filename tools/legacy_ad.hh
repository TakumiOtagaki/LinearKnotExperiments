// Historical BPP diagnostic only; not used in E0 PF parity reports.
// Reverse differentiation of the oracle's own PF arithmetic, in an external
// observation build. Numeric branch decisions and local weights are unchanged.
#pragma once
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

namespace oracle_ad {
constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
struct Node { std::size_t a, b; double da, db, adj; };
struct Tag { std::size_t node; int i, j; double value; };
inline std::vector<Node> tape;
inline std::vector<Tag> tags;
struct Real {
  double value = 0;
  std::size_t id = none;
  Real() = default;
  Real(double v) : value(v) {}
  operator double() const { return value; }
  Real &operator+=(Real b);
  Real &operator-=(Real b);
  Real &operator*=(Real b);
  Real &operator/=(Real b);
};
inline Real operation(double value, Real a, Real b, double da, double db) {
  Real result(value);
  if (a.id != none || b.id != none) {
    result.id = tape.size(); tape.push_back({a.id, b.id, da, db, 0});
  }
  return result;
}
inline Real operator+(Real a, Real b) { return operation(a.value+b.value,a,b,1,1); }
inline Real operator-(Real a, Real b) { return operation(a.value-b.value,a,b,1,-1); }
inline Real operator*(Real a, Real b) { return operation(a.value*b.value,a,b,b.value,a.value); }
inline Real operator/(Real a, Real b) { return operation(a.value/b.value,a,b,1/b.value,-a.value/(b.value*b.value)); }
inline Real operator-(Real a) { return Real(0)-a; }
#define ORACLE_MIXED(OP) \
template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
inline Real operator OP(Real a, T b) { return a OP Real(double(b)); } \
template<class T, std::enable_if_t<std::is_arithmetic_v<T>, int> = 0> \
inline Real operator OP(T a, Real b) { return Real(double(a)) OP b; }
ORACLE_MIXED(+)
ORACLE_MIXED(-)
ORACLE_MIXED(*)
ORACLE_MIXED(/)
#undef ORACLE_MIXED
inline Real &Real::operator+=(Real b) { return *this=*this+b; }
inline Real &Real::operator-=(Real b) { return *this=*this-b; }
inline Real &Real::operator*=(Real b) { return *this=*this*b; }
inline Real &Real::operator/=(Real b) { return *this=*this/b; }
inline Real log(Real a) { return operation(std::log(a.value),a,Real(),1/a.value,0); }
inline Real exp(Real a) { double v=std::exp(a.value); return operation(v,a,Real(),v,0); }
inline Real pow(Real a, double b) { return operation(std::pow(a.value,b),a,Real(),b*std::pow(a.value,b-1),0); }
inline Real tag(Real a, int i, int j) {
  // Independent log multiplier at this pair-emitting chart assignment.
  Real r(a.value); r.id=tape.size(); tape.push_back({a.id,none,1,0,0});
  tags.push_back({r.id,i,j,a.value}); return r;
}
inline void reverse(Real energy) {
  if (energy.id == none) return;
  tape[energy.id].adj = 1;
  for (std::size_t n=tape.size(); n-- > 0;) {
    const auto v=tape[n];
    if (v.adj == 0) continue; // Avoid 0*NaN on unused zero-weight cells.
    if (v.a != none) tape[v.a].adj += v.adj*v.da;
    if (v.b != none) tape[v.b].adj += v.adj*v.db;
  }
}
}
