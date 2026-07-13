// GL degree-one vortex interval certificate
//
// Compact CAPD adaptation of the AYMH Part I v6 certificate style to the
// non-magnetic Ginzburg-Landau vortex used in Section 6 of GL-arxiv-final.

#include <capd/capdlib.h>
#ifdef GL_USE_MP
#include <capd/mpcapdlib.h>
#endif

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace capd;

namespace glcert {

#ifdef GL_USE_MP
using R = MpFloat;
using I = MpInterval;
using Vec = MpIVector;
using Map = MpIMap;
using Solver = MpIOdeSolver;
using TimeMap = MpITimeMap;
using Set = MpC0TripletonSet;
#else
using R = double;
using I = DInterval;
using Vec = IVector;
using Map = IMap;
using Solver = IOdeSolver;
using TimeMap = ITimeMap;
using Set = C0TripletonSet;
#endif

#ifdef GL_USE_MP
static R RF(const char* s) { return R(s); }
#else
static R RF(const char* s) { return std::stod(s); }
#endif

static I IV(const char* s) { return I(RF(s)); }
static I IV(const char* a, const char* b) {
#ifdef GL_USE_MP
  return I(RF(a), RF(b));
#else
  double left = std::nextafter(RF(a), -std::numeric_limits<double>::infinity());
  double right = std::nextafter(RF(b), std::numeric_limits<double>::infinity());
  return I(left, right);
#endif
}

static I zero() { return IV("0"); }
static I one() { return IV("1"); }
static I two() { return IV("2"); }
static I three() { return IV("3"); }
static I isqr(const I& x) { return x * x; }
static I icube(const I& x) { return x * x * x; }

static I ipow(I x, unsigned n) {
  I y = one();
  while (n) {
    if (n & 1U) y *= x;
    x *= x;
    n >>= 1U;
  }
  return y;
}

static bool positive(const I& x) { return x.leftBound() > RF("0"); }
static bool negative(const I& x) { return x.rightBound() < RF("0"); }
static bool upperLess(const I& x, const char* b) { return x.rightBound() < RF(b); }
static std::string pass(bool ok) { return ok ? "PASS" : "FAIL"; }

static I symmetricError(const I& eps) {
  R e = eps.rightBound();
  return I(-e, e);
}
static I widen(const I& x, const I& eps) { return x + symmetricError(eps); }
static I upperSingleton(const I& x) { return I(x.rightBound()); }
static I absUpperAsInterval(const I& x) {
  R a = std::abs(x.leftBound());
  R b = std::abs(x.rightBound());
  return I(a > b ? a : b);
}
static I nonnegativePartUpper(const I& x) {
  if (x.rightBound() <= RF("0")) return zero();
  if (x.leftBound() >= RF("0")) return x;
  return I(RF("0"), x.rightBound());
}
static I intervalWidth(const I& x) { return I(x.rightBound() - x.leftBound()); }
static void requireTrue(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
static void printWidth(const std::string& label, const I& x) {
  std::cout << "  width(" << label << ") <= " << intervalWidth(x) << "\n";
}

struct Settings {
  static constexpr int default_precision_bits = 220;
  int precision_bits = default_precision_bits;
  int taylor_order = 40;
  int series_order = 90;

  I c_box = IV("0.583189495860", "0.583189495862");
  I c_lower = IV("0.583189495860");
  I c_upper = IV("0.583189495862");

  I r_origin = IV("0.1");
  I r_shoot = IV("20");
  I r_turn = IV("0.6140");
  I r_tail = IV("7");

  int turn_cells = 6000;
  int lt_cells = 50000;
};

// U(r)=sum a_n r^(2n+1), a_0=c.
// 4(n+1)(n+2)a_{n+1} + a_n - sum_{i+j+k=n-1} a_i a_j a_k = 0.
static std::vector<I> glCoeffs(const I& c, int nMax) {
  std::vector<I> a(nMax + 1, zero());
  a[0] = c;
  for (int n = 0; n < nMax; ++n) {
    I cubic = zero();
    if (n >= 1) {
      for (int i = 0; i <= n - 1; ++i) {
        for (int j = 0; j <= n - 1 - i; ++j) {
          int k = n - 1 - i - j;
          cubic += a[i] * a[j] * a[k];
        }
      }
    }
    I den = IV(std::to_string(4 * (n + 1) * (n + 2)).c_str());
    a[n + 1] = (cubic - a[n]) / den;
  }
  return a;
}

struct SeriesTail { I valueTail, derivativeTail, observedRatio, lastTerm; };

static SeriesTail seriesTailCertificate(const std::vector<I>& coeffs, const I& x,
                                        const char* ratioMax) {
  int M = static_cast<int>(coeffs.size()) - 1;
  requireTrue(M >= 10, "series tail needs at least 10 coefficients");
  I rho = IV(ratioMax);
  I ratio = zero();
  for (int j = M - 8; j <= M; ++j) {
    I prev = absUpperAsInterval(coeffs[j - 1]) * ipow(x, static_cast<unsigned>(j - 1));
    I curr = absUpperAsInterval(coeffs[j]) * ipow(x, static_cast<unsigned>(j));
    if (prev.rightBound() > RF("0")) {
      I r = curr / prev;
      if (r.rightBound() > ratio.rightBound()) ratio = r;
    }
  }
  requireTrue(upperLess(ratio, ratioMax), "origin series geometric tail ratio failed");
  I last = absUpperAsInterval(coeffs[M]) * ipow(x, static_cast<unsigned>(M));
  I oneMinus = one() - rho;
  I valueTail = last * rho / oneMinus;
  I m = IV(std::to_string(M).c_str());
  I weighted = m * rho / oneMinus + rho / isqr(oneMinus);
  I derivativeTail = last * weighted / x;
  return {valueTail, derivativeTail, ratio, last};
}

struct GLState { I u, up; };

static GLState glSeries(const I& r, const I& c, int nMax, bool verbose = false) {
  std::vector<I> a = glCoeffs(c, nMax);
  I x = isqr(r);
  I p = zero();
  I px = zero();
  for (int n = static_cast<int>(a.size()) - 1; n >= 0; --n) p = p * x + a[n];
  for (int n = static_cast<int>(a.size()) - 1; n >= 1; --n) {
    px = px * x + IV(std::to_string(n).c_str()) * a[n];
  }
  SeriesTail tail = seriesTailCertificate(a, x, "0.02");
  I u = widen(r * p, r * tail.valueTail);
  I up = widen(p + two() * x * px, tail.valueTail + two() * x * tail.derivativeTail);
  if (verbose) {
    std::cout << "origin Frobenius tail at r=" << r << "\n";
    std::cout << "  observed ratio = " << tail.observedRatio << "\n";
    std::cout << "  last term <= " << tail.lastTerm << "\n";
    std::cout << "  U tail <= " << r * tail.valueTail << "\n";
    std::cout << "  U' tail <= " << tail.valueTail + two() * x * tail.derivativeTail << "\n";
  }
  return {u, up};
}

static const std::string VF_GL_FWD =
    "var:r,u,v,c;"
    "fun:1,v,-v/r+u/(r*r)-u+u*u*u,0;";

static Vec integrate(const std::string& vf, const Vec& x0, const I& T, int order) {
  Map map(vf);
  Solver solver(map, order);
  TimeMap tm(solver);
  Set set(x0);
  return tm(T, set);
}

struct StepResult { Vec endpoint, tube; };

static StepResult integrateTube(const std::string& vf, const Vec& x0, const I& T, int order) {
  Map map(vf);
  Solver solver(map, order);
  TimeMap tm(solver);
  TimeMap::SolutionCurve curve(I(RF("0")));
  Set set(x0);
  Vec endpoint = tm(T, set, curve);
  I domain(RF("0"), T.rightBound());
  Vec tube = curve(domain);
  return {endpoint, tube};
}

static Vec initialGL(const Settings& S, const I& c, bool verbose = false) {
  GLState g = glSeries(S.r_origin, c, S.series_order, verbose);
  Vec x(4);
  x[0] = S.r_origin; x[1] = g.u; x[2] = g.up; x[3] = c;
  return x;
}

static I potentialV(const I& r, const I& u) {
  return one() / isqr(r) - three() * (one() - isqr(u));
}
static I ltIntegrandPositivePart(const I& r, const I& u) {
  I g = three() * (one() - isqr(u)) - one() / isqr(r);
  I gp = nonnegativePartUpper(g);
  return icube(gp) * upperSingleton(r);
}

static bool checkShootingBracket(const Settings& S) {
  std::cout << "\n[1] GL vortex shooting bracket\n";
  std::cout << "c_box = " << S.c_box << "\n";
  Vec xl = initialGL(S, S.c_lower, true);
  Vec xu = initialGL(S, S.c_upper, false);
  Vec yl = integrate(VF_GL_FWD, xl, S.r_shoot - S.r_origin, S.taylor_order);
  Vec yu = integrate(VF_GL_FWD, xu, S.r_shoot - S.r_origin, S.taylor_order);
  I lowerUminus1 = yl[1] - one();
  I lowerUp = yl[2];
  I upperUminus1 = yu[1] - one();
  std::cout << "lower endpoint U(20)-1 = " << lowerUminus1 << "\n";
  std::cout << "lower endpoint U'(20)  = " << lowerUp << "\n";
  std::cout << "upper endpoint U(20)-1 = " << upperUminus1 << "\n";
  std::cout << "Analytic input: standard shooting dichotomy/uniqueness for the degree-one GL vortex.\n";
  bool lowerBad = negative(lowerUminus1) && negative(lowerUp);
  bool upperBad = positive(upperUminus1);
  bool ok = lowerBad && upperBad;
  std::cout << "lower endpoint turn-over check: " << pass(lowerBad) << "\n";
  std::cout << "upper endpoint overshoot check:  " << pass(upperBad) << "\n";
  std::cout << pass(ok) << "\n";
  return ok;
}

static bool checkTurningPositivity(const Settings& S) {
  std::cout << "\n[2] Positivity before the turning radius\n";
  std::cout << "r_turn = " << S.r_turn << "\n";
  std::cout << "For 0<r<=0.1, U(r)<=r gives V(r)>=r^-2-3>0 analytically.\n";
  Vec x = initialGL(S, S.c_box, false);
  I h = (S.r_turn - S.r_origin) / IV(std::to_string(S.turn_cells).c_str());
  I minV = IV("1e9");
  for (int j = 0; j < S.turn_cells; ++j) {
    StepResult st = integrateTube(VF_GL_FWD, x, h, S.taylor_order);
    I V = potentialV(st.tube[0], st.tube[1]);
    if (V.leftBound() < minV.leftBound()) minV = V;
    x = st.endpoint;
  }
  std::cout << "minimum lower enclosure for V = r^-2 - 3(1-U^2): " << minV << "\n";
  printWidth("turning V enclosure", minV);
  bool ok = positive(minV);
  std::cout << pass(ok) << "\n";
  return ok;
}

static bool checkTailEntry(const Settings& S) {
  std::cout << "\n[3] Tail-entry inequality at r=7\n";
  Vec x = initialGL(S, S.c_box, false);
  Vec y = integrate(VF_GL_FWD, x, S.r_tail - S.r_origin, S.taylor_order);
  I boundary = one() - IV("0.55") / isqr(S.r_tail);
  I margin = y[1] - boundary;
  std::cout << "U(7) interval = " << y[1] << "\n";
  std::cout << "U'(7) interval = " << y[2] << "\n";
  std::cout << "U(7) - (1 - 0.55/7^2) = " << margin << "\n";
  bool ok = positive(margin);
  std::cout << pass(ok) << "\n";
  return ok;
}

static I tailAllowanceRaw(const Settings& S) {
  return icube(IV("2.3")) / (IV("4") * ipow(S.r_tail, 4));
}

static bool checkLTIntegral(const Settings& S) {
  std::cout << "\n[4] Lieb-Thirring compact integral\n";
  Vec x = initialGL(S, S.c_box, false);
  x = integrate(VF_GL_FWD, x, S.r_turn - S.r_origin, S.taylor_order);
  I h = (S.r_tail - S.r_turn) / IV(std::to_string(S.lt_cells).c_str());
  I compact = zero();
  for (int j = 0; j < S.lt_cells; ++j) {
    StepResult st = integrateTube(VF_GL_FWD, x, h, S.taylor_order);
    compact += ltIntegrandPositivePart(st.tube[0], st.tube[1]) * h;
    x = st.endpoint;
  }
  I compactOver6 = compact / IV("6");
  I tailRaw = tailAllowanceRaw(S);
  I tailOver6 = tailRaw / IV("6");
  I total = compactOver6 + tailOver6;
  I lambdaLower = IV("2") - sqrt(total);
  std::cout << "compact raw integral upper = " << compact << "\n";
  std::cout << "compact / 6 upper = " << compactOver6 << "\n";
  std::cout << "tail raw upper = " << tailRaw << "\n";
  std::cout << "tail / 6 upper = " << tailOver6 << "\n";
  std::cout << "LT total upper = " << total << "\n";
  std::cout << "lambda lower interval 2 - sqrt(total) = " << lambdaLower << "\n";
  printWidth("compact/6", compactOver6);
  bool compactOk = upperLess(compactOver6, "0.4454");
  bool totalOk = upperLess(total, "0.4456");
  std::cout << "compact target (<0.4454): " << pass(compactOk) << "\n";
  std::cout << "total target (<0.4456):   " << pass(totalOk) << "\n";
  std::cout << pass(compactOk && totalOk) << "\n";
  return compactOk && totalOk;
}

static void writeIntervalCsv(std::ostream& os, const I& x) {
  os << x.leftBound() << "," << x.rightBound();
}
static void emitProfileTubeCsv(const Settings& S, const std::string& path, int cells) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open profile CSV: " + path);
  out << std::setprecision(40);
  out << "cell,r_left,r_right,U_endpoint_lo,U_endpoint_hi,Up_endpoint_lo,Up_endpoint_hi,"
         "U_tube_lo,U_tube_hi,Up_tube_lo,Up_tube_hi\n";
  Vec x = initialGL(S, S.c_box, false);
  I h = (S.r_tail - S.r_origin) / IV(std::to_string(cells).c_str());
  for (int j = 0; j < cells; ++j) {
    I rLeft = x[0];
    StepResult st = integrateTube(VF_GL_FWD, x, h, S.taylor_order);
    I rRight = st.endpoint[0];
    out << j << "," << rLeft.leftBound() << "," << rRight.rightBound() << ",";
    writeIntervalCsv(out, st.endpoint[1]); out << ",";
    writeIntervalCsv(out, st.endpoint[2]); out << ",";
    writeIntervalCsv(out, st.tube[1]); out << ",";
    writeIntervalCsv(out, st.tube[2]); out << "\n";
    x = st.endpoint;
  }
}

} // namespace glcert

int main() {
  using namespace glcert;
  try {
#ifdef GL_USE_MP
    R::setDefaultPrecision(Settings::default_precision_bits);
#endif
    Settings S;
    if (const char* csvPath = std::getenv("GL_EXPORT_PROFILE_TUBE")) {
      int cells = 6900;
      if (const char* rawCells = std::getenv("GL_EXPORT_PROFILE_TUBE_CELLS")) cells = std::atoi(rawCells);
      emitProfileTubeCsv(S, csvPath, cells);
      std::cout << "wrote profile tube CSV: " << csvPath << "\n";
      return EXIT_SUCCESS;
    }
    std::cout.precision(30);
    std::cout << "GL vortex CAPD interval certificate\n";
#ifdef GL_USE_MP
    std::cout << "interval backend = MpInterval/MpFloat, requested precision bits = " << S.precision_bits << "\n";
#else
    std::cout << "interval backend = DInterval/double\n";
#endif
    std::cout << "Taylor order = " << S.taylor_order << "\n";
    bool ok = true;
    ok = checkShootingBracket(S) && ok;
    ok = checkTurningPositivity(S) && ok;
    ok = checkTailEntry(S) && ok;
    ok = checkLTIntegral(S) && ok;
    if (!ok) {
      std::cout << "\nOVERALL: FAIL\n";
      return EXIT_FAILURE;
    }
    std::cout << "\nOVERALL: PASS FOR THE GL SECTION 6 NUMERICAL CERTIFICATE\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "\nException: " << e.what() << "\n";
    std::cerr << "OVERALL: FAIL\n";
    return EXIT_FAILURE;
  }
}
