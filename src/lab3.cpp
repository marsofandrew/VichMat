//
// Created by marsofandrew on 18.02.18.
//
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include "libs/cmath.h"
#include <stdexcept>
#include <array>

namespace detail
{
  void lab3UseRkf45();

  void countAdams(const double h, const double xStart, const double xEnd);

  template<const std::size_t N>
  std::array<double, N> adamsIteration(const double h, const std::array<double, N> xn, const double (&fn)[2][N]);
}

int function(int n, double t, double x[], double dxdt[])
{
  dxdt[0] = -73 * x[0] - 210 * x[1] + std::log(1 + t * t);
  dxdt[1] = x[0] + std::exp(-t) + t * t + 1;
  return 0;
}

void lab3()
{
  detail::lab3UseRkf45();
  std::cout << "\nh = 0.025\n";
  detail::countAdams(0.025, 0, 1);

  std::cout << "\nh = 0.0\n";
  detail::countAdams(0.01, 0, 1);
  
}

void detail::lab3UseRkf45()
{
  int n = 2;
  double x[n] = {-3, 1};
  double out[n];
  double tStart = 0;
  double tEnd = 1;
  const double printStep = 0.05;
  double rErr = 0.0001;
  double aErr = 0.0001;
  const int maxFuncEval = 5000;
  int funcEvalCounter = 0;
  int flag = 0;

  rkfinit(n, &flag);
  if (flag != 0) {
    throw std::runtime_error("some problem in rkfinit");
  }

  flag = 1;
  double nextStepSizeSuggested = 0;
  for (double tCurrent = tStart; tCurrent < tEnd;) {
    rkf45(function, n, x, out, &tCurrent, tCurrent + printStep, &rErr, aErr, &nextStepSizeSuggested, &funcEvalCounter,
          maxFuncEval, &flag);

    std::cout << "tcurrent = " << tCurrent << "          x[0] = " << x[0] << "            x[1] = " << x[1] << "\n";
  }

}

void detail::countAdams(const double h, const double xStart, const double xEnd)
{

  const int n = 2;
  double z = 0;
  std::array<double, 2> x = {-3, 1};
  double t = xStart;
  double out[2];
  double f[2][n];
  double rErr = 0.0001;
  int flag = 1;
  int maxFuncEval = 5000;
  int funcEvalCounter = 0;
  double nextStepSizeSuggested = 0;

  for (int i = 0; i < 2; ++i) {
    rkf45(function, n, x.data(), out, &t, t + h, &rErr, 0.0001,
          &nextStepSizeSuggested, &funcEvalCounter, maxFuncEval, &flag);
    std::cout << "t = " << t << "     x0 = " << x[0] << "   x1 = " << x[1] << "\n";
    function(n, t, x.data(), f[i]);
  }

  t += h;
  for (; t <= xEnd + h / 2; t += h) {

    x = detail::adamsIteration(h, x, f);
    std::cout << "t = " << t << "     x0 = " << x[0] << "   x1 = " << x[1] << "\n";

/*<--------------------------- Correct---------------------------->*/

    for (int i = 0; i < 1; i++) {
      for (int j = 0; j < n; j++) {
        f[i][j] = f[i + 1][j];
      }
    }
    function(n, t, x.data(), f[1]);

/*<--------------------------- Correct---------------------------->*/
  }
}

template<const std::size_t N>
std::array<double, N> detail::adamsIteration(const double h, const std::array<double, N> xn, const double (&fn)[2][N])
{
  std::array<double, N> array = {};
  for (std::size_t i = 0; i < N; i++) {
    array[i] = xn[i] + h * (3 * fn[1][i] - fn[0][i]) / 2;
  }

  return array;
}


