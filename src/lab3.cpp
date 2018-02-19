//
// Created by marsofandrew on 18.02.18.
//
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>
#include "libs/cmath.h"
#include <stdexcept>

namespace detail
{
  void lab3UseRkf45();
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

    std::cout << "tcurrent = " << tCurrent << "          x[0] = " << x[0] << "            x[1] = " << x[1]<< "\n";
  }

}


