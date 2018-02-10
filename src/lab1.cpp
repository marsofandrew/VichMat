//
// Created by marsofandrew on 10.02.18.
//
#include "libs/cmath.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include "Common.hpp"

double x = 0;

double function(const double z)
{
  return 1 / (std::pow(M_E, z) * (z + x));
}

double a = 0;
double b = 20;
const double xCountStart = 1.1875;

void lab1(const double xStart, const double xEnd, const double h)
{
  std::vector<std::pair<double, double>> points;

  for (x = xStart; x <= xEnd; x += h) {
    double epsA = 0;
    double epsR = 0.0001;
    double result = 0;
    double error = 0;
    int nofun = 0;
    double posn = 0;
    int flag = 0;

    quanc8(&function, a, b, epsA, epsR, &result, &error, &nofun, &posn, &flag);

    if (flag == 0) {
      points.emplace_back(x, result);
    }
  }
  double xArray[9];
  double yArray[9];
  for (int i = 0; i < points.size(); ++i) {
    xArray[i] = points[i].first;
    yArray[i] = points[i].second;
  }
  double b[9];
  double c[9];
  double d[9];
  int flag = 0;
  spline(9, 0, 0, 0, 0, xArray, yArray, b, c, d, &flag);//TODO: ask about unnecessary arguments
  for (int k = 0; k <= 7; k++) {
    double xp = xCountStart + k * h;
    int last = 0;
    double lagrange = countLagrange(points, 8, xp);
    double spline_ = seval(9, xp, xArray, yArray, b, c, d, &last);
    std::cout << k << ")  " << "x = " << xp << "      " << "Lagrange = " << lagrange << "       Spline = "
              << spline_ << "   abs(L-S) = " << std::abs(lagrange - spline_) << "\n";
  }
}
