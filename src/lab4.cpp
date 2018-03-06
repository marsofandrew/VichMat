#include <cmath>
#include "libs/cmath.h"
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "libs/FORSYTHE.H"
//
// Created by marsofandrew on 05.03.18.
//

double integr(double x)
{
  return 1 / (std::sqrt(x) * (1 + std::pow(x, (1 / 3))));
}

REAL fx(REAL x)
{
  return std::exp(-x) + std::pow(x, 2) - 2;
}

double countM()
{

  REAL x1 = ZEROIN(0, 10000, fx, 0.000001);
  std::cout << "x1 = " << x1 << "\n";

  return 0.7598945 * x1;
}

double countL()
{
  double epsA = 0;
  double epsR = 0.0001;
  double result = 0;
  double error = 0;
  int nofun = 0;
  double posn = 0;
  int flag = 0;
  double down = 1;
  double up = 5;
  quanc8(&integr, down, up, epsA, epsR, &result, &error, &nofun, &posn, &flag);
  if (flag != 0) {
    throw std::runtime_error("problem in quanc8");
  }
  return 1 + std::pow((result - 1.0507242), 4);
}

constexpr double g = 9.81;
double K = 39.24;
double M = countM();
double L = countL();

int functionLab4(int n, double t, double u[], double dudt[])
{
  dudt[0] = u[1];
  dudt[2] = u[3];
  dudt[1] = -(K / M) * u[0] - g * (1 - std::cos(u[2])) + (L + u[0]) * std::pow(u[3], 2);
  dudt[3] = -(g / (L + u[0])) * std::sin(u[2]) - (2 * u[1] * u[3] / (L + u[0]));
  return 0;
}

void countStartExpressions(double array[4])
{
  const unsigned int size = 4;
  double matrix[size][size] = {{10, 1,  4,  0},
                               {1,  10, 5,  -1},
                               {4,  5,  10, 7},
                               {0,  -1, 7,  9}};

  double vector[size] = {16, 20, 40, 28};

  double cond = 1;
  int pivot[size];
  int flag = 0;
  decomp(size, size, &matrix[0][0], &cond, pivot, &flag);
  std::cout << "cond = " << cond << "\n";
  solve(size, size, &matrix[0][0], vector, pivot);
  for (int i = 0; i < size; ++i) {
    array[i] = vector[i];
  }

}

void lab4(const double tStart, const double tEnd, const std::string &&fileNameX, const std::string &&fileNameTeta)
{
  std::cout << "L = " << L << "    K = " << K << "     M = " << M << "\n";
  int n = 4;
  double array[4];
  countStartExpressions(array);
  std::cout << array[0] << " " << array[1] << " " << array[2] << " " << array[3] << "\n";

  double u[n] = {array[0], array[1], array[3], array[2]};

  double out[n];
  const double printStep = 0.035;
  double rErr = 0.0001;
  double aErr = 0.0001;
  const int maxFuncEval = 500000;
  int funcEvalCounter = 0;
  int flag = 0;

  rkfinit(n, &flag);
  if (flag != 0) {
    throw std::runtime_error("some problem in rkfinit");
  }

  flag = 1;
  double nextStepSizeSuggested = 0;
  std::ofstream outX(fileNameX);
  std::ofstream outTeta(fileNameTeta);
  outX << "# out of X\n";
  outTeta << "# out of Teta\n";

  outX << tStart << "\t" << u[0] << "\n";
  outTeta << tStart << "\t" << u[2] << "\n";
  std::cout << "tcurrent = " << tStart << "              x = " << u[0] << "               Tetta = " << u[2] << "\n";
  for (double tCurrent = tStart; tCurrent < tEnd;) {
    rkf45(functionLab4, n, u, out, &tCurrent, tCurrent + printStep, &rErr, aErr, &nextStepSizeSuggested,
          &funcEvalCounter,
          maxFuncEval, &flag);
    outX << tCurrent << "\t" << u[0] << "\n";
    outTeta << tCurrent << "\t" << u[2] << "\n";
    std::cout << "tcurrent = " << tCurrent << "          x = " << u[0] << "            Tetta = " << u[2] << "\n";
  }

}



