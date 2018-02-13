#include <initializer_list>
#include <numeric>
#include <cmath>
#include <iostream>
#include "Common.hpp"
#include "libs/cmath.h"
//
// Created by marsofandrew on 10.02.18.
//

const double E[8][8] = {{1, 0, 0, 0, 0, 0, 0, 0},
                        {0, 1, 0, 0, 0, 0, 0, 0},
                        {0, 0, 1, 0, 0, 0, 0, 0},
                        {0, 0, 0, 1, 0, 0, 0, 0},
                        {0, 0, 0, 0, 1, 0, 0, 0},
                        {0, 0, 0, 0, 0, 1, 0, 0},
                        {0, 0, 0, 0, 0, 0, 1, 0},
                        {0, 0, 0, 0, 0, 0, 0, 1}};

template<typename T, const std::size_t size>
double countCond(const T (&matrix)[size][size])
{
  T matrix1[size][size];
  double cond = 1;
  int pivot[size];
  int flag = 0;
  decomp(size, size, &matrix1[0][0], &cond, pivot, &flag);

  return cond;
}

void lab2(const std::initializer_list<double> p)
{
  for (const double p1:p) {
    double array1[8][8] = {{p1 + 6, 2,   6,   8,  -2, 1,   8,   -5},
                           {6,      -22, -2,  -1, 0,  5,   -6,  4},
                           {-2,     -3,  -16, 0,  0,  -4,  2,   -5},
                           {1,      1,   4,   9,  1,  0,   0,   -6},
                           {0,      2,   0,   2,  -3, -5,  7,   5},
                           {6,      -2,  -4,  2,  -8, -12, 3,   -3},
                           {-6,     -6,  0,   -8, 0,  5,   -15, 0},
                           {0,      7,   6,   0,  -5, -8,  -5,  -3}};
    double cond = 1;
    int pivot[8];
    int flag = 0;
    decomp(8, 8, &array1[0][0], &cond, pivot, &flag);
    std::cout << "p = " << p1 << " cond = " << cond << "  ";
    //<COUNT INVERSE MATRIX>
    double inverse[8][8] = {};
    for (int i = 0; i < 8; ++i) {
      double vector[8] = {};
      vector[i] = 1;
      solve(8, 8, &array1[0][0], vector, pivot);
      for (int j = 0; j < 8; ++j) {
        inverse[j][i] = vector[j];
      }
    }
    //<COUNT INVERSE MATRIX>

    double matrixEn[8][8] = {};
// <--------------------multiple matrix----------------->
    for (int i = 0; i < 8; i++) {
      for (int k = 0; k < 8; k++) {
        for (int j = 0; j < 8; j++) {
          matrixEn[i][j] += array1[i][j] * inverse[k][j];
        }
      }
    }
//<-------------------------------multiple matrix---------------->
    double R[8][8] = {};
    for (int l = 0; l < 8; ++l) {
      for (int i = 0; i < 8; ++i) {
        R[l][i] = E[l][i] - matrixEn[l][i];
      }
    }
    std::cout << " Norm R = " << countNorm(R) << "\n";
  }
}

