//
// Created by marsofandrew on 10.02.18.
//

#ifndef LABS_COMMON_HPP
#define LABS_COMMON_HPP

#include <cmath>
#include <vector>

double countLagrange(const std::vector<std::pair<double, double>> &points, const std::size_t n, const double x);

namespace detail
{
  template<typename T>
  T abs(T &value)
  {
    return (value >= 0 ? value : -value);
  }
}

template<const std::size_t rows, const std::size_t columns>
double countNorm(const double (&matrix)[rows][columns])
{
  double result = 0;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      result += std::abs(matrix[i][j]) * std::abs(matrix[i][j]);
    }
  }
  result = std::sqrt(result);
  return result;
}

#endif //LABS_COMMON_HPP
