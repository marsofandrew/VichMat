//
// Created by marsofandrew on 10.02.18.
//

#ifndef LABS_COMMON_HPP
#define LABS_COMMON_HPP

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

template<typename T, const std::size_t rows, const std::size_t columns>
double countNorm(const T (&matrix)[rows][columns])
{
  T result = 0;
  for (std::size_t i = 0; i < rows; i++) {
    T columnSum = 0;
    for (const T t:matrix[i]) {
      columnSum += detail::abs(t);
    }
    result = std::max(result, columnSum);
  }
  return result;
}

#endif //LABS_COMMON_HPP
