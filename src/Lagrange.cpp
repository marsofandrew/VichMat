//
// Created by marsofandrew on 10.02.18.
//
#include <vector>
#include <stdexcept>

double l(const std::size_t i, const double x, const std::vector<std::pair<double, double> > &points)
{
  double result = 1;
  double xi = points[i].first;
  for (std::size_t k = 0; k < points.size(); k++) {
    double xk = points[k].first;
    if (k != i) {
      result *= (x - points[k].first) / (xi - xk);
    }
  }
  return result;
}

//@param points - vector<x,y>,
//@param n - order of polinomium
//@param x - counting point

double countLagrange(const std::vector<std::pair<double, double>> &points, const std::size_t n, const double x)
{

  if (n >= points.size()) {
    throw std::invalid_argument("n must be < points.size()");
  }
  double result = 0;
  for (std::size_t i = 0; i <= n; i++) {
    result += points[i].second * l(i, x, points);
  }

  return result;
}





