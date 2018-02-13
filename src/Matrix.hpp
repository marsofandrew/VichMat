//
// Created by marsofandrew on 12.02.18.
//

#ifndef LABS_MATRIX_HPP
#define LABS_MATRIX_HPP

#include <cstddef>
#include <vector>
#include <algorithm>

namespace detail1
{
  template<typename T>
  T abs(T &value)
  {
    return (value >= 0 ? value : -value);
  }
}
//@param T = numeric type
template<typename T, const unsigned long rows, const unsigned long columns>
class Matrix
{
public:
  Matrix();

  Matrix(const T (&matrix)[rows][columns]);

  Matrix(const Matrix<T, rows, columns> &lhs);

  virtual ~Matrix() = default;

  Matrix<T, rows, columns> &operator=(const Matrix<T, rows, columns> &lhs);

  Matrix<T, rows, columns> &operator+=(const Matrix<T, rows, columns> &lhs);

  Matrix<T, rows, columns> &operator-=(const Matrix<T, rows, columns> &lhs);

  Matrix<T, rows, columns> &operator-(const Matrix<T, rows, columns> &lhs, const Matrix<T, rows, columns> &rhs);

  Matrix<T, rows, columns> &operator+(const Matrix<T, rows, columns> &lhs, const Matrix<T, rows, columns> &rhs);

  /*
  Matrix<T, rows, rows> operator*(const Matrix<T, rows, columns> &lhs, const Matrix<T, columns, rows> &rhs);
  */
  double countNorm() const;

protected:
  std::vector<std::vector<T>> matrix_;
};

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns>::Matrix()
{
  std::vector<T> vector = {};
  vector.assign(columns, 0);
  matrix_.assign(rows, vector);
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns>::Matrix(const T (&matrix)[rows][columns])
{
  for (std::size_t i = 0; i < rows; ++i) {
    matrix_.emplace_back(&matrix[i][0], &matrix[i][columns - 1] + 1);
  }
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns>::Matrix(const Matrix<T, rows, columns> &lhs):
  matrix_(lhs.matrix_) {}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns> &Matrix<T, rows, columns>::operator=(const Matrix<T, rows, columns> &lhs)
{
  matrix_ = lhs.matrix_;
  return *this;
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns> &Matrix<T, rows, columns>::operator+=(const Matrix<T, rows, columns> &lhs)
{
  for (std::size_t l = 0; l < rows; ++l) {
    for (std::size_t i = 0; i < columns; ++i) {
      matrix_[l][i] = matrix_[l][i] + lhs.matrix_[l][i];
    }
  }
  return *this;
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns> &Matrix<T, rows, columns>::operator-=(const Matrix<T, rows, columns> &lhs)
{
  for (std::size_t l = 0; l < rows; ++l) {
    for (std::size_t i = 0; i < columns; ++i) {
      matrix_[l][i] = matrix_[l][i] - lhs.matrix_[l][i];
    }
  }
  return *this;
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns> &Matrix<T, rows, columns>::operator+(const Matrix<T, rows, columns> &lhs,
                                                              const Matrix<T, rows, columns> &rhs)
{
  Matrix<T, rows, columns> tmp(lhs);
  tmp += rhs;
  return tmp;
}

template<typename T, const unsigned long rows, const unsigned long columns>
Matrix<T, rows, columns> &Matrix<T, rows, columns>::operator-(const Matrix<T, rows, columns> &lhs,
                                                              const Matrix<T, rows, columns> &rhs)
{
  Matrix<T, rows, columns> tmp(lhs);
  tmp -= rhs;
  return tmp;
}

template<typename T, const unsigned long rows, const unsigned long columns>
double Matrix<T, rows, columns>::countNorm() const
{
  double result = 0;
  for (std::size_t i = 0; i < rows; i++) {
    T columnSum = 0;
    for (const T t:matrix_[i]) {
      columnSum += detail1::abs(t);
    }
    result = std::max(result, columnSum);
  }
  return result;
}

#endif //LABS_MATRIX_HPP
