/*
 * matrix.h
 *
 *  Created on: Sep 15, 2015
 *      Author: hossein
 */

#ifndef MATRIX_H_
#define MATRIX_H_


#include <vector>

template <typename T> class Matrix {
 private:
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

 public:
  Matrix(unsigned _rows, unsigned _cols, const T& _initial);
  Matrix(unsigned _rows, unsigned _cols);
  Matrix(const Matrix<T>& rhs);
  virtual ~Matrix();

  // Operator overloading, for "standard" mathematical matrix operations
  Matrix<T>& operator=(const Matrix<T>& rhs);

  // Matrix mathematical operations
  Matrix<T> operator+(const Matrix<T>& rhs);
  Matrix<T>& operator+=(const Matrix<T>& rhs);
  Matrix<T> operator-(const Matrix<T>& rhs);
  Matrix<T>& operator-=(const Matrix<T>& rhs);
  Matrix<T> operator*(const Matrix<T>& rhs);
  Matrix<T>& operator*=(const Matrix<T>& rhs);
  Matrix<T> transpose();

  // Matrix/scalar operations
  Matrix<T> operator+(const T& rhs);
  Matrix<T> operator-(const T& rhs);
  Matrix<T> operator*(const T& rhs);
  Matrix<T> operator/(const T& rhs);

  // Matrix/vector operations
  std::vector<T> operator*(const std::vector<T>& rhs);
  std::vector<T> diag_vec();

  // Access the individual elements
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;

  // Access the row and column sizes
  unsigned get_rows() const;
  unsigned get_cols() const;

};


#endif /* MATRIX_H_ */
