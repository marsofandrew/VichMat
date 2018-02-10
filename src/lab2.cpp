#include <initializer_list>

//
// Created by marsofandrew on 10.02.18.
//

double cond(const std::size_t row, const double (*array)[8])
{

}

void lab2(const double ...p)
{
  for (const double p1:std::initializer_list{p...}) {
    const double array[8][8] = {{p + 6, 2,   6,   8,  -2, 1,   8,   -5},
                                {6,     -22, -2,  -1, 0,  5,   -6,  4},
                                {-2,    -3,  -16, 0,  0,  -4,  2,   -5},
                                {1,     1,   4,   9,  1,  0,   0,   -6},
                                {0,     2,   0,   2,  -3, -5,  7,   5},
                                {6,     -2,  -4,  2,  -8, -12, 3,   -3},
                                {-6,    -6,  0,   -8, 0,  5,   -15, 0},
                                {0,     7,   6,   0,  -5, -8,  -5,  -3}};
    cond(8,array);
  }
}

