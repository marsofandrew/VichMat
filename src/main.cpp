#include <iostream>
#include "impl.hpp"
#include "Common.hpp"

int main()
{
  std::cout << "Lab 1, Var 14 \n \n";
  lab1(1, 4, 0.375, 1.1875);
  std::cout << "\nLab 2, Var 14 \n \n";

  lab2({1.0, 0.1, 0.01, 0.0001, 0.000001});
  return 0;
}