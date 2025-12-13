#include "AgentModel.h"
#include "CompModel.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
  int N = 90;
  int T = 80;
  int series = 50;
  std::vector<double> medium(T * 3, 0.0);
  //for (int i = 0; i < series; i++) {
  //  AgentModel3 A(N, 30, 30, 30, T,
  //    0.1, 0.2, 0.1, 0.3, 0.2);
  //  A.Execution();

  //  std::vector<int> curres = A.OutResult();
  //  for (int j = 0; j < T * 3; j++)
  //    medium[j] += curres[j];
  //}

  //for (int i = 0; i < T * 3; i++)
  //  medium[i] /= series;

  //// ---- formatted console output ----
  //std::cout << std::setw(8) << "z"
  //  << std::setw(8) << "y"
  //  << std::setw(8) << "v" << std::endl;

  //for (int i = 0; i < T; i++) {
  //  std::cout << std::setw(8) << std::fixed << std::setprecision(2) << medium[i * 3]
  //    << std::setw(8) << medium[i * 3 + 1]
  //    << std::setw(8) << medium[i * 3 + 2]
  //    << std::endl;
  //}

  CompModel3 C(N, 30, 30, 30, T, 0.01, 0.2, 0.1, 0.3, 0.2);
  C.rungeKutta(0.01);
  C.Print();

  //std::ofstream fout("resultAgent.csv");
  //fout << "z;y;v\n";

  //for (int i = 0; i < T; i++) {
  //  fout << medium[i * 3] << ";"
  //       << medium[i * 3 + 1] << ";"
  //       << medium[i * 3 + 2] << "\n";
  //}

  //fout.close();

  return 0;
}