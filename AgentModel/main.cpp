#include "AgentModel.h"
#include "CompModel.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
  int N = 100;
  int T = 100;
  int series = 100;
  std::vector<double> medium(T * 3, 0.0);
  for (int i = 0; i < series; i++) {
    AgentModel3 A(N, 10, 10, 80, T, 0.001, 0.5, 0.0, 0.1, 0.1);
    A.Execution();

    std::vector<int> curres = A.OutResult();
    for (int j = 0; j < T * 3; j++)
      medium[j] += curres[j];
  }

  for (int i = 0; i < T * 3; i++)
    medium[i] /= series;

  // ---- formatted console output ----
  std::cout << std::setw(8) << "z"
    << std::setw(8) << "y"
    << std::setw(8) << "v" << std::endl;

  for (int i = 0; i < T; i++) {
    std::cout << std::setw(8) << std::fixed << std::setprecision(2) << medium[i * 3]
      << std::setw(8) << medium[i * 3 + 1]
      << std::setw(8) << medium[i * 3 + 2]
      << std::endl;
  }

  //CompModel3 C(N, 10, 10, 80, T, 0.001, 0.5, 0.0, 0.1, 0.1);
  //C.rungeKutta(0.1);
  //C.Print();
  //C.SaveToExcel("resultComp.csv");

  std::ofstream fout("resultAgent.csv");
  fout << "z;y;v\n";
  for (int i = 0; i < T; i++) {
    fout << medium[i * 3] << ";"
         << medium[i * 3 + 1] << ";"
         << medium[i * 3 + 2] << "\n";
  }

  fout.close();

  return 0;
}