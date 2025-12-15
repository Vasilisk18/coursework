#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>

class AgentModel3 {
  int N, T;
  double c, r, p, s, q, dt;

  std::vector<int> agents;   // 0 = Z, 1 = Y, 2 = V
  std::vector<int> result;   // z,y,v per step

  std::mt19937 gen;
  std::uniform_real_distribution<> U;

public:
  AgentModel3(int N_, int z0, int y0, int v0,
    int T_,
    double c_, double r_, double p_,
    double s_, double q_,
    double dt_ = 0.1)
    : N(N_), T(T_), c(c_), r(r_), p(p_), s(s_), q(q_), dt(dt_),
    gen(std::random_device{}()), U(0.0, 1.0)
  {
    agents.resize(N, 2);
    result.resize(T * 3);

    int k = 0;
    while (z0-- > 0) agents[k++] = 0;
    while (y0-- > 0) agents[k++] = 1;
    while (k < N)   agents[k++] = 2;

    save_counts(0);
  }

  void Execution() {
    for (int t = 1; t < T; ++t) {

      // текущие доли
      int cz = result[(t - 1) * 3];
      int cy = result[(t - 1) * 3 + 1];
      int cv = N - cz - cy;

      double z = (double)cz / N;
      double y = (double)cy / N;
      double v = (double)cv / N;

      for (int i = 0; i < N; ++i) {
        double u = U(gen);

        // ---------- спонтанные переходы ----------
        if (agents[i] == 0) {              // Z
          if (u < r * dt) { agents[i] = 1; continue; }
          if (u < (r + s) * dt) { agents[i] = 2; continue; }
        }
        else if (agents[i] == 1) {         // Y
          if (u < s * dt) { agents[i] = 2; continue; }
        }

        // ---------- парное взаимодействие ----------
        int j = gen() % N;
        if (j == i) continue;

        int a = agents[i];
        int b = agents[j];

        // Z + Y ? Y
        if (a == 0 && b == 1 && U(gen) < p * dt) {
          agents[i] = 1;
        }

        // Y + V ? V
        else if (a == 1 && b == 2 && U(gen) < q * dt) {
          agents[i] = 2;
        }

        // V + Z ? Z   (c z? v)
        else if (a == 2 && b == 0 && U(gen) < c * z * dt) {
          agents[i] = 0;
        }

        // V + Y ? Y   (c y z v)
        else if (a == 2 && b == 1 && U(gen) < c * z * dt) {
          agents[i] = 1;
        }
      }

      save_counts(t);
    }
  }

  void save_counts(int t) {
    int z = 0, y = 0;
    for (int a : agents) {
      if (a == 0) ++z;
      else if (a == 1) ++y;
    }
    result[t * 3] = z;
    result[t * 3 + 1] = y;
    result[t * 3 + 2] = N - z - y;
  }

  void Print() const {
    std::cout << std::setw(6) << "z"
      << std::setw(6) << "y"
      << std::setw(6) << "v" << "\n";
    for (int t = 0; t < T; ++t)
      std::cout << std::setw(6) << result[t * 3]
      << std::setw(6) << result[t * 3 + 1]
      << std::setw(6) << result[t * 3 + 2] << "\n";
  }

  std::vector<int> OutResult() const {
    return result;
  }
};
