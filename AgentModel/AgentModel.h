#pragma once

#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <cmath>

// Agent-based approximation of the reduced 2D system:
// z' = c z^2 v - r z - p z y - s z
// y' = c y z v + r z + p z y - q v y - s y
// v = N - z - y

class AgentModel3 {
  int N;
  int T;
  double c, r, p, s, q;
  double dt;

  // agent states: 0 = Z, 1 = Y, 2 = V
  std::vector<int> agents;
  std::vector<int> result; // (z,y,v) per timestep

public:
  AgentModel3() {
    N = 90;
    T = 100;
    c = r = p = s = q = 0.0;
    dt = 0.1;
    agents.resize(N, 2);
    result.resize(T * 3, 0);
  }

  AgentModel3(int N_, int z_, int y_, int v_, int T_,
    double c_, double r_, double p_, double s_, double q_,
    double dt_ = 0.1) {

    N = N_;
    T = T_;
    c = c_;
    r = r_;
    p = p_;
    s = s_;
    q = q_;
    dt = dt_;

    agents.resize(N, 2);
    result.resize(T * 3, 0);

    int idx = 0;
    while (z_-- > 0 && idx < N) agents[idx++] = 0;
    while (y_-- > 0 && idx < N) agents[idx++] = 1;
    while (idx < N) agents[idx++] = 2;

    save_counts(0);
  }

  void Execution() {
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> U(0.0, 1.0);

    for (int t = 1; t < T; ++t) {

      int cz = 0, cy = 0, cv = 0;
      for (int a : agents) {
        if (a == 0) ++cz;
        else if (a == 1) ++cy;
        else ++cv;
      }

      double z = (double)cz / N;
      double y = (double)cy / N;
      double v = (double)cv / N;

      for (int i = 0; i < N; ++i) {
        int state = agents[i];
        double u = U(gen);

        // ---- spontaneous transitions ----
        if (state == 0) {          // Z
          if (u < r * dt) { agents[i] = 1; continue; } // Z ? Y
          if (u < (r + s) * dt) { agents[i] = 2; continue; } // Z ? V
        }
        else if (state == 1) {     // Y
          if (u < s * dt) { agents[i] = 2; continue; } // Y ? V
        }

        // ---- pairwise interactions ----
        int j = gen() % N;
        if (j == i) continue;

        int partner = agents[j];

        // Z + Y ? Y + Y   (p z y)
        if (state == 0 && partner == 1) {
          if (U(gen) < p * dt) {
            agents[i] = 1;
            continue;
          }
        }

        // Y + V ? V + V   (q v y)
        if (state == 1 && partner == 2) {
          if (U(gen) < q * dt) {
            agents[i] = 2;
            continue;
          }
        }

        // V ? Z   (c z^2 v)
        if (state == 2) {
          double prob = c * z * z * dt;
          if (prob > 1.0) prob = 1.0;
          if (U(gen) < prob) {
            agents[i] = 0;
            continue;
          }
        }

        // Z + Y + V ? Y   (c y z v)
        if (state == 2 && partner == 1) {
          double prob = c * y * z * dt;
          if (prob > 1.0) prob = 1.0;
          if (U(gen) < prob) {
            agents[i] = 1;
            continue;
          }
        }
      }

      save_counts(t);
    }
  }

  void save_counts(int t) {
    int cz = 0, cy = 0, cv = 0;
    for (int a : agents) {
      if (a == 0) ++cz;
      else if (a == 1) ++cy;
      else ++cv;
    }
    result[t * 3] = cz;
    result[t * 3 + 1] = cy;
    result[t * 3 + 2] = cv;
  }

  void Print() {
    std::cout << std::setw(5) << "z"
      << std::setw(5) << "y"
      << std::setw(5) << "v" << std::endl;

    for (int t = 0; t < T; ++t) {
      std::cout << std::setw(5) << result[t * 3]
        << std::setw(5) << result[t * 3 + 1]
        << std::setw(5) << result[t * 3 + 2] << std::endl;
    }
  }

  std::vector<int> OutResult() const {
    return result;
  }
};
