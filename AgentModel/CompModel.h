#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

class CompModel3 {
  int N;                // total population (constant)
  double z0, y0;        // initial z, y
  int T;
  double c, r, p, s, q;
  std::vector<int> result; // store rounded results (z,y,v)

public:
  CompModel3() {
    N = 90;
    z0 = 30;
    y0 = 30.0;
    T = 100;
    c = r = p = s = q = 0.0;

    result.resize(T * 3);
    result[0] = static_cast<int>(round(z0));
    result[1] = static_cast<int>(round(y0));
    result[2] = N - result[0] - result[1];
  }

  CompModel3(int N_, double z1_, double y1_, double v1_,
    int T_, double c_, double r_, double p_, double s_, double q_) {
    N = N_;
    z0 = z1_;
    y0 = y1_;
    T = T_;
    c = c_;
    r = r_;
    p = p_;
    s = s_;
    q = q_;

    result.resize(T * 3);
    result[0] = static_cast<int>(round(z0));
    result[1] = static_cast<int>(round(y0));
    result[2] = N - result[0] - result[1];
  }

  // RHS of reduced 2D system
  std::vector<double> derivatives(const std::vector<double>& cur) {
    double z = cur[0];
    double y = cur[1];
    double v = N - z - y;

    std::vector<double> dz(2, 0.0);

    dz[0] = c * z * z * v - r * z - p * z * y - s * z;
    dz[1] = c * y * z * v + r * z + p * z * y - q * v * y - s * y;

    return dz;
  }

  // RK4 step (2D)
  std::vector<double> rk4_step(const std::vector<double>& z, double dt) {
    int n = 2;
    std::vector<double> k1 = derivatives(z);
    std::vector<double> z_temp(n, 0.0);

    for (int i = 0; i < n; ++i)
      z_temp[i] = z[i] + dt * k1[i] / 2.0;
    std::vector<double> k2 = derivatives(z_temp);

    for (int i = 0; i < n; ++i)
      z_temp[i] = z[i] + dt * k2[i] / 2.0;
    std::vector<double> k3 = derivatives(z_temp);

    for (int i = 0; i < n; ++i)
      z_temp[i] = z[i] + dt * k3[i];
    std::vector<double> k4 = derivatives(z_temp);

    std::vector<double> z_next(n, 0.0);
    for (int i = 0; i < n; ++i) {
      z_next[i] = z[i] + dt / 6.0 *
        (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);

      if (z_next[i] < 0.0)
        z_next[i] = 0.0;
    }

    // enforce invariant: z + y <= N
    double sum = z_next[0] + z_next[1];
    if (sum > N) {
      z_next[0] *= N / sum;
      z_next[1] *= N / sum;
    }

    return z_next;
  }

  void rungeKutta(double dt) {
    std::vector<double> z = { z0, y0 };

    for (int t = 1; t < T; ++t) {
      z = rk4_step(z, dt);

      int zi = static_cast<int>(round(z[0]));
      int yi = static_cast<int>(round(z[1]));
      int vi = N - zi - yi;

      if (vi < 0) vi = 0;

      result[t * 3 + 0] = zi;
      result[t * 3 + 1] = yi;
      result[t * 3 + 2] = vi;
    }
  }

  void Print() {
    std::cout << std::setw(5) << "z"
      << std::setw(5) << "y"
      << std::setw(5) << "v" << std::endl;

    for (int i = 0; i < T; ++i) {
      std::cout << std::setw(5) << result[i * 3]
        << std::setw(5) << result[i * 3 + 1]
        << std::setw(5) << result[i * 3 + 2] << std::endl;
    }
  }

  std::vector<int> OutResult() const {
    return result;
  }
};
