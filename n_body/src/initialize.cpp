#include "../include/initialize.h"
#include <cmath>

void initialize(Stars &stars, int N, double sigma_pos, Random& rng)
{
    stars.resize(N);

    const double m = 1.0 / N;
    const double two_sigma_sq = 2.0 * sigma_pos * sigma_pos;

    for (int i = 0; i < N; ++i) {
        stars.m[i]  = m;
        stars.x[i]  = sigma_pos * rng.normal();
        stars.y[i]  = sigma_pos * rng.normal();
        stars.ax[i] = 0.0;
        stars.ay[i] = 0.0;
    }

    double cx = 0.0, cy = 0.0;
    for (int i = 0; i < N; ++i) {
        cx += stars.m[i] * stars.x[i];
        cy += stars.m[i] * stars.y[i];
    }
    for (int i = 0; i < N; ++i) {
        stars.x[i] -= cx;
        stars.y[i] -= cy;
    }

    for (int i = 0; i < N; ++i) {
        double x = stars.x[i];
        double y = stars.y[i];
        double r = std::sqrt(x * x + y * y);
        double M_enc = 1.0 - std::exp(-r * r / two_sigma_sq);
        double v_circ = (r > 0.0) ? std::sqrt(G * M_enc / r) : 0.0;
        stars.vx[i] = -v_circ * y / r;
        stars.vy[i] =  v_circ * x / r;
    }
}

void virial_rescale(Stars& s, double eps)
{
    const std::size_t N = s.size();
    const double eps_sq = eps * eps;

    double K = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        K += 0.5 * s.m[i] * (s.vx[i] * s.vx[i] + s.vy[i] * s.vy[i]);
    }

    double U = 0.0;
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j) {
            double dx = s.x[i] - s.x[j];
            double dy = s.y[i] - s.y[j];
            double r  = std::sqrt(dx * dx + dy * dy + eps_sq);
            U -= G * s.m[i] * s.m[j] / r;
        }

    const double factor = std::sqrt(-U / (2.0 * K));
    for (std::size_t i = 0; i < N; ++i) {
        s.vx[i] *= factor;
        s.vy[i] *= factor;
    }
}
