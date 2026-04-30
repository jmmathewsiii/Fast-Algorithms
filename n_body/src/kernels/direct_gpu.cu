#include "../../include/force_gpu.h"
#include "../../include/quadtree_gpu.h"
#include <chrono>
#include <cstdlib>
#include <cstdio>

#define CUDA_CHECK(expr_to_check) do {            \
    cudaError_t result  = expr_to_check;          \
    if(result != cudaSuccess)                     \
    {                                             \
        fprintf(stderr,                           \
                "CUDA Runtime Error: %s:%i:%d = %s\n", \
                __FILE__,                         \
                __LINE__,                         \
                result,\
                cudaGetErrorString(result));      \
                exit(EXIT_FAILURE);               \
    }                                             \
} while(0)


__global__ void direct_kernel(DeviceStars s, double eps)
{
    int star_idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (star_idx >= s.N) return;
    double ax_sum = 0.;
    double ay_sum = 0.;

    double eps_sq = eps * eps;

    for (int i = 0; i < s.N; ++i)
    {
        if (star_idx == i) continue;
        double dx = s.x[star_idx] - s.x[i];
        double dy = s.y[star_idx] - s.y[i];
        double r_sq = dx * dx + dy * dy + eps_sq;
        double r_inv = rsqrt(r_sq);
        double r3_inv = r_inv * r_inv * r_inv;

        ax_sum -= s.m[i] * dx * r3_inv;
        ay_sum -= s.m[i] * dy * r3_inv;
    }
    ax_sum *= G;
    ay_sum *= G;

    s.ax[star_idx] = ax_sum;
    s.ay[star_idx] = ay_sum;
}

namespace {
    double g_force_seconds = 0.0;
}

void Direct::reset_timers() { g_force_seconds = 0.0; }
double Direct::force_compute_seconds() { return g_force_seconds; }

void Direct::calculate_accelerations(DeviceStars s, double eps)
{
    const static int block_size = 32;
    const int num_blocks = (s.N + block_size - 1) / block_size;

    auto t0 = std::chrono::steady_clock::now();
    direct_kernel<<<num_blocks, block_size>>>(s, eps);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    g_force_seconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
}
