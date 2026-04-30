#include "../include/integrator.h"
#include "../include/force_gpu.h"
#include "../include/plotter.h"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <cmath>

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

namespace {

__global__ void kick_drift(DeviceStars s, double dt)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i >= s.N) return;

    s.vx[i] += 0.5 * s.ax[i] * dt;
    s.vy[i] += 0.5 * s.ay[i] * dt;
    s.x[i]  += s.vx[i] * dt;
    s.y[i]  += s.vy[i] * dt;
}

__global__ void kick(DeviceStars s, double dt)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;

    if (i >= s.N) return;

    s.vx[i] += 0.5 * s.ax[i] * dt;
    s.vy[i] += 0.5 * s.ay[i] * dt;
}

}  // namespace

void Direct::simulate(Stars &s, std::size_t n_iter, double dt, double eps,
                      const std::string& plotname)
{
    auto t_start = std::chrono::steady_clock::now();

    DeviceStars d;
    uint32_t N = s.size();
    d.N = N;

    d.m = nullptr;
    d.x = nullptr; d.y = nullptr;
    d.vx = nullptr; d.vy = nullptr;
    d.ax = nullptr; d.ay = nullptr;

    CUDA_CHECK(cudaMalloc(&d.m, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.x, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.y, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.vx, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.vy, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.ax, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d.ay, N * sizeof(double)));

    CUDA_CHECK(cudaMemcpy(d.m, s.m.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.x, s.x.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.y, s.y.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.vx, s.vx.data(), N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d.vy, s.vy.data(), N * sizeof(double), cudaMemcpyHostToDevice));

    Direct::calculate_accelerations(d, eps);

    Plotter::Animator anim = Plotter::animator_begin(plotname);
    Plotter::animator_add_frame(anim, s);

    constexpr int block_size = 32;
    const int num_blocks = (N + block_size - 1) / block_size;

    for (std::size_t i = 0; i < n_iter; ++i)
    {
        kick_drift<<<num_blocks, block_size>>>(d, dt);
        Direct::calculate_accelerations(d, eps);
        kick<<<num_blocks, block_size>>>(d, dt);

        if (i % 10 == 9)
        {
            CUDA_CHECK(cudaMemcpy(s.x.data(), d.x, N * sizeof(double), cudaMemcpyDeviceToHost));
            CUDA_CHECK(cudaMemcpy(s.y.data(), d.y, N * sizeof(double), cudaMemcpyDeviceToHost));
            Plotter::animator_add_frame(anim, s);
            if (i % 200 == 199)
            {
                CUDA_CHECK(cudaMemcpy(s.vx.data(), d.vx, N * sizeof(double), cudaMemcpyDeviceToHost));
                CUDA_CHECK(cudaMemcpy(s.vy.data(), d.vy, N * sizeof(double), cudaMemcpyDeviceToHost));
                double energy = total_energy(s, eps);
                std::cout << "[direct] Total Energy at step " << i + 1 << ": " << energy << "\n";
            }
        }
    }

    Plotter::animator_end(anim, /*half_range=*/4.0, /*pause_sec=*/0.005);

    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFree(d.m));
    CUDA_CHECK(cudaFree(d.x));  CUDA_CHECK(cudaFree(d.y));
    CUDA_CHECK(cudaFree(d.vx)); CUDA_CHECK(cudaFree(d.vy));
    CUDA_CHECK(cudaFree(d.ax)); CUDA_CHECK(cudaFree(d.ay));

    double total = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count();
    std::cout << "[direct] Total run time: " << total << " s\n";
}

void BH::simulate(Stars &s, std::size_t n_iter, double dt, double eps,
                  const std::string& plotname, const std::string& treename)
{std::cerr << "Not implemented yet.\n";}
