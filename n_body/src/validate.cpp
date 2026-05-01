#include "../include/validate.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

static std::string final_path(const std::string& tag)
{
    return "./output/" + tag + "_final.dat";
}

static std::string snapshot_path(const std::string& tag, std::size_t step)
{
    std::ostringstream oss;
    oss << "./output/lyapunov/" << tag << "_" << std::setw(7) << std::setfill('0')
        << step << ".dat";
    return oss.str();
}

void save_snapshot(const Stars& s, const std::string& tag, std::size_t step)
{
    std::string path = snapshot_path(tag, step);
    std::ofstream f(path.c_str(), std::ios::out);
    f << std::setprecision(17);
    const std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i) {
        f << s.x[i] << " " << s.y[i] << "\n";
    }
}

static bool load_snapshot(Stars& s, const std::string& tag, std::size_t step)
{
    std::string path = snapshot_path(tag, step);
    std::ifstream f(path.c_str());
    if (!f) return false;
    s.x.clear();
    s.y.clear();
    double x, y;
    while (f >> x >> y) {
        s.x.push_back(x);
        s.y.push_back(y);
    }
    return true;
}

void compare_snapshots(const std::string& baseline_tag,
                       const std::string& other_tag,
                       std::size_t n_iter, std::size_t stride,
                       const std::string& label)
{
    std::string out_path = "./output/lyapunov/error_" + label + ".csv";
    std::ofstream out(out_path.c_str());
    out << std::setprecision(17);
    out << "step,abs_max,abs_sum,rel_sum\n";

    std::size_t n_written = 0;
    for (std::size_t step = 0; step <= n_iter; step += stride) {
        Stars baseline, other;
        if (!load_snapshot(baseline, baseline_tag, step)) continue;
        if (!load_snapshot(other,    other_tag,    step)) continue;
        const std::size_t N = baseline.x.size();
        if (other.x.size() != N) continue;

        double max_err = 0.0, abs_err = 0.0, rel_err = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            double dx = other.x[i] - baseline.x[i];
            double dy = other.y[i] - baseline.y[i];
            double err = std::sqrt(dx * dx + dy * dy);
            if (err > max_err) max_err = err;
            abs_err += err;
            double mag = std::sqrt(baseline.x[i] * baseline.x[i] +
                                   baseline.y[i] * baseline.y[i]);
            if (mag > 0.0) rel_err += err / mag;
        }
        out << step << "," << max_err << "," << abs_err << "," << rel_err << "\n";
        ++n_written;
    }
    std::cout << "[validate] wrote " << n_written << " rows to " << out_path << "\n";

    // Emit a gnuplot script: log(error) vs step, semilog y.
    std::string script_path = "./scripts/lyapunov/error_" + label + ".gp";
    std::ofstream gp(script_path.c_str());
    gp << "set datafile separator ','\n"
       << "set key autotitle columnhead\n"
       << "datafile = '../../output/lyapunov/error_" << label << ".csv'\n"
       << "\n"
       << "# Time per step. EDIT if dt changes in main.cpp.\n"
       << "dt = 0.001\n"
       << "\n"
       << "# --- Lyapunov fit window (in STEPS). Tune by eye after first viewing. ---\n"
       << "# Exclude the initial transient and the saturation tail; keep the\n"
       << "# straight-line region in log(error) vs time.\n"
       << "fit_step_lo = 200\n"
       << "fit_step_hi = 1500\n"
       << "\n"
       << "f(x) = a + b*x\n"
       << "set fit quiet\n"
       << "fit [fit_step_lo*dt:fit_step_hi*dt] f(x) datafile \\\n"
       << "    using ($1*dt):(log($3)) via a, b\n"
       << "lambda = b\n"
       << "print sprintf('Lyapunov exponent lambda = %.5f  (1/time units)', lambda)\n"
       << "print sprintf('Fit window: steps %d-%d  (t = %.3f-%.3f)', \\\n"
       << "    fit_step_lo, fit_step_hi, fit_step_lo*dt, fit_step_hi*dt)\n"
       << "\n"
       << "set logscale y\n"
       << "set xlabel 'time (step * dt)'\n"
       << "set ylabel 'position error (log scale)'\n"
       << "set title sprintf('BH error vs direct baseline  -  lambda = %.4f', lambda)\n"
       << "set grid\n"
       << "plot datafile using ($1*dt):2 with lines title 'abs max', \\\n"
       << "     datafile using ($1*dt):3 with lines title 'abs sum', \\\n"
       << "     datafile using ($1*dt):4 with lines title 'rel sum', \\\n"
       << "     [fit_step_lo*dt:fit_step_hi*dt] exp(f(x)) with lines lw 2 lc rgb 'black' \\\n"
       << "         title sprintf('fit: exp(%.3f + %.4f t)', a, lambda)\n"
       << "pause -1 'press enter to exit'\n";
    std::cout << "[validate] wrote gnuplot script: " << script_path
              << " (run: gnuplot " << script_path << ")\n";
}

void save_final_positions(const Stars& s, const std::string& tag)
{
    std::string path = final_path(tag);
    std::ofstream f(path.c_str(), std::ios::out);
    f << std::setprecision(17);
    const std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i) {
        f << s.x[i] << " " << s.y[i] << "\n";
    }
    std::cout << "Saved final positions: " << path << "\n";
}

bool load_final_positions(Stars& s, const std::string& tag)
{
    std::string path = final_path(tag);
    std::ifstream f(path.c_str());
    if (!f) return false;

    s.x.clear();
    s.y.clear();
    double x, y;
    while (f >> x >> y) {
        s.x.push_back(x);
        s.y.push_back(y);
    }
    return true;
}

void compare_positions(const Stars& baseline, const Stars& other,
                       const std::string& label)
{
    const std::size_t N = baseline.x.size();
    if (other.x.size() != N) {
        std::cerr << "[validate] " << label
                  << ": size mismatch (" << other.x.size()
                  << " vs baseline " << N << "), skipping.\n";
        return;
    }

    double max_err = 0.0;
    double abs_err = 0.0;
    double rel_err = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        double dx = other.x[i] - baseline.x[i];
        double dy = other.y[i] - baseline.y[i];
        double err = std::sqrt(dx * dx + dy * dy);

        if (err > max_err) max_err = err;
        abs_err += err;

        double mag = std::sqrt(baseline.x[i] * baseline.x[i] +
                               baseline.y[i] * baseline.y[i]);
        if (mag > 0.0) rel_err += err / mag;
    }

    std::cout << "[validate] " << label << " vs direct-serial baseline:\n"
              << "    abs max error : " << max_err << "\n"
              << "    abs sum error : " << abs_err << "\n"
              << "    rel sum error : " << rel_err << "\n";
}
