#include "../include/validate.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

static std::string final_path(const std::string& tag)
{
    return "./output/" + tag + "_final.dat";
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
