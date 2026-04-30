#include "../include/plotter.h"
#include "../include/quadtree_cpu.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {

// Zero-pad a frame index to 5 digits.
std::string pad(std::size_t frame)
{
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << frame;
    return oss.str();
}

std::string frame_path(const std::string& name, std::size_t frame, const std::string& tag)
{
    // tag is "" for plain frames, or "_tree" / "_particles" for tree-animator frames.
    return "./output/" + name + "_" + pad(frame) + tag + ".plt";
}

void write_particles(std::ofstream& f, const Stars& s)
{
    const std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i) {
        f << s.x[i] << " " << s.y[i] << "\n";
    }
}

}  // namespace

void Plotter::plot_positions(const Stars& stars, const std::string& name)
{
    const std::size_t N = stars.size();

    std::string data_filename = name + ".plt";
    std::string data_filepath = "./output/" + data_filename;

    std::ofstream f(data_filepath.c_str(), std::ios::out);

    double x_min =  1e300, x_max = -1e300;
    double y_min =  1e300, y_max = -1e300;

    for (std::size_t i = 0; i < N; ++i) {
        double x = stars.x[i];
        double y = stars.y[i];
        f << x << " " << y << "\n";

        if (x < x_min) x_min = x;
        if (x > x_max) x_max = x;
        if (y < y_min) y_min = y;
        if (y > y_max) y_max = y;
    }
    f.close();

    std::cout << "Data file written: " << data_filepath << "\n";

    double x_pad = 0.1 * (x_max - x_min);
    double y_pad = 0.1 * (y_max - y_min);

    // Script lives in ./scripts/; references the data file via ../output/.
    std::string cmd_filepath = "./scripts/" + name;
    f.open(cmd_filepath.c_str(), std::ios::out);
    f << "set title '" << name << "' font ', 14'\n";
    f << "set tics font ', 14'\n";
    f << "set size square\n";
    f << "set xrange [" << x_min - x_pad << ":" << x_max + x_pad << "]\n";
    f << "set yrange [" << y_min - y_pad << ":" << y_max + y_pad << "]\n";
    f << "set xlabel 'x'\nset ylabel 'y'\n";
    f << "unset key\n";
    f << "plot '../output/" << data_filename << "' with points pt 7 ps 0.5 lc rgb 'black'\n";
    f.close();

    std::cout << "Command file written: " << cmd_filepath << "\n";
}

Plotter::Animator Plotter::animator_begin(const std::string& name)
{
    Animator a;
    a.name = name;
    return a;
}

void Plotter::animator_add_frame(Animator& a, const Stars& s)
{
    std::ofstream f(frame_path(a.name, a.n_frames, "").c_str(), std::ios::out);
    write_particles(f, s);
    f.close();
    ++a.n_frames;
}

void Plotter::animator_end(Animator& a, double half_range, double pause_sec)
{
    std::string cmd_filepath = "./scripts/" + a.name;
    std::ofstream f(cmd_filepath.c_str(), std::ios::out);
    f << "set size square\n";
    f << "set xrange [" << -half_range << ":" << half_range << "]\n";
    f << "set yrange [" << -half_range << ":" << half_range << "]\n";
    f << "set xlabel 'x'\nset ylabel 'y'\n";
    f << "unset key\n";
    f << "do for [i=0:" << a.n_frames - 1 << "] {\n";
    f << "    set title sprintf('" << a.name << "  frame %d', i)\n";
    f << "    plot sprintf('../output/" << a.name << "_%05d.plt', i) "
         "with points pt 7 ps 0.5 lc rgb 'black'\n";
    f << "    pause " << pause_sec << "\n";
    f << "}\n";
    f.close();

    std::cout << "Animation: " << a.n_frames << " frames at ./output/"
              << a.name << "_*.plt; script: " << cmd_filepath << "\n";
}

Plotter::TreeAnimator Plotter::tree_animator_begin(const std::string& name)
{
    TreeAnimator a;
    a.name = name;
    return a;
}

void Plotter::tree_animator_add_frame(TreeAnimator& a, const QuadTree& tree, const Stars& s)
{
    {
        std::ofstream f(frame_path(a.name, a.n_frames, "_tree").c_str(), std::ios::out);
        tree.write_divisions(f);
    }
    {
        std::ofstream f(frame_path(a.name, a.n_frames, "_particles").c_str(), std::ios::out);
        write_particles(f, s);
    }
    ++a.n_frames;
}

void Plotter::tree_animator_end(TreeAnimator& a, double half_range, double pause_sec)
{
    std::string cmd_filepath = "./scripts/" + a.name;
    std::ofstream f(cmd_filepath.c_str(), std::ios::out);
    f << "set size square\n";
    f << "set xrange [" << -half_range << ":" << half_range << "]\n";
    f << "set yrange [" << -half_range << ":" << half_range << "]\n";
    f << "set xlabel 'x'\nset ylabel 'y'\n";
    f << "unset key\n";
    f << "do for [i=0:" << a.n_frames - 1 << "] {\n";
    f << "    set title sprintf('" << a.name << "  frame %d', i)\n";
    f << "    plot sprintf('../output/" << a.name << "_%05d_tree.plt', i) "
         "with lines lc rgb 'gray60', \\\n";
    f << "         sprintf('../output/" << a.name << "_%05d_particles.plt', i) "
         "with points pt 7 ps 0.5 lc rgb 'black'\n";
    f << "    pause " << pause_sec << "\n";
    f << "}\n";
    f.close();

    std::cout << "Tree animation: " << a.n_frames << " frames at ./output/"
              << a.name << "_*_{tree,particles}.plt; script: " << cmd_filepath << "\n";
}
