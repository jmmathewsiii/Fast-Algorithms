#include "../include/plotter.h"
#include <fstream>
#include <iostream>

void Plotter::plot_positions(const Stars& stars, const std::string& name)
{
    const std::size_t N = stars.size();

    std::string data_filename = name + ".plt";
    std::string data_filepath = "./output/" + data_filename;

    std::ofstream f;
    f.open(data_filepath.c_str(), std::ios::out);

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

    std::string cmd_filepath = "./output/" + name;

    f.open(cmd_filepath.c_str(), std::ios::out);
    f << "set title '" << name << "' font ', 14'\n";
    f << "set tics font ', 14'\n";
    f << "set size square\n";
    f << "set xrange [" << x_min - x_pad << ":" << x_max + x_pad << "]\n";
    f << "set yrange [" << y_min - y_pad << ":" << y_max + y_pad << "]\n";
    f << "set xlabel 'x'\nset ylabel 'y'\n";
    f << "unset key\n";
    f << "plot '" << data_filename << "' with points pt 7 ps 0.5 lc rgb 'black'\n";

    f.close();

    std::cout << "Command file written: " << cmd_filepath << "\n";
}

Plotter::Animator Plotter::animator_begin(const std::string& name)
{
    Animator a;
    a.name = name;
    a.data.open(("./output/" + name + ".plt").c_str(), std::ios::out);
    return a;
}

void Plotter::animator_add_frame(Animator& a, const Stars& s)
{
    const std::size_t N = s.size();
    for (std::size_t i = 0; i < N; ++i) {
        a.data << s.x[i] << " " << s.y[i] << "\n";
    }
    a.data << "\n\n";
    ++a.n_frames;
}

void Plotter::animator_end(Animator& a, double half_range, double pause_sec)
{
    a.data.close();

    std::string cmd_filepath = "./output/" + a.name;
    std::ofstream f(cmd_filepath.c_str(), std::ios::out);
    f << "set size square\n";
    f << "set xrange [" << -half_range << ":" << half_range << "]\n";
    f << "set yrange [" << -half_range << ":" << half_range << "]\n";
    f << "set xlabel 'x'\nset ylabel 'y'\n";
    f << "unset key\n";
    f << "do for [i=0:" << a.n_frames - 1 << "] {\n";
    f << "    set title sprintf('" << a.name << "  frame %d', i)\n";
    f << "    plot '" << a.name << ".plt' index i "
         "with points pt 7 ps 0.5 lc rgb 'black'\n";
    f << "    pause " << pause_sec << "\n";
    f << "}\n";
    f.close();

    std::cout << "Animation data: ./output/" << a.name << ".plt"
              << " (" << a.n_frames << " frames)\n";
    std::cout << "Animation script: " << cmd_filepath << "\n";
}
