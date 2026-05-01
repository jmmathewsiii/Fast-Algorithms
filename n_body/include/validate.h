#ifndef VALIDATE_H
#define VALIDATE_H

#include "star.h"
#include <string>

void save_final_positions(const Stars& s, const std::string& tag);
bool load_final_positions(Stars& s, const std::string& tag);

// Prints absolute max error, total absolute error, and total relative error
// of `other` against `baseline`. `label` is a human-readable name for stdout.
void compare_positions(const Stars& baseline, const Stars& other,
                       const std::string& label);

// Per-step snapshot dump: writes output/<tag>_<step>.dat (zero-padded step).
void save_snapshot(const Stars& s, const std::string& tag, std::size_t step);

// Walks paired snapshots (baseline_tag vs other_tag) for steps
// {0, stride, 2*stride, ..., n_iter} and writes a CSV
// output/error_<label>.csv with columns: step,abs_max,abs_sum,rel_sum.
void compare_snapshots(const std::string& baseline_tag,
                       const std::string& other_tag,
                       std::size_t n_iter, std::size_t stride,
                       const std::string& label);

#endif // !VALIDATE_H
