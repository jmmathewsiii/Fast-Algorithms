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

#endif // !VALIDATE_H
