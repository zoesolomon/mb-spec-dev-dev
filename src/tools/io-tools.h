#ifndef IO_TOOLS_H
#define IO_TOOLS_H

#include <vector>
#include <string>
#include <iostream>

#include "macros.h"

#include "struct-lammpstrj-frame.h"

namespace opt {

void print_error(size_t lineno, size_t error);

void calculate_rcell(size_t, const double*, double*);

bool match_name(const char*, const char*);

void print_atom(std::ostream&, const char*);

} // namespace opt

#endif // IO_TOOLS_H
