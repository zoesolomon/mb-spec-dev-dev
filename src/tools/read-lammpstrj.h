#ifndef READ_LAMMPSTRJ_H
#define READ_LAMMPSTRJ_H

#include <vector>
#include <string>
#include <iostream>

#include "macros.h"

#include "struct-lammpstrj-frame.h"

namespace opt {

void read_lammpstrj_frame(size_t&, size_t&, std::istream&, opt::lammpstrj_frame&, double);

size_t read_lammpstrj(const char *, std::vector<opt::lammpstrj_frame>&, size_t, double);

void read_xyz_frame(size_t&, size_t&, std::istream&, opt::lammpstrj_frame&, double, std::vector<double>);

} // namespace opt

#endif // READ_HIST_H
