#ifndef UTILITY_H
#define UTILITY_H

#include "QSGS.h"

void get_D(double an_x, double an_y, double an_z, double (&D)[26]);

void run_once(const std::string &output_dir, int idx, double cd, double P2, const double (&aniso)[3], const int (&resolution)[3]);

#endif