#include <iostream>
#include <string>
#include "QSGS.h"

using namespace std;

int main()
{
    // std_map(100, 0.5, 100);
    // std_curve(0.9, 0.1);
    // std_curve(0.9, 0.05);
    // std_curve(0.9, 0.02);

    // std_map(100, 0.5, 0.01, 10);
    std_curve(0.9, 0.01, 0.1);
    std_curve(0.9, 0.01, 0.05);
    std_curve(0.9, 0.01, 0.02);

    // save_structure(1, 20, 0.1);
    // save_structure(1, 20, 0.4);
    // save_structure(1, 100, 0.1);
    // save_structure(1, 100, 0.4);

    // save_structure(1, 20, 0.02, 0.1);
    // save_structure(1, 20, 0.02, 0.4);
    // save_structure(1, 100, 0.01, 0.25);
    // save_structure(1, 100, 0.02, 0.4);

    // save_parallel(20, 0.2);
    // save_serial(20, 0.2);

    // save_structure(1, 100, 0.01, 0.25, 10, 0.5, 0.5, "min");
}