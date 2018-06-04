#include <iostream>
#include <string>
#include "QSGS.h"

using namespace std;

int main()
{
    // std_map(100, 0.5, 100);
    // std_map(100, 0.5, 0.01, 10);
    // std_map(100, 0.5, 0.02, 10);

    // std_curve(0.9, 0.1);
    // std_curve(0.9, 0.05);
    // std_curve(0.9, 0.02);

    // std_curve(0.9, 0.01, 0.1);
    // std_curve(0.9, 0.01, 0.05);
    // std_curve(0.9, 0.01, 0.02);

    // std_curve(0.9, 0.02, 0.1);
    // std_curve(0.9, 0.02, 0.05);
    // std_curve(0.9, 0.02, 0.02);

    // save_structure("core", 50, 0.25, 1);
    // save_structure("parallel", 20, 0.2, 1);
    // save_structure("serial", 20, 0.2, 1);

    // save_structure("iso", 50, 0.25, 0.01, 10);
    // save_structure("aniso", 50, 0.25, 0.01, 10, 5, 0.5, 0.5, "min");

    RandParam rp{0.005, 0.02, 15, 0.2};
    save_structure("rand", 50, 0.25, 100, rp, "min");
}
