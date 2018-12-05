#include <iostream>
#include <string>
#include <fstream>
#include <time.h>

#include "QSGS.h"
#include "utility.h"

using namespace std;

double rand_range(double r1, double r2)
{
    double tmp = rand() * 1.0 / RAND_MAX;
    return r1 + tmp * (r2 -r1);
}

int main()
{
    string root = "../..";
    string output_dir  = root + SPT + "qsgs_3d" + SPT + "run_1_valid/output" + SPT;
    make_directory(output_dir);

    int resolution[3] = {100, 100, 100};
    double P2 = 0.35;
    double cd = 0.01;

    ofstream summary("./output/summary.txt");
    summary << "index  cd   an_x  an_y  an_z" << std::endl;
    srand(time(NULL));
    for (int idx = 0; idx < 500; ++idx)
    {
        double aniso[3] = {rand_range(1, 20), rand_range(1, 20), rand_range(1, 20)};
        aniso[rand() % 3] = 1;
        run_once(output_dir, idx, cd, P2, aniso, resolution);
        summary << idx << " " << cd << " " <<  aniso[0] << " " << aniso[1] << " " << aniso[2] << std::endl;
    }
    summary.close();

    return 0;
}
