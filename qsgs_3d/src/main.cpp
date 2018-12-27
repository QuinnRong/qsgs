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

void dataset_generate()
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
}

void post_process()
{
    
    int idx = 381;
    string input_file = "./output/run_1_valid/output/3D_";
    input_file += to_string(idx) + ".dat";
    QSGS q(input_file, 100);

    string output_dir = "./output/fig/";
    q.get_structure("", idx);
    q.get_section(output_dir, idx, "x", 0);
    q.get_section(output_dir, idx, "y", 0);
    q.get_section(output_dir, idx, "z", 0);
}

int main()
{
    dataset_generate();
    // post_process();
}
