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

void dataset_generate(const string &output_dir, int num, int res)
{   /*
    generate num samples in output_dir
    */
    make_directory(output_dir);

    int resolution[3] = {res, res, res};
    double P2 = 0.35;
    double cd = 0.01;

    string log_file = output_dir + "/summary.txt";
    ofstream summary(log_file);
    summary << "index  cd   an_x  an_y  an_z" << std::endl;
    srand(time(NULL));
    for (int idx = 0; idx < num; ++idx)
    {
        double aniso[3] = {rand_range(1, 20), rand_range(1, 20), rand_range(1, 20)};
        // double aniso[3] = {1, 1, 1}; // for isotropic composites
        aniso[rand() % 3] = 1;
        string output_file = output_dir + "/3D_" + to_string(idx);
        run_once(output_file, cd, P2, aniso, resolution);
        summary << idx << " " << cd << " " <<  aniso[0] << " " << aniso[1] << " " << aniso[2] << std::endl;
    }
    summary.close();
}

int main(int argc, char *argv[])
{
    string out_path = "./output/";
    dataset_generate(out_path + argv[1], atoi(argv[2]), atoi(argv[3]));
    // post_process("./output/run_1_valid/3D_0.dat", 20);
    // cout << get_std(0.01, 0.35, 10, 100) << endl;   // cd, P2, resolution, iteration
}
