#include <cstdio>
#include <vector>
#include <string>
#include "math.h"

#include "utility.h"

double mix_prob(const std::vector<double> &vd, const std::string &str)
{   /*
    return average or min of the vector of double
    */
    if (str == "ave")
    {
        double res = 0;
        for (auto val : vd) res += val;
        return res / vd.size();
    }
    else
    {
        double res = vd[0];
        for (auto val : vd) if (val < res) res = val;
        return res;
    }
}

void get_D(double an_x, double an_y, double an_z, double (&D)[26])
{   /*
    get directional growth probability
    it is better that min(an_x, an_y, an_z) == 1
    */
    double p_centre = 0.001;
    double p_edge   = 0.001 / 2;
    double p_corner = 0.001 / 8;

    for (int i = 0; i < 26; ++i)
    {
        D[i] = 0;
        std::vector<double> vd;
        int n = 0;
        if (offsets[i].x != 0)
        {
            ++n;
            vd.push_back(an_x);
        }
        if (offsets[i].y != 0)
        {
            ++n;
            vd.push_back(an_y);
        }
        if (offsets[i].z != 0)
        {
            ++n;
            vd.push_back(an_z);
        }

        D[i] = mix_prob(vd, "min");

        if (n == 1)
            D[i] *= p_centre;
        else if (n == 2)
            D[i] *= p_edge;
        else if (n == 3)
            D[i] *= p_corner;

        // printf("%2d,%2d,%2d -> %f\n", offsets[i].x, offsets[i].y, offsets[i].z, D[i]);
    }
}

void run_once(const std::string &output_file, double cd, double P2, const double (&aniso)[3], const int (&resolution)[3])
{
    QuartetParams params;
    params.cd = cd;
    params.P2 = P2;
    get_D(aniso[0], aniso[1], aniso[2], params.D);

    QSGS q(resolution[0], resolution[1], resolution[2]);
    q.generation(params);
    q.growth(params);

    q.get_fenics_input(output_file + ".dat");
    q.get_structure(output_file + ".txt");
    q.get_section(output_file, "x", 0);
    q.get_section(output_file, "y", 0);
    q.get_section(output_file, "z", 0);
}

void post_process(const std::string &input_file, int resolution)
{
    QSGS q(input_file, resolution);

    std::string output_file = "./output/post";

    q.get_structure(output_file + ".txt");
    q.get_section(output_file, "x", 0);
    q.get_section(output_file, "y", 0);
    q.get_section(output_file, "z", 0);
}

double get_std(double cd, double P2, int resolution, int iter)
{
    QuartetParams params;
    params.cd = cd;
    params.P2 = P2;
    get_D(1, 1, 1, params.D);
    
    double std = 0;
    QSGS q(resolution, resolution, resolution);
    for (int i = 0; i < iter; ++i)
    {
        q.generation(params);
        q.growth(params);
        double section_std = q.get_section_std(params);
        std += section_std * section_std;
    }
    return sqrt(std / iter);
}
