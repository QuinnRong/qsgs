#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <unistd.h>
#include "QSGS.h"

#define SIMBLE "\\"
// #define SIMBLE "/"

QSGS::QSGS(const int &dim): NX(dim), NY(dim), NZ(dim)
{
    std::vector<int> vi (NX, 0);
    std::vector<std::vector<int>> vvi (NY, vi);
    arrgrid = std::vector<std::vector<std::vector<int>>> (NZ, vvi);

    vi = std::vector<int> (3, 0);
    soild = std::vector<std::vector<int>> (NX * NY * NZ, vi);
}

void std_map(const int &max_dim, const double &max_frac, const int &num)
{
    std::string filename = "core_only_map.txt";
    std::string path = "statistic";
    path += SIMBLE + filename;
    delete_file(path);
    for (int dim = 10; dim <= max_dim; dim += 10)
    {
        std::cout << "dim = " << dim << std::endl;
        QSGS myq(dim);
        for (double frac = 0.01; frac <= max_frac; frac += 0.01)
        {
            myq.repeat(frac, num);
            myq.dump_statistic(filename);
        }
    }   
}

void std_map(const int &max_dim, const double &max_frac, const double &cdd, const int &num)
{
    std::string filename = "core_grow_map_" + std::to_string(cdd) + ".txt";
    std::string path = "statistic";
    path += SIMBLE + filename;
    delete_file(path);
    for (int dim = 10; dim <= max_dim; dim += 10)
    {
        std::cout << "dim = " << dim << std::endl;
        QSGS myq(dim);
        for (double frac = 0.02; frac <= max_frac; frac += 0.01)
        {
            myq.repeat(frac, cdd, num);
            myq.dump_statistic(filename);
        }
    }   
}

void std_curve(const double &max_frac, const double &threshold)
{
    std::string filename = "core_only_curve.txt";
    int min_dim = 1;
    for (double frac = max_frac; frac >= 0.049; frac -= 0.05)
    {
        std::cout << "frac = " << frac <<  ", ";
        for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
        {
            std::cout << dim << " ";
            QSGS myq(dim);
            myq.repeat(frac, (dim > 80) ? 10 : 50);
            if (myq.get_norm_std() < threshold)
            {
                std::cout << "dim = " << dim << std::endl;
                myq.dump_statistic(filename);
                min_dim = dim;
                break;
            }
        }
    }
}

void std_curve(const double &max_frac, const double &cdd, const double &threshold)
{
    std::string filename = "core_grow_curve_" + std::to_string(cdd) + ".txt";
    int min_dim = 10;
    for (double frac = max_frac; frac >= 0.049; frac -= 0.05)
    {
        std::cout << "frac = " << frac <<  ", ";
        for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
        {
            std::cout << dim << " ";
            QSGS myq(dim);
            myq.repeat(frac, cdd, (dim > 80) ? 10 : 50);
            if (myq.get_norm_std() < threshold)
            {
                std::cout << "dim = " << dim << std::endl;
                myq.dump_statistic(filename);
                min_dim = dim;
                break;
            }
        }
    }
}

void save_structure(const std::string &mode, const int &dim, const double &frac, const int &num)
{
    QSGS myq(dim);
    for (int i = 0; i < num; ++i)
    {
        if (mode == "parallel")
            myq.special_parallel(frac);
        else if (mode == "serial")
            myq.special_serial(frac);
        else if (mode == "core")
            myq.QuartetStructureGenerationSet(frac);
        else
        {
            std::cout << "Error in structure mode!!!\n";
            return;
        }
        std::cout << i << ": " << myq.volume_farction() << std::endl;
        myq.dump_structure(i, mode + "-" + std::to_string(dim) + "-" + std::to_string(frac));
    }
}

void save_structure(const std::string &mode, const int &dim, const double &frac,  const double &cdd, const int &num,
    const double &px, const double &py, const double &pz, const std::string &str)
{
    QSGS myq(dim);
    for (int i = 0; i < num; ++i)
    {
        if (mode == "iso")
            myq.QuartetStructureGenerationSet(frac, cdd);
        else if (mode == "aniso")
            myq.QuartetStructureGenerationSet(frac, cdd, px, py, pz, str);
        else
        {
            std::cout << "Error in structure mode!!!\n";
            return;
        }
        std::cout << i << ": " << myq.volume_farction() << std::endl;
        myq.dump_structure(i, mode + "-" + std::to_string(dim) + "-" + std::to_string(frac) + "-" + std::to_string(cdd));
    }
}

double rand_range(double r1, double r2)
{
    double tmp = rand() * 1.0 / RAND_MAX;
    return r1 + tmp * (r2 -r1);
}

void get_aniso(const double &ani, const double &rat, double &px, double &py, double &pz)
{
    double res[3] = {1, 1, 1};
    res[rand() % 3] = ani;
    if (rand() * 1.0 / RAND_MAX > 0.5)
        res[rand() % 3] = ani;
    px = rat * res[0];
    py = rat * res[1];
    pz = rat * res[2];
}

void save_structure(const std::string &mode, const int &dim, const double &frac, const int &num,
    const RandParam &rp, const std::string &str)
{
    QSGS myq(dim);
    for (int i = 0; i < num; ++i)
    {
        if (mode == "rand")
        {
            double cdd = rand_range(rp.cdd1, rp.cdd2);
            std::cout << "cdd = " << cdd << std::endl;
            double ani = rand_range(1, rp.anis);
            std::cout << "ani = " << ani << std::endl;
            double rat = rand_range(rp.rate, 1);
            std::cout << "rat = " << rat << std::endl;
            double px, py, pz;
            get_aniso(ani, rat, px, py, pz);
            std::cout << "px = " << px << ", py = " << py << ", pz = " << pz << std::endl;
            myq.QuartetStructureGenerationSet(frac, cdd, px, py, pz, str);
        }
        else
        {
            std::cout << "Error in structure mode!!!\n";
            return;
        }
        std::cout << i << ": " << myq.volume_farction() << std::endl;
        myq.dump_structure(i, mode + "-" + std::to_string(dim) + "-" + std::to_string(frac));
    }
}

void QSGS::generate_core(const double &cdd)
{
    numsoild = 0;
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            for (int k = 0; k < NZ; ++k)
            {
                if ((rand() * 1.0 / RAND_MAX) <= cdd)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
                else
                {
                    arrgrid[i][j][k] = 0;
                }
            }
        }
    }
}

void QSGS::QuartetStructureGrow(const double &frac)
{
    int numtotal_need = frac * NX * NY * NZ;
    int num_acc = numsoild;   //第一步生长出的核的数目
    // srand((unsigned)time(NULL));

    while (num_acc < numtotal_need)
    {
        for (int index_soild = 0; index_soild < num_acc; index_soild++)
        {
            if (numsoild < numtotal_need)
            {
                Axis axis = {soild[index_soild][0], soild[index_soild][1], soild[index_soild][2]};
                for (int i = 0; i < Delt.size(); ++i)
                    QuartetStructureSingle(axis, Delt[i], prob[i]);
            }
        }
        num_acc = numsoild;
    }
}

void QSGS::QuartetStructureSingle(const Axis &a, const Axis &delt, const double &p)
{
    Axis b{a.x + delt.x, a.y + delt.y, a.z + delt.z};
    RoundBoundary(b);
    if (arrgrid[b.x][b.y][b.z] == 0 && (rand() * 1.0 / RAND_MAX) < p)
    {
        arrgrid[b.x][b.y][b.z] = 1;
        soild[numsoild][0] = b.x;
        soild[numsoild][1] = b.y;
        soild[numsoild][2] = b.z;
        ++numsoild;
    }
}

double mix_prob(const std::vector<double> &vd, const std::string &str)
{
    if (str == "ave")
    {
        double res = 0;
        for (auto val : vd) res += val;
        return res / vd.size();
    }
    else if (str == "min")
    {
        double res = vd[0];
        for (auto val : vd) if (val < res) res = val;
        return res;
    }
    else
    {
        std::cout << "Error in mix_prob!!!\n";
        return 0;
    }
}

double QSGS::get_prob(const Axis &delt, const std::string &str)
{
    double p = 0;
    std::vector<double> vd;
    int n = 0;
    if (delt.x != 0)
    {
        n += 1;
        vd.push_back(aniso.p1);
    }
    if (delt.y != 0)
    {
        n += 1;
        vd.push_back(aniso.p2);
    }
    if (delt.z != 0)
    {
        n += 1;
        vd.push_back(aniso.p3);
    }
    p = mix_prob(vd, str);
    if (n == 1)
        p *= params.p1;
    else if (n == 2)
        p *= params.p2;
    else if (n == 3)
        p *= params.p3;
    else
        std::cout << "Error in get_prob, n = " << n << std::endl;
    return p;
}

void QSGS::get_prob(const std::string &str)
{
    for (int i = 0; i < 26; ++i)
    {
        prob[i] = get_prob(Delt[i], str);
        // std::cout << prob[i] << std::endl;
    }
}

void QSGS::RoundBoundary(Axis &a)
{
    if (a.x < 0) a.x += NX;
    if (a.y < 0) a.y += NY;
    if (a.z < 0) a.z += NZ;
    if (a.x >= NX) a.x -= NX;
    if (a.y >= NY) a.y -= NY;
    if (a.z >= NZ) a.z -= NZ;
    if (!WithinCell(a))
    {
        std::cout << "Error in RoundBoundary!!!" << std::endl;
        printf("x = %d, y = %d, z = %d\n", a.x, a.y, a.z);
    }
}

bool QSGS::WithinCell(const Axis &a)
{
    if (a.x >= 0 && a.x <NX)
    {
        if (a.y >= 0 && a.y <NY)
        {
            if (a.z >= 0 && a.z <NZ)
            {
                return true;
            }
        }
    }
    return false;
}

void QSGS::repeat(const double &frac, const int &iter)
{
    statistic = std::vector<double>();
    statistic.push_back(frac);
    statistic.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        QuartetStructureGenerationSet(frac);
        standard_dev();
        mean_sum += mean;
        std_sum += std;
    }
    statistic.push_back(mean_sum / iter);
    statistic.push_back(std_sum / iter);
    statistic.push_back(std_sum / mean_sum);
}

void QSGS::repeat(const double &frac, const double &cdd, const int &iter)
{
    statistic = std::vector<double>();
    statistic.push_back(frac);
    statistic.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        QuartetStructureGenerationSet(frac, cdd);
        standard_dev();
        mean_sum += mean;
        std_sum += std;
    }
    statistic.push_back(mean_sum / iter);
    statistic.push_back(std_sum / iter);
    statistic.push_back(std_sum / mean_sum);
}

void QSGS::dump_statistic(const std::string &str)
{
    // mkdir
    std::string path = "statistic";
    if (access(path.c_str(), 0) == -1)
    {
        std::string cmd = "mkdir " + path;
        system(cmd.c_str());
    }

    FILE *out;
    out = fopen((path + SIMBLE + str).c_str(), "a");
    fprintf(out, "%.3f\t%.0f\t%.6f\t%.6f\t%.6f\n", statistic[0], statistic[1],
        statistic[2], statistic[3], statistic[4]);
    fclose(out);
}

void QSGS::dump_structure(const int &n, const std::string &p)
{
    // mkdir
    std::string root = "structure", path, filename;
    path = root;
    make_directory(path);
    path += SIMBLE + p;
    make_directory(path);
    path += SIMBLE + std::to_string(n);
    delete_directory(path);
    make_directory(path);
    // save to file
    for (int k = 0; k < NZ; ++k)
    {
        filename = path + SIMBLE +
            "3D_" + std::to_string(n) + "_" + std::to_string(k) + ".dat";
        std::ofstream out(filename);
        for (int j = 0; j < NY; ++j)
        {
            for (int i = 0; i < NX; ++i)
            {
                out << arrgrid[i][j][k] << " ";
            }
            out << std::endl;
        }
        out.close();
    }
    // mkdir
    path = root + SIMBLE + p + SIMBLE + "character";
    make_directory(path);
    // save to file
    filename = path + SIMBLE +
        "character_" + std::to_string(n) + ".txt";
    FILE *out;
    out = fopen(filename.c_str(),"w");
    standard_dev();
    fprintf(out, "%-8s\n%-8.3f%-8.3f%-8.3f\n", "aniso:", aniso.p1, aniso.p2, aniso.p3);
    fprintf(out, "%-8s%-8s%-8s\n", "mean:", "std:", "normstd:");
    fprintf(out, "%-8.3f%-8.3f%-8.3f\n\n", mean, std, std / mean);
    for (int k = 0; k < NZ; ++k)
        fprintf(out, "%f\n", vf_layer[k]);
    fclose(out);
}

double QSGS::volume_farction()
{
    double sum = 0;
    for (int k = 0; k < NZ; ++k)
    {
        for (int j = 0; j < NY; ++j)
        {
            for (int i = 0; i < NX; ++i)
            { 
                if (arrgrid[i][j][k] == 1)
                {
                   sum += 1;
                }
            }
        }
    }
    vf = sum / (NX * NY * NZ);
    return vf;
}

void QSGS::standard_dev()
{
    vf_layer = std::vector<double>();
    for (int k = 0; k < NZ; ++k)
    {
        double sum = 0;
        for (int j = 0; j < NY; ++j)
        {
            for (int i = 0; i < NX; ++i)
            {
                if ( arrgrid[i][j][k] == 1)
                {
                   sum += 1;
                }
            }
        }
        vf_layer.push_back(sum / (NX * NY));
    }
    mean = std::accumulate(std::begin(vf_layer), std::end(vf_layer), 0.0) / vf_layer.size();
    double accum = 0;
    std::for_each (std::begin(vf_layer), std::end(vf_layer), [&](const double d) { accum += (d-mean)*(d-mean); });
    std = sqrt(accum/(vf_layer.size()-1));
}

void QSGS::special_serial(const double &ratio)
{
    numsoild = 0;
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            for (int k = 0; k < NZ; ++k)
            {
                if (k < round(NZ * ratio))
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
                else
                {
                    arrgrid[i][j][k] = 0;
                }
            }
        }
    }
}

void QSGS::special_parallel(const double &ratio)
{
    numsoild = 0;
    for (int i = 0; i < NX; ++i)
    {
        for (int j = 0; j < NY; ++j)
        {
            for (int k = 0; k < NZ; ++k)
            {
                if (i < round(NX * ratio))
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
                else
                {
                    arrgrid[i][j][k] = 0;
                }
            }
        }
    }
}

void make_directory(const std::string &path)
{
    // path not exist
    if (access(path.c_str(), 0) == -1)
    {
        std::string cmd = "mkdir " + path;
        system(cmd.c_str());
    }
}

void delete_directory(const std::string &path)
{
    // path exist
    if (access(path.c_str(), 0) == 0)
    {
        std::string cmd = "rm -rf " + path;
        system(cmd.c_str());
    }
}

void delete_file(const std::string &path)
{
    // file exist
    if (access(path.c_str(), 0) == 0)
    {
        std::string cmd = "rm " + path;
        system(cmd.c_str());
    }
}
