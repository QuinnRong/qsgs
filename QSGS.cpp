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
    for (int dim = 10; dim <= max_dim; dim += 10)
    {
        std::cout << "dim = " << dim << std::endl;
        QSGS myq(dim);
        for (double frac = 0.01; frac <= max_frac; frac += 0.01)
            myq.core_only(frac, num);
        myq.dump_statistic("core_only_map.txt");
    }   
}

void std_map(const int &max_dim, const double &max_frac, const double &cdd, const int &num)
{
    for (int dim = 10; dim <= max_dim; dim += 10)
    {
        std::cout << "dim = " << dim << std::endl;
        QSGS myq(dim);
        for (double frac = 0.02; frac <= max_frac; frac += 0.01)
            myq.core_grow(frac, cdd, num);
        myq.dump_statistic("core_grow_map.txt");
    }   
}

void QSGS::core_only(const double &frac, const int &iter)
{
    std::vector<double> result;
    result.push_back(frac);
    result.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        generate_core(frac);
        volume_farction();
        // std::cout << "volume fraction = " << vf << std::endl;
        standard_dev();
        // std::cout << "mean = " << mean << ", std = " << std << std::endl;
        mean_sum += mean;
        std_sum += std;
    }
    result.push_back(mean_sum / iter);
    result.push_back(std_sum / iter);
    result.push_back(std_sum / mean_sum);
    statistic.push_back(result);
    ++count;
}

void QSGS::core_grow(const double &frac, const double &cdd, const int &iter)
{
    std::vector<double> result;
    result.push_back(frac);
    result.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        QuartetStructureGenerationSet(frac, cdd);
        volume_farction();
        // std::cout << "volume fraction = " << vf << std::endl;
        standard_dev();
        // std::cout << "mean = " << mean << ", std = " << std << std::endl;
        mean_sum += mean;
        std_sum += std;
    }
    result.push_back(mean_sum / iter);
    result.push_back(std_sum / iter);
    result.push_back(std_sum / mean_sum);
    statistic.push_back(result);
    ++count;
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
    out = fopen((path + SIMBLE + str).c_str(), "w");
    for (int i = 0; i < statistic.size(); ++i)
    {
        fprintf(out, "%.3f\t%.0f\t%.6f\t%.6f\t%.6f\n", statistic[i][0], statistic[i][1],
            statistic[i][2], statistic[i][3], statistic[i][4]);
    }
    fclose(out);
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

void make_directory(const std::string &path)
{
    // file not exist
    if (access(path.c_str(), 0) == -1)
    {
        std::string cmd = "mkdir " + path;
        system(cmd.c_str());
    }
}

void delete_directory(const std::string &path)
{
    // file exist
    if (access(path.c_str(), 0) == 0)
    {
        std::string cmd = "rm -rf " + path;
        system(cmd.c_str());
    }
}

void QSGS::output(const int &n, const std::string &p)
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

void QSGS::QuartetStructureSingle(const Axis &a, const Axis &delt, const double &p)
{
    Axis b{a.x + delt.x, a.y + delt.y, a.z + delt.z};
    // if (WithinCell(b))
    {
        RoundBoundary(b);
        if (arrgrid[b.x][b.y][b.z] == 0 && ((rand() % 1000) / 1000.0) < p)
        {
            arrgrid[b.x][b.y][b.z] = 1;
            soild[numsoild][0] = b.x;
            soild[numsoild][1] = b.y;
            soild[numsoild][2] = b.z;
            ++numsoild;
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
