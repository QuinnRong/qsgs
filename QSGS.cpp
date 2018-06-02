#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <time.h>
#include "QSGS.h"

QSGS::QSGS(const int &dim): NX(dim), NY(dim), NZ(dim)
{
    std::vector<int> vi (NX, 0);
    std::vector<std::vector<int>> vvi (NY, vi);
    arrgrid = std::vector<std::vector<std::vector<int>>> (NZ, vvi);

    vi = std::vector<int> (3, 0);
    soild = std::vector<std::vector<int>> (NX * NY * NZ, vi);
}

void QSGS::core_only_test(const double &cdd, const int &iter)
{
    std::vector<double> result;
    result.push_back(cdd);
    result.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        generate_core(cdd);
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

void QSGS::core_grow_test(const double &cdd, const double &frac, const int &iter)
{
    std::vector<double> result;
    result.push_back(frac);
    result.push_back(NX);
    double mean_sum = 0, std_sum = 0;
    for (int i = 0; i < iter; ++i)
    {
        QuartetStructureGenerationSet(cdd, frac);
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
    FILE *out;
    out = fopen(str.c_str(),"a");
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

void QSGS::output(const int &n, const std::string &p)
{
    // mkdir
    std::string path = p + "\\" + std::to_string(n);
    std::string cmd = "md " + path;
    std::cout << cmd << std::endl;
    system(cmd.c_str());

    // save to file
    for (int k = 0; k < NZ; ++k)
    {
        std::string filename = path + "\\"+
            "3D_" + std::to_string(n) + "_" + std::to_string(k) + ".dat";
        // std::cout << filename << std::endl;
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
    path = p + "\\" + "character";
    cmd = "md " + path;
    std::cout << cmd << std::endl;
    system(cmd.c_str());

    // save to file
    std::string filename = path + "\\"+
        "character_" + std::to_string(n) + ".txt";
    std::ofstream out(filename);
    mean = 0;
    std = 0;
    standard_dev();
    out << "mean: " << mean << std::endl;
    out << "std:  " << std << std::endl;
    for (int k = 0; k < NZ; ++k)
    {
        out << vf_layer[k] << std::endl;
    }
    out.close();
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

double QSGS::get_prob(const Axis &delt)
{
    double p = 0;
    int n = 0;
    if (delt.x != 0)
    {
        n += 1;
        p += aniso.p1;
    }
    if (delt.y != 0)
    {
        n += 1;
        p += aniso.p2;
    }
    if (delt.z != 0)
    {
        n += 1;
        p += aniso.p3;
    }
    p /= n;
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

void QSGS::get_prob()
{
    for (int i = 0; i < 26; ++i)
    {
        prob[i] = get_prob(Delt[i]);
        std::cout << prob[i] << std::endl;
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
