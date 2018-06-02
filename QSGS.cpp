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
}

void QSGS::QuartetStructureGrow(const double &frac)
{
    int i, j, k;
    double d1_6 = 0.02;
    double d7_18 = d1_6 / 2.;
    double d19_26 = d1_6 / 8.;
    double numtotal_need = frac * NX * NY * NZ;
    int Tnumsoild;
    srand((unsigned)time(NULL));

    Tnumsoild = numsoild;   //第一步生长出的核的数目

    //第二步，从固相核向周围26个方向生长
    while (Tnumsoild < numtotal_need)
    {
        for (int index_soild = 0; index_soild < Tnumsoild; index_soild++)
        {
            int index_i = soild[index_soild][0];
            int index_j = soild[index_soild][1];
            int index_k = soild[index_soild][2];

            //1向右方向生长
            if (index_j < NY - 1)
            {
                i = index_i;
                j = index_j + 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }


            //2向后方向生长
            if (index_i > 0)
            {
                i = index_i - 1;
                j = index_j;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }


            //3向左方向生长
            if (index_j > 0)
            {
                i = index_i;
                j = index_j - 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }


            //4向前方向生长
            if (index_i < NX - 1)
            {
                i = index_i + 1;
                j = index_j;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //5向上方向生长
            if (index_k < NZ - 1)
            {
                i = index_i;
                j = index_j;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //6向下方向生长       
            if (index_k > 0)
            {
                i = index_i;
                j = index_j;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d1_6)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //7向水平方向右前生长
            if (index_i < NX - 1 && index_j < NY - 1)
            {
                i = index_i + 1;
                j = index_j + 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //8向水平方向左前生长
            if (index_i < NX - 1 && index_j > 0)
            {
                i = index_i + 1;
                j = index_j - 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //9向水平方向右后生长
            if (index_i > 0 && index_j < NY - 1)
            {
                i = index_i - 1;
                j = index_j + 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //10向水平方向左后生长
            if (index_i > 0 && index_j > 0)
            {
                i = index_i - 1;
                j = index_j - 1;
                k = index_k;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //11向右上方向生长
            if (index_j < NY - 1 && index_k < NZ - 1)
            {
                i = index_i;
                j = index_j + 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //12向右下方向生长
            if (index_j <NY - 1 && index_k >0)
            {
                i = index_i;
                j = index_j + 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //13向左上方向生长
            if (index_j > 0 && index_k < NZ - 1)
            {
                i = index_i;
                j = index_j - 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //14向左下方向生长
            if (index_j > 0 && index_k > 0)
            {
                i = index_i;
                j = index_j - 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //15向前上方向生长
            if (index_i < NX - 1 && index_k < NZ - 1)
            {
                i = index_i + 1;
                j = index_j;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //16向前下方向生长
            if (index_i <NX - 1 && index_k >0)
            {
                i = index_i + 1;
                j = index_j;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //17向后上方向生长
            if (index_i > 0 && index_k < NZ - 1)
            {
                i = index_i - 1;
                j = index_j;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //18向后下方向生长
            if (index_i > 0 && index_k > 0)
            {
                i = index_i - 1;
                j = index_j;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d7_18)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //19向右前上对角线方向生长
            if (index_i < NX - 1 && index_j < NY - 1 && index_k < NZ - 1)
            {
                i = index_i + 1;
                j = index_j + 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //20向右后上对角线方向生长
            if (index_i > 0 && index_j < NY - 1 && index_k < NZ - 1)
            {
                i = index_i - 1;
                j = index_j + 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //21向左后上对角线方向生长
            if (index_i > 0 && index_j > 0 && index_k < NZ - 1)
            {
                i = index_i - 1;
                j = index_j - 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //22向左前上对角线方向生长
            if (index_i<NX - 1 && index_j>0 && index_k < NZ - 1)
            {
                i = index_i + 1;
                j = index_j - 1;
                k = index_k + 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //23向右前下对角线方向生长
            if (index_i < NX - 1 && index_j<NY - 1 && index_k>0)
            {
                i = index_i + 1;
                j = index_j + 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //24向右后下对角线方向生长
            if (index_i > 0 && index_j<NY - 1 && index_k>0)
            {
                i = index_i - 1;
                j = index_j + 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //25向左后下对角线方向生长
            if (index_i > 0 && index_j > 0 && index_k > 0)
            {
                i = index_i - 1;
                j = index_j - 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }

            //26向左前下对角线方向生长
            if (index_i<NX - 1 && index_j>0 && index_k > 0)
            {
                i = index_i + 1;
                j = index_j - 1;
                k = index_k - 1;
                if (arrgrid[i][j][k] == 0 && ((rand() % 1000) / 1000.0) < d19_26)
                {
                    arrgrid[i][j][k] = 1;
                    soild[numsoild][0] = i;
                    soild[numsoild][1] = j;
                    soild[numsoild][2] = k;
                    ++numsoild;
                }
            }
        }
        Tnumsoild = numsoild;
    }
}
