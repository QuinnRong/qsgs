#ifndef QSGS_H
#define QSGS_H

#include <string>
#include <vector>

struct Axis
{
    int x;
    int y;
    int z;
};

struct Params
{
    double p1;
    double p2;
    double p3;
};

class QSGS
{
public:
    // constructors
    QSGS(const int &dim);

    // member functions
    void core_only(const double &frac, const int &iter);
    /* get statistic for iter structures with given cdd */
    void core_grow(const double &frac, const double &cdd, const int &iter);
    /* get statistic for iter structures with given cdd, frac */
    void dump_statistic(const std::string&);
    /* dump statistic for all testes with different cdd, frac */
    void output(const int&, const std::string&);
    /* dump one structure */
    double get_norm_std() { return statistic[count - 1][4]; }
    double volume_farction();

    void QuartetStructureGenerationSet(const double &frac)
    {
        generate_core(frac);
    }

    void QuartetStructureGenerationSet(const double &frac, const double &cdd,
        const double &px=1, const double &py=1, const double &pz=1, const std::string &str="min")
    {
        aniso = {px, py, pz};
        get_prob(str);
        generate_core(cdd);
        QuartetStructureGrow(frac);
    }
    
    void special_serial(const double&);
    void special_parallel(const double&);

    // void character(int tot);

    // data members
private:
    const int NX, NY, NZ;
    double d1_6 = 0.02;
    double d7_18 = d1_6 / 2.;
    double d19_26 = d1_6 / 8.;
    const Params params{d1_6, d7_18, d19_26};
    Params aniso{1, 1, 1};
    std::vector<Axis> Delt = {
        { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}, {0, 0, -1},
        { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0},
        { 0, 1, 1}, { 0, 1,-1}, { 0,-1, 1}, { 0,-1,-1},
        { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
        { 1, 1, 1}, {-1, 1, 1}, {-1,-1, 1}, { 1,-1, 1}, { 1, 1,-1}, {-1, 1,-1}, {-1,-1,-1}, { 1,-1,-1}};
    double prob[26];

    void generate_core(const double&);
    std::vector<std::vector<std::vector<int>>> arrgrid;
    int numsoild;
    std::vector<std::vector<int>> soild;

    void QuartetStructureGrow(const double&);
    void QuartetStructureSingle(const Axis&, const Axis&, const double&);
    bool WithinCell(const Axis&);
    void RoundBoundary(Axis&);
    void get_prob(const std::string &);
    double get_prob(const Axis&, const std::string &);

    double vf;

    void standard_dev();
    std::vector<double> vf_layer;
    double mean, std;

    // different parameters with the same dim;
    int count = 0;
    std::vector<std::vector<double>> statistic;
    /* frac, dim, mean, std, std/mean */
};

void std_map(const int &max_dim, const double &max_frac, const int &num=100);

void std_map(const int &max_dim, const double &max_frac, const double &cdd, const int &num=10);

// void std_curve(const double &max_frac, const double &threshold)
// {
//     int min_dim = 1;
//     for (double cdd = max_frac; cdd >= 0.049; cdd -= 0.05)
//     {
//         std::cout << "cdd = " << cdd <<  ", ";
//         for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
//         {
//             std::cout << dim << " ";
//             QSGS myq(dim);
//             myq.core_only(cdd, (dim > 80) ? 10 : 50);
//             if (myq.get_norm_std() < threshold)
//             {
//                 std::cout << "dim = " << dim << std::endl;
//                 myq.dump_statistic("core_only_curve.txt");
//                 min_dim = dim;
//                 break;
//             }
//         }
//     }
// }

// void std_curve(const double &cdd, const double &max_frac, const double &threshold)
// {
//     int min_dim = 10;
//     for (double frac = max_frac; frac >= 0.049; frac -= 0.05)
//     {
//         std::cout << "frac = " << frac <<  ", ";
//         for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
//         {
//             std::cout << dim << " ";
//             QSGS myq(dim);
//             myq.core_grow(cdd, frac, (dim > 80) ? 2 : 10);
//             if (myq.get_norm_std() < threshold)
//             {
//                 std::cout << "dim = " << dim << std::endl;
//                 myq.dump_statistic("core_grow_curve.txt");
//                 min_dim = dim;
//                 break;
//             }
//         }
//     }
// }

// void save_structure(const std::string &mode, const int &dim, const double &frac, const int &num)
// {
//     QSGS myq(dim);
//     for (int i = 0; i < num; ++i)
//     {
//         if (mode == "parallel")
//             myq.special_parallel(frac);
//         else if (mode == "serial")
//             myq.special_serial(frac);
//         else if (mode == "core")
//             myq.QuartetStructureGenerationSet(frac);
//         else
//         {
//             std::cout << "Error in structure mode!!!\n";
//             return;
//         }
//         std::cout << i << ": " << myq.volume_farction() << std::endl;
//         myq.output(i, mode + "-" + std::to_string(dim) + "-" + std::to_string(frac));
//     }
// }

// void save_structure(const std::string &mode, const int &dim, const double &frac,  const double &cdd, const int &num,
//     const double &px=1, const double &py=1, const double &pz=1, const std::string &str="min")
// {
//     QSGS myq(dim);
//     for (int i = 0; i < num; ++i)
//     {
//         if (mode == "iso")
//             myq.QuartetStructureGenerationSet(frac, cdd);
//         else if (mode == "aniso")
//             myq.QuartetStructureGenerationSet(frac, cdd, px, py, pz, str);
//         else
//         {
//             std::cout << "Error in structure mode!!!\n";
//             return;
//         }
//         std::cout << i << ": " << myq.volume_farction() << std::endl;
//         myq.output(i, mode + "-" + std::to_string(dim) + "-" + std::to_string(frac) + "-" + std::to_string(cdd));
//     }
// }

#endif