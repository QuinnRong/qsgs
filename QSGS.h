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
    /* generate one structure */
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
    /* generate special structure */
    void special_serial(const double&);
    void special_parallel(const double&);
    /* dump one structure */
    void dump_structure(const int&, const std::string&);
    /* get statistic for iter structures */
    void repeat(const double &frac, const int &iter);
    void repeat(const double &frac, const double &cdd, const int &iter);
    /* dump statistic for those structures */
    void dump_statistic(const std::string&);
    /* calculate mean and std */
    void standard_dev();
    /* calculate vf and return */
    double volume_farction();
    /* return normalized std */
    double get_norm_std() { return statistic[4]; }
    // data members
private:
    const int NX, NY, NZ;
    int numsoild;
    double d1_6 = 0.02;
    double d7_18 = d1_6 / 2.;
    double d19_26 = d1_6 / 8.;
    Params params{d1_6, d7_18, d19_26}, aniso{1, 1, 1};
    double prob[26];
    std::vector<Axis> Delt = {
        { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}, {0, 0, -1},
        { 1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0},
        { 0, 1, 1}, { 0, 1,-1}, { 0,-1, 1}, { 0,-1,-1},
        { 1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},
        { 1, 1, 1}, {-1, 1, 1}, {-1,-1, 1}, { 1,-1, 1}, { 1, 1,-1}, {-1, 1,-1}, {-1,-1,-1}, { 1,-1,-1}};

    std::vector<std::vector<std::vector<int>>> arrgrid;
    std::vector<std::vector<int>> soild;
    std::vector<double> vf_layer;
    double vf, mean, std;
    /* frac, dim, mean, std, std/mean (about iter structures) */
    std::vector<double> statistic;
    
    void generate_core(const double&);
    void QuartetStructureGrow(const double&);
    void QuartetStructureSingle(const Axis&, const Axis&, const double&);
    bool WithinCell(const Axis&);
    void RoundBoundary(Axis&);
    void get_prob(const std::string &);
    double get_prob(const Axis&, const std::string &);
};

/* use repeat and dump_statistic */
void std_map(const int &max_dim, const double &max_frac, const int &num=100);
void std_map(const int &max_dim, const double &max_frac, const double &cdd, const int &num=10);
void std_curve(const double &max_frac, const double &threshold);
void std_curve(const double &max_frac, const double &cdd, const double &threshold);
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