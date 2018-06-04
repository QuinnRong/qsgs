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

struct RandParam
{
    double cdd1, cdd2;
    double anis, rate;
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
    /* dump one structure */
    void dump_structure(const int&, const std::string&);

    /* get statistic for iter structures */
    void repeat(const double &frac, const int &iter);
    void repeat(const double &frac, const double &cdd, const int &iter);
    /* dump statistic for those structures */
    void dump_statistic(const std::string&);
    /* return normalized std */
    double get_norm_std() { return statistic[4]; }

    /* calculate mean and std */
    void standard_dev();
    /* calculate vf and return */
    double volume_farction();

    /* generate special structure */
    void special_serial(const double&);
    void special_parallel(const double&);
    
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
/* use QuartetStructureGenerationSet and dump_structure */
void save_structure(const std::string &mode, const int &dim, const double &frac, const int &num);
void save_structure(const std::string &mode, const int &dim, const double &frac,  const double &cdd, const int &num,
    const double &px=1, const double &py=1, const double &pz=1, const std::string &str="min");
void save_structure(const std::string &mode, const int &dim, const double &frac, const int &num,
    const RandParam &rp, const std::string &str="min");

double mix_prob(const std::vector<double> &vd, const std::string &str);
void make_directory(const std::string &path);
void delete_directory(const std::string &path);
void delete_file(const std::string &path);

#endif