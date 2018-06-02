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
    void core_only_test(const double &cdd, const int &iter);
    /* get statistic for iter structures with given cdd */
    void core_grow_test(const double &cdd, const double &frac, const int &iter);
    /* get statistic for iter structures with given cdd, frac */
    void dump_statistic(const std::string&);
    /* dump statistic for all testes with different cdd, frac */
    void output(const int&, const std::string&);
    /* dump one structure */
    double get_norm_std() { return statistic[count - 1][4]; }
    double volume_farction();

    void QuartetStructureGenerationSet(const double &cdd)
    {
        generate_core(cdd);
    }

    void QuartetStructureGenerationSet(const double &cdd, const double &frac, const double &px=1, const double &py=1, const double &pz=1)
    {
        aniso = {px, py, pz};
        get_prob();
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
    Params aniso;
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
    void get_prob();
    double get_prob(const Axis&);

    double vf;

    void standard_dev();
    std::vector<double> vf_layer;
    double mean, std;


    // different frac;
    int count = 0;
    std::vector<std::vector<double>> statistic;
    /* frac, dim, mean, std, std/mean */
};

#endif