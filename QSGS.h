#ifndef QSGS_H
#define QSGS_H

#include <string>
#include <vector>

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

    void QuartetStructureGenerationSet(const double &cdd, const double &frac)
    {
        generate_core(cdd);
        QuartetStructureGrow(frac);
    }
    
    void special_serial(const double&);
    void special_parallel(const double&);

    // void character(int tot);

    // data members
private:
    const int NX, NY, NZ;

    void generate_core(const double&);
    std::vector<std::vector<std::vector<int>>> arrgrid;
    int numsoild;
    std::vector<std::vector<int>> soild;

    void QuartetStructureGrow(const double&);

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