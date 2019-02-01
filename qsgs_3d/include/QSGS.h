#ifndef QSGS_H
#define QSGS_H

#include <string>

struct QuartetParams
{
    double cd;          // core distribution probability
    double P2;          // volume fraction of the second phase
    double D[26];       // directional growth probability
    // double Im_n;     // phase interaction growth probability, not used
};

struct Coords
{
    int x, y, z;
};

extern const Coords offsets[26];

class QSGS
{
public:
    // constructor and destructor
    QSGS(int dim): QSGS(dim, dim, dim) {};
    QSGS(int dim_x, int dim_y, int dim_z);
    QSGS(const std::string &file, int dim);
    ~QSGS();

    // public member functions
    void generation(const QuartetParams &params);
    void growth(const QuartetParams &params);

    // save as ovito form
    void get_section(const std::string &filename, const std::string &dire, int dist);
    void get_structure(const std::string &filename);
    // save as fenics form
    void get_fenics_input(const std::string &filename);
    // compute std
    double get_section_std(const QuartetParams &params);

private:
    // data members
    const int Nx, Ny, Nz;   // mesh size in three directions
    int ***grid3D;          // Nz*Ny*Nx, cell flags are 0 or 1
    int g_count;            // number of cells with flag 1
    Coords *coords;         // maximum Nx*Ny*Nz, coordinates of cells with flag 1

    // private member functions
    void growth_neighbor(const Coords &center, const double (&D)[26]);
    void round_boundary(Coords &a);
    bool in_cell(const Coords &a);
};

void make_directory(const std::string &path);
void delete_directory(const std::string &path);

#endif