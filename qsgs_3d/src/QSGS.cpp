#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <time.h>

#include "QSGS.h"

const Coords offsets[26] =                              // offsets of 26 neighbors
{
    { 0, 0, 1},                                         // centre block (top)
    { 0, 0,-1},                                         // centre block (bottom)
    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0},     // centre blocks (middle)
    { 1, 0, 1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1, 1},     // edge blocks (top)
    { 1, 0,-1}, {-1, 0,-1}, { 0, 1,-1}, { 0,-1,-1},     // edge blocks (bottom)
    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0},     // edge blocks (middle)
    { 1, 1, 1}, {-1,-1, 1}, { 1,-1, 1}, {-1, 1, 1},     // corner blocks (top)
    { 1, 1,-1}, {-1,-1,-1}, { 1,-1,-1}, {-1, 1,-1}      // corner blocks (bottom)
};

//==============================================================
// Constructor
//==============================================================
QSGS::QSGS(int dim_x, int dim_y, int dim_z): Nx(dim_x), Ny(dim_y), Nz(dim_z)
{
    grid3D = new int**[Nz];
    for (int z = 0; z < Nz; ++z)
    {
        grid3D[z] = new int*[Ny];
        for (int y = 0; y < Ny; ++y)
        {
            grid3D[z][y] = new int[Nx];
        }
    }

    coords = new Coords[Nx * Ny * Nz];
}

QSGS::QSGS(const std::string &file, int dim): QSGS(dim)
{
    std::ifstream in(file);
    int val;

    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                in >> val;
                // std::cout << val << " ";
                grid3D[z][y][x] = val;
            }
        }
    }

    in.close();
}

//==============================================================
// Destructor
//==============================================================
QSGS::~QSGS()
{
    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            delete[] grid3D[z][y];
        }
        delete[] grid3D[z];
    }
    delete[] grid3D;

    delete[] coords;
}

//==============================================================
// generation process
//==============================================================
void QSGS::generation(const QuartetParams &params)
{
    srand(time(NULL));
    g_count = 0;            // reset g_count of 1-cells
    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                if ((rand() * 1.0 / RAND_MAX) <= params.cd)
                {
                    grid3D[z][y][x] = 1;        // set flags to 1 with probability cd
                    coords[g_count].x = x;      // record 1-cells' coordinates
                    coords[g_count].y = y;
                    coords[g_count].z = z;
                    ++g_count;                  // increase number of 1-cells
                }
                else
                {
                    grid3D[z][y][x] = 0;        // set flags to 0 with probability 1 - cd
                }
            }
        }
    }
}

//==============================================================
// growth process
//==============================================================
void QSGS::growth(const QuartetParams &params)
{
    srand(time(NULL));
    int total = params.P2 * Nx * Ny * Nz;   // total number of 1-cells needed
    int prev = g_count;
    while (prev < total)
    {
        for (int id = 0; id < prev; ++id)
        {
            if (g_count < total)
            {
                growth_neighbor(coords[id], params.D);
            }
        }
        prev = g_count;
    }
}

void QSGS::growth_neighbor(const Coords &center, const double (&D)[26])
{
    for (int i = 0; i < 26; ++i)
    {
        Coords neighbor{center.x + offsets[i].x, center.y + offsets[i].y, center.z + offsets[i].z};
        round_boundary(neighbor);
        if (grid3D[neighbor.z][neighbor.y][neighbor.x] == 0 && (rand() * 1.0 / RAND_MAX) < D[i])
        {
            grid3D[neighbor.z][neighbor.y][neighbor.x] = 1;     // set flags to 1 with probability D[i]
            coords[g_count].x = neighbor.x;                     // record 1-cells' coordinates
            coords[g_count].y = neighbor.y;
            coords[g_count].z = neighbor.z;
            ++g_count;                                          // increase number of 1-cells
        }
    }
}

void QSGS::round_boundary(Coords &a)
{
    if (a.x < 0) a.x += Nx;
    if (a.y < 0) a.y += Ny;
    if (a.z < 0) a.z += Nz;
    if (a.x >= Nx) a.x -= Nx;
    if (a.y >= Ny) a.y -= Ny;
    if (a.z >= Nz) a.z -= Nz;
    if (!in_cell(a))
    {
        std::cout << "Error in round_boundary" << std::endl;
        printf("x = %d, y = %d, z = %d\n", a.x, a.y, a.z);
    }
}

bool QSGS::in_cell(const Coords &a)
{
    if (a.x >= 0 && a.x <Nx)
    {
        if (a.y >= 0 && a.y <Ny)
        {
            if (a.z >= 0 && a.z <Nz)
            {
                return true;
            }
        }
    }
    return false;
}

//==============================================================
// save to file
//==============================================================
void SaveFile(const std::string& file, int* type, int dimX, int dimY, int dimZ)
{
    FILE *out = fopen(file.c_str(), "w");
    // std::cout << file << std::endl;

    fprintf(out, "\n%d atoms\n\n", dimX*dimY*dimZ);
    fprintf(out, "2 atom types\n\n");
    fprintf(out, "0 %d xlo xhi\n", dimX);
    fprintf(out, "0 %d ylo yhi\n", dimY);
    fprintf(out, "0 %d zlo zhi\n\n", dimZ);
    fprintf(out, "Masses\n\n");
    fprintf(out, "1 1\n");
    fprintf(out, "2 2\n\n");
    fprintf(out, "Atoms\n\n");

    int atomID = 0;
    for (int z = 0; z < dimZ; ++z)
    {
        for (int y = 0; y < dimY; ++y)
        {
            for (int x = 0; x < dimX; ++x)
            {
                fprintf(out, "%d %d %.1f %.1f %.1f\n", atomID + 1, type[atomID] + 1, x + 0.5, y + 0.5, z + 0.5);
                ++atomID;
            }
        }
    }

    fclose(out);
}

void QSGS::get_section(const std::string &root, int idx, const std::string &dire, int dist)
{
    std::string file = root;
    file += std::to_string(idx) + "_" + dire + "_" + std::to_string(dist) + ".txt";

    if (dire == "x")
    {
        if (dist < 0 || dist >= Nx)
        {
            std::cout << "invalid x offset" << std::endl;
            return;
        }
        int* type = new int[Ny*Nz];
        int atomID = 0;
        for (int z = 0; z < Nz; ++z)
        {
            for (int y = 0; y < Ny; ++y)
            {
                type[atomID] = grid3D[z][y][dist];
                ++atomID;
            }
        }
        SaveFile(file, type, 1, Ny, Nz);
        delete[] type;
    }
    else if (dire == "y")
    {
        if (dist < 0 || dist >= Ny)
        {
            std::cout << "invalid y offset" << std::endl;
            return;
        }
        int* type = new int[Nx*Nz];
        int atomID = 0;
        for (int z = 0; z < Nz; ++z)
        {
            for (int x = 0; x < Nx; ++x)
            {
                type[atomID] = grid3D[z][dist][x];
                ++atomID;
            }
        }
        SaveFile(file, type, Nx, 1, Nz);
        delete[] type;
    }
    else if (dire == "z")
    {
        if (dist < 0 || dist >= Nz)
        {
            std::cout << "invalid z offset" << std::endl;
            return;
        }
        int* type = new int[Nx*Ny];
        int atomID = 0;
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                type[atomID] = grid3D[dist][y][x];
                ++atomID;
            }
        }
        SaveFile(file, type, Nx, Ny, 1);
        delete[] type;
    }
}

void QSGS::get_structure(const std::string &root, int idx)
{
    std::string file = root + std::to_string(idx) + ".txt";;
    int* type = new int[Nx*Ny*Nz];

    int atomID = 0;
    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                type[atomID] = grid3D[z][y][x];
                // std::cout << type[atomID];
                ++atomID;
            }
            // std::cout << std::endl;
        }
        // std::cout << std::endl;
    }
    SaveFile(file, type, Nx, Ny, Nz);

    delete[] type;
}

void QSGS::get_fenics_input(const std::string& root, int idx)
{
    // save to file
    std::string filename = root + "3D_" + std::to_string(idx) + ".dat";
    std::ofstream out(filename);
    for (int z = 0; z < Nz; ++z)
    {
        for (int y = 0; y < Ny; ++y)
        {
            for (int x = 0; x < Nx; ++x)
            {
                out << grid3D[z][y][x] << " ";
            }
            out << std::endl;
        }
    }
    out.close();
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
