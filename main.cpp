#include "QSGS.h"
#include <iostream>
#include <string>
#include <io.h>

using namespace std;

void core_only_map(const int &max_dim, const double &max_cdd)
{
	for (int dim = 10; dim <= max_dim; dim += 10)
	{
		cout << "dim = " << dim << endl;
		QSGS myq(dim);
		for (double cdd = 0.01; cdd <= max_cdd; cdd += 0.01)
			myq.core_only_test(cdd, 100);
		myq.dump_statistic("core_only_map.txt");
	}	
}

void core_grow_map(const int &max_dim, const double &cdd, const double &max_frac)
{
	for (int dim = 10; dim <= max_dim; dim += 10)
	{
		cout << "dim = " << dim << endl;
		QSGS myq(dim);
		for (double frac = 0.02; frac <= max_frac; frac += 0.01)
			myq.core_grow_test(cdd, frac, 10);
		myq.dump_statistic("core_grow_map.txt");
	}	
}

void core_only_curve(const double &max_cdd, const double &threshold)
{
	int min_dim = 1;
	for (double cdd = max_cdd; cdd >= 0.049; cdd -= 0.05)
	{
		cout << "cdd = " << cdd <<  ", ";
		for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
		{
			cout << dim << " ";
			QSGS myq(dim);
			myq.core_only_test(cdd, (dim > 80) ? 10 : 50);
			if (myq.get_norm_std() < threshold)
			{
				cout << "dim = " << dim << endl;
				myq.dump_statistic("core_only_curve.txt");
				min_dim = dim;
				break;
			}
		}
	}
}

void core_grow_curve(const double &cdd, const double &max_frac, const double &threshold)
{
	int min_dim = 10;
	for (double frac = max_frac; frac >= 0.049; frac -= 0.05)
	{
		cout << "frac = " << frac <<  ", ";
		for (int dim = min_dim; dim <= 200; dim += (min_dim / 20 > 1) ? min_dim / 20 : 1)
		{
			cout << dim << " ";
			QSGS myq(dim);
			myq.core_grow_test(cdd, frac, (dim > 80) ? 2 : 10);
			if (myq.get_norm_std() < threshold)
			{
				cout << "dim = " << dim << endl;
				myq.dump_statistic("core_grow_curve.txt");
				min_dim = dim;
				break;
			}
		}
	}
}

void save_structure(const int &num, const int &dim, const double &cdd)
{
	QSGS myq(dim);
	for (int i = 0; i < num; ++i)
	{
		myq.QuartetStructureGenerationSet(cdd);
		cout << i << ": " << myq.volume_farction() << endl;
		myq.output(i, to_string(dim) + "-" + to_string(cdd));
	}
}

void save_structure(const int &num, const int &dim, const double &cdd, const double &frac)
{
	QSGS myq(dim);
	for (int i = 0; i < num; ++i)
	{
		myq.QuartetStructureGenerationSet(cdd, frac);
		cout << i << ": " << myq.volume_farction() << endl;
		myq.output(i, to_string(dim) + "-" + to_string(frac) + "-" + to_string(cdd));
	}
}

void save_parallel(const int &dim, const double &ratio)
{
	QSGS myq(dim);
	myq.special_parallel(ratio);
	cout << 0 << ": " << myq.volume_farction() << endl;
	myq.output(0, to_string(dim) + "-" + to_string(ratio) + "-parallel");
}

void save_serial(const int &dim, const double &ratio)
{
	QSGS myq(dim);
	myq.special_serial(ratio);
	cout << 0 << ": " << myq.volume_farction() << endl;
	myq.output(0, to_string(dim) + "-" + to_string(ratio) + "-serial");
}

void test_rand(int n)
{
    float sum = 0, max = 0, min = 1;
    for (int i = 0; i < n; ++i)
    {
        float x = rand() / float(RAND_MAX);
        sum += x;
        if (x > max) max = x;
        if (x < min) min = x;
        cout << x << endl;
    }
    cout << "average = " << sum / n << endl;
    cout << "max = " << max << ", min = " << min << endl;
}

void save_aniso(const int &num, const int &dim, const double &cdd, const double &frac, const double &px=1, const double &py=1, const double &pz=1)
{
    // mkdir
    string path = "structure\\" + to_string(dim) + "-" + to_string(frac) + "-" + to_string(cdd) + "-aniso";
    if (_access(path.c_str(), 0) == -1)
    {
        string cmd = "md " + path;
        cout << cmd << endl;
        system(cmd.c_str());
    }

    // QSGS myq(dim);
    // for (int i = 0; i < num; ++i)
    // {
    //     myq.QuartetStructureGenerationSet(cdd, frac);
    //     cout << i << ": " << myq.volume_farction() << endl;
    //     myq.output(i, to_string(dim) + "-" + to_string(frac) + "-" + to_string(cdd));
    // }
}

int main()
{
	// core_only_map(100, 0.5);
	// core_only_curve(0.9, 0.1);
	// core_only_curve(0.9, 0.05);
	// core_only_curve(0.9, 0.02);

	// core_grow_map(100, 0.02, 0.5);
	// core_grow_curve(0.02, 0.9, 0.1);
	// core_grow_curve(0.02, 0.9, 0.05);
	// core_grow_curve(0.02, 0.9, 0.02);

	// save_structure(1, 20, 0.1);
	// save_structure(1, 20, 0.4);
	// save_structure(1, 100, 0.1);
	// save_structure(1, 100, 0.4);

	// save_structure(1, 20, 0.02, 0.1);
	// save_structure(1, 20, 0.02, 0.4);
	save_structure(1, 100, 0.01, 0.25);
	// save_structure(1, 100, 0.02, 0.4);

	// save_parallel(20, 0.2);
	// save_serial(20, 0.2);
}