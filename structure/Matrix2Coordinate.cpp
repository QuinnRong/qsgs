#include <iostream>
#include <cstdio>
#include <string>

using namespace std;

void Maxtrix2Coordinate(const string& name, const string& idx, const int& dim)
{
	string saveFile = name + "_" + idx + ".txt";
	
	freopen(saveFile.c_str(), "w", stdout);
	printf("\n%d atoms\n\n", dim*dim*dim);
	printf("2 atom types\n\n");
	printf("0 %d xlo xhi\n", dim);
	printf("0 %d ylo yhi\n", dim);
	printf("0 %d zlo zhi\n\n", dim);
	printf("Masses\n\n");
	printf("1 1\n");
	printf("2 2\n\n");
	printf("Atoms\n\n");
	fclose(stdout);

	int atomID = 0;
	for (int z = 0; z < dim; ++z)
	{
		string path = name + "/" + idx + "/" + "3D_" + idx + "_" + to_string(z) + ".dat";
		// cout << path << endl;
		freopen(path.c_str(), "r", stdin);
		freopen(saveFile.c_str(), "a", stdout);

		int type;
		for (int y = 0; y < dim; ++y)
		{
			for (int x = 0; x < dim; ++x)
			{
				cin >> type;
				printf("%d %d %d %d %d\n", ++atomID, type + 1, x, y, z);
			}
		}

		fclose(stdout);
		fclose(stdin);
	}
}

int main()
{
	string name = "rand-50-0.250000_batch0";
	int dim = 50;
	string idx;
	idx = "23";
	idx = "59";
	Maxtrix2Coordinate(name, idx, dim);
}