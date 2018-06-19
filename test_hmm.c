#include "hmm.h"
#include "model.h"
#include <vector>
#include <math.h>
#include <fstream>
#include <string>

int main(int argc, char const *argv[])
{
	vector<pair<int, double> > result;
	HMM hmms[5];
	load_models( "modellist.txt", hmms, 5);
	ModelManager test_mgr;
	test_mgr.loadSeq(argv[2]);
	test_mgr.evaluate(hmms, result);

	ofstream outfile;
	string filename = argv[3];
	outfile.open(filename.c_str(), ios::out);
	for (int i = 0; i < result.size(); ++i) {
		outfile << "model_0" << result[i].first+1 << ".txt " << result[i].second << endl;
	}

	printf("%f\n", log(1.5) ); // make sure the math library is included
	return 0;
}
