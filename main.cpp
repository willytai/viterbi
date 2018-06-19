#include <iostream>
#include "hmm.h"
#include "model.h"
using namespace std;

inline static bool myStr2Int(const string& str, int& num) {
   num = 0;
   size_t i = 0;
   int sign = 1;
   if (str[0] == '-') { sign = -1; i = 1; }
   bool valid = false;
   for (; i < str.size(); ++i) {
      if (isdigit(str[i])) {
         num *= 10;
         num += int(str[i] - '0');
         valid = true;
      }
      else return false;
   }
   num *= sign;
   return valid;
}

int main(int argc, char const *argv[])
{
	if (argc != 5) {
		cout << argc << endl;
		cerr << "Wrong input format!!" << endl;
		cout << "USAGE: ";
		cout << "./train <iteration> model_init.txt <seq_model_0X.txt> <model_0X.txt>" << endl;
		return 0;
	}

	HMM hmm_initial;
	loadHMM( &hmm_initial, "model_init.txt" );
	int itr;
	myStr2Int(argv[1], itr);
	ModelManager mgr(&hmm_initial, size_t(itr));
	mgr.loadSeq(argv[3]);
	mgr.train();

	FILE* fp = fopen(argv[4], "w");
	dumpHMM(fp, mgr.get_hmm());
	fclose(fp);

	return 0;
}