#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;

int main() {
    ifstream file_ref("../testing_answer.txt");
    ifstream file_test("result1.txt");

    string buf_ref;
    string buf_test;
    vector<string> ref, result;
    while (getline(file_ref, buf_ref) && getline(file_test, buf_test)) {
        if (buf_ref == "") {
            assert(!(buf_ref.compare(buf_test)));
            break;
        }
        ref.push_back(buf_ref);
        result.push_back(buf_test.substr(0, 12));
    }

    float err = 0;
    assert(ref.size() == result.size());
    for (int i = 0; i < ref.size(); ++i) {
        if (ref[i].compare(result[i])) ++err;
    }
    int num = ref.size() - err;
    err /= ref.size();
    float acc = 1 - err;
    cout << "accuracy: " << acc << endl;
    cout << "Number of correct answers: " << num << endl;

    ofstream outfile;
    outfile.open("acc.txt", ios::out);
    outfile << acc << endl;
}