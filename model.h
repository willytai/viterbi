#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include "hmm.h"
using namespace std;

class ModelManager{
public:
    ModelManager() {};
    ModelManager(HMM* h, size_t it) { _hmm = h; _itr = it; }
    ~ModelManager() {}

    void loadSeq(string filename) {
        ifstream file(filename.c_str());
        string buf;
        while(getline(file, buf)) {
            if (buf == "") break;
            vector<int> tmp;
            for (int i = 0; i < buf.length(); ++i) {
                switch(buf[i]) {
                    case 'A':
                        tmp.push_back(0);
                        break;
                    case 'B':
                        tmp.push_back(1);
                        break;
                    case 'C':
                        tmp.push_back(2);
                        break;
                    case 'D':
                        tmp.push_back(3);
                        break;
                    case 'E':
                        tmp.push_back(4);
                        break;
                    case 'F':
                        tmp.push_back(5);
                        break;
                    default:
                        cerr << "Sequence file contains illegal characters" << endl;
                        assert(false);
                        break;
                }
            }
            seq.push_back(tmp);
        }
        _modelName = "model_0";
        _modelName += filename[14];
    }

    void cal_alpha(const vector<int>& seq) { // call this every time you get a new sequence
        // initialize alpha
        for (int state = 0; state < 6; ++state) {
            _alpha[state][0] = _hmm->initial[state] * obsProb(state, seq[0]);
        }
        
        // recursive
        int t = 1; // ith character in the sequence
        while (t < 50) {
            for (int state = 0; state < 6; ++state) {
                _alpha[state][t] = 0;
                for (int l = 0; l < 6; ++l)
                    _alpha[state][t] += _alpha[l][t-1] * transition(l, state);
                _alpha[state][t] *= obsProb(state, seq[t]);
            }
            ++t;
        }
    }

    void cal_beta(const vector<int>& seq) { // call this every time you get a new sequence
        // initialize beta
        for (int state = 0; state < 6; ++state)
            _beta[state][49] = 1;

        // recursive
        int t = 48; // ith character in the sequence
        while (t >= 0) {
            for (int state = 0; state < 6; ++state) {
                _beta[state][t] = 0;
                for (int N = 0; N < 6; ++N) // N is the number of states
                    _beta[state][t] += transition(state, N) * obsProb(N, seq[t+1]) * _beta[N][t+1];
            }
            --t;
        }
    }

    void cal_gamma(const vector<int>& seq) { // accumulate throught all of the sequences
        // given model, the probability of observing sequence
        double denominator = 0;
        for (int state = 0; state < 6; ++state) {
            denominator += _alpha[state][49];
        }
        assert (denominator);
        
        for (int state = 0; state < 6; ++state) {
            for (int t = 0; t < 50; ++t) {
                _gamma[state][t] += _alpha[state][t]*_beta[state][t] / denominator;
                // accumulate _update_obs
                _update_obs[seq[t]][state] += _alpha[state][t]*_beta[state][t] / denominator;
            }
        }
    }

    void cal_epsilon(const vector<int>& seq) { // accumulate through all of the sequences
        // given model, the probability of observing sequence
        double denominator = 0;
        for (int state = 0; state < 6; ++state)
            denominator += _alpha[state][49];

        assert (denominator);

        for (int t = 0; t < 49; ++t) {
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j)
                    _epsilon[t][i][j] += _alpha[i][t]*transition(i, j)*obsProb(j ,seq[t+1])*_beta[j][t+1] / denominator;
            }
        }
    }

    void train() { cout << "training " << _modelName << "......" << endl;
        // cout << endl << "initial" << endl;
        // dumpHMM(stderr, _hmm);
        for (int iteration = 0; iteration < _itr; ++iteration) { cout << '\r' << "epoch " << iteration+1 << '/' << _itr << flush;
            init_param(); // only gamma, epsilon, _update_obs need initialization
            for (int seq_num = 0; seq_num < seq.size(); ++seq_num) {
                assert (seq[seq_num].size() == 50);
                cal_alpha(seq[seq_num]);
                cal_beta(seq[seq_num]);
                cal_gamma(seq[seq_num]);
                cal_epsilon(seq[seq_num]);
            }
            update();
            // dumpHMM(stderr, _hmm);
        }
        cout << endl << endl;
        // dumpHMM(stderr, _hmm);
    }

    void test() {
        init_param();
        cal_alpha(seq[0]);
        cal_beta(seq[0]);
        cal_gamma(seq[0]);
        cal_epsilon(seq[0]);
        cout << "===============alpha===============" << endl;
        for (int i = 0; i < 6; ++i) {
            for (int k = 0; k < 10; ++k) {
                cout << _alpha[i][k] << ' ';
            } cout << endl;
        }
        
        cout << "===============beta===============" << endl;
        for (int i = 0; i < 6; ++i) {
            for (int k = 0; k < 10; ++k) {
                cout << _beta[i][k] << ' ';
            } cout << endl;
        }

        cout << "===============gamma===============" << endl;
        for (int i = 0; i < 6; ++i) {
            for (int k = 0; k < 10; ++k) {
                cout << _gamma[i][k] << ' ';
            } cout << endl;
        }
    }

    void evaluate(const HMM* hmms, vector<pair<int, double> >& result) {
        cout << "evaluating" << endl;

        // evaluate every seqs
        for (int i = 0; i < seq.size(); ++i) { cout << "evaluating seq "  << i+1 << endl;
            double score[5];
            // caculate the max prob for the five hmms
            for (int j = 0; j < 5; ++j) {
                score[j] = viterbi(&hmms[j], seq[i]);
            }
            // get the max value in socore
            pair<int, double> max(0,score[0]);
            for (int i = 1; i < 5; ++i)
                if (score[max.first] < score[i]) {
                    max.first = i;          // model num
                    max.second = score[i];  // the according probability
                }
            result.push_back(max);
        }
    }

    HMM* get_hmm() { return _hmm; }

private:
    HMM*                     _hmm; // this is only for training
    vector<vector<int> >     seq;
    double                   _alpha[6][50], _beta[6][50], _gamma[6][50]; // six states; length of sequence if 50
    double                   _epsilon[50][6][6]; // first index is depth, second idx transition to third index
    double                   _update_obs[6][6];  // accumulates the probabitlity of observing a paritcular charater given a state through all of the sequences
    double                   _delta[6][50];      // for viterbi
    size_t                   _itr;
    string                   _modelName;

    double obsProb(int state, int obs) {
        assert(state < 6);
        assert(obs < 6);
        return _hmm->observation[obs][state];
    }

    void check_alpha_beta() {
        double ref = 0;
        for (int i = 0; i < 6; ++i)
            ref += _alpha[i][49];
        
        int test_idx =20;
        double test = 0;
        for (int i = 0; i < 6; ++i)
            test += _alpha[i][test_idx] * _beta[i][test_idx];

        cout << "ref:  " << ref << endl;
        cout << "test: " << test << endl;
    }

    double transition(int i, int j) {
        return _hmm->transition[i][j];
    }

    void init_param() {
        for (int t = 0; t < 50; ++t) {
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    _epsilon[t][i][j] = 0;
                }
            }
        }
        for (int t = 0; t < 50; ++t) {
            for (int i = 0; i < 6; ++i) {
                _gamma[i][t] = 0;
            }
        }
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j)
                _update_obs[i][j] = 0;
        }
    }

    void update() {
        // update init_prob
        for (int state = 0; state < 6; ++state) {
            _hmm->initial[state] = _gamma[state][0] / seq.size();
        }

        // update transition matrix
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                double denominator = 0; // sum over _gamma
                _hmm->transition[i][j] = 0;
                for (int t = 0; t < 49; ++t) {
                    denominator            += _gamma[i][t];
                    _hmm->transition[i][j] += _epsilon[t][i][j];
                }
                _hmm->transition[i][j] /= denominator;
            }
        }

        // update observation matrix
        for (int state = 0; state < 6; ++state) {
            double denominator = 0;
            for (int t = 0; t < 50; ++t)
                denominator += _gamma[state][t];
            for (int obs = 0; obs < 6; ++obs)
                _hmm->observation[obs][state] = _update_obs[obs][state] / denominator;
        }
    }

    double viterbi(const HMM* hmm, const vector<int>& seq) {
        // initialize
        for (int state = 0; state < 6; ++ state) {
            _delta[state][0] = hmm->initial[state] * hmm->observation[seq[0]][state];
        }
        // recursive
        int t = 0;
        while (t < 49) {
            for (int state = 0; state < 6; ++ state) {
                _delta[state][t+1] = hmm->observation[seq[t+1]][state] * Max_transition(state, t, hmm);
            }
            ++t;
        }
        // get the max value
        double max = -1;
        for (int state = 0; state < 6; ++state) {
            if (max < _delta[state][49]) max = _delta[state][49];
        }
        return max;
    }

    double Max_transition(const int& state, const int& t, const HMM* hmm) {
        double max = -1;
        for (int prestate = 0; prestate < 6; ++prestate) {
            double tmp = hmm->transition[prestate][state] * _delta[prestate][t];
            if (max < tmp) max = tmp;
        }
        return max;
    }
};