#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <io.h>

#include "hmm.h"

using namespace std;

void initDict();
void loadFilePath();
void readUnigram();
void readBigram();
void setTransProb();
void inputData(string path);
void mapping();
double calObservationProb(int time, int state);
void viterbi();
vector<string> backtracking();
void outputData(string path, vector<string> result);

class Data {
	public :	
		int index;
		int start;
		int end;

		vector<int> dict;

		Data(int _index, int _start, int _end, vector<int> _dict)
		{
			index = _index;
			start = _start;
			end = _end;

			dict = _dict;
		}
};

enum MonoPhone {
	m_f, m_k, m_n, m_r, m_s, m_t, m_v, m_w, m_z, m_ah, m_ao, m_ay, m_eh, m_ey, m_ih, m_iy, m_ow, m_th, m_uw, m_sil, m_sp
};

double trans[127][127];
double unigram[12];
double bigram[12][12];
double **input;
double m[500][127];
int s[500][127];
int t, d;
vector<string> filePath;
map<string, Data*> word;
map<int, vector<int>> stateToindex;
ofstream ofs("recognized2.txt");
ofstream chk("check.txt");

int main() {
	loadFilePath();

	initDict();

	readUnigram();
	readBigram();

	mapping();
	setTransProb();
	
	ofs << "#!MLF!#" << endl;

	for (int i = 0; i < filePath.size(); ++i) {
		inputData(filePath[i]);
		viterbi();
		vector<string> result = backtracking();
		outputData(filePath[i], result);
	}

	ofs.close();
	chk.close();

	return 0;
}

void loadFilePath() {
	_finddata_t fd;
	string path = "tst";
	long handle = _findfirst((path + "\\*.*").c_str(), &fd);

	int result = 0;

	do
	{
		if (!strcmp(fd.name, ".") || !strcmp(fd.name, "..")) {
			result = _findnext(handle, &fd);
			continue;
		}

		_finddata_t fd2;
		int result2 = 0;
		long handle2 = _findfirst((path + "\\" + fd.name + "\\*.*").c_str(), &fd2);

		do
		{
			if (!strcmp(fd2.name, ".") || !strcmp(fd2.name, "..")) {
				result2 = _findnext(handle2, &fd2);
				continue;
			}

			_finddata_t fd3;
			int result3 = 0;
			long handle3 = _findfirst((path + "\\" + fd.name + "\\" + fd2.name + "\\*.txt").c_str(), &fd3);

			do
			{
				if (!strcmp(fd3.name, ".") || !strcmp(fd3.name, "..")) {
					result3 = _findnext(handle3, &fd3);
					continue;
				}

				filePath.push_back(path + "\\" + fd.name + "\\" + fd2.name + "\\" + fd3.name);

				result3 = _findnext(handle3, &fd3);
			} while (result3 != -1);

			_findclose(handle3);

			result2 = _findnext(handle2, &fd2);
		} while (result2 != -1);

		_findclose(handle2);

		result = _findnext(handle, &fd);
	} while (result != -1);

	_findclose(handle);
}

void initDict() {
	word["<s>"] = new Data(0, 1, 3, vector<int>{(MonoPhone::m_sil)});
	word["eight"] = new Data(1, 4, 10, vector<int>{MonoPhone::m_ey, MonoPhone::m_t, MonoPhone::m_sp});
	word["five"] = new Data(2, 11, 20, vector<int>{MonoPhone::m_f, MonoPhone::m_ay, MonoPhone::m_v, MonoPhone::m_sp});
	word["four"] = new Data(3, 21, 30, vector<int>{MonoPhone::m_f, MonoPhone::m_ao, MonoPhone::m_r, MonoPhone::m_sp});
	word["nine"] = new Data(4, 31, 40, vector<int>{MonoPhone::m_n, MonoPhone::m_ay, MonoPhone::m_n, MonoPhone::m_sp});
	word["oh"] = new Data(5, 41, 44, vector<int>{MonoPhone::m_ow, MonoPhone::m_sp});
	word["one"] = new Data(6, 45, 54, vector<int>{MonoPhone::m_w, MonoPhone::m_ah, MonoPhone::m_n, MonoPhone::m_sp});
	word["seven"] = new Data(7, 55, 70, vector<int>{MonoPhone::m_s, MonoPhone::m_eh, MonoPhone::m_v, MonoPhone::m_ah, MonoPhone::m_n, MonoPhone::m_sp});
	word["six"] = new Data(8, 71, 83, vector<int>{MonoPhone::m_s, MonoPhone::m_ih, MonoPhone::m_k, MonoPhone::m_s, MonoPhone::m_sp});
	word["three"] = new Data(9, 84, 93, vector<int>{MonoPhone::m_th, MonoPhone::m_r, MonoPhone::m_iy, MonoPhone::m_sp});
	word["two"] = new Data(10, 94, 100, vector<int>{MonoPhone::m_t, MonoPhone::m_uw, MonoPhone::m_sp});
	word["zero"] = new Data(11, 101, 113, vector<int>{MonoPhone::m_z, MonoPhone::m_ih, MonoPhone::m_r, MonoPhone::m_ow, MonoPhone::m_sp});
	word["zero2"] = new Data(11, 114, 126, vector<int>{MonoPhone::m_z, MonoPhone::m_iy, MonoPhone::m_r, MonoPhone::m_ow, MonoPhone::m_sp});
}

void readUnigram() {
	ifstream ifs("unigram.txt");

	string key;
	double data;

	while (!ifs.eof()) {
		ifs >> key >> data;

		unigram[word[key]->index] = data;
	}

	ifs.close();
}

void readBigram() {
	ifstream ifs("bigram.txt");

	string first_key;
	string second_key;
	double data;

	while (!ifs.eof()) {
		ifs >> first_key >> second_key >> data;

		bigram[word[first_key]->index][word[second_key]->index] = data;
	}

	ifs.close();
}

void mapping() {
	for (int i = 1; i < 127; ++i) {
		map<string, Data*>::iterator iter;

		for (iter = word.begin(); iter != word.end(); ++iter) {
			if (word[iter->first]->start <= i && i <= word[iter->first]->end) {
				break;
			}
		}

		int idx = i - word[iter->first]->start;

		int phoneIdx = word[iter->first]->dict[idx / N_STATE];
		int stateIdx = idx % N_STATE;

		stateToindex[i] = vector<int>{ phoneIdx, stateIdx };
	}
}

void setTransProb() {
	map<string, Data*>::iterator iter;

	//Initial Prob
	for (iter = word.begin(); iter != word.end(); ++iter) {
		trans[0][word[iter->first]->start] = unigram[word[iter->first]->index];
	}

	//MonoPhone Interal Prob
	for (iter = word.begin(); iter != word.end(); ++iter) {
		int index = 0;

		for (int j = 0; j<word[iter->first]->dict.size(); ++j) {
			if (word[iter->first]->dict[j] == MonoPhone::m_sp)
				trans[word[iter->first]->end][word[iter->first]->end] = phones[MonoPhone::m_sp].tp[1][1];
			else {
				for (int a = 0; a < N_STATE; ++a) {
					for (int b = 0; b < N_STATE + 1; ++b) {
						trans[word[iter->first]->start + index + a][word[iter->first]->start + index + b] = phones[word[iter->first]->dict[j]].tp[a + 1][b + 1];
					}
				}

				index += 3;
			}
		}
	}

	//Word to Word Prob
	for (iter = word.begin(); iter != word.end(); ++iter) {
		map<string, Data*>::iterator iter2;

		for (iter2 = word.begin(); iter2 != word.end(); ++iter2) {
			if (iter->first == "<s>") {
				if(iter2->first == "<s>")
					trans[word[iter->first]->end][word[iter2->first]->start] = phones[word[iter->first]->dict[word[iter->first]->dict.size() - 1]].tp[3][1];
				else
					trans[word[iter->first]->end][word[iter2->first]->start] = phones[word[iter->first]->dict[word[iter->first]->dict.size() - 1]].tp[3][4] * bigram[word[iter->first]->index][word[iter2->first]->index];
			}
			else {
				trans[word[iter->first]->end - 1][word[iter2->first]->start] = phones[word[iter->first]->dict[word[iter->first]->dict.size() - 2]].tp[3][4] * phones[MonoPhone::m_sp].tp[0][2] * bigram[word[iter->first]->index][word[iter2->first]->index];
				trans[word[iter->first]->end][word[iter2->first]->start] = phones[MonoPhone::m_sp].tp[1][2] * bigram[word[iter->first]->index][word[iter2->first]->index];
			}

			trans[word[iter->first]->end - 1][word[iter2->first]->end] *= phones[MonoPhone::m_sp].tp[0][1];
		}
	}
}

double calObservationProb(int time, int state) {
	vector<int> idx = stateToindex[state];
	
	int phoneIdx = idx[0];
	int stateIdx = idx[1];

	double p[N_PDF] = { 0, };

	int maxIndex = -1;
	double maxValue = -1 * DBL_MAX;

	for(int i=0;i<N_PDF;++i) {
		double weight = phones[phoneIdx].state[stateIdx].pdf[i].weight;

		p[i] = log(weight);

		for (int j = 0; j < N_DIMENSION; ++j) {
			double data = input[time][j];

			double mean = phones[phoneIdx].state[stateIdx].pdf[i].mean[j];
			double var = phones[phoneIdx].state[stateIdx].pdf[i].var[j];

			double a = log(1 / sqrt(2 * M_PI * var));
			double b = -0.5 * pow(data - mean, 2) / var;

			p[i] = p[i] + a + b;
		}

		if (maxValue < p[i]) {
			maxIndex = i;
			maxValue = p[i];
		}
	}

	double temp = p[0];
	p[0] = maxValue;
	p[maxIndex] = temp;

	double result = 1;

	for (int i = 1; i < N_PDF; ++i)
		result = result + pow(M_E, p[i] - p[0]);

	result = p[0] + log(result);

	return result;
}

void viterbi() {
	for (int i = 0; i < 500; ++i) {
		for (int j = 0; j < 127; ++j) {
			m[i][j] = -1 * DBL_MAX;
			s[i][j] = 0;
		}
	}

	for (int i = 1; i < 127; ++i) {
		if (trans[0][i] != 0) {
			m[0][i] = log(trans[0][i]) + calObservationProb(0, i);
		}
	}

	for (int time = 1; time < t; ++time) {
		for (int i = 1; i < 127; ++i) {
			double op = calObservationProb(time, i);

			double maxValue = -1 * DBL_MAX;
			int maxState = -1;

			for (int j = 1; j < 127; ++j) {
				if (trans[j][i] != 0) {
					double temp = op + log(trans[j][i]) + m[time - 1][j];

					if (temp > maxValue) {
						maxValue = temp;
						maxState = j;
					}
				}
			}

			m[time][i] = maxValue;
			s[time][i] = maxState;
		}
	}
}

vector<string> backtracking() {
	double maxValue = -1 * DBL_MAX;
	int maxState;

	for (int i = 1; i < 127; ++i) {
		if (maxValue < m[t - 1][i]) {
			maxValue = m[t - 1][i];
			maxState = i;
		}
	}

	int stateSeq[500];

	stateSeq[t - 1] = maxState;

	for (int i = t - 2; i >= 0; --i)
		stateSeq[i] = s[i + 1][stateSeq[i + 1]];

	vector<string> result;

	for (int i = 0; i < t; ++i) {
		map<string, Data*>::iterator iter;

		for (iter = word.begin(); iter != word.end(); ++iter) {
			if (iter->first == "<s>")
				continue;
			else {
				if (stateSeq[i] == word[iter->first]->start) {
					result.push_back(iter->first);

					while (word[iter->first]->end - 1 != stateSeq[i]) {
						i++;
						if (i == t) {
							break;
						}
					}
					break;
				}
			}
		}
	}

	return result;
}

void inputData(string path) {
	FILE *fp = NULL;

	fp = fopen(path.c_str(), "rt");

	if (fp == NULL) {
		cout << "Fail " << endl;
	}
	else {
		fscanf(fp, "%d %d\n", &t, &d);
	}

	input = new double*[t];

	for (int i = 0; i < t; ++i)
		input[i] = new double[d];

	for (int i = 0; i < t; ++i)
		for (int j = 0; j < d; ++j)
			fscanf(fp, "%lf", &input[i][j]);

	fclose(fp);
}

void outputData(string path, vector<string> result) {
	string name = "";

	for (int i = 0; i < path.size(); ++i) {
		if (path[i] == '\\')
			name = name + "/";
		else if (path[i] == '.')
			break;
		else
			name = name + path[i];
	}

	ofs << "\"" << name << ".rec\"" << endl;

	if (result.size() > 7)
		chk << path << endl;

	if (result.size() == 9) {
		for (int i = 1; i < result.size() - 1; ++i)
			ofs << result[i] << endl;
	}
	else if (result.size() == 8) {
		for (int i = 0; i < result.size() - 2; ++i)
			ofs << result[i] << endl;
	}
	else {
		int count = 0;
		bool flag = false;

		for (int i = 0; i < result.size(); ++i) {
			if (count == 7)
				break;

			if (result[i] == "oh" && flag == false) {
				flag = true;
				continue;
			}
			else if (result[i] == "oh" && flag == true) {
				ofs << result[i] << endl;
				count++;
				flag = false;
			}
			else {
				ofs << result[i] << endl;
				count++;
			}
		}
	}

	ofs << "." << endl;

	for(int i=0;i<t;++i)
		delete input[i];

	delete input;
}