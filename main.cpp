#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "Cephes.h"

using namespace std;

int LINE_SIZE = 128;
int BLOCK_SIZE = 128;

bool getFileContent(string fileName, vector<string> & vecOfStrs)
{
    ifstream in(fileName.c_str());
    if(!in)
    {
        cerr << "Cannot open the File : "<<fileName<<endl;
        return false;
    }
    string str;
    int file_size = 0;
    while (getline(in, str))
    {
        if(str.size() > 0)
            vecOfStrs.push_back(str);
        file_size++;
    }
    in.close();
    return true;
}

int longestRun(string input) {
    int maxLength = 0;
    int currentLength = 0;

    for (char c : input) {
        if (c == '1') currentLength++;
        else {
            if (currentLength > maxLength) maxLength = currentLength;
            currentLength = 1;
        }
    }

    if (currentLength > maxLength) maxLength = currentLength;
    return maxLength;
}

void frequencyTest(vector<string> vecOfStr) {
    int verdict = 0;
    for(string & line : vecOfStr) {
        int s = 0;
        float s_obs = 0.0;
        float p_val = 0.0;
            
        for(int i=0;i<LINE_SIZE;i++) {
            if(line[i]=='0') s -= 1;
            else s++;
        }
        s_obs = abs(s) / sqrt(LINE_SIZE);
        p_val = erfc(s_obs / sqrt(2));
        if(p_val > 0.01) verdict++;
    }
    printf("Frequency Test: %d keys are random.\n", verdict);
}

void blockFrequencyTest(vector<string> vecOfStr) {
    string all;
    for (auto c : vecOfStr)
    {
        all += c;
    }
    int verdict = 0;
    for(int i=0;i<all.length();i+=BLOCK_SIZE) {
        int s = 0;
        float s_obs = 0.0;
        float p_val = 0.0;
        for(int j=i;j<BLOCK_SIZE;j++) {
            if(all[j]=='0') s -= 1;
            else s++;
        }
        s_obs = abs(s) / sqrt(LINE_SIZE);
        p_val = erfc(s_obs / sqrt(2));
        if(p_val > 0.01) verdict++;
    }
    printf("Block Frequency Test: %d blocks are random.\n", verdict);
}

void runsTest(vector<string> vecOfStr) {
    float tao = 2 / sqrt(LINE_SIZE);
    int verdict = 0;
    for(string & line : vecOfStr) {
        int s = 0;
        float p_val = 0.0;
        float count = 0.0;
        for (char c : line) {
            if (c == '1') count+=1.0;
        }
        float pi = 0.0;
        pi = count / LINE_SIZE;
        if(abs(pi - 0.5) >= tao) continue;
        else {
            for(int i=0;i<LINE_SIZE;i++) {
                if(line[i]=='1') s += 1;
            }
        }
        p_val = erfc((s - 2 * pi * LINE_SIZE * (1 - pi)) / 2 * sqrt(2 * LINE_SIZE) * pi * (1 - pi));
        if(p_val > 0.01) verdict++;
        }
        printf("Runs Test: %d keys are random.\n", verdict);
    }

void longestRunOfOnesTest(vector<string> vecOfStr) {
    string all;
    for (auto c : vecOfStr)
    {
        all += c;
    }
    int verdict = 0;
    double M=8.0, K=3.0, N=16.0;
    int v[4] = {0,0,0,0};
    double pi_levels[4] = {0.2148, 0.3672, 0.2305, 0.1875};
    for(int i=0;i<all.length();i+=1600) {
        int gowno = 8;
        double chi_obs = 0.0;
        double p_val = 0.0;
        while(gowno > 0) {
            int x = longestRun(all.substr(i, 1600));
            switch(x) {
                case 0 ... 1:
                    v[0]++;
                    break;
                case 2:
                    v[1]++;
                    break;
                case 3:
                    v[2]++;
                    break;
                case 4 ... 8:
                    v[3]++;
                    break;
            }
            gowno--;
        }
        for (i = 0; i <= K; ++i){
            double td = v[i] - (N * pi_levels[i]);
            chi_obs += (td * td) / (N * pi_levels[i]);
        }
        p_val = Cephes::cephes_igamc(K/2,chi_obs/2);
        if(p_val > 0.01) verdict++;
    }
    printf("Longest Run of Ones Test: %d keys are random.\n", verdict);
}
    

int main()
{
    vector<string> vecOfStr;
    bool result = getFileContent("data128.txt", vecOfStr);
    
    if(result) {
        // frequencyTest(vecOfStr);
        // blockFrequencyTest(vecOfStr);
        // runsTest(vecOfStr);
        longestRunOfOnesTest(vecOfStr); // not gonna work due to underflow; g++ -o myfile cephes.cpp main.cpp

    }

    return 0;
}