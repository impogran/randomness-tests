#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex>

#include "Cephes.h"

using namespace std;

double LINE_SIZE = 128.0;
double BLOCK_SIZE = 128.0;
double PI = 3.131592;

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

int longestRun(string text) {
    int ctr = 0;
    int result = 0;  
    for (int i = 0; i < text.size(); i++) 
    { 
    
        if (text[i] == '0') 
        ctr = 0; 
    
        else
        { 
        ctr++;
        result = max(result, ctr); 
        } 
    } 
    string result1(result, '1');
    return result1.length(); 
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

// void longestRunOfOnesTest(vector<string> vecOfStr) {
//     int verdict = 0;
//     double M=8.0, K=3.0, N=16.0;
//     int v[4] = {0,0,0,0};
//     double pi_levels[4] = {0.2148, 0.3672, 0.2305, 0.1875};
//     for(auto line : vecOfStr) {
//         double chi_obs = 0.0;
//         double p_val = 0.0;
//         for(int i=0;i<LINE_SIZE;i+=N) {
//             int x = longestRun(line);
//             switch(x) {
//                 case 0 ... 1:
//                     v[0]++;
//                     break;
//                 case 2:
//                     v[1]++;
//                     break;
//                 case 3:
//                     v[2]++;
//                     break;
//                 case 4 ... 8:
//                     v[3]++;
//                     break;
//             }
//             for (int o = 0; o <= K; ++o){
//                 double td = v[i] - (N * pi_levels[i]);
//                 chi_obs += (td * td) / (N * pi_levels[i]);
//             }
//             p_val = Cephes::cephes_igamc(K/2,chi_obs/2);
//             printf("%f", p_val);
//             if(p_val > 0.01) verdict++;
//         }
//     }
//     printf("Longest Run of Ones Test: %d keys are random.\n", verdict);
// }

void discreteFourierTransformTest(string input) {
    int len = 100;
    double m[len/2+1];
    float Xr[len];
    float Xn[len];
    int xn[len];
    int N = len;
    for(int i=0;i<len;i++) {
        if(input[i] == '1') xn[i] = 1;
        else xn[i] = -1;
    }
    m[0] = sqrt(xn[0]*xn[0]);
    for (int k = 0; k < len/2; k++) {
        Xr[k] = 0;
        Xn[k] = 0;
        for (int n = 0; n < len; n++) {
            Xr[k]
                = (Xr[k]
                   + xn[n] * cos(2 * 3.141592 * k * n / N));
            Xn[k]
                = (Xn[k]
                   - xn[n] * sin(2 * 3.141592 * k * n / N));
        }
    }
    float T = sqrt(2.995732274*(double)(len));
    float N0 = len * 0.95 / 2.0;
    double count = 0;
    for (int z=0; z<len/2; z++ )
        if ( Xr[z]-Xn[z] < T )
            count++;
    float N1 = count;
    float d = (N1 - N0)/sqrt((double)(len)/4.0*0.95*0.05);
    float p_value = erfc(fabs(d)/sqrt(2.0));
    printf("Discreet Fourier Transform Test; p-value: %f.\n", p_value);
}

double min(double x, double y)
{
    return (x > y) ? y : x;
}

double max(double x, double y)
{
    return (x < y) ? y : x;
}

void RandomExcursionsVariantTest(string all) {
    int len = all.length()+1;
    int J = len/2;
    int constraint = (int)max(0.005*pow(len, 0.5), 500.0);
    int norm[len];
    for(int i=0;i<len;i++) {
        if(all[i] == '1') norm[i] = 1;
        else norm[i] = -1;
    }
    int partial_sums[len];
    for(int j=0;j<len;j++) {
        int partial_sum = 0;
        for(int k=0;k<j;k++) {
            partial_sum += norm[k];
        }
            partial_sums[j] = partial_sum;
    }
    int sums[len+2];
    sums[0] = 0;
    sums[len+1] = 0;
    for(int g=1;g<len+2;g++) 
    {
        sums[g] = partial_sums[g-1];
    }
    int	stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double p_value, min_p_value;
    int x, count;
    if(constraint > J) {
        printf("%s %f %f", "Not enough cycles...", constraint, J);
    } else {
        for (int p=0; p<=17; p++ ) {
            x = stateX[p];
            count = 0;
            for (int r=0; r<len; r++ )
                if ( all[p] == x )
                    count++;
            p_value = erfc(fabs(count-J)/(sqrt(2.0*J*(4.0*fabs(x)-2))));
            if (p == 0)
                min_p_value = p_value;
            else
                min_p_value = min(min_p_value,p_value);
            printf("%f \n", min_p_value);
            }
            printf("Smallest p in line: %f", min_p_value);
    } 
}


int main()
{
    vector<string> vecOfStr;
    bool result = getFileContent("C:\\Users\\gerar\\Desktop\\randomness-tests\\data128.txt", vecOfStr);
    
    if(result) {
        string all;
        for (auto c : vecOfStr)
        {
            all += c;
        }
        // frequencyTest(vecOfStr);
        // blockFrequencyTest(vecOfStr);
        // runsTest(vecOfStr);
        // discreteFourierTransformTest(all); //this one needs at least 1,000 bits to work
        //RandomExcursionsVariantTest(all); //this one needs at least 1,000,000 bits to work; make a lot of keys

    }

    return 0;
}