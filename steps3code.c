#include "iostream"
#include "algorithm"
#include <time.h>
#include <cmath>
using namespace std;
#define LINE_NUM 3
#define BLOCK_NUM (LINE_NUM * LINE_NUM * LINE_NUM)
#define  LevelA   1
#define  LevelB   2
#define  LevelC   3
#define  MaxLevel  LevelC
static double _blockAlpha[BLOCK_NUM];
static double  _blockarragement[BLOCK_NUM];
static double levels[MaxLevel][LINE_NUM * LINE_NUM ];
static double  maxvalue;
static double  minvalue;
static int levelspattern[LINE_NUM][LINE_NUM][LINE_NUM] = {
    // Bottom
    {{LevelA,LevelC,LevelB},
        {LevelB,LevelA,LevelC},
        {LevelC,LevelB,LevelA}},
    // Middle
    {{LevelB,LevelA,LevelC},
        {LevelC,LevelB,LevelA},
        {LevelA,LevelC,LevelB}},
    // Top
    {{LevelC,LevelB,LevelA},
        {LevelA,LevelC,LevelB},
        {LevelB,LevelA,LevelC}}
};


//do the init T(x)
void init1() {
    
    double alpha[BLOCK_NUM] = {0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79};
    
    
    for (int i = 0; i < BLOCK_NUM - 1; i++) {
        for (int j = i + 1; j < BLOCK_NUM; j++) {
            if (alpha[i] > alpha[j]) {
                double temp = alpha[i];
                alpha[i] = alpha[j];
                alpha[j] = temp;
            }
        }
    }
    
    memcpy(_blockAlpha, alpha, sizeof(double) * BLOCK_NUM);
    maxvalue = _blockAlpha[26];
    minvalue = _blockAlpha[0];
}

void divide_level(){
    double levelAmaxvalue  = minvalue + (maxvalue - minvalue)/3;
    double levelBmaxvalue = minvalue + 2*(maxvalue - minvalue)/3;
    
    int  j =0;
    int k = 0;
    int m = 0;
    for (int i = 0; i < BLOCK_NUM; i++) {
        if (_blockAlpha[i] < levelAmaxvalue ) {
            levels[0][j++] = _blockAlpha[i];
        } else if(_blockAlpha[i] < levelBmaxvalue ){
            levels[1][k++] = _blockAlpha[i];
        } else {
            levels[2][m++] =  _blockAlpha[i];
        }
    }
}

void arrange_pattern(){
    int aindex = 0;
    int bindex = 0;
    int cindex = 0;
    
    for (int i = 0; i < LINE_NUM; i++) {
        for (int j = 0; j < LINE_NUM;j++){
            for (int k = 0; k < LINE_NUM;k++){
                switch(levelspattern[i][j][k]){
                    case LevelA:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelA-1][aindex++];
                        break;
                    case LevelB:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelB-1][bindex++];
                        break;
                        
                    case LevelC:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelC-1][cindex++];
                        break;
                }
                
            }
            
        }
    }
    
    
}
double calculateSD(double data[])
{
    double avgAlpha = 0.0,mean, standardDeviation = 0.0;
    double sumxyz[BLOCK_NUM];
    
    int i;
    for (int i = 0; i < BLOCK_NUM; i++) {
        avgAlpha += _blockarragement[i];
    }
    mean = avgAlpha / BLOCK_NUM * LINE_NUM;
    
    for(i = 0; i < 9 ; i++){
        //for x
        sumxyz[i] = _blockarragement[3*i] + _blockarragement[3*i+1] + _blockarragement[3*i+2];
        
    }
    for(i = 0; i < 3 ; i++){
        
        // for y
        sumxyz[9+i] = _blockarragement[i] + _blockarragement[i+3] + _blockarragement[i+6];
        sumxyz[9+3+i] = _blockarragement[i+9] + _blockarragement[i+12] + _blockarragement[i+15];
        sumxyz[9+3+3+i] = _blockarragement[i+18] + _blockarragement[i+21] + _blockarragement[i+24];
    }
    for(i = 0; i < 9 ; i++){
        //for z
        sumxyz[18+i] = _blockarragement[i] + _blockarragement[i+9] + _blockarragement[i+18];
    }
    
    for(i = 0; i < BLOCK_NUM; ++i)
        standardDeviation += pow(sumxyz[i] - mean, 2);
    
    return sqrt(standardDeviation / BLOCK_NUM);
}


//output
void outResult() {
    double avgAlpha = 0;
    for (int i = 0; i < BLOCK_NUM; i++) {
        avgAlpha += _blockAlpha[i];
    }
    avgAlpha = avgAlpha / BLOCK_NUM * LINE_NUM;
    cout << "the  T(average) is " << avgAlpha << endl;
    
    
    cout <<"the arrangement of the Cube:\n";
    
    cout << "---------------\n";
    for (int i = 0; i < 3; i++) {
        switch (i){
            case 0:
                cout << "\nBottom:\n";
                break;
            case 1:
                cout << "\nMiddle:\n";
                break;
                
            case 2:
                cout << "\nTop:\n";
        }
        
        
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                cout << _blockarragement[i * 9 + j * 3 + k] << " ";
                
            }
            cout << "\n";
        }
        
    }
    cout << "---------------\n";
    cout << "the SD : ";
    cout << calculateSD(_blockarragement) << "\n";
    
}


int main(int argc, const char * argv[]) {
    clock_t start;
    start = clock();
    init1();
    divide_level();
    arrange_pattern();
    cout << "\ntime consume（ticks）: " << (clock() - start) << endl;
    
    outResult();
    
    return 0;
}
