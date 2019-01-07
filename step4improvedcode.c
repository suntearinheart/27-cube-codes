#include "iostream"
#include "algorithm"
#include <time.h>
#include <cmath>
using namespace std;
#define LINE_NUM 3
#define BLOCK_NUM (LINE_NUM * LINE_NUM * LINE_NUM)
enum {
    LevelAM = 1,
    LevelAS = 2,
    LevelAL = 3,
    LevelBM = 4,
    LevelBS = 5,
    LevelBL = 6,
    LevelCM = 7,
    LevelCS = 8,
    LevelCL = 9,
};
#define LevelA  1
#define LevelB  2
#define LevelC  3

#define  MaxLevel  3
static double _blockAlpha[BLOCK_NUM];
static double  _blockarragement[BLOCK_NUM];
static double levels[MaxLevel][LINE_NUM * LINE_NUM ];
static double  maxvalue;
static double  minvalue;
static int levelspattern[LINE_NUM][LINE_NUM][LINE_NUM] = {
    // Bottom
    {{LevelAL,LevelCM,LevelBS},
     {LevelBM,LevelAS,LevelCL},
     {LevelCS,LevelBL,LevelAM}},
    // Middle
    {{LevelBS,LevelAL,LevelCM},
     {LevelCL,LevelBM,LevelAS},
     {LevelAM,LevelCS,LevelBL}},
    // Top
    {{LevelCM,LevelBS,LevelAL},
     {LevelAS,LevelCL,LevelBM},
     {LevelBL,LevelAM,LevelCS}}
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
    int asindex = 0;
    int amindex = 3;
    int alindex = 6;
    int bsindex = 0;
    int bmindex = 3;
    int blindex = 6;
    int csindex = 0;
    int cmindex = 0;
    int clindex = 0;
    
    
    for (int i = 0; i < LINE_NUM; i++) {
        for (int j = 0; j < LINE_NUM;j++){
            for (int k = 0; k < LINE_NUM;k++){
                switch(levelspattern[i][j][k]){
                    case    LevelAM:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelA-1][amindex++];
                        break;
                    case    LevelAS:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelA-1][asindex++];
                        break;
                    case    LevelAL:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelA-1][alindex++];
                        break;
                    case    LevelBM :
                        _blockarragement[i*9 + 3*j +k] = levels[LevelB-1][bmindex++];
                        break;
                    case   LevelBS :
                        _blockarragement[i*9 + 3*j +k] = levels[LevelB-1][bsindex++];
                        break;
                    case    LevelBL:
                        _blockarragement[i*9 + 3*j +k] = levels[LevelB-1][blindex++];
                        break;
                    case    LevelCM :
                        _blockarragement[i*9 + 3*j +k] = levels[LevelC-1][cmindex++];
                        break;
                    case     LevelCS :
                        _blockarragement[i*9 + 3*j +k] = levels[LevelC-1][csindex++];
                        break;
                    case     LevelCL :
                        _blockarragement[i*9 + 3*j +k] = levels[LevelC-1][clindex++];
                        break;
                        
                }
                
            }
            
        }
    }
    
    
}

double calculateSD(double data[])
{
    double avgAlpha = 0.0,mean, standardDeviation = 0.0;
    int i;
    for (int i = 0; i < BLOCK_NUM; i++) {
        avgAlpha += _blockAlpha[i];
    }
    mean = avgAlpha / BLOCK_NUM * LINE_NUM;

    for(i = 0; i < BLOCK_NUM; ++i)
        standardDeviation += pow(data[i] - mean, 2);
    
    return sqrt(standardDeviation / 10);
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

