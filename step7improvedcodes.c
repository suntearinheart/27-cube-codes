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
static double levelabest[3][3];
static double levelbbest[3][3];
static double levelcbest[3][3];
static  double best_sd;

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
    int cmindex = 3;
    int clindex = 6;
    
    
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
    double sumxyz[BLOCK_NUM];
    
    int i;
    for (int i = 0; i < BLOCK_NUM; i++) {
        avgAlpha += data[i];
    }
    mean = avgAlpha / BLOCK_NUM * LINE_NUM;
    
    for(i = 0; i < 9 ; i++){
        //for x
        sumxyz[i] = data[3*i] + data[3*i+1] + data[3*i+2];
        
    }
    for(i = 0; i < 3 ; i++){
        
        // for y
        sumxyz[9+i] = data[i] + data[i+3] + data[i+6];
        sumxyz[9+3+i] = data[i+9] + data[i+12] + data[i+15];
        sumxyz[9+3+3+i] = data[i+18] + data[i+21] + data[i+24];
    }
    for(i = 0; i < 9 ; i++){
        //for z
        sumxyz[18+i] = data[i] + data[i+9] + data[i+18];
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
    cout << "the best SD: " << " ";
    cout << best_sd<< "\n";
    
}


bool  next_loop(double * levell_p,double * levelm_p,double * levels_p){
    bool levels_loop = true;
    bool levelm_loop = true;

    do {
        if(!levelm_loop){
       
            return true;
        }
        
        do {
            
            if (levels_loop&&next_permutation(levels_p,levels_p+3)){
         
                return true;
            } else {
                if(!levels_loop){
       
                    return true;
                }
                sort(levels_p,levels_p+3);
                levels_loop = false;
            }
        } while (next_permutation(levelm_p,levelm_p+3));
        
         levelm_loop = false;
         sort(levels_p,levels_p+3);
         sort(levelm_p,levelm_p+3);
    } while (next_permutation(levell_p,levell_p+3) );
    
    return false;
}

bool checkbest_arrangement(double  block[BLOCK_NUM]){

    double  mean;
  
    mean = calculateSD(block);
  
    if(best_sd >  mean){
        best_sd =  mean;
        return true;
    }
    
    return false;
 
}

void findout_thebest(){
    // 0~2 : large   3~5: Medium 6~8: Small
    int i = 0;
    double levela[3][3];
    double levelb[3][3];
    double levelc[3][3];
    double  block[BLOCK_NUM];
    bool levels_loop = true;
    bool levelm_loop = true;
    
    double * levelal_p = &levela[0][0];
    double * levelam_p = &levela[1][0];
    double * levelas_p = &levela[2][0];
    
    double * levelbl_p = &levelb[0][0];
    double * levelbm_p = &levelb[1][0];
    double * levelbs_p = &levelb[2][0];
    
    double * levelcl_p = &levelc[0][0];
    double * levelcm_p = &levelc[1][0];
    double * levelcs_p = &levelc[2][0];

    levela[0][0] = _blockarragement[0];
    levela[0][1] = _blockarragement[10];
    levela[0][2] = _blockarragement[20];
    levela[1][0] = _blockarragement[8];
    levela[1][1] = _blockarragement[15];
    levela[1][2] = _blockarragement[25];
    levela[2][0] = _blockarragement[4];
    levela[2][1] = _blockarragement[14];
    levela[2][2] = _blockarragement[21];
    
    levelb[0][0] = _blockarragement[7];
    levelb[0][1] = _blockarragement[17];
    levelb[0][2] = _blockarragement[24];
    levelb[1][0] = _blockarragement[3];
    levelb[1][1] = _blockarragement[13];
    levelb[1][2] = _blockarragement[23];
    levelb[2][0] = _blockarragement[2];
    levelb[2][1] = _blockarragement[9];
    levelb[2][2] = _blockarragement[19];
    
    levelc[0][0] = _blockarragement[5];
    levelc[0][1] = _blockarragement[12];
    levelc[0][2] = _blockarragement[22];
    levelc[1][0] = _blockarragement[1];
    levelc[1][1] = _blockarragement[11];
    levelc[1][2] = _blockarragement[18];
    levelc[2][0] = _blockarragement[6];
    levelc[2][1] = _blockarragement[16];
    levelc[2][2] = _blockarragement[26];
    
    /*init the best arrangement*/
    memcpy(levelabest, levela, sizeof(double)*9);
    memcpy(levelbbest, levelb, sizeof(double)*9);
    memcpy(levelcbest, levelc, sizeof(double)*9);
    
    best_sd = calculateSD(_blockarragement);
    
    do {
        if(!levelm_loop){
            i++;
            levels_loop = true;
            levelm_loop = true;
            block[0] = levela[0][0];
            block[10] = levela[0][1];
            block[20] = levela[0][2];
            block[8] = levela[1][0] ;
            block[15] = levela[1][1] ;
            block[25]  = levela[1][2] ;
            block[4] = levela[2][0];
            block[14] = levela[2][1] ;
            block[21] = levela[2][2];
            
            block[7] = levelb[0][0]  ;
            block[17] = levelb[0][1] ;
            block[24] = levelb[0][2];
            block[3] = levelb[1][0] ;
            block[13] = levelb[1][1] ;
            block[23] = levelb[1][2]  ;
            block[2] = levelb[2][0] ;
            block[9] = levelb[2][1] ;
            block[19] = levelb[2][2] ;
            
            block[5] = levelc[0][0];
            block[12] = levelc[0][1] ;
            block[22] = levelc[0][2];
            block[1] = levelc[1][0] ;
            block[11] = levelc[1][1] ;
            block[18] = levelc[1][2];
            block[6] = levelc[2][0] ;
            block[16] = levelc[2][1];
            block[26] = levelc[2][2] ;
            
            if( true == checkbest_arrangement(block)){
                memcpy(levelabest, levela, sizeof(double)*9);
                memcpy(levelbbest, levelb, sizeof(double)*9);
                memcpy(levelcbest, levelc, sizeof(double)*9);
            }
        }
        
        do {
            if(!levels_loop){
                i++;
                block[0] = levela[0][0];
                block[10] = levela[0][1];
                block[20] = levela[0][2];
                block[8] = levela[1][0] ;
                block[15] = levela[1][1] ;
                block[25]  = levela[1][2] ;
                block[4] = levela[2][0];
                block[14] = levela[2][1] ;
                block[21] = levela[2][2];
                
                block[7] = levelb[0][0]  ;
                block[17] = levelb[0][1] ;
                block[24] = levelb[0][2];
                block[3] = levelb[1][0] ;
                block[13] = levelb[1][1] ;
                block[23] = levelb[1][2]  ;
                block[2] = levelb[2][0] ;
                block[9] = levelb[2][1] ;
                block[19] = levelb[2][2] ;
                
                block[5] = levelc[0][0];
                block[12] = levelc[0][1] ;
                block[22] = levelc[0][2];
                block[1] = levelc[1][0] ;
                block[11] = levelc[1][1] ;
                block[18] = levelc[1][2];
                block[6] = levelc[2][0] ;
                block[16] = levelc[2][1];
                block[26] = levelc[2][2] ;
                
                if( true == checkbest_arrangement(block)){
                    memcpy(levelabest, levela, sizeof(double)*9);
                    memcpy(levelbbest, levelb, sizeof(double)*9);
                    memcpy(levelcbest, levelc, sizeof(double)*9);
                }

                levels_loop = true;
            }
            
            while (levels_loop&&next_loop(levelal_p,levelam_p,levelas_p))
            {
          
                i++;
                block[0] = levela[0][0];
                block[10] = levela[0][1];
                block[20] = levela[0][2];
                block[8] = levela[1][0] ;
                block[15] = levela[1][1] ;
                block[25]  = levela[1][2] ;
                block[4] = levela[2][0];
                block[14] = levela[2][1] ;
                block[21] = levela[2][2];
                
                block[7] = levelb[0][0]  ;
                block[17] = levelb[0][1] ;
                block[24] = levelb[0][2];
                block[3] = levelb[1][0] ;
                block[13] = levelb[1][1] ;
                block[23] = levelb[1][2]  ;
                block[2] = levelb[2][0] ;
                block[9] = levelb[2][1] ;
                block[19] = levelb[2][2] ;
                
                block[5] = levelc[0][0];
                block[12] = levelc[0][1] ;
                block[22] = levelc[0][2];
                block[1] = levelc[1][0] ;
                block[11] = levelc[1][1] ;
                block[18] = levelc[1][2];
                block[6] = levelc[2][0] ;
                block[16] = levelc[2][1];
                block[26] = levelc[2][2] ;
                
                if( true == checkbest_arrangement(block)){
                    memcpy(levelabest, levela, sizeof(double)*9);
                    memcpy(levelbbest, levelb, sizeof(double)*9);
                    memcpy(levelcbest, levelc, sizeof(double)*9);
                }
            };
                

            sort(levelal_p,levelal_p+3);
            sort(levelam_p,levelam_p+3);
            sort(levelas_p,levelas_p+3);
            levels_loop = false;

        } while (next_loop(levelbl_p,levelbm_p,levelbs_p));
        
        levelm_loop = false;
        sort(levelal_p,levelal_p+3);
        sort(levelam_p,levelam_p+3);
        sort(levelas_p,levelas_p+3);
        
        sort(levelbl_p,levelbl_p+3);
        sort(levelbm_p,levelbm_p+3);
        sort(levelbs_p,levelbs_p+3);
    } while (next_loop(levelcl_p,levelcm_p,levelcs_p) );

    _blockarragement[0] = levelabest[0][0];
    _blockarragement[10] = levelabest[0][1];
    _blockarragement[20] = levelabest[0][2];
    _blockarragement[8] = levelabest[1][0] ;
    _blockarragement[15] = levelabest[1][1] ;
    _blockarragement[25]  = levelabest[1][2] ;
    _blockarragement[4] = levelabest[2][0];
    _blockarragement[14] = levelabest[2][1] ;
    _blockarragement[21] = levelabest[2][2];
    
    _blockarragement[7] = levelbbest[0][0]  ;
    _blockarragement[17] = levelbbest[0][1] ;
    _blockarragement[24] = levelbbest[0][2];
    _blockarragement[3] = levelbbest[1][0] ;
    _blockarragement[13] = levelbbest[1][1] ;
    _blockarragement[23] = levelbbest[1][2]  ;
    _blockarragement[2] = levelbbest[2][0] ;
    _blockarragement[9] = levelbbest[2][1] ;
    _blockarragement[19] = levelbbest[2][2] ;
    
    _blockarragement[5] = levelcbest[0][0];
    _blockarragement[12] = levelcbest[0][1] ;
    _blockarragement[22] = levelcbest[0][2];
    _blockarragement[1] = levelcbest[1][0] ;
    _blockarragement[11] = levelcbest[1][1] ;
    _blockarragement[18] = levelcbest[1][2];
    _blockarragement[6] = levelcbest[2][0] ;
    _blockarragement[16] = levelcbest[2][1];
    _blockarragement[26] = levelcbest[2][2] ;
    cout <<  "the total loop:" << i;
}

int main(int argc, const char * argv[]) {
    clock_t start;
    start = clock();
    init1();
    divide_level();
    arrange_pattern();
    findout_thebest();
    cout << "\ntime consume（ticks）: " << (clock() - start) << endl;
    cout << "time consume（sec）: " << (clock() - start)/CLOCKS_PER_SEC << endl;

    outResult();
    
    return 0;
}

