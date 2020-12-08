#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void E_update(
    double ***Er, double ***Eth, double ***Eph,
    double ***nDr, double ***nDth, double ***nDph,
    double ***oDr, double ***oDth, double ***oDph
){
    for( int i = 0; i < Nr; i++ ){
        for( int j = 1; j < Nth; j++ ){
            for( int k = 1; k < Nph; k++ ){
                Er[i][j][k] = Er[i][j][k] +
                    (nDr[i][j][k] - oDr[i][j][k])/EPS0;
                oDr[i][j][k] = nDr[i][j][k];
            }
        }
    }

    for( int i = 1; i < Nr; i++ ){
        for( int j = 0; j < Nth; j++ ){
            for( int k = 1; k < Nph; k++ ){
                Eth[i][j][k] = Eth[i][j][k] + 
                    (nDth[i][j][k] - oDth[i][j][k])/EPS0;
                oDth[i][j][k] = nDth[i][j][k];
            }
        }
    }

    for( int i = 1; i < Nr; i++ ){
        for( int j = 1; j < Nth; j++ ){
            for( int k = 0; k < Nph; k++ ){
                Eph[i][j][k] = Eph[i][j][k] + 
                        (nDph[i][j][k] - oDph[i][j][k])/EPS0;
                oDph[i][j][k] = nDph[i][j][k];
            }
        }
    }

}