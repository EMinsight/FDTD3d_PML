#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void H_update(
    double ***Er, double ***Eth, double ***Eph,
    double ***Hr, double ***Hth, double ***Hph
){
    double CHr1, CHr2;
    double CHth1, CHth2;
    double CHph1, CHph2;

    double r_i1, r_i2, r_i3;
    double si_th1, si_th2, si_th3;

    for( int i = 1; i < Nr + 1; i++ ){
        r_i1 = dist(i);
        for( int j = 0; j < Nth; j++ ){
            si_th1 = std::sin(th(j));
            si_th2 = std::sin(th(j+0.5));
            si_th3 = std::sin(th(j+1.0));

            CHr1 = Dt/MU0/r_i1/si_th2/delta_th;
            CHr2 = Dt/MU0/r_i1/si_th2/delta_ph;

            for( int k = 0; k < Nph; k++ ){

                Hr[i][j][k] = Hr[i][j][k] - CHr1*(si_th3*Eph[i][j+1][k] - si_th1*Eph[i][j][k])
                    + CHr2*(Eth[i][j][k+1] - Eth[i][j][k]);
            }
        }
    }

    for( int i = 0; i < Nr; i++ ){
        r_i1 = dist(i);
        r_i2 = dist(i+0.5);
        r_i3 = dist(i+1.0);
        for( int j = 1; j < Nth + 1; j++ ){
            si_th1 = std::sin(th(j));
            
            CHth1 = Dt/MU0/r_i2/si_th1/delta_ph;
            CHth2 = Dt/MU0/r_i2/delta_r;
            for( int k = 0; k < Nph; k++ ){

                Hth[i][j][k] = Hth[i][j][k] - CHth1*(Er[i][j][k+1] - Er[i][j][k])
                        + CHth2*(r_i3*Eph[i+1][j][k] - r_i1*Eph[i][j][k]);

            }
        }
    }

    for( int i = 0; i < Nr; i++ ){
        r_i1 = dist(i);
        r_i2 = dist(i+0.5);
        r_i3 = dist(i+1.0);
        for( int j = 0; j < Nth; j++ ){
            CHph1 = Dt/MU0/r_i2/delta_r;
            CHph2 = Dt/MU0/r_i2/delta_th;
            for( int k = 1; k < Nph + 1; k++ ){

                Hph[i][j][k] = Hph[i][j][k] - CHph1*(r_i3*Eth[i+1][j][k] - r_i1*Eth[i][j][k])
                        + CHph2*(Er[i][j+1][k] - Er[i][j][k]);
            }
        }
    }

}