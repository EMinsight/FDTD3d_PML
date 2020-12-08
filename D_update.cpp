#define _USE_MATH_DEFINES
#include <cmath>
#include "fdtd3d.h"

void D_update( 
    double ***nDr, double ***nDth, double***nDph,
     double ***oDr, double ***oDth, double ***oDph,
    double ***Hr, double ***Hth, double ***Hph )
    {
        double CDr1, CDr2;
        double CDth1, CDth2;
        double CDph1, CDph2;
        
        double r_i1, r_i2, r_i3;
        double si_th1, si_th2, si_th3;

        for( int i = 0; i < Nr; i++ ){
            r_i2 = dist(i+0.5);
            for( int j = 1; j < Nth; j++ ){
                si_th1 = std::sin(th(j-0.5));
                si_th2 = std::sin(th(j));
                si_th3 = std::sin(th(j+0.5));
                
                CDr1 = Dt/r_i2/si_th2/delta_th;
                CDr2 = Dt/r_i2/si_th2/delta_ph;
                
                for( int k = 1; k < Nph; k++ ){

                    nDr[i][j][k] = oDr[i][j][k] + CDr1*( si_th3*Hph[i][j][k] - si_th1*Hph[i][j-1][k] )
                                - CDr2*( Hth[i][j][k] - Hth[i][j][k-1] );
                }
            }
        }

        for( int i = 1; i < Nr; i++ ){
            r_i1 = dist(i-0.5);
            r_i2 = dist(i);
            r_i3 = dist(i+0.5);
            for( int j = 0; j < Nth; j++ ){
                si_th3 = std::sin(th(j+0.5));
                for( int k = 1; k < Nph; k++ ){
                    CDth1 = Dt/r_i2/si_th3/delta_ph;
                    CDth2 = Dt/r_i2/delta_r;

                    nDth[i][j][k] = oDth[i][j][k] + CDth1*( Hr[i][j][k] - Hr[i][j][k-1] )
                                - CDth2*( r_i3*Hph[i][j][k] - r_i1*Hph[i-1][j][k] );
                }
            }
        }

        for( int i = 1; i < Nr; i++ ){
            r_i1 = dist(i-0.5);
            r_i2 = dist(i);
            r_i3 = dist(i+0.5);
            for( int j = 1; j < Nth; j++ ){
                for( int k = 0; k < Nph; k++ ){
                    CDph1 = Dt/r_i2/delta_r;
                    CDph2 = Dt/r_i2/delta_th;

                    nDph[i][j][k] = oDph[i][j][k] + CDph1*( r_i3*Hth[i][j][k] - r_i1*Hth[i-1][j][k] )
                                - CDph2*( Hr[i][j][k] - Hr[i][j-1][k] );
                }
            }
        }
        
    }