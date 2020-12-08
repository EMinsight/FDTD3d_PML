#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include "fdtd3d.h"

const int Nr{100};
const int Nth{100};
const int Nph{1000};

constexpr double R_r{100.0e3};

const double delta_r{ R_r/(double)Nr };
const double delta_th{ 1.0e3/(double)R0 };
const double delta_ph{ 1.0e3/(double)R0 };
const double Dt { double( 0.99/C0/std::sqrt(1.0/delta_r/delta_r
 + 1.0/R0/R0/delta_th/delta_th
 + 1.0/R0/R0/std::sin(THETA0)/std::sin(THETA0)/delta_ph/delta_ph) ) };
const double inv_Dt{ 1.0/Dt };
const double sigma_t{ 7.0*Dt };
const double t0{ 6.0*sigma_t };

// center point //
const int i_0{ Nr/2 };
const int j_0{ Nth/2 };
const int k_0{ Nph/2 };

// Source point //
const int i_s{1};
const int j_s{50};
const int k_s{100};

// Receive Point //
const int i_r{1};
const int j_r{50};
const int k_r{ Nph - 50 };

void calc_fdtd(void)
{   
    double ***Hr, ***Hth, ***Hph;
    Hr = memory_allocate3d( Nr+1, Nth, Nph, 0.0);
    Hth = memory_allocate3d( Nr, Nth+1, Nph, 0.0 );
    Hph = memory_allocate3d( Nr, Nth, Nph+1, 0.0 );

    double ***Er, ***Eth, ***Eph;
    Er = memory_allocate3d( Nr, Nth+1, Nph+1, 0.0 );
    Eth = memory_allocate3d( Nr+1, Nth, Nph+1, 0.0 );
    Eph = memory_allocate3d( Nr+1, Nth+1, Nph, 0.0 );

    double ***nDr, ***nDth, ***nDph;
    nDr = memory_allocate3d( Nr, Nth+1, Nph+1, 0.0 );
    nDth = memory_allocate3d( Nr+1, Nth, Nph+1, 0.0 );
    nDph = memory_allocate3d( Nr+1, Nth+1, Nph, 0.0 );

    double ***oDr, ***oDth, ***oDph;
    oDr = memory_allocate3d( Nr, Nth+1, Nph+1, 0.0 );
    oDth = memory_allocate3d( Nr+1, Nth, Nph+1, 0.0 );
    oDph = memory_allocate3d( Nr+1, Nth+1, Nph, 0.0 );

    int time_step = 1700;

    std::chrono::system_clock::time_point start
        = std::chrono::system_clock::now();

    for( int n = 0; n < time_step; n++ ){

        //std::cout  << n << " / " << time_step << "\n";

        double t = double(n - 0.5)*Dt;

        double J = -((t - t0)/sigma_t/sigma_t/delta_r/(dist(i_s + 0.5)*delta_th)/(dist(i_s + 0.5)*delta_ph))
            *std::exp(-std::pow(t - t0, 2.0)/2.0/std::pow(sigma_t, 2.0));
        
        Eth[i_s][j_s][k_s] = Eth[i_s][j_s][k_s] + J;

        //std::cout << Eth[i_s][j_s][k_s + 1] << "\n";

        D_update(
            nDr, nDth, nDph,
            oDr, oDth, oDph,
            Hr, Hth, Hph
        );

        E_update(
            Er, Eth, Eph,
            nDr, nDth, nDph,
            oDr, oDth, oDph
        );

        H_update(
            Er, Eth, Eph,
            Hr, Hth, Hph
        );

        //std::cout << n << " : " << Eth[i_s][j_s][k_s] << "\n";

        /*std::string filename = "./result/ez_" + std::to_string(n) + ".dat";
        std::ofstream ofs(filename.c_str());
        for( int i = 0; i < Nr; i++ ){
            for( int k = 1; k < Nph; k++ ){
                ofs << k << " " << i << " " << Eth[i][j_s][k] << "\n";
            }
            ofs << "\n";
        }
        ofs.close();*/

    }

    std::chrono::system_clock::time_point end
            = std::chrono::system_clock::now();

    double total_time;
    total_time = std::chrono::duration_cast <std::chrono::milliseconds>
    (end - start).count();

    std::cout << "elapsed time : " << total_time*1.e-3 << "\n";
    std::cout << "Source point : " << Er[i_s][j_s][k_s] << "\n";
    std::cout << "Receive point : " << Er[i_r][j_r][k_r] << "\n";

}