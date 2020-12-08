#include "fdtd3d.h"

double*** memory_allocate3d( int d_1, int d_2, int d_3, double ini)
{
    double ***array;
    array = new double** [d_1];

    for( int i = 0; i < d_1; i++ ){
        array[i] = new double* [d_2];
        for( int j = 0; j < d_2; j++){
            array[i][j] = new double[d_3];
            for( int k = 0; k < d_3; k++){
                array[i][j][k] = ini;
            }
        }
    }

    return array;
}