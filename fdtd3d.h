#define _USE_MATH_DEFINES
#include <cmath>

#define C0 (3.0e8)
#define MU0 (4.0*M_PI*1.0e-7)
#define EPS0 (1.0/MU0/C0/C0)
#define R0 (6370e3)
#define THETA0 (M_PI*0.5 - std::atan(50e3/R0))
#define E_Q (1.6e-19)
#define E_M (9.11e-31)

extern const int Nr;
extern const int Nth;
extern const int Nph;

extern const double delta_r;
extern const double delta_th;
extern const double delta_ph;
extern const double Dt;
extern const double inv_Dt;

double*** memory_allocate3d(
    int, int, int, double
);

void D_update(
    double ***nDr, double ***nDth, double ***nDph,
    double ***oDr, double ***oDth, double ***oDph,
    double ***Hr, double ***Hth, double ***Hph
);

void E_update(
    double ***Er, double ***Eth, double ***Eph,
    double ***nDr, double ***nDth, double ***nDph,
    double ***oDr, double ***oDth, double ***oDph
);

void H_update(
    double ***Er, double ***Eth, double ***Eph,
    double ***Hr, double ***Hth, double ***Hph
);

inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_th;};
inline double ph(double k){return k*delta_ph;};