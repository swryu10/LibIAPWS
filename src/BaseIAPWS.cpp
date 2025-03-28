#include<math.h>
#include"BaseIAPWS.h"

namespace IAPWS {

bool flag_verbose_ = false;

int n_iter_max_ = 1000000;
double eps_precision_ = 5. * 1.0e-8;

void get_complex_product(double *z_1in,
                         double *z_2in,
                         double *z_out) {
    z_out[0] = z_1in[0] * z_2in[0] -
               z_1in[1] * z_2in[1];

    z_out[1] = z_1in[0] * z_2in[1] +
               z_1in[1] * z_2in[0];

    return;
}

void get_complex_inverse(double *z_in,
                         double *z_out) {
    double z_sqr = z_in[0] * z_in[0] +
                   z_in[1] * z_in[1];

    z_out[0] = z_in[0] / z_sqr;
    z_out[1] = -z_in[1] / z_sqr;

    return;
}

void get_complex_log(double *z_in,
                     double *z_out) {

    double z_abs = sqrt(z_in[0] * z_in[0] +
                        z_in[1] * z_in[1]);

    z_out[0] = log(z_abs);
    z_out[1] = acos(z_in[0] / z_abs);
    if (z_in[1] < 0.) {
        z_out[1] = -z_out[1];
    }

    return;
}

double get_pow_int(double x, int n) {
    int n_abs = n;
    double fac_mult = x;
    if (n < 0) {
        n_abs = -n_abs;
        fac_mult = 1. / x;
    }

    double x_pow_n = 1.;
    for (int i = 0; i < n_abs; i++) {
        x_pow_n *= fac_mult;
    }

    return x_pow_n;
}

} // end namespace IAPWS
