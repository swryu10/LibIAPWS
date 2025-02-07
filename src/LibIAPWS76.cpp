#include<math.h>
#include"LibIAPWS76.h"

double LibIAPWS76::get_tension_surf(double temperature_in) {
    double tau = 1. - temperature_in / temperature_crit_;

    double sigma =
        coeff_B_ * pow(tau, coeff_mu_) * (1 + coeff_b_ * tau);

    return sigma;
}
