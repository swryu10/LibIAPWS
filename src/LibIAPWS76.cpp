#include<math.h>
#include"LibIAPWS76.h"

namespace IAPWS {

void Lib76::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R1-76 (2014)\n");
    fprintf(ptr_fout, "Surface Tension of Ordinary Water Substance\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib76::get_tension_surf(double temperature_in) {
    double tau = 1. - temperature_in / temperature_crit_;

    double sigma =
        coeff_B_ * pow(tau, coeff_mu_) * (1 + coeff_b_ * tau);

    return sigma;
}

} // end namespace IAPWS
