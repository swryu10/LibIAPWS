#include<math.h>
#include"LibIAPWS08V.h"

namespace IAPWS {

void Lib08V::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R12-08 (2008)\n");
    fprintf(ptr_fout, "Release on the IAPWS Formulation 2008\n");
    fprintf(ptr_fout,
            "for the Viscosity of Ordinary Water Substance\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib08V::get_param_viscosity(double mdensity_in,
                                   double tempearture_in) {
    double viscosity_out = viscosity_ref_ *
        get_param_visc0_dimless(tempearture_in) *
        get_param_visc1_dimless(mdensity_in, tempearture_in);

    return viscosity_out;
}

double Lib08V::get_param_visc0_dimless(double tempearture_in) {
    double temp_dimless = tempearture_in / temperature_ref_;

    double nomin = 100. * sqrt(temp_dimless);

    double denom = 0.;
    for (int i = 0; i <= 3; i++) {
        double temp_pow_i = 1.;
        for (int j = 0; j < i; j++) {
            temp_pow_i *= temp_dimless;
        }

        denom += coeff0_H_[i] / temp_pow_i;
    }

    double visc0_dimless = nomin / denom;

    if (flag_verbose_) {
        fprintf(stderr,
                "        visc0_dimless = %e\n", visc0_dimless);
    }

    return visc0_dimless;
}

double Lib08V::get_param_visc1_dimless(double mdensity_in,
                                       double tempearture_in) {
    double mden_dimless = mdensity_in / mdensity_ref_;
    double temp_dimless = tempearture_in / temperature_ref_;

    double *list_dtemp_pow_i = new double[6];
    list_dtemp_pow_i[0] = 1.;
    for (int i = 0; i < 5; i++) {
        list_dtemp_pow_i[i + 1] =
            list_dtemp_pow_i[i] * (1. / temp_dimless - 1.);
    }

    double *list_dmden_pow_j = new double[7];
    list_dmden_pow_j[0] = 1.;
    for (int j = 0; j < 6; j++) {
        list_dmden_pow_j[j + 1] =
            list_dmden_pow_j[j] * (mden_dimless - 1.);
    }

    double exponent = 0.;
    for (int i = 0; i <= 5; i++) {
        for (int j = 0; j <= 6; j++) {
            exponent += coeff1_H_[i][j] *
                list_dtemp_pow_i[i] * list_dmden_pow_j[j];
        }
    }

    double visc1_dimless = exp(mden_dimless * exponent);

    delete [] list_dtemp_pow_i;
    delete [] list_dmden_pow_j;

    if(flag_verbose_) {
        fprintf(stderr,
                "        visc1_dimless = %e\n", visc1_dimless);
    }

    return visc1_dimless;
}

} // end namespace IAPWS
