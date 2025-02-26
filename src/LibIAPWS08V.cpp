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
        get_param_visc1_dimless(mdensity_in, tempearture_in) *
        get_param_visc2_dimless(mdensity_in, tempearture_in);

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

double Lib08V::get_param_visc2_dimless(double mdensity_in,
                                       double tempearture_in) {
    double visc2_dimless = 1.;

    bool in_range_crit =
        mdensity_in > range_visc2_mdensity_min_ &&
        mdensity_in < range_visc2_mdensity_max_ &&
        tempearture_in > range_visc2_temperature_min_ &&
        tempearture_in < range_visc2_temperature_max_;
    if (!in_range_crit) {
        return visc2_dimless;
    }

    if (ptr_lib95_ == NULL) {
        return visc2_dimless;
    }

    double mden_dimless = mdensity_in / mdensity_ref_;
    double temp_dimless = tempearture_in / temperature_ref_;

    double temp_R = coeff2_TR_ * temperature_ref_;
    double dchi =
        mden_dimless * (pressure_ref_ / mdensity_ref_) *
        (1. / ptr_lib95_->get_param_dpress_drho(mdensity_in,
                                                tempearture_in) -
         (coeff2_TR_ / temp_dimless) /
         ptr_lib95_->get_param_dpress_drho(mdensity_in, temp_R));
    if (dchi < 0.) {
        dchi = 0.;
    }

    double xi = coeff2_xi0_ *
        pow(dchi / coeff2_Gamma0_,
            coeff2_nu_ / coeff2_gamma_);
    double fac_q_C_xi = coeff2_q_C_ * xi;
    double fac_q_D_xi = coeff2_q_D_ * xi;

    double psi_D =
        acos(1. / sqrt(1. + fac_q_D_xi * fac_q_D_xi));

    double func_w =
        pow(fabs((fac_q_C_xi - 1.) / (fac_q_C_xi + 1.)), 0.5) *
        tan(0.5 * psi_D);
    double func_L = 0.;
    if (fac_q_C_xi > 1.) {
        func_L = log((1. + func_w) / (1. - func_w));
    } else {
        func_L = 2. * atan(fabs(func_w));
    }

    double func_Y = 0.;
    if (xi < 0.3817016416) {
        func_Y = fac_q_C_xi * pow(fac_q_D_xi, 5.) *
                 (1. - fac_q_C_xi + fac_q_C_xi * fac_q_C_xi -
                  765. * fac_q_D_xi * fac_q_D_xi / 504.) / 5.;
    } else {
        func_Y =
            sin(3. * psi_D) / 12. -
            sin(2. * psi_D) / (4. * fac_q_C_xi) +
            sin(psi_D) * (1. - 5. * fac_q_C_xi * fac_q_C_xi / 4.) /
            (fac_q_C_xi * fac_q_C_xi) -
            ((1. - 3. * fac_q_C_xi * fac_q_C_xi / 2.) * psi_D -
             pow(fabs(fac_q_C_xi * fac_q_C_xi - 1.), 1.5) * func_L) /
            (fac_q_C_xi * fac_q_C_xi * fac_q_C_xi);
    }

    visc2_dimless = exp(coeff2_x_mu_ * func_Y);

    if(flag_verbose_) {
        fprintf(stderr,
                "        xi = %.8e nm\n", xi);
        fprintf(stderr,
                "        visc2_dimless = %.8e\n", visc2_dimless);
    }

    return visc2_dimless;
}

} // end namespace IAPWS
