#include<math.h>
#include"LibIAPWS06.h"

namespace IAPWS {

void Lib06::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R10-06 (2009)\n");
    fprintf(ptr_fout, "Revised Release on the Equation of State ");
    fprintf(ptr_fout, "2006 for H2O Ice Ih\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib06::get_param_g0(double pressure) {
    double press_dpi =
        (pressure - pressure_norm_) / pressure_trip_;

    double g0 = 0.;
    for (int l = 0; l <= 4; l++) {
        double pow_dpi_l = 1.;
        for (int ll = 0; ll < l; ll++) {
            pow_dpi_l *= press_dpi;
        }

        g0 += coeff_g0_[l] * pow_dpi_l;
    }

    return g0;
}

double Lib06::get_param_dg0_dp(double pressure) {
    double press_dpi =
        (pressure - pressure_norm_) / pressure_trip_;

    double dg0_dp = 0.;
    for (int l = 1; l <= 4; l++) {
        double pow_dpi_l = 1.;
        for (int ll = 1; ll < l; ll++) {
            pow_dpi_l *= press_dpi;
        }

        dg0_dp += coeff_g0_[l] * pow_dpi_l *
                  static_cast<double>(l) / pressure_trip_;
    }

    return dg0_dp;
}

double Lib06::get_param_d2g0_dp_dp(double pressure) {
    double press_dpi =
        (pressure - pressure_norm_) / pressure_trip_;

    double d2g0_dp_dp = 0.;
    for (int l = 2; l <= 4; l++) {
        double pow_dpi_l = 1.;
        for (int ll = 2; ll < l; ll++) {
            pow_dpi_l *= press_dpi;
        }

        d2g0_dp_dp += coeff_g0_[l] * pow_dpi_l *
                      static_cast<double>(l * (l - 1)) /
                      (pressure_trip_ * pressure_trip_);
    }

    return d2g0_dp_dp;
}

void Lib06::get_param_r(double pressure, int k,
                        double *r_out) {
    if (k == 1) {
        r_out[0] = coeff_r1_[0];
        r_out[1] = coeff_r1_[1];
    } else if (k == 2) {
        double press_dpi =
            (pressure - pressure_norm_) / pressure_trip_;

        r_out[0] = 0.;
        r_out[1] = 0.;

        for (int l = 0; l <= 2; l++) {
            double pow_dpi_l = 1.;
            for (int ll = 0; ll < l; ll++) {
                pow_dpi_l *= press_dpi;
            }

            r_out[0] += coeff_r2_[l][0] * pow_dpi_l;
            r_out[1] += coeff_r2_[l][1] * pow_dpi_l;
        }
    } else {
        r_out[0] = 0.;
        r_out[1] = 0.;
    }

    return;
}

void Lib06::get_param_dr_dp(double pressure, int k,
                            double *dr_dp_out) {
    if (k == 1) {
        dr_dp_out[0] = 0.;
        dr_dp_out[1] = 0.;
    } else if (k == 2) {
        double press_dpi =
            (pressure - pressure_norm_) / pressure_trip_;

        dr_dp_out[0] = 0.;
        dr_dp_out[1] = 0.;

        for (int l = 1; l <= 2; l++) {
            double pow_dpi_l = 1.;
            for (int ll = 1; ll < l; ll++) {
                pow_dpi_l *= press_dpi;
            }

            double fac_pow_dpi =
                pow_dpi_l * static_cast<double>(l) /
                pressure_trip_;

            dr_dp_out[0] += coeff_r2_[l][0] * fac_pow_dpi;
            dr_dp_out[1] += coeff_r2_[l][1] * fac_pow_dpi;
        }
    } else {
        dr_dp_out[0] = 0.;
        dr_dp_out[1] = 0.;
    }

    return;
}

void Lib06::get_param_d2r_dp_dp(double pressure, int k,
                                double *d2r_dp_dp_out) {
    if (k == 1) {
        d2r_dp_dp_out[0] = 0.;
        d2r_dp_dp_out[1] = 0.;
    } else if (k == 2) {
        /*
        double press_dpi =
            (pressure - pressure_norm_) / pressure_trip_;

        d2r_dp_dp_out[0] = 0.;
        d2r_dp_dp_out[1] = 0.;

        for (int l = 2; l <= 2; l++) {
            double pow_dpi_l = 1.;
            for (int ll = 2; ll < l; ll++) {
                pow_dpi_l *= press_dpi;
            }

            double fac_pow_dpi =
                pow_dpi_l * static_cast<double>(l * (l - 1)) /
                (pressure_trip_ * pressure_trip_);

            d2r_dp_dp_out[0] += coeff_r2_[l][0] * fac_pow_dpi;
            d2r_dp_dp_out[1] += coeff_r2_[l][1] * fac_pow_dpi;
        }
        */

        d2r_dp_dp_out[0] = 2. * coeff_r2_[2][0] /
                           (pressure_trip_ * pressure_trip_);
        d2r_dp_dp_out[1] = 2. * coeff_r2_[2][1] /
                           (pressure_trip_ * pressure_trip_);
    } else {
        d2r_dp_dp_out[0] = 0.;
        d2r_dp_dp_out[1] = 0.;
    }

    return;
}

void Lib06::get_param_q(double temperature, int k,
                        double *q_out) {
    q_out[0] = 0.;
    q_out[1] = 0.;

    if (k < 1 || k > 2) {
        return;
    }

    double tau = temperature / temperature_trip_;

    double *t_pos = new double[2];
    double *t_neg = new double[2];

    double *ln_t_pos = new double[2];
    double *ln_t_neg = new double[2];
    double *ln_t_org = new double[2];

    double *t_ln_t_pos = new double[2];
    double *t_ln_t_neg = new double[2];
    double *t_ln_t_org = new double[2];

    double *t_inv = new double[2];

    t_pos[0] = coeff_t_[k][0] + tau;
    t_pos[1] = coeff_t_[k][1];
    get_complex_log(t_pos, ln_t_pos);
    get_complex_product(t_pos, ln_t_pos, t_ln_t_pos);

    t_neg[0] = coeff_t_[k][0] - tau;
    t_neg[1] = coeff_t_[k][1];
    get_complex_log(t_neg, ln_t_neg);
    get_complex_product(t_neg, ln_t_neg, t_ln_t_neg);

    get_complex_log(coeff_t_[k], ln_t_org);
    get_complex_product(coeff_t_[k], ln_t_org, t_ln_t_org);

    get_complex_inverse(coeff_t_[k], t_inv);

    q_out[0] =
        t_ln_t_pos[0] + t_ln_t_neg[0] - 2. * t_ln_t_org[0] -
        tau * tau * t_inv[0];
    q_out[1] =
        t_ln_t_pos[1] + t_ln_t_neg[1] - 2. * t_ln_t_org[1] -
        tau * tau * t_inv[1];

    delete [] t_pos;
    delete [] t_neg;

    delete [] ln_t_pos;
    delete [] ln_t_neg;
    delete [] ln_t_org;

    delete [] t_ln_t_pos;
    delete [] t_ln_t_neg;
    delete [] t_ln_t_org;

    delete [] t_inv;

    return;
}

void Lib06::get_param_dq_dT(double temperature, int k,
                            double *dq_dT_out) {
    dq_dT_out[0] = 0.;
    dq_dT_out[1] = 0.;

    if (k < 1 || k > 2) {
        return;
    }

    double tau = temperature / temperature_trip_;

    double *t_pos = new double[2];
    double *t_neg = new double[2];

    double *ln_t_pos = new double[2];
    double *ln_t_neg = new double[2];

    double *t_inv = new double[2];

    t_pos[0] = coeff_t_[k][0] + tau;
    t_pos[1] = coeff_t_[k][1];
    get_complex_log(t_pos, ln_t_pos);

    t_neg[0] = coeff_t_[k][0] - tau;
    t_neg[1] = coeff_t_[k][1];
    get_complex_log(t_neg, ln_t_neg);

    get_complex_inverse(coeff_t_[k], t_inv);

    dq_dT_out[0] =
        (ln_t_pos[0] - ln_t_neg[0] - 2. * tau * t_inv[0]) /
        temperature_trip_;
    dq_dT_out[1] =
        (ln_t_pos[1] - ln_t_neg[1] - 2. * tau * t_inv[1]) /
        temperature_trip_;

    delete [] t_pos;
    delete [] t_neg;

    delete [] ln_t_pos;
    delete [] ln_t_neg;

    delete [] t_inv;

    return;
}

void Lib06::get_param_d2q_dT_dT(double temperature, int k,
                                double *d2q_dT_dT_out) {
    d2q_dT_dT_out[0] = 0.;
    d2q_dT_dT_out[1] = 0.;

    if (k < 1 || k > 2) {
        return;
    }

    double tau = temperature / temperature_trip_;

    double *t_pos = new double[2];
    double *t_neg = new double[2];

    double *t_pos_inv = new double[2];
    double *t_neg_inv = new double[2];

    double *t_inv = new double[2];

    t_pos[0] = coeff_t_[k][0] + tau;
    t_pos[1] = coeff_t_[k][1];
    get_complex_inverse(t_pos, t_pos_inv);

    t_neg[0] = coeff_t_[k][0] - tau;
    t_neg[1] = coeff_t_[k][1];
    get_complex_inverse(t_neg, t_neg_inv);

    get_complex_inverse(coeff_t_[k], t_inv);

    d2q_dT_dT_out[0] =
        (t_pos_inv[0] + t_neg_inv[0] - 2. * t_inv[0]) /
        (temperature_trip_ * temperature_trip_);
    d2q_dT_dT_out[1] =
        (t_pos_inv[1] + t_neg_inv[1] - 2. * t_inv[1]) /
        (temperature_trip_ * temperature_trip_);

    delete [] t_pos;
    delete [] t_neg;

    delete [] t_pos_inv;
    delete [] t_neg_inv;

    delete [] t_inv;

    return;
}

double Lib06::get_param_g(double temperature,
                          double pressure) {
    double tau = temperature / temperature_trip_;

    double g = 0.;

    g += get_param_g0(pressure);

    double entropy_res_now = 0.;
    if (flag_s_absolute_) {
        entropy_res_now = entropy_res_absol_;
    } else {
        entropy_res_now = entropy_res_Lib95_;
    }

    g -= entropy_res_now * temperature_trip_ * tau;

    double *func_r = new double[2];
    double *func_q = new double[2];

    for (int k = 1; k <= 2; k++) {
        get_param_r(pressure, k, func_r);
        get_param_q(temperature, k, func_q);

        g += temperature_trip_ *
             (func_r[0] * func_q[0] -
              func_r[1] * func_q[1]);
    }

    delete [] func_r;
    delete [] func_q;

    return g;
}

double Lib06::get_param_dg_dT(double temperature,
                              double pressure) {
    double tau = temperature / temperature_trip_;

    double dg_dT = 0.;

    double entropy_res_now = 0.;
    if (flag_s_absolute_) {
        entropy_res_now = entropy_res_absol_;
    } else {
        entropy_res_now = entropy_res_Lib95_;
    }

    dg_dT -= entropy_res_now;

    double *func_r = new double[2];
    double *func_q = new double[2];

    for (int k = 1; k <= 2; k++) {
        get_param_r(pressure, k, func_r);
        get_param_dq_dT(temperature, k, func_q);

        dg_dT += temperature_trip_ *
                 (func_r[0] * func_q[0] -
                  func_r[1] * func_q[1]);
    }

    delete [] func_r;
    delete [] func_q;

    return dg_dT;
}

double Lib06::get_param_dg_dp(double temperature,
                              double pressure) {
    double tau = temperature / temperature_trip_;

    double dg_dp = 0.;

    dg_dp += get_param_dg0_dp(pressure);

    double *func_r = new double[2];
    double *func_q = new double[2];

    get_param_dr_dp(pressure, 2, func_r);
    get_param_q(temperature, 2, func_q);

    dg_dp += temperature_trip_ *
             (func_r[0] * func_q[0] -
              func_r[1] * func_q[1]);

    delete [] func_r;
    delete [] func_q;

    return dg_dp;
}

double Lib06::get_param_d2g_dT_dT(double temperature,
                                  double pressure) {
    double tau = temperature / temperature_trip_;

    double d2g_dT_dT = 0.;

    double *func_r = new double[2];
    double *func_q = new double[2];

    for (int k = 1; k <= 2; k++) {
        get_param_r(pressure, k, func_r);
        get_param_d2q_dT_dT(temperature, k, func_q);

        d2g_dT_dT += temperature_trip_ *
                     (func_r[0] * func_q[0] -
                      func_r[1] * func_q[1]);
    }

    delete [] func_r;
    delete [] func_q;

    return d2g_dT_dT;
}

double Lib06::get_param_d2g_dT_dp(double temperature,
                                  double pressure) {
    double tau = temperature / temperature_trip_;

    double d2g_dT_dp = 0.;

    double *func_r = new double[2];
    double *func_q = new double[2];

    get_param_dr_dp(pressure, 2, func_r);
    get_param_dq_dT(temperature, 2, func_q);

    d2g_dT_dp += temperature_trip_ *
                 (func_r[0] * func_q[0] -
                  func_r[1] * func_q[1]);

    delete [] func_r;
    delete [] func_q;

    return d2g_dT_dp;
}

double Lib06::get_param_d2g_dp_dp(double temperature,
                                  double pressure) {
    double tau = temperature / temperature_trip_;

    double d2g_dp_dp = 0.;

    d2g_dp_dp += get_param_d2g0_dp_dp(pressure);

    double *func_r = new double[2];
    double *func_q = new double[2];

    get_param_d2r_dp_dp(pressure, 2, func_r);
    get_param_q(temperature, 2, func_q);

    d2g_dp_dp += temperature_trip_ *
                 (func_r[0] * func_q[0] -
                  func_r[1] * func_q[1]);

    delete [] func_r;
    delete [] func_q;

    return d2g_dp_dp;
}

double Lib06::get_param_mdensity(double temperature_in,
                                 double pressure_in) {
    double mden =
        1. / get_param_dg_dp(temperature_in,
                             pressure_in);

    return mden;
}

double Lib06::get_param_entropy(double temperature_in,
                                double pressure_in) {
    double entropy =
        -get_param_dg_dT(temperature_in,
                         pressure_in);

    return entropy;
}

double Lib06::get_param_heat_c_p(double temperature_in,
                                 double pressure_in) {
    double c_p =
        -temperature_in *
         get_param_d2g_dT_dT(temperature_in,
                             pressure_in);

    return c_p;
}

double Lib06::get_param_enthalpy(double temperature_in,
                                 double pressure_in) {
    double enthalpy =
        get_param_g(temperature_in, pressure_in) -
        temperature_in * get_param_dg_dT(temperature_in,
                                         pressure_in);

    return enthalpy;
}

double Lib06::get_param_erg_int(double temperature_in,
                                double pressure_in) {
    double erg_int =
        get_param_g(temperature_in, pressure_in) -
        temperature_in * get_param_dg_dT(temperature_in,
                                         pressure_in) -
        pressure_in * get_param_dg_dp(temperature_in,
                                      pressure_in);

    return erg_int;
}

double Lib06::get_param_f(double temperature_in,
                          double pressure_in) {
    double f =
        get_param_g(temperature_in, pressure_in) -
        pressure_in * get_param_dg_dp(temperature_in,
                                      pressure_in);

    return f;
}

double Lib06::get_param_coeff_alpha(double temperature_in,
                                    double pressure_in) {
    double alpha =
        get_param_d2g_dT_dp(temperature_in,
                            pressure_in) /
        get_param_dg_dp(temperature_in,
                        pressure_in);

    return alpha;
}

double Lib06::get_param_coeff_beta(double temperature_in,
                                   double pressure_in) {
    double beta =
        -get_param_d2g_dT_dp(temperature_in,
                             pressure_in) /
         get_param_d2g_dp_dp(temperature_in,
                             pressure_in);;

    return beta;
}

double Lib06::get_param_comp_kappa_T(double temperature_in,
                                     double pressure_in) {
    double kappa_T =
        -get_param_d2g_dp_dp(temperature_in,
                             pressure_in) /
         get_param_dg_dp(temperature_in,
                         pressure_in);

    return kappa_T;
}

double Lib06::get_param_comp_kappa_s(double temperature_in,
                                     double pressure_in) {
    double g_p = get_param_dg_dp(temperature_in,
                                 pressure_in);
    double g_TT = get_param_d2g_dT_dT(temperature_in,
                                      pressure_in);
    double g_Tp = get_param_d2g_dT_dp(temperature_in,
                                      pressure_in);
    double g_pp = get_param_d2g_dp_dp(temperature_in,
                                      pressure_in);

    double kappa_s =
        (g_Tp * g_Tp - g_TT * g_pp) / (g_p * g_TT);

    return kappa_s;
}

} // end namespace IAPWS
