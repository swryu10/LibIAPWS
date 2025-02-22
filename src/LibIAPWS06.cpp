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

bool Lib06::find_state_coex(Lib95 *ptr_lib95eos,
                            double temperature_in,
                            double &pressure_out,
                            double &mden_vap_out,
                            double &mden_ice_out) {
    double press_now = pressure_out;
    double mden_vap = mden_vap_out;
    double mden_ice = mden_ice_out;

    int i_iter = 0;
    bool found_state = false;
    while (!found_state) {
        i_iter += 1;

        double press_prev = press_now;

        bool found_root_vap =
            ptr_lib95eos->find_root_mdensity(temperature_in,
                                             press_prev,
                                             mden_vap);

        if (!found_root_vap) {
            break;
        }

        mden_ice =
            get_param_mdensity(temperature_in,
                               press_prev);

        double g_vap =
            ptr_lib95eos->get_param_g(mden_vap,
                                      temperature_in);
        double g_ice =
            get_param_g(temperature_in, press_prev);

        if (fabs(g_vap - g_ice) <
            0.5 * eps_precision_ * fabs(g_vap + g_ice)) {
            pressure_out = press_prev;
            mden_vap_out = mden_vap;
            mden_ice_out = mden_ice;

            found_state = true;
        }

        if (i_iter > n_iter_max_) {
            break;
        }

        press_now = press_prev +
            0.1 * (g_ice - g_vap) * mden_vap;
    }

    /*
    fprintf(stderr, "      temperature = %.8e degK, ", temperature_in);
    fprintf(stderr, "      pressure = %.8e Pa, ", pressure_out);
    fprintf(stderr, "      mden_vap = %.8e kg / m^3, ", mden_vap_out);
    fprintf(stderr, "      mden_ice = %.8e kg / m^3\n", mden_ice_out);
    */

    return found_state;
}

void Lib06::make_tab_coex(Lib95 *ptr_lib95eos,
                          int nbin_in,
                          double temperature_min) {
    reset_tab_coex();

    if (temperature_min > temperature_trip_) {
        return;
    }

    nbin_coex_ = nbin_in;
    temperature_coex_min_ = temperature_min;
    temperature_coex_max_ = temperature_trip_;

    double delta_temp =
        (temperature_coex_max_ - temperature_coex_min_) /
        static_cast<double>(nbin_coex_);

    have_tab_coex_ = true;
    tab_coex_temperature_ = new double[nbin_coex_ + 1];
    tab_coex_pressure_ = new double[nbin_coex_ + 1];
    tab_coex_mden_vap_ = new double[nbin_coex_ + 1];
    tab_coex_mden_ice_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_ice_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_ice_ = new double[nbin_coex_ + 1];

    bool made_tab = true;
    if (flag_verbose_) {
        fprintf(stderr, "    temperature (degK), ");
        fprintf(stderr, "    pressure (Pa), ");
        fprintf(stderr, "    mden_vap (kg / m^3), ");
        fprintf(stderr, "    mden_ice (kg / m^3)\n");
    }
    for (int it = nbin_coex_; it >= 0; it--) {
        tab_coex_temperature_[it] =
            temperature_coex_min_ +
            delta_temp * static_cast<double>(it);

        if (it == nbin_coex_) {
            tab_coex_pressure_[it] = pressure_trip_num_;
            tab_coex_mden_vap_[it] = 4.85457548e-03;
            tab_coex_mden_ice_[it] =
                get_param_mdensity(tab_coex_temperature_[it],
                                   tab_coex_pressure_[it]);
        } else {
            tab_coex_pressure_[it] = tab_coex_pressure_[it + 1];
            tab_coex_mden_vap_[it] = tab_coex_mden_vap_[it + 1];
            tab_coex_mden_ice_[it] = tab_coex_mden_ice_[it + 1];
        }

        bool found_state =
            find_state_coex(ptr_lib95eos,
                            tab_coex_temperature_[it],
                            tab_coex_pressure_[it],
                            tab_coex_mden_vap_[it],
                            tab_coex_mden_ice_[it]);

        if (!found_state) {
            made_tab = false;
            break;
        }

        tab_coex_enthalpy_vap_[it] =
            ptr_lib95eos->get_param_enthalpy(tab_coex_mden_vap_[it],
                                             tab_coex_temperature_[it]);
        tab_coex_enthalpy_ice_[it] =
            get_param_enthalpy(tab_coex_temperature_[it],
                               tab_coex_pressure_[it]);

        tab_coex_entropy_vap_[it] =
            ptr_lib95eos->get_param_entropy(tab_coex_mden_vap_[it],
                                            tab_coex_temperature_[it]);
        tab_coex_entropy_ice_[it] =
            get_param_entropy(tab_coex_temperature_[it],
                              tab_coex_pressure_[it]);

        if (flag_verbose_) {
            fprintf(stderr, "      %.8e",
                    tab_coex_temperature_[it]);
            fprintf(stderr, "      %.8e",
                    tab_coex_pressure_[it]);
            fprintf(stderr, "      %.8e",
                    tab_coex_mden_vap_[it]);
            fprintf(stderr, "      %.8e\n",
                    tab_coex_mden_ice_[it]);
        }
    }

    if (!made_tab) {
        reset_tab_coex();
    }

    set_cspline_coex();

    return;
}

void Lib06::export_tab_coex(char *filename) {
    if (!have_tab_coex_) {
        return;
    }

    FILE *ptr_fout;
    ptr_fout = fopen(filename, "w");
    if (ptr_fout == NULL) {
        return;
    }

    fprintf(ptr_fout, "# table for coexisting phases in IAPWS06\n");
    fprintf(ptr_fout, "nbin_coex    %d\n", nbin_coex_);
    fprintf(ptr_fout, "temperature_coex_min    %e\n",
            temperature_coex_min_);
    fprintf(ptr_fout, "temperature_coex_max    %e\n",
            temperature_coex_max_);
    fprintf(ptr_fout, "# temperature (degK)");
    fprintf(ptr_fout, "    pressure (Pa)");
    fprintf(ptr_fout, "    mden_vap (kg / m^3)");
    fprintf(ptr_fout, "    mden_ice (kg / m^3)");
    fprintf(ptr_fout, "    enthalpy_vap (J / kg)");
    fprintf(ptr_fout, "    enthalpy_ice (J / kg)");
    fprintf(ptr_fout, "    entropy_vap (J / kg / degK)");
    fprintf(ptr_fout, "    entropy_ice (J / kg / degK)");
    fprintf(ptr_fout, "\n");

    for (int it = 0; it <= nbin_coex_; it++) {
        fprintf(ptr_fout, "    %.8e",
                tab_coex_temperature_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_pressure_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_mden_vap_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_mden_ice_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_enthalpy_vap_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_enthalpy_ice_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_entropy_vap_[it]);
        fprintf(ptr_fout, "    %.8e\n",
                tab_coex_entropy_ice_[it]);
    }

    fclose(ptr_fout);

    return;
}

void Lib06::import_tab_coex(char *filename) {
    reset_tab_coex();

    FILE *ptr_fin;
    ptr_fin = fopen(filename, "r");
    if (ptr_fin == NULL) {
        return;
    }

    char line_current[1000];

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        if (flag_verbose_) {
            fprintf(stderr, "%s", line_current);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %d", buffer, &nbin_coex_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %d\n", buffer, nbin_coex_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %lf", buffer, &temperature_coex_min_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %f\n", buffer, temperature_coex_min_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %lf", buffer, &temperature_coex_max_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %f\n", buffer, temperature_coex_max_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        if (flag_verbose_) {
            fprintf(stderr, "%s", line_current);
        }
    } else {
        return;
    }

    have_tab_coex_ = true;
    tab_coex_temperature_ = new double[nbin_coex_ + 1];
    tab_coex_pressure_ = new double[nbin_coex_ + 1];
    tab_coex_mden_vap_ = new double[nbin_coex_ + 1];
    tab_coex_mden_ice_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_ice_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_ice_ = new double[nbin_coex_ + 1];

    bool made_tab = true;
    for (int it = 0; it <= nbin_coex_; it++) {
        if (fgets(line_current, sizeof(line_current), ptr_fin) == NULL) {
            made_tab = false;
            break;
        }

        sscanf(line_current, "%lf %lf %lf %lf %lf %lf %lf %lf",
               &tab_coex_temperature_[it],
               &tab_coex_pressure_[it],
               &tab_coex_mden_vap_[it],
               &tab_coex_mden_ice_[it],
               &tab_coex_enthalpy_vap_[it],
               &tab_coex_enthalpy_ice_[it],
               &tab_coex_entropy_vap_[it],
               &tab_coex_entropy_ice_[it]);

        if (flag_verbose_) {
            fprintf(stderr, "    %.8e",
                    tab_coex_temperature_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_pressure_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_mden_vap_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_mden_ice_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_enthalpy_vap_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_enthalpy_ice_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_entropy_vap_[it]);
            fprintf(stderr, "    %.8e\n",
                    tab_coex_entropy_ice_[it]);
        }
    }

    fclose(ptr_fin);

    if (!made_tab) {
        reset_tab_coex();
    }

    set_cspline_coex();

    return;
}

void Lib06::set_cspline_coex() {
    if (!have_tab_coex_) {
        return;
    }

    csp_coex_pressure_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_pressure_);
    csp_coex_mden_vap_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_mden_vap_);
    csp_coex_mden_ice_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_mden_ice_);
    csp_coex_enthalpy_vap_.init(nbin_coex_,
                                tab_coex_temperature_,
                                tab_coex_enthalpy_vap_);
    csp_coex_enthalpy_ice_.init(nbin_coex_,
                                tab_coex_temperature_,
                                tab_coex_enthalpy_ice_);
    csp_coex_entropy_vap_.init(nbin_coex_,
                               tab_coex_temperature_,
                               tab_coex_entropy_vap_);
    csp_coex_entropy_ice_.init(nbin_coex_,
                               tab_coex_temperature_,
                               tab_coex_entropy_ice_);

    return;
}

bool Lib06::find_state_melt(Lib95 *ptr_lib95eos,
                            double pressure_in,
                            double &temperature_out,
                            double &mden_liq_out,
                            double &mden_ice_out) {
    double temp_now = temperature_out;
    double mden_liq = mden_liq_out;
    double mden_ice = mden_ice_out;

    int i_iter = 0;
    bool found_state = false;
    while (!found_state) {
        i_iter += 1;

        double temp_prev = temp_now;

        bool found_root_vap =
            ptr_lib95eos->find_root_mdensity(temp_prev,
                                             pressure_in,
                                             mden_liq);

        if (!found_root_vap) {
            break;
        }

        mden_ice =
            get_param_mdensity(temp_prev,
                               pressure_in);

        double g_liq =
            ptr_lib95eos->get_param_g(mden_liq,
                                      temp_prev);
        double g_ice =
            get_param_g(temp_prev, pressure_in);

        if (fabs(g_liq - g_ice) <
            0.5 * eps_precision_ * fabs(g_liq + g_ice)) {
            temperature_out = temp_prev;
            mden_liq_out = mden_liq;
            mden_ice_out = mden_ice;

            found_state = true;
        }

        if (i_iter > n_iter_max_) {
            break;
        }

        temp_now = temp_prev +
            0.1 * (g_ice - g_liq) /
            get_param_entropy(temp_prev, pressure_in);
    }

    /*
    fprintf(stderr, "      pressure = %.8e Pa, ", pressure_in);
    fprintf(stderr, "      temperature = %.8e degK, ", temperature_out);
    fprintf(stderr, "      mden_liq = %.8e kg / m^3, ", mden_liq_out);
    fprintf(stderr, "      mden_ice = %.8e kg / m^3\n", mden_ice_out);
    */

    return found_state;
}

void Lib06::make_tab_melt(Lib95 *ptr_lib95eos,
                          int nbin_in,
                          double pressure_max) {
    reset_tab_melt();

    if (pressure_max < pressure_trip_num_) {
        return;
    }

    nbin_melt_ = nbin_in;
    pressure_melt_min_ = pressure_trip_num_;
    pressure_melt_max_ = pressure_max;

    double d_log_press =
        log(pressure_melt_max_ / pressure_melt_min_) /
        static_cast<double>(nbin_melt_);

    have_tab_melt_ = true;
    tab_melt_pressure_ = new double[nbin_melt_ + 1];
    tab_melt_temperature_ = new double[nbin_melt_ + 1];
    tab_melt_mden_liq_ = new double[nbin_melt_ + 1];
    tab_melt_mden_ice_ = new double[nbin_melt_ + 1];
    tab_melt_enthalpy_liq_ = new double[nbin_melt_ + 1];
    tab_melt_enthalpy_ice_ = new double[nbin_melt_ + 1];
    tab_melt_entropy_liq_ = new double[nbin_melt_ + 1];
    tab_melt_entropy_ice_ = new double[nbin_melt_ + 1];

    bool made_tab = true;
    if (flag_verbose_) {
        fprintf(stderr, "    pressure (Pa), ");
        fprintf(stderr, "    temperature (degK), ");
        fprintf(stderr, "    mden_vap (kg / m^3), ");
        fprintf(stderr, "    mden_ice (kg / m^3)\n");
    }
    for (int ip = 0; ip <= nbin_melt_; ip++) {
        tab_melt_pressure_[ip] =
            pressure_melt_min_ *
            exp(d_log_press * static_cast<double>(ip));

        if (ip == 0) {
            tab_melt_temperature_[ip] = temperature_trip_;
            tab_melt_mden_liq_[ip] = 9.99792520e+02;
            tab_melt_mden_ice_[ip] =
                get_param_mdensity(tab_melt_temperature_[ip],
                                   tab_melt_pressure_[ip]);
        } else {
            tab_melt_temperature_[ip] = tab_melt_temperature_[ip - 1];
            tab_melt_mden_liq_[ip] = tab_melt_mden_liq_[ip - 1];
            tab_melt_mden_ice_[ip] = tab_melt_mden_ice_[ip - 1];
        }

        bool found_state =
            find_state_melt(ptr_lib95eos,
                            tab_melt_pressure_[ip],
                            tab_melt_temperature_[ip],
                            tab_melt_mden_liq_[ip],
                            tab_melt_mden_ice_[ip]);

        if (!found_state) {
            made_tab = false;
            break;
        }

        tab_melt_enthalpy_liq_[ip] =
            ptr_lib95eos->get_param_enthalpy(tab_melt_mden_liq_[ip],
                                             tab_melt_temperature_[ip]);
        tab_melt_enthalpy_ice_[ip] =
            get_param_enthalpy(tab_melt_temperature_[ip],
                               tab_melt_pressure_[ip]);

        tab_melt_entropy_liq_[ip] =
            ptr_lib95eos->get_param_entropy(tab_melt_mden_liq_[ip],
                                            tab_melt_temperature_[ip]);
        tab_melt_entropy_ice_[ip] =
            get_param_entropy(tab_melt_temperature_[ip],
                              tab_melt_pressure_[ip]);

        if (flag_verbose_) {
            fprintf(stderr, "      %.8e",
                    tab_melt_pressure_[ip]);
            fprintf(stderr, "      %.8e",
                    tab_melt_temperature_[ip]);
            fprintf(stderr, "      %.8e",
                    tab_melt_mden_liq_[ip]);
            fprintf(stderr, "      %.8e\n",
                    tab_melt_mden_ice_[ip]);
        }
    }

    if (!made_tab) {
        reset_tab_melt();
    }

    set_cspline_melt();

    return;
}

void Lib06::export_tab_melt(char *filename) {
    if (!have_tab_melt_) {
        return;
    }

    FILE *ptr_fout;
    ptr_fout = fopen(filename, "w");
    if (ptr_fout == NULL) {
        return;
    }

    fprintf(ptr_fout, "# table for melting phases in IAPWS06\n");
    fprintf(ptr_fout, "nbin_melt    %d\n", nbin_melt_);
    fprintf(ptr_fout, "pressure_melt_min    %e\n",
            pressure_melt_min_);
    fprintf(ptr_fout, "pressure_melt_max    %e\n",
            pressure_melt_max_);
    fprintf(ptr_fout, "# pressure (Pa)");
    fprintf(ptr_fout, "    temperature (degK)");
    fprintf(ptr_fout, "    mden_liq (kg / m^3)");
    fprintf(ptr_fout, "    mden_ice (kg / m^3)");
    fprintf(ptr_fout, "    enthalpy_liq (J / kg)");
    fprintf(ptr_fout, "    enthalpy_ice (J / kg)");
    fprintf(ptr_fout, "    entropy_liq (J / kg / degK)");
    fprintf(ptr_fout, "    entropy_ice (J / kg / degK)");
    fprintf(ptr_fout, "\n");

    for (int ip = 0; ip <= nbin_melt_; ip++) {
        fprintf(ptr_fout, "    %.8e",
                tab_melt_pressure_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_temperature_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_mden_liq_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_mden_ice_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_enthalpy_liq_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_enthalpy_ice_[ip]);
        fprintf(ptr_fout, "    %.8e",
                tab_melt_entropy_liq_[ip]);
        fprintf(ptr_fout, "    %.8e\n",
                tab_melt_entropy_ice_[ip]);
    }

    fclose(ptr_fout);

    return;
}

void Lib06::import_tab_melt(char *filename) {
    reset_tab_melt();

    FILE *ptr_fin;
    ptr_fin = fopen(filename, "r");
    if (ptr_fin == NULL) {
        return;
    }

    char line_current[1000];

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        if (flag_verbose_) {
            fprintf(stderr, "%s", line_current);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %d", buffer, &nbin_melt_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %d\n", buffer, nbin_melt_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %lf", buffer, &pressure_melt_min_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %f\n", buffer, pressure_melt_min_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        char buffer[100];
        sscanf(line_current, "%s %lf", buffer, &pressure_melt_max_);
        if (flag_verbose_) {
            fprintf(stderr, "%s    %f\n", buffer, pressure_melt_max_);
        }
    } else {
        return;
    }

    if (fgets(line_current, sizeof(line_current), ptr_fin) != NULL) {
        if (flag_verbose_) {
            fprintf(stderr, "%s", line_current);
        }
    } else {
        return;
    }

    have_tab_melt_ = true;
    tab_melt_pressure_ = new double[nbin_melt_ + 1];
    tab_melt_temperature_ = new double[nbin_melt_ + 1];
    tab_melt_mden_liq_ = new double[nbin_melt_ + 1];
    tab_melt_mden_ice_ = new double[nbin_melt_ + 1];
    tab_melt_enthalpy_liq_ = new double[nbin_melt_ + 1];
    tab_melt_enthalpy_ice_ = new double[nbin_melt_ + 1];
    tab_melt_entropy_liq_ = new double[nbin_melt_ + 1];
    tab_melt_entropy_ice_ = new double[nbin_melt_ + 1];

    bool made_tab = true;
    for (int ip = 0; ip <= nbin_melt_; ip++) {
        if (fgets(line_current, sizeof(line_current), ptr_fin) == NULL) {
            made_tab = false;
            break;
        }

        sscanf(line_current, "%lf %lf %lf %lf %lf %lf %lf %lf",
               &tab_melt_pressure_[ip],
               &tab_melt_temperature_[ip],
               &tab_melt_mden_liq_[ip],
               &tab_melt_mden_ice_[ip],
               &tab_melt_enthalpy_liq_[ip],
               &tab_melt_enthalpy_ice_[ip],
               &tab_melt_entropy_liq_[ip],
               &tab_melt_entropy_ice_[ip]);

        if (flag_verbose_) {
            fprintf(stderr, "    %.8e",
                    tab_melt_pressure_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_temperature_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_mden_liq_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_mden_ice_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_enthalpy_liq_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_enthalpy_ice_[ip]);
            fprintf(stderr, "    %.8e",
                    tab_melt_entropy_liq_[ip]);
            fprintf(stderr, "    %.8e\n",
                    tab_melt_entropy_ice_[ip]);
        }
    }

    fclose(ptr_fin);

    if (!made_tab) {
        reset_tab_melt();
    }

    set_cspline_melt();

    return;
}

void Lib06::set_cspline_melt() {
    if (!have_tab_melt_) {
        return;
    }

    csp_melt_temperature_.init(nbin_melt_,
                               tab_melt_pressure_,
                               tab_melt_temperature_);
    csp_melt_mden_liq_.init(nbin_melt_,
                            tab_melt_pressure_,
                            tab_melt_mden_liq_);
    csp_melt_mden_ice_.init(nbin_melt_,
                            tab_melt_pressure_,
                            tab_melt_mden_ice_);
    csp_melt_enthalpy_liq_.init(nbin_melt_,
                                tab_melt_pressure_,
                                tab_melt_enthalpy_liq_);
    csp_melt_enthalpy_ice_.init(nbin_melt_,
                                tab_melt_pressure_,
                                tab_melt_enthalpy_ice_);
    csp_melt_entropy_liq_.init(nbin_melt_,
                               tab_melt_pressure_,
                               tab_melt_entropy_liq_);
    csp_melt_entropy_ice_.init(nbin_melt_,
                               tab_melt_pressure_,
                               tab_melt_entropy_ice_);

    return;
}

} // end namespace IAPWS
