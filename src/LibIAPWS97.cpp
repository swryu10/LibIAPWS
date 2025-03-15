#include<math.h>
#include"LibIAPWS97.h"

namespace IAPWS {

void Lib97::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R7-97(2012)\n");
    fprintf(ptr_fout, "Revised Release on the IAPWS");
    fprintf(ptr_fout, " Industrial Formulation 1997\n");
    fprintf(ptr_fout, "for the Thermodynamic Properties");
    fprintf(ptr_fout, " of Water and Steam\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib97::get_param_g(double temperature_in,
                          double pressure_in,
                          bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double fn_gamma;

    switch(n_reg) {
      case 1 :
        fn_gamma =
            get_param1_gamma(temperature_in,
                             pressure_in);
        break;
      case 2 :
        fn_gamma =
            get_param2_gamma_ide(temperature_in,
                                 pressure_in,
                                 flag_metastable) +
            get_param2_gamma_res(temperature_in,
                                 pressure_in,
                                 flag_metastable);
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        fn_gamma = 0.;
    }

    double g =
        const_R_spec_ * temperature_in * fn_gamma;

    return g;
}

double Lib97::get_param_vol_spec(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double ppi = 0.;
    double fn_dgamma_dppi;

    double vol_spec;

    switch(n_reg) {
      case 1 :
        ppi = pressure_in / pressure_ref1_;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        break;
      case 2 :
        if (flag_metastable) {
            ppi = pressure_in / pressure_ref2mst_;
        } else {
            ppi = pressure_in / pressure_ref2_;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        vol_spec = 0.;
    }

    return vol_spec;
}

double Lib97::get_param_mdensity(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    double mdensity =
        1. / get_param_vol_spec(temperature_in,
                                pressure_in,
                                flag_metastable);

    return mdensity;
}

double Lib97::get_param_erg_int(double temperature_in,
                                double pressure_in,
                                bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double ppi = 0.;
    double fn_dgamma_dtau;
    double fn_dgamma_dppi;

    double erg_int;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        ppi = pressure_in / pressure_ref1_;
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * fn_dgamma_dtau -
             ppi * fn_dgamma_dppi);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
            ppi = pressure_in / pressure_ref2mst_;
        } else {
            tau = temperature_ref2_ / temperature_in;
            ppi = pressure_in / pressure_ref2_;
        }
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * fn_dgamma_dtau -
             ppi * fn_dgamma_dppi);
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        erg_int = 0.;
    }

    return erg_int;
}

double Lib97::get_param_entropy(double temperature_in,
                                double pressure_in,
                                bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_gamma;
    double fn_dgamma_dtau;

    double entropy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_gamma =
            get_param1_gamma(temperature_in,
                             pressure_in);
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        entropy =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_gamma =
            get_param2_gamma_ide(temperature_in,
                                 pressure_in,
                                 flag_metastable) +
            get_param2_gamma_res(temperature_in,
                                 pressure_in,
                                 flag_metastable);
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        entropy =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        entropy = 0.;
    }

    return entropy;
}

double Lib97::get_param_enthalpy(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dtau;

    double enthalpy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        enthalpy = 0.;
    }

    return enthalpy;
}

double Lib97::get_param_heat_c_p(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_d2gamma_dtau_dtau;

    double heat_c_p;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        heat_c_p =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        heat_c_p =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau;
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        heat_c_p = 0.;
    }

    return heat_c_p;
}

double Lib97::get_param_heat_c_v(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi;
    double fn_d2gamma_dppi_dtau;
    double fn_d2gamma_dtau_dtau;

    double heat_c_v;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param1_d2gamma_dppi_dppi(temperature_in,
                                         pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param1_d2gamma_dppi_dtau(temperature_in,
                                         pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        heat_c_v =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau +
            const_R_spec_ *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
            fn_d2gamma_dppi_dppi;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_d2gamma_dppi_dppi =
            get_param2_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dppi_dtau =
            get_param2_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        heat_c_v =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau +
            const_R_spec_ *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
            fn_d2gamma_dppi_dppi;
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        heat_c_v = 0.;
    }

    return heat_c_v;
}

double Lib97::get_param_speed_sound(double temperature_in,
                                    double pressure_in,
                                    bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi;
    double fn_d2gamma_dppi_dtau;
    double fn_d2gamma_dtau_dtau;

    double v2_s;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param1_d2gamma_dppi_dppi(temperature_in,
                                         pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param1_d2gamma_dppi_dtau(temperature_in,
                                         pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        v2_s = const_R_spec_ * temperature_in *
            fn_dgamma_dppi * fn_dgamma_dppi /
            ((fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
             (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
             (tau * tau * fn_d2gamma_dtau_dtau) -
             fn_d2gamma_dppi_dppi);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_d2gamma_dppi_dppi =
            get_param2_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dppi_dtau =
            get_param2_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        v2_s = const_R_spec_ * temperature_in *
            fn_dgamma_dppi * fn_dgamma_dppi /
            ((fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
             (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
             (tau * tau * fn_d2gamma_dtau_dtau) -
             fn_d2gamma_dppi_dppi);
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        v2_s = 0.;
    }

    return sqrt(v2_s);
}

int Lib97::get_region(double temperature_in,
                      double pressure_in,
                      bool flag_metastable) {
    int n_reg = 0;

    if (flag_metastable) {
        n_reg = 2;
    } else if (temperature_in >= 273.15 &&
               temperature_in <= 623.15) {
        double press_sat =
            get_param4_sat_pressure(temperature_in);

        if (pressure_in <= press_sat) {
            n_reg = 2;
        } else if (pressure_in > press_sat &&
                   pressure_in < 1.0e+8) {
            n_reg = 1;
        }
    } else if (temperature_in > 623.15 &&
               temperature_in <= 863.15) {
        double press_B23 =
            get_paramB23_pressure(temperature_in);

        if (pressure_in <= press_B23) {
            n_reg = 2;
        } else if (pressure_in > press_B23 &&
                   pressure_in < 1.0e+8) {
            n_reg = 3;
        }
    } else if (temperature_in > 863.15 &&
               temperature_in <= 1073.15) {
        if (pressure_in > 0. &&
            pressure_in < 1.0e+8) {
            n_reg = 2;
        }
    }

    return n_reg;
}

double Lib97::get_paramB23_pressure(double temperature_in) {
    double tau = temperature_in / temperature_refB23_;

    double ppi =
        coeffB23_n_[1] +
        coeffB23_n_[2] * tau +
        coeffB23_n_[3] * tau * tau;

    return pressure_refB23_ * ppi;
}

double Lib97::get_paramB23_temperature(double pressure_in) {
    double ppi = pressure_in / pressure_refB23_;

    double tau =
        coeffB23_n_[4] +
        sqrt((ppi - coeffB23_n_[5]) / coeffB23_n_[3]);

    return temperature_refB23_ * tau;
}

double Lib97::get_param1_gamma(double temperature_in,
                               double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double gamma = 0.;
    for (int i = 1; i <= 34; i++) {
        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        gamma += coeff1_n_[i] * dppi_pow_I * dtau_pow_J;
    }

    return gamma;
}

double Lib97::get_param1_dgamma_dppi(double temperature_in,
                                     double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double dgamma_dppi = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 1);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        dgamma_dppi -=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i]);
    }

    return dgamma_dppi;
}

double Lib97::get_param1_dgamma_dtau(double temperature_in,
                                     double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double dgamma_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_J_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 1);

        dgamma_dtau +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param1_d2gamma_dppi_dppi(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dppi_dppi = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0 ||
            coeff1_I_[i] == 1) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 2);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        d2gamma_dppi_dppi +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i] * (coeff1_I_[i] - 1));
    }

    return d2gamma_dppi_dppi;
}

double Lib97::get_param1_d2gamma_dppi_dtau(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dppi_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0 ||
            coeff1_J_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 1);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 1);

        d2gamma_dppi_dtau -=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i] * coeff1_J_[i]);
    }

    return d2gamma_dppi_dtau;
}

double Lib97::get_param1_d2gamma_dtau_dtau(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dtau_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_J_[i] == 0 ||
            coeff1_J_[i] == 1) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 2);

        d2gamma_dtau_dtau +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_J_[i] * (coeff1_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param1_temperature_ph(double pressure_in,
                                        double enthalpy_in) {
    double ppi = pressure_in / pressure_ref1Tph_;
    double eta = enthalpy_in / enthalpy_ref1Tph_;

    double fn_theta = 0.;
    for (int i = 1; i <= 20; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff1Tph_I_[i]);
        double eta_pow_J =
            get_pow_int(eta + 1., coeff1Tph_J_[i]);

        fn_theta +=
            coeff1Tph_n_[i] * ppi_pow_I * eta_pow_J;
    }

    return temperature_ref1Tph_ * fn_theta;
}

double Lib97::get_param1_temperature_ps(double pressure_in,
                                        double entropy_in) {
    double ppi = pressure_in / pressure_ref1Tps_;
    double sig = entropy_in / entropy_ref1Tps_;

    double fn_theta = 0.;
    for (int i = 1; i <= 20; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff1Tps_I_[i]);
        double sig_pow_J =
            get_pow_int(sig + 2., coeff1Tps_J_[i]);

        fn_theta +=
            coeff1Tps_n_[i] * ppi_pow_I * sig_pow_J;
    }

    return temperature_ref1Tps_ * fn_theta;
}

double Lib97::get_param2_gamma_ide(double temperature_in,
                                   double pressure_in,
                                   bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double gamma = log(ppi);

    for (int i = 1; i <= 9; i++) {
        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i]);
        gamma +=
            coeff_n * dtau_pow_J;
    }

    return gamma;
}

double Lib97::get_param2_dgamma_ide_dppi(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dgamma_dppi = 1. / ppi;

    return dgamma_dppi;
}

double Lib97::get_param2_dgamma_ide_dtau(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    double dgamma_dtau = 0.;

    for (int i = 1; i <= 9; i++) {
        if (coeff2_ide_J_[i] == 0) {
            continue;
        }

        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i] - 1);
        dgamma_dtau +=
            coeff_n * dtau_pow_J *
            static_cast<double>(coeff2_ide_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param2_d2gamma_ide_dppi_dppi(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double d2gamma_dppi_dppi = -1. / (ppi * ppi);

    return d2gamma_dppi_dppi;
}

double Lib97::get_param2_d2gamma_ide_dppi_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    return 0.;
}

double Lib97::get_param2_d2gamma_ide_dtau_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    double d2gamma_dtau_dtau = 0.;

    for (int i = 1; i <= 9; i++) {
        if (coeff2_ide_J_[i] == 0 ||
            coeff2_ide_J_[i] == 1) {
            continue;
        }

        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i] - 2);
        d2gamma_dtau_dtau +=
            coeff_n * dtau_pow_J *
            static_cast<double>(coeff2_ide_J_[i] *
                                (coeff2_ide_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param2_gamma_res(double temperature_in,
                                   double pressure_in,
                                   bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double gamma = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            gamma +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J;
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            gamma +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J;
        }
    }

    return gamma;
}

double Lib97::get_param2_dgamma_res_dppi(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double dgamma_dppi = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            dgamma_dppi +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            dgamma_dppi +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i]);
        }
    }

    return dgamma_dppi;
}

double Lib97::get_param2_dgamma_res_dtau(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double dgamma_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 1);

            dgamma_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_J_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 1);

            dgamma_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_J_[i]);
        }
    }

    return dgamma_dtau;
}

double Lib97::get_param2_d2gamma_res_dppi_dppi(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dppi_dppi = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0 ||
                coeff2mst_res_I_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 2);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            d2gamma_dppi_dppi +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i] *
                                    (coeff2mst_res_I_[i] - 1));
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0 ||
                coeff2_res_I_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 2);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            d2gamma_dppi_dppi +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i] *
                                    (coeff2_res_I_[i] - 1));
        }
    }

    return d2gamma_dppi_dppi;
}

double Lib97::get_param2_d2gamma_res_dppi_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dppi_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0 ||
                coeff2mst_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 1);

            d2gamma_dppi_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i] *
                                    coeff2mst_res_J_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0 ||
                coeff2_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 1);

            d2gamma_dppi_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i] *
                                    coeff2_res_J_[i]);
        }
    }

    return d2gamma_dppi_dtau;
}

double Lib97::get_param2_d2gamma_res_dtau_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dtau_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_J_[i] == 0 ||
                coeff2mst_res_J_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 2);

            d2gamma_dtau_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_J_[i] *
                                    (coeff2mst_res_J_[i] - 1));
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_J_[i] == 0 ||
                coeff2_res_J_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 2);

            d2gamma_dtau_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_J_[i] *
                                    (coeff2_res_J_[i] - 1));
        }
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param4_sat_pressure(double temperature_in) {
    double fn_tau = temperature_in / temperature_ref4_;
    double fn_theta =
        fn_tau + coeff4_n_[9] / (fn_tau - coeff4_n_[10]);

    double fn_A =
        fn_theta * fn_theta +
        fn_theta * coeff4_n_[1] + coeff4_n_[2];
    double fn_B =
        fn_theta * fn_theta * coeff4_n_[3] +
        fn_theta * coeff4_n_[4] + coeff4_n_[5];
    double fn_C =
        fn_theta * fn_theta * coeff4_n_[6] +
        fn_theta * coeff4_n_[7] + coeff4_n_[8];

    double fn_beta =
        2. * fn_C / (sqrt(fn_B * fn_B - 4. * fn_A * fn_C) - fn_B);

    double press_sat =
        pressure_ref4_ * get_pow_int(fn_beta, 4);

    return press_sat;
}

double Lib97::get_param4_sat_temperature(double pressure_in) {
    double fn_beta = pow(pressure_in / pressure_ref4_, 0.25);

    double fn_E =
        fn_beta * fn_beta +
        fn_beta * coeff4_n_[3] + coeff4_n_[6];
    double fn_F =
        fn_beta * fn_beta * coeff4_n_[1] +
        fn_beta * coeff4_n_[4] + coeff4_n_[7];
    double fn_G =
        fn_beta * fn_beta * coeff4_n_[2] +
        fn_beta * coeff4_n_[5] + coeff4_n_[8];

    double fn_D =
        -2. * fn_G / (sqrt(fn_F * fn_F - 4. * fn_E * fn_G) + fn_F);
    double fn_tau =
        (coeff4_n_[10] + fn_D -
         sqrt((coeff4_n_[10] + fn_D) * (coeff4_n_[10] + fn_D) -
              4. * (coeff4_n_[9] + coeff4_n_[10] * fn_D))) / 2.;

    double temp_sat = temperature_ref4_ * fn_tau;

    return temp_sat;
}

void Lib97::set_coefficients() {
    temperature_crit_ = 647.096;
    mdensity_crit_ = 322.;
    pressure_crit_ = 22.064 * 1.0e+6;
    const_R_spec_ = 461.526;

    temperature_refB23_ = 1.;
    pressure_refB23_ = 1.0e+6;

    coeffB23_n_[0] = 0.;
    // coefficient boundary between regions 2 and 3 n_i
    coeffB23_n_[1] = 0.34805185628969 * 1.0e+3;
    coeffB23_n_[2] = -0.11671859879975 * 1.0e+1;
    coeffB23_n_[3] = 0.10192970039326 * 1.0e-2;
    coeffB23_n_[4] = 0.57254459862746 * 1.0e+3;
    coeffB23_n_[5] = 0.13918839778870 * 1.0e+2;

    temperature_ref1_ = 1386.;
    pressure_ref1_ = 16.53 * 1.0e+6;
    tau_ref1_ = 1.222;
    ppi_ref1_ = 7.1;

    coeff1_I_[0] = 0;
    // coefficient region 1 I_i
    coeff1_I_[1] = 0;
    coeff1_I_[2] = 0;
    coeff1_I_[3] = 0;
    coeff1_I_[4] = 0;
    coeff1_I_[5] = 0;
    coeff1_I_[6] = 0;
    coeff1_I_[7] = 0;
    coeff1_I_[8] = 0;
    coeff1_I_[9] = 1;
    coeff1_I_[10] = 1;
    coeff1_I_[11] = 1;
    coeff1_I_[12] = 1;
    coeff1_I_[13] = 1;
    coeff1_I_[14] = 1;
    coeff1_I_[15] = 2;
    coeff1_I_[16] = 2;
    coeff1_I_[17] = 2;
    coeff1_I_[18] = 2;
    coeff1_I_[19] = 2;
    coeff1_I_[20] = 3;
    coeff1_I_[21] = 3;
    coeff1_I_[22] = 3;
    coeff1_I_[23] = 4;
    coeff1_I_[24] = 4;
    coeff1_I_[25] = 4;
    coeff1_I_[26] = 5;
    coeff1_I_[27] = 8;
    coeff1_I_[28] = 8;
    coeff1_I_[29] = 21;
    coeff1_I_[30] = 23;
    coeff1_I_[31] = 29;
    coeff1_I_[32] = 30;
    coeff1_I_[33] = 31;
    coeff1_I_[34] = 32;

    coeff1_J_[0] = 0;
    // coefficient region 1 J_i
    coeff1_J_[1] = -2;
    coeff1_J_[2] = -1;
    coeff1_J_[3] = 0;
    coeff1_J_[4] = 1;
    coeff1_J_[5] = 2;
    coeff1_J_[6] = 3;
    coeff1_J_[7] = 4;
    coeff1_J_[8] = 5;
    coeff1_J_[9] = -9;
    coeff1_J_[10] = -7;
    coeff1_J_[11] = -1;
    coeff1_J_[12] = 0;
    coeff1_J_[13] = 1;
    coeff1_J_[14] = 3;
    coeff1_J_[15] = -3;
    coeff1_J_[16] = 0;
    coeff1_J_[17] = 1;
    coeff1_J_[18] = 3;
    coeff1_J_[19] = 17;
    coeff1_J_[20] = -4;
    coeff1_J_[21] = 0;
    coeff1_J_[22] = 6;
    coeff1_J_[23] = -5;
    coeff1_J_[24] = -2;
    coeff1_J_[25] = 10;
    coeff1_J_[26] = -8;
    coeff1_J_[27] = -11;
    coeff1_J_[28] = -6;
    coeff1_J_[29] = -29;
    coeff1_J_[30] = -31;
    coeff1_J_[31] = -38;
    coeff1_J_[32] = -39;
    coeff1_J_[33] = -40;
    coeff1_J_[34] = -41;

    coeff1_n_[0] = 0.;
    // coefficient region 1 n_i
    coeff1_n_[1] = 0.14632971213167;
    coeff1_n_[2] = -0.84548187169114;
    coeff1_n_[3] = -0.37563603672040 * 1.0e+1;
    coeff1_n_[4] = 0.33855169168385 * 1.0e+1;
    coeff1_n_[5] = -0.95791963387872;
    coeff1_n_[6] = 0.15772038513228;
    coeff1_n_[7] = -0.16616417199501 * 1.0e-1;
    coeff1_n_[8] = 0.81214629983568 * 1.0e-3;
    coeff1_n_[9] = 0.28319080123804 * 1.0e-3;
    coeff1_n_[10] = -0.60706301565874 * 1.0e-3;
    coeff1_n_[11] = -0.18990068218419 * 1.0e-1;
    coeff1_n_[12] = -0.32529748770505 * 1.0e-1;
    coeff1_n_[13] = -0.21841717175414 * 1.0e-1;
    coeff1_n_[14] = -0.52838357969930 * 1.0e-4;
    coeff1_n_[15] = -0.47184321073267 * 1.0e-3;
    coeff1_n_[16] = -0.30001780793026 * 1.0e-3;
    coeff1_n_[17] = 0.47661393906987 * 1.0e-4;
    coeff1_n_[18] = -0.44141845330846 * 1.0e-5;
    coeff1_n_[19] = -0.72694996297594 * 1.0e-15;
    coeff1_n_[20] = -0.31679644845054 * 1.0e-4;
    coeff1_n_[21] = -0.28270797985312 * 1.0e-5;
    coeff1_n_[22] = -0.85205128120103 * 1.0e-9;
    coeff1_n_[23] = -0.22425281908000 * 1.0e-5;
    coeff1_n_[24] = -0.65171222895601 * 1.0e-6;
    coeff1_n_[25] = -0.14341729937924 * 1.0e-12;
    coeff1_n_[26] = -0.40516996860117 * 1.0e-6;
    coeff1_n_[27] = -0.12734301741641 * 1.0e-8;
    coeff1_n_[28] = -0.17424871230634 * 1.0e-9;
    coeff1_n_[29] = -0.68762131295531 * 1.0e-18;
    coeff1_n_[30] = 0.14478307828521 * 1.0e-19;
    coeff1_n_[31] = 0.26335781662795 * 1.0e-22;
    coeff1_n_[32] = -0.11947622640071 * 1.0e-22;
    coeff1_n_[33] = 0.18228094581404 * 1.0e-23;
    coeff1_n_[34] = -0.93537087292458 * 1.0e-25;

    temperature_ref1Tph_ = 1.;
    pressure_ref1Tph_ = 1.0e+6;
    enthalpy_ref1Tph_ = 2500. * 1.0e+3;

    coeff1Tph_I_[0] = 0;
    // backward T(p,h) coefficient region 1 I_i
    coeff1Tph_I_[1] = 0;
    coeff1Tph_I_[2] = 0;
    coeff1Tph_I_[3] = 0;
    coeff1Tph_I_[4] = 0;
    coeff1Tph_I_[5] = 0;
    coeff1Tph_I_[6] = 0;
    coeff1Tph_I_[7] = 1;
    coeff1Tph_I_[8] = 1;
    coeff1Tph_I_[9] = 1;
    coeff1Tph_I_[10] = 1;
    coeff1Tph_I_[11] = 1;
    coeff1Tph_I_[12] = 1;
    coeff1Tph_I_[13] = 1;
    coeff1Tph_I_[14] = 2;
    coeff1Tph_I_[15] = 2;
    coeff1Tph_I_[16] = 3;
    coeff1Tph_I_[17] = 3;
    coeff1Tph_I_[18] = 4;
    coeff1Tph_I_[19] = 5;
    coeff1Tph_I_[20] = 6;

    coeff1Tph_J_[0] = 0;
    // backward T(p,h) coefficient region 1 J_i
    coeff1Tph_J_[1] = 0;
    coeff1Tph_J_[2] = 1;
    coeff1Tph_J_[3] = 2;
    coeff1Tph_J_[4] = 6;
    coeff1Tph_J_[5] = 22;
    coeff1Tph_J_[6] = 32;
    coeff1Tph_J_[7] = 0;
    coeff1Tph_J_[8] = 1;
    coeff1Tph_J_[9] = 2;
    coeff1Tph_J_[10] = 3;
    coeff1Tph_J_[11] = 4;
    coeff1Tph_J_[12] = 10;
    coeff1Tph_J_[13] = 32;
    coeff1Tph_J_[14] = 10;
    coeff1Tph_J_[15] = 32;
    coeff1Tph_J_[16] = 10;
    coeff1Tph_J_[17] = 32;
    coeff1Tph_J_[18] = 32;
    coeff1Tph_J_[19] = 32;
    coeff1Tph_J_[20] = 32;

    coeff1Tph_n_[0] = 0.;
    // backward T(p,h) coefficient region 1 n_i
    coeff1Tph_n_[1] = -0.23872489924521 * 1.0e+3;
    coeff1Tph_n_[2] = 0.40421188637945 * 1.0e+3;
    coeff1Tph_n_[3] = 0.11349746881718 * 1.0e+3;
    coeff1Tph_n_[4] = -0.58457616048039 * 1.0e+1;
    coeff1Tph_n_[5] = -0.15285482413140 * 1.0e-3;
    coeff1Tph_n_[6] = -0.10866707695377 * 1.0e-5;
    coeff1Tph_n_[7] = -0.13391744872602 * 1.0e+2;
    coeff1Tph_n_[8] = 0.43211039183559 * 1.0e+2;
    coeff1Tph_n_[9] = -0.54010067170506 * 1.0e+2;
    coeff1Tph_n_[10] = 0.30535892203916 * 1.0e+2;
    coeff1Tph_n_[11] = -0.65964749423638 * 1.0e+1;
    coeff1Tph_n_[12] = 0.93965400878363 * 1.0e-2;
    coeff1Tph_n_[13] = 0.11573647505340 * 1.0e-6;
    coeff1Tph_n_[14] = -0.25858641282073 * 1.0e-4;
    coeff1Tph_n_[15] = -0.40644363084799 * 1.0e-8;
    coeff1Tph_n_[16] = 0.66456186191635 * 1.0e-7;
    coeff1Tph_n_[17] = 0.80670734103027 * 1.0e-10;
    coeff1Tph_n_[18] = -0.93477771213947 * 1.0e-12;
    coeff1Tph_n_[19] = 0.58265442020601 * 1.0e-14;
    coeff1Tph_n_[20] = -0.15020185953503 * 1.0e-16;

    temperature_ref1Tps_ = 1.;
    pressure_ref1Tps_ = 1.0e+6;
    entropy_ref1Tps_ = 1.0e+3;

    coeff1Tps_I_[0] = 0;
    // backward T(p,s) coefficient region 1 I_i
    coeff1Tps_I_[1] = 0;
    coeff1Tps_I_[2] = 0;
    coeff1Tps_I_[3] = 0;
    coeff1Tps_I_[4] = 0;
    coeff1Tps_I_[5] = 0;
    coeff1Tps_I_[6] = 0;
    coeff1Tps_I_[7] = 1;
    coeff1Tps_I_[8] = 1;
    coeff1Tps_I_[9] = 1;
    coeff1Tps_I_[10] = 1;
    coeff1Tps_I_[11] = 1;
    coeff1Tps_I_[12] = 1;
    coeff1Tps_I_[13] = 2;
    coeff1Tps_I_[14] = 2;
    coeff1Tps_I_[15] = 2;
    coeff1Tps_I_[16] = 2;
    coeff1Tps_I_[17] = 2;
    coeff1Tps_I_[18] = 3;
    coeff1Tps_I_[19] = 3;
    coeff1Tps_I_[20] = 4;

    coeff1Tps_J_[0] = 0;
    // backward T(p,s) coefficient region 1 J_i
    coeff1Tps_J_[1] = 0;
    coeff1Tps_J_[2] = 1;
    coeff1Tps_J_[3] = 2;
    coeff1Tps_J_[4] = 3;
    coeff1Tps_J_[5] = 11;
    coeff1Tps_J_[6] = 31;
    coeff1Tps_J_[7] = 0;
    coeff1Tps_J_[8] = 1;
    coeff1Tps_J_[9] = 2;
    coeff1Tps_J_[10] = 3;
    coeff1Tps_J_[11] = 12;
    coeff1Tps_J_[12] = 31;
    coeff1Tps_J_[13] = 0;
    coeff1Tps_J_[14] = 1;
    coeff1Tps_J_[15] = 2;
    coeff1Tps_J_[16] = 9;
    coeff1Tps_J_[17] = 31;
    coeff1Tps_J_[18] = 10;
    coeff1Tps_J_[19] = 32;
    coeff1Tps_J_[20] = 32;

    coeff1Tps_n_[0] = 0.;
    // backward T(p,s) coefficient region 1 n_i
    coeff1Tps_n_[1] = 0.17478268058307 * 1.0e+3;
    coeff1Tps_n_[2] = 0.34806930892873 * 1.0e+2;
    coeff1Tps_n_[3] = 0.65292584978455 * 1.0e+1;
    coeff1Tps_n_[4] = 0.33039981775489;
    coeff1Tps_n_[5] = -0.19281382923196 * 1.0e-6;
    coeff1Tps_n_[6] = -0.24909197244573 * 1.0e-22;
    coeff1Tps_n_[7] = -0.26107636489332;
    coeff1Tps_n_[8] = 0.22592965981586;
    coeff1Tps_n_[9] = -0.64256463395226 * 1.0e-1;
    coeff1Tps_n_[10] = 0.78876289270526 * 1.0e-2;
    coeff1Tps_n_[11] = 0.35672110607366 * 1.0e-9;
    coeff1Tps_n_[12] = 0.17332496994895 * 1.0e-23;
    coeff1Tps_n_[13] = 0.56608900654837 * 1.0e-3;
    coeff1Tps_n_[14] = -0.32635483139717 * 1.0e-3;
    coeff1Tps_n_[15] = 0.44778286690632 * 1.0e-4;
    coeff1Tps_n_[16] = -0.51322156908507 * 1.0e-9;
    coeff1Tps_n_[17] = -0.42522657042207 * 1.0e-25;
    coeff1Tps_n_[18] = 0.26400441360689 * 1.0e-12;
    coeff1Tps_n_[19] = 0.78124600459723 * 1.0e-28;
    coeff1Tps_n_[20] = -0.30732199903668 * 1.0e-30;

    temperature_ref2_ = 540.;
    pressure_ref2_ = 1.0e+6;
    tau_ref2_res_ = 0.5;

    coeff2_ide_J_[0] = 0;
    // coefficient region 2 ideal J_i
    coeff2_ide_J_[1] = 0;
    coeff2_ide_J_[2] = 1;
    coeff2_ide_J_[3] = -5;
    coeff2_ide_J_[4] = -4;
    coeff2_ide_J_[5] = -3;
    coeff2_ide_J_[6] = -2;
    coeff2_ide_J_[7] = -1;
    coeff2_ide_J_[8] = 2;
    coeff2_ide_J_[9] = 3;

    coeff2_ide_n_[0] = 0.;
    // coefficient region 2 ideal n_i
    coeff2_ide_n_[1] = -0.96927686500217 * 1.0e+1;
    coeff2_ide_n_[2] = 0.10086655968018 * 1.0e+2;
    coeff2_ide_n_[3] = -0.56087911283020 * 1.0e-2;
    coeff2_ide_n_[4] = 0.71452738081455 * 1.0e-1;
    coeff2_ide_n_[5] = -0.40710498223928;
    coeff2_ide_n_[6] = 0.14240819171444 * 1.0e+1;
    coeff2_ide_n_[7] = -0.43839511319450 * 1.0e+1;
    coeff2_ide_n_[8] = -0.28408632460772;
    coeff2_ide_n_[9] = 0.21268463753307 * 1.0e-1;

    // coefficient region 2 ideal n_i metastable case
    coeff2mst_ide_n1_ = -0.96937268393049 * 1.0e+1;
    coeff2mst_ide_n2_ = 0.10087275970006 * 1.0e+2;

    coeff2_res_I_[0] = 0;
    // coefficient region 2 residual I_i
    coeff2_res_I_[1] = 1;
    coeff2_res_I_[2] = 1;
    coeff2_res_I_[3] = 1;
    coeff2_res_I_[4] = 1;
    coeff2_res_I_[5] = 1;
    coeff2_res_I_[6] = 2;
    coeff2_res_I_[7] = 2;
    coeff2_res_I_[8] = 2;
    coeff2_res_I_[9] = 2;
    coeff2_res_I_[10] = 2;
    coeff2_res_I_[11] = 3;
    coeff2_res_I_[12] = 3;
    coeff2_res_I_[13] = 3;
    coeff2_res_I_[14] = 3;
    coeff2_res_I_[15] = 3;
    coeff2_res_I_[16] = 4;
    coeff2_res_I_[17] = 4;
    coeff2_res_I_[18] = 4;
    coeff2_res_I_[19] = 5;
    coeff2_res_I_[20] = 6;
    coeff2_res_I_[21] = 6;
    coeff2_res_I_[22] = 6;
    coeff2_res_I_[23] = 7;
    coeff2_res_I_[24] = 7;
    coeff2_res_I_[25] = 7;
    coeff2_res_I_[26] = 8;
    coeff2_res_I_[27] = 8;
    coeff2_res_I_[28] = 9;
    coeff2_res_I_[29] = 10;
    coeff2_res_I_[30] = 10;
    coeff2_res_I_[31] = 10;
    coeff2_res_I_[32] = 16;
    coeff2_res_I_[33] = 16;
    coeff2_res_I_[34] = 18;
    coeff2_res_I_[35] = 20;
    coeff2_res_I_[36] = 20;
    coeff2_res_I_[37] = 20;
    coeff2_res_I_[38] = 21;
    coeff2_res_I_[39] = 22;
    coeff2_res_I_[40] = 23;
    coeff2_res_I_[41] = 24;
    coeff2_res_I_[42] = 24;
    coeff2_res_I_[43] = 24;

    coeff2_res_J_[0] = 0;
    // coefficient region 2 residual J_i
    coeff2_res_J_[1] = 0;
    coeff2_res_J_[2] = 1;
    coeff2_res_J_[3] = 2;
    coeff2_res_J_[4] = 3;
    coeff2_res_J_[5] = 6;
    coeff2_res_J_[6] = 1;
    coeff2_res_J_[7] = 2;
    coeff2_res_J_[8] = 4;
    coeff2_res_J_[9] = 7;
    coeff2_res_J_[10] = 36;
    coeff2_res_J_[11] = 0;
    coeff2_res_J_[12] = 1;
    coeff2_res_J_[13] = 3;
    coeff2_res_J_[14] = 6;
    coeff2_res_J_[15] = 35;
    coeff2_res_J_[16] = 1;
    coeff2_res_J_[17] = 2;
    coeff2_res_J_[18] = 3;
    coeff2_res_J_[19] = 7;
    coeff2_res_J_[20] = 3;
    coeff2_res_J_[21] = 16;
    coeff2_res_J_[22] = 35;
    coeff2_res_J_[23] = 0;
    coeff2_res_J_[24] = 11;
    coeff2_res_J_[25] = 25;
    coeff2_res_J_[26] = 8;
    coeff2_res_J_[27] = 36;
    coeff2_res_J_[28] = 13;
    coeff2_res_J_[29] = 4;
    coeff2_res_J_[30] = 10;
    coeff2_res_J_[31] = 14;
    coeff2_res_J_[32] = 29;
    coeff2_res_J_[33] = 50;
    coeff2_res_J_[34] = 57;
    coeff2_res_J_[35] = 20;
    coeff2_res_J_[36] = 35;
    coeff2_res_J_[37] = 48;
    coeff2_res_J_[38] = 21;
    coeff2_res_J_[39] = 53;
    coeff2_res_J_[40] = 39;
    coeff2_res_J_[41] = 26;
    coeff2_res_J_[42] = 40;
    coeff2_res_J_[43] = 58;

    coeff2_res_n_[0] = 0.;
    // coefficient region 2 residual n_i
    coeff2_res_n_[1] = -0.17731742473213 * 1.0e-2;
    coeff2_res_n_[2] = -0.17834862292358 * 1.0e-1;
    coeff2_res_n_[3] = -0.45996013696365 * 1.0e-1;
    coeff2_res_n_[4] = -0.57581259083432 * 1.0e-1;
    coeff2_res_n_[5] = -0.50325278727930 * 1.0e-1;
    coeff2_res_n_[6] = -0.33032641670203 * 1.0e-4;
    coeff2_res_n_[7] = -0.18948987516315 * 1.0e-3;
    coeff2_res_n_[8] = -0.39392777243355 * 1.0e-2;
    coeff2_res_n_[9] = -0.43797295650573 * 1.0e-1;
    coeff2_res_n_[10] = -0.26674547914087 * 1.0e-4;
    coeff2_res_n_[11] = 0.20481737692309 * 1.0e-7;
    coeff2_res_n_[12] = 0.43870667284435 * 1.0e-6;
    coeff2_res_n_[13] = -0.32277677238570 * 1.0e-4;
    coeff2_res_n_[14] = -0.15033924542148 * 1.0e-2;
    coeff2_res_n_[15] = -0.40668253562649 * 1.0e-1;
    coeff2_res_n_[16] = -0.78847309559367 * 1.0e-9;
    coeff2_res_n_[17] = 0.12790717852285 * 1.0e-7;
    coeff2_res_n_[18] = 0.48225372718507 * 1.0e-6;
    coeff2_res_n_[19] = 0.22922076337661 * 1.0e-5;
    coeff2_res_n_[20] = -0.16714766451061 * 1.0e-10;
    coeff2_res_n_[21] = -0.21171472321355 * 1.0e-2;
    coeff2_res_n_[22] = -0.23895741934104 * 1.0e+2;
    coeff2_res_n_[23] = -0.59059564324270 * 1.0e-17;
    coeff2_res_n_[24] = -0.12621808899101 * 1.0e-5;
    coeff2_res_n_[25] = -0.38946842435739 * 1.0e-1;
    coeff2_res_n_[26] = 0.11256211360459 * 1.0e-10;
    coeff2_res_n_[27] = -0.82311340897998 * 1.0e+1;
    coeff2_res_n_[28] = 0.19809712802088 * 1.0e-7;
    coeff2_res_n_[29] = 0.10406965210174 * 1.0e-18;
    coeff2_res_n_[30] = -0.10234747095929 * 1.0e-12;
    coeff2_res_n_[31] = -0.10018179379511 * 1.0e-8;
    coeff2_res_n_[32] = -0.80882908646985 * 1.0e-10;
    coeff2_res_n_[33] = 0.10693031879409;
    coeff2_res_n_[34] = -0.33662250574171;
    coeff2_res_n_[35] = 0.89185845355421 * 1.0e-24;
    coeff2_res_n_[36] = 0.30629316876232 * 1.0e-12;
    coeff2_res_n_[37] = -0.42002467698208 * 1.0e-5;
    coeff2_res_n_[38] = -0.59056029685639 * 1.0e-25;
    coeff2_res_n_[39] = 0.37826947613457 * 1.0e-5;
    coeff2_res_n_[40] = -0.12768608934681 * 1.0e-14;
    coeff2_res_n_[41] = 0.73087610595061 * 1.0e-28;
    coeff2_res_n_[42] = 0.55414715350778 * 1.0e-16;
    coeff2_res_n_[43] = -0.94369707241210 * 1.0e-6;

    temperature_ref2mst_ = 540.;
    pressure_ref2mst_ = 1.0e+6;
    tau_ref2mst_res_ = 0.5;

    coeff2mst_res_I_[0] = 0;
    // coefficient region 2 metastable residual I_i
    coeff2mst_res_I_[1] = 1;
    coeff2mst_res_I_[2] = 1;
    coeff2mst_res_I_[3] = 1;
    coeff2mst_res_I_[4] = 1;
    coeff2mst_res_I_[5] = 2;
    coeff2mst_res_I_[6] = 2;
    coeff2mst_res_I_[7] = 2;
    coeff2mst_res_I_[8] = 3;
    coeff2mst_res_I_[9] = 3;
    coeff2mst_res_I_[10] = 4;
    coeff2mst_res_I_[11] = 4;
    coeff2mst_res_I_[12] = 5;
    coeff2mst_res_I_[13] = 5;

    coeff2mst_res_J_[0] = 0;
    // coefficient region 2 metastable residual J_i
    coeff2mst_res_J_[1] = 0;
    coeff2mst_res_J_[2] = 2;
    coeff2mst_res_J_[3] = 5;
    coeff2mst_res_J_[4] = 11;
    coeff2mst_res_J_[5] = 1;
    coeff2mst_res_J_[6] = 7;
    coeff2mst_res_J_[7] = 16;
    coeff2mst_res_J_[8] = 4;
    coeff2mst_res_J_[9] = 16;
    coeff2mst_res_J_[10] = 7;
    coeff2mst_res_J_[11] = 10;
    coeff2mst_res_J_[12] = 9;
    coeff2mst_res_J_[13] = 10;

    coeff2mst_res_n_[0] = 0.;
    // coefficient region 2 metastable residual n_i
    coeff2mst_res_n_[1] = -0.73362260186506 * 1.0e-2;
    coeff2mst_res_n_[2] = -0.88223831943146 * 1.0e-1;
    coeff2mst_res_n_[3] = -0.72334555213245 * 1.0e-1;
    coeff2mst_res_n_[4] = -0.40813178534455 * 1.0e-2;
    coeff2mst_res_n_[5] = 0.20097803380207 * 1.0e-2;
    coeff2mst_res_n_[6] = -0.53045921898642 * 1.0e-1;
    coeff2mst_res_n_[7] = -0.76190409086970 * 1.0e-2;
    coeff2mst_res_n_[8] = -0.63498037657313 * 1.0e-2;
    coeff2mst_res_n_[9] = -0.86043093028588 * 1.0e-1;
    coeff2mst_res_n_[10] = 0.75321581522770 * 1.0e-2;
    coeff2mst_res_n_[11] = -0.79238375446139 * 1.0e-2;
    coeff2mst_res_n_[12] = -0.22888160778447 * 1.0e-3;
    coeff2mst_res_n_[13] = -0.26456501482810 * 1.0e-2;

    temperature_ref4_ = 1.;
    pressure_ref4_ = 1.0e+6;

    coeff4_n_[0] = 0.;
    // coefficient region 4 n_i
    coeff4_n_[1] = 0.11670521452767 * 1.0e+4;
    coeff4_n_[2] = -0.72421316703206 * 1.0e+6;
    coeff4_n_[3] = -0.17073846940092 * 1.0e+2;
    coeff4_n_[4] = 0.12020824702470 * 1.0e+5;
    coeff4_n_[5] = -0.32325550322333 * 1.0e+7;
    coeff4_n_[6] = 0.14915108613530 * 1.0e+2;
    coeff4_n_[7] = -0.48232657361591 * 1.0e+4;
    coeff4_n_[8] = 0.40511340542057 * 1.0e+6;
    coeff4_n_[9] = -0.23855557567849;
    coeff4_n_[10] = 0.65017534844798 * 1.0e+3;

    return;
}

} // end namespace IAPWS
