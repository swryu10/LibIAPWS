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
                          double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double gamma;

    switch(n_reg) {
      case 1 :
        gamma = get_param1_gamma(temperature_in,
                                 pressure_in);
        break;
      case 2 :
        break;
      case 3 :
        break;
      case 5 :
        break;
      default :
        gamma = 0.;
    }

    double g =
        const_R_spec_ * temperature_in * gamma;

    return g;
}

double Lib97::get_param_vol_spec(double temperature_in,
                                 double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double ppi = 0.;

    double vol_spec;

    switch(n_reg) {
      case 1 :
        ppi = pressure_in / pressure_ref1_;
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * get_param1_dgamma_dppi(temperature_in,
                                         pressure_in);
        break;
      case 2 :
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
                                 double pressure_in) {
    double mdensity =
        1. / get_param_vol_spec(temperature_in,
                                pressure_in);

    return mdensity;
}

double Lib97::get_param_erg_int(double temperature_in,
                                double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;
    double ppi = 0.;

    double erg_int;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        ppi = pressure_in / pressure_ref1_;
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * get_param1_dgamma_dtau(temperature_in,
                                          pressure_in) -
             ppi * get_param1_dgamma_dppi(temperature_in,
                                          pressure_in));
        break;
      case 2 :
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
                                double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;

    double entropy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        entropy =
            const_R_spec_ *
            (tau * get_param1_dgamma_dtau(temperature_in,
                                          pressure_in) -
             get_param1_gamma(temperature_in,
                              pressure_in));
        break;
      case 2 :
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
                                 double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;

    double enthalpy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * get_param1_dgamma_dtau(temperature_in,
                                          pressure_in);
        break;
      case 2 :
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
                                 double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;
    double fn_d2gamma_dtau_dtau = 0.;

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
                                 double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi = 0.;
    double fn_d2gamma_dppi_dtau = 0.;
    double fn_d2gamma_dtau_dtau = 0.;

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
                                    double pressure_in) {
    int n_reg = get_region(temperature_in, pressure_in);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi = 0.;
    double fn_d2gamma_dppi_dtau = 0.;
    double fn_d2gamma_dtau_dtau = 0.;

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
                      double pressure_in) {
    int n_reg = 1;

    return n_reg;
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
        if (coeff1_I_[i] == 0 || coeff1_I_[i] == 1) {
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
        if (coeff1_I_[i] == 0 || coeff1_J_[i] == 0) {
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
        if (coeff1_J_[i] == 0 || coeff1_J_[i] == 1) {
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

void Lib97::set_coefficients() {
    temperature_crit_ = 647.096;
    mdensity_crit_ = 322.;
    pressure_crit_ = 22.064 * 1.0e+6;
    const_R_spec_ = 461.526;

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

    return;
}

} // end namespace IAPWS
