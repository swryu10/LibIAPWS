#include<stdio.h>
#include<math.h>
#include"LibIAPWS95.h"

namespace IAPWS {

void Lib95::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R6-95 (2018)\n");
    fprintf(ptr_fout, "Thermodynamic Properties of Ordinary Water Substance ");
    fprintf(ptr_fout, "for General and Scientific Use\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib95::get_param_phi_ide(double mdensity_in,
                                double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double phi =
        log(delta) +
        coeff_ide_n_[1] +
        coeff_ide_n_[2] * tau +
        coeff_ide_n_[3] * log(tau);
    for (int i = 4; i <= 8; i++) {
        double exp_gamma_tau = exp(-coeff_ide_gamma_[i] * tau);

        phi += coeff_ide_n_[i] * log(1. - exp_gamma_tau);
    }

    return phi;
}

double Lib95::get_param_dphi_ide_ddelta(double mdensity_in,
                                        double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double dphi_ddelta = 1. / delta;

    return dphi_ddelta;
}

double Lib95::get_param_dphi_ide_dtau(double mdensity_in,
                                      double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double dphi_dtau =
        coeff_ide_n_[2] +
        coeff_ide_n_[3] / tau;
    for (int i = 4; i <= 8; i++) {
        double exp_gamma_tau = exp(-coeff_ide_gamma_[i] * tau);

        dphi_dtau +=
            coeff_ide_n_[i] * coeff_ide_gamma_[i] *
            (1. / (1. - exp_gamma_tau) - 1.);
    }

    return dphi_dtau;
}

double Lib95::get_param_d2phi_ide_ddelta_ddelta(double mdensity_in,
                                                double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_ddelta = -1. / (delta * delta);

    return d2phi_ddelta_ddelta;
}

double Lib95::get_param_d2phi_ide_ddelta_dtau(double mdensity_in,
                                              double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_dtau = 0.;

    return d2phi_ddelta_dtau;
}

double Lib95::get_param_d2phi_ide_dtau_dtau(double mdensity_in,
                                            double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double d2phi_dtau_dtau = -coeff_ide_n_[3] / (tau * tau);
    for (int i = 4; i <= 8; i++) {
        double exp_gamma_tau = exp(-coeff_ide_gamma_[i] * tau);

        d2phi_dtau_dtau -=
            coeff_ide_n_[i] *
            coeff_ide_gamma_[i] * coeff_ide_gamma_[i] *
            exp_gamma_tau /
            ((1. - exp_gamma_tau) * (1. - exp_gamma_tau));
    }

    return d2phi_dtau_dtau;
}

double Lib95::get_param_phi_res(double mdensity_in,
                                double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double phi = 0.;
    for (int i = 1; i <= 7; i++) {
        phi +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]);
    }
    for (int i = 8; i <= 51; i++) {
        phi +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        phi +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt);
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        phi += coeff_res_n_[i] * fn_Delta_powb * delta * fn_psi;
    }

    return phi;
}

double Lib95::get_param_dphi_res_ddelta(double mdensity_in,
                                        double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double dphi_ddelta = 0.;
    for (int i = 1; i <= 7; i++) {
        dphi_ddelta +=
            coeff_res_n_[i] *
            coeff_res_d_[i] *
            pow(delta, coeff_res_d_[i] - 1.) *
            pow(tau, coeff_res_t_[i]);
    }
    for (int i = 8; i <= 51; i++) {
        dphi_ddelta +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i] - 1.) *
            pow(tau, coeff_res_t_[i]) *
            exp(-pow(delta, coeff_res_c_[i])) *
            (coeff_res_d_[i] - coeff_res_c_[i] * pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        dphi_ddelta +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt) *
            (coeff_res_d_[i] / delta - 2. * coeff_res_alpha_[i - 52] * dd);
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        double fn_dDelta_ddelta =
            (delta - 1.) * (coeff_res_A_[i - 55] * fn_Theta *
                            (2. / coeff_res_beta_[i - 52]) *
                            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.) +
                            2. * coeff_res_t_[i] * coeff_res_c_[i] *
                            pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 2.));
        double fn_dDelta_powb_ddelta =
            coeff_res_d_[i] * pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_ddelta;
        double fn_dpsi_ddelta =
            -2. * coeff_res_C_[i - 55] * (delta - 1.) * fn_psi;

        dphi_ddelta +=
            coeff_res_n_[i] *
            (fn_Delta_powb * (fn_psi + delta * fn_dpsi_ddelta) +
             fn_dDelta_powb_ddelta * delta * fn_psi);
    }

    return dphi_ddelta;
}

double Lib95::get_param_dphi_res_dtau(double mdensity_in,
                                      double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double dphi_dtau = 0.;
    for (int i = 1; i <= 7; i++) {
        dphi_dtau +=
            coeff_res_n_[i] *
            coeff_res_t_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i] - 1.);
    }
    for (int i = 8; i <= 51; i++) {
        dphi_dtau +=
            coeff_res_n_[i] *
            coeff_res_t_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i] - 1.) *
            exp(-pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        dphi_dtau +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt) *
            (coeff_res_t_[i] / tau - 2. * coeff_res_beta_[i - 52] * dt);
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        double fn_dDelta_dtau = -2. * fn_Theta;
        double fn_dDelta_powb_dtau =
            coeff_res_d_[i] * pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_dtau;
        double fn_dpsi_dtau =
            -2. * coeff_res_D_[i - 55] * (tau - 1.) * fn_psi;

        dphi_dtau +=
            coeff_res_n_[i] * delta *
            (fn_dDelta_powb_dtau * fn_psi + fn_Delta_powb * fn_dpsi_dtau);
    }

    return dphi_dtau;
}

double Lib95::get_param_d2phi_res_ddelta_ddelta(double mdensity_in,
                                                double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_ddelta = 0.;
    for (int i = 1; i <= 7; i++) {
        d2phi_ddelta_ddelta +=
            coeff_res_n_[i] *
            coeff_res_d_[i] * (coeff_res_d_[i] - 1.) *
            pow(delta, coeff_res_d_[i] - 2.) *
            pow(tau, coeff_res_t_[i]);
    }
    for (int i = 8; i <= 51; i++) {
        d2phi_ddelta_ddelta +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i] - 2.) *
            pow(tau, coeff_res_t_[i]) *
            exp(-pow(delta, coeff_res_c_[i])) *
            ((coeff_res_d_[i] -
                coeff_res_c_[i] * pow(delta, coeff_res_c_[i])) *
             (coeff_res_d_[i] - 1. -
                coeff_res_c_[i] * pow(delta, coeff_res_c_[i])) -
             coeff_res_c_[i] * coeff_res_c_[i] * pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        d2phi_ddelta_ddelta +=
            coeff_res_n_[i] *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt) *
            (-2. * coeff_res_alpha_[i - 52] * pow(delta, coeff_res_d_[i]) +
             4. * coeff_res_alpha_[i - 52] * coeff_res_alpha_[i - 52] *
                dd * dd * pow(delta, coeff_res_d_[i]) -
             4. * coeff_res_d_[i] * coeff_res_alpha_[i - 52] *
                dd * pow(delta, coeff_res_d_[i] - 1.) +
             coeff_res_d_[i] * (coeff_res_d_[i] - 1.) *
                pow(delta, coeff_res_d_[i] - 2.));
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        double fn_dDelta_ddelta =
            (delta - 1.) * (coeff_res_A_[i - 55] * fn_Theta *
                            (2. / coeff_res_beta_[i - 52]) *
                            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.) +
                            2. * coeff_res_t_[i] * coeff_res_c_[i] *
                            pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 2.));
        double fn_dDelta_powb_ddelta =
            coeff_res_d_[i] * pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_ddelta;
        double fn_dpsi_ddelta =
            -2. * coeff_res_C_[i - 55] * (delta - 1.) * fn_psi;

        double fn_d2Delta_ddelta_ddelta =
            fn_dDelta_ddelta / (delta - 1.) +
            (delta - 1.) * (delta - 1.) *
            (4. * coeff_res_t_[i] * coeff_res_c_[i] *
                (coeff_res_c_[i] - 1.) *
                pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 4.) +
             2. * (coeff_res_A_[i - 55] / coeff_res_beta_[i - 52]) *
                  (coeff_res_A_[i - 55] / coeff_res_beta_[i - 52]) *
                pow(fabs(delta - 1.), 2. / coeff_res_beta_[i - 52] - 4.) +
             coeff_res_A_[i - 55] * fn_Theta * (4. / coeff_res_beta_[i - 52]) *
                (1. / (2. * coeff_res_beta_[i - 52]) - 1.) *
                pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 4.));
        double fn_d2Delta_powb_ddelta_ddelta =
            coeff_res_d_[i] *
            (pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_d2Delta_ddelta_ddelta +
             (coeff_res_d_[i] - 1.) * pow(fn_Delta, coeff_res_d_[i] - 2.) *
                fn_dDelta_ddelta * fn_dDelta_ddelta);
        double fn_d2psi_ddelta_ddelta =
            (2. * coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) - 1.) *
             2. * coeff_res_C_[i - 55] * fn_psi;

        d2phi_ddelta_ddelta +=
            coeff_res_n_[i] *
            (fn_Delta_powb *
                (2. * fn_dpsi_ddelta + delta * fn_d2psi_ddelta_ddelta) +
             2. * fn_dDelta_powb_ddelta *
                (fn_psi + delta * fn_dpsi_ddelta) +
             fn_d2Delta_powb_ddelta_ddelta * delta * fn_psi);
    }

    return d2phi_ddelta_ddelta;
}

double Lib95::get_param_d2phi_res_ddelta_dtau(double mdensity_in,
                                              double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_dtau = 0.;
    for (int i = 1; i <= 7; i++) {
        d2phi_ddelta_dtau +=
            coeff_res_n_[i] *
            coeff_res_d_[i] * coeff_res_t_[i] *
            pow(delta, coeff_res_d_[i] - 1.) *
            pow(tau, coeff_res_t_[i] - 1.);
    }
    for (int i = 8; i <= 51; i++) {
        d2phi_ddelta_dtau +=
            coeff_res_n_[i] *
            coeff_res_t_[i] *
            pow(delta, coeff_res_d_[i] - 1.) *
            pow(tau, coeff_res_t_[i] - 1.) *
            exp(-pow(delta, coeff_res_c_[i])) *
            (coeff_res_d_[i] - coeff_res_c_[i] * pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        d2phi_ddelta_dtau +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt) *
            (coeff_res_d_[i] / delta - 2. * coeff_res_alpha_[i - 52] * dd) *
            (coeff_res_t_[i] / tau - 2. * coeff_res_beta_[i - 52] * dt);
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        double fn_dDelta_ddelta =
            (delta - 1.) * (coeff_res_A_[i - 55] * fn_Theta *
                            (2. / coeff_res_beta_[i - 52]) *
                            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.) +
                            2. * coeff_res_t_[i] * coeff_res_c_[i] *
                            pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 2.));
        double fn_dDelta_powb_ddelta =
            coeff_res_d_[i] * pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_ddelta;
        double fn_dpsi_ddelta =
            -2. * coeff_res_C_[i - 55] * (delta - 1.) * fn_psi;
        double fn_dDelta_dtau = -2. * fn_Theta;
        double fn_dDelta_powb_dtau =
            coeff_res_d_[i] * pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_dtau;
        double fn_dpsi_dtau =
            -2. * coeff_res_D_[i - 55] * (tau - 1.) * fn_psi;

        double fn_d2Delta_powb_ddelta_dtau =
            -coeff_res_A_[i - 55] * coeff_res_d_[i] *
                (2. / coeff_res_beta_[i - 52]) *
                pow(fn_Delta, coeff_res_d_[i] - 1.) * (delta - 1.) *
                pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.) -
            2. * fn_Theta * coeff_res_d_[i] * (coeff_res_d_[i] - 1.) *
                pow(fn_Delta, coeff_res_d_[i] - 2.) * fn_dDelta_ddelta;
        double fn_d2psi_ddelta_dtau =
            4. * coeff_res_C_[i - 55] * (delta - 1.) *
                 coeff_res_D_[i - 55] * (tau - 1.) * fn_psi;

        d2phi_ddelta_dtau +=
            coeff_res_n_[i] *
            (fn_Delta_powb *
                (fn_dpsi_dtau + delta * fn_d2psi_ddelta_dtau) +
             delta * fn_dDelta_powb_ddelta * fn_dpsi_dtau +
             fn_dDelta_powb_dtau *
                (fn_psi + delta * fn_dpsi_ddelta) +
             fn_d2Delta_powb_ddelta_dtau * delta * fn_psi);
    }

    return d2phi_ddelta_dtau;
}

double Lib95::get_param_d2phi_res_dtau_dtau(double mdensity_in,
                                            double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double d2phi_dtau_dtau = 0.;
    for (int i = 1; i <= 7; i++) {
        d2phi_dtau_dtau +=
            coeff_res_n_[i] *
            coeff_res_t_[i] * (coeff_res_t_[i] - 1.) *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i] - 2.);
    }
    for (int i = 8; i <= 51; i++) {
        d2phi_dtau_dtau +=
            coeff_res_n_[i] *
            coeff_res_t_[i] * (coeff_res_t_[i] - 1.) *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i] - 2.) *
            exp(-pow(delta, coeff_res_c_[i]));
    }
    for (int i = 52; i <= 54; i++) {
        double dd = delta - coeff_res_epsilon_[i - 52];
        double dt = tau - coeff_res_gamma_[i - 52];

        d2phi_dtau_dtau +=
            coeff_res_n_[i] *
            pow(delta, coeff_res_d_[i]) *
            pow(tau, coeff_res_t_[i]) *
            exp(-coeff_res_alpha_[i - 52] * dd * dd) *
            exp(-coeff_res_beta_[i - 52] * dt * dt) *
            ((coeff_res_t_[i] / tau - 2. * coeff_res_beta_[i - 52] * dt) *
             (coeff_res_t_[i] / tau - 2. * coeff_res_beta_[i - 52] * dt) -
             coeff_res_t_[i] / (tau * tau) -
             2. * coeff_res_beta_[i - 52]);
    }
    for (int i = 55; i <= 56; i++) {
        double fn_Theta =
            1. - tau +
            coeff_res_A_[i - 55] *
            pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52]);
        double fn_Delta =
            fn_Theta * fn_Theta +
            coeff_res_t_[i] *
            pow(fabs(delta - 1.), 2. * coeff_res_c_[i]);
        double fn_Delta_powb = pow(fn_Delta, coeff_res_d_[i]);
        double fn_psi =
            exp(-coeff_res_C_[i - 55] * (delta - 1.) * (delta - 1.) -
                 coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.));

        double fn_dDelta_dtau = -2. * fn_Theta;
        double fn_dDelta_powb_dtau =
            coeff_res_d_[i] *
            pow(fn_Delta, coeff_res_d_[i] - 1.) * fn_dDelta_dtau;
        double fn_dpsi_dtau =
            -2. * coeff_res_D_[i - 55] * (tau - 1.) * fn_psi;

        double fn_d2Delta_powb_dtau_dtau =
            2. * coeff_res_d_[i] *
                pow(fn_Delta, coeff_res_d_[i] - 1.) +
            4. * fn_Theta * fn_Theta *
                coeff_res_d_[i] * (coeff_res_d_[i] - 1.) *
                pow(fn_Delta, coeff_res_d_[i] - 2.);
        double fn_d2psi_dtau_dtau =
            (2. * coeff_res_D_[i - 55] * (tau - 1.) * (tau - 1.) - 1.) *
             2. * coeff_res_D_[i - 55] * fn_psi;

        d2phi_dtau_dtau +=
            coeff_res_n_[i] * delta *
            (fn_d2Delta_powb_dtau_dtau * fn_psi +
             2. * fn_dDelta_powb_dtau * fn_dpsi_dtau +
             fn_Delta_powb * fn_d2psi_dtau_dtau);
    }

    return d2phi_dtau_dtau;
}

double Lib95::get_param_f(double mdensity_in,
                          double temperature_in) {
    double f =
        const_R_spec_ * temperature_in *
        (get_param_phi_ide(mdensity_in,
                           temperature_in) +
         get_param_phi_res(mdensity_in,
                           temperature_in));

    return f;
}

double Lib95::get_param_g(double mdensity_in,
                          double temperature_in) {
    double g =
        get_param_f(mdensity_in,
                    temperature_in) +
        get_param_pressure(mdensity_in,
                           temperature_in) / mdensity_in;

    return g;
}

double Lib95::get_param_pressure(double mdensity_in,
                                 double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double press =
        mdensity_in * const_R_spec_ * temperature_in *
        (1. + delta * get_param_dphi_res_ddelta(mdensity_in,
                                                temperature_in));

    return press;
}

double Lib95::get_param_dpress_drho(double mdensity_in,
                                    double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double dpress_drho =
        const_R_spec_ * temperature_in *
        (1. + 2. * delta * get_param_dphi_res_ddelta(mdensity_in,
                                                     temperature_in) +
         delta * delta * get_param_d2phi_res_ddelta_ddelta(mdensity_in,
                                                           temperature_in));

    return dpress_drho;
}

double Lib95::get_param_erg_int(double mdensity_in,
                                double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double erg_int =
        const_R_spec_ * temperature_in * tau *
        (get_param_dphi_ide_dtau(mdensity_in, temperature_in) +
         get_param_dphi_res_dtau(mdensity_in, temperature_in));

    return erg_int;
}

double Lib95::get_param_entropy(double mdensity_in,
                                double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double entropy =
        const_R_spec_ *
        (tau * (get_param_dphi_ide_dtau(mdensity_in,
                                        temperature_in) +
                get_param_dphi_res_dtau(mdensity_in,
                                        temperature_in)) -
         get_param_phi_ide(mdensity_in, temperature_in) -
         get_param_phi_res(mdensity_in, temperature_in));

    return entropy;
}

double Lib95::get_param_enthalpy(double mdensity_in,
                                 double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double enthalpy =
        const_R_spec_ * temperature_in *
        (1. + tau *
         (get_param_dphi_ide_dtau(mdensity_in, temperature_in) +
          get_param_dphi_res_dtau(mdensity_in, temperature_in)) +
         delta * get_param_dphi_res_ddelta(mdensity_in, temperature_in));

    return enthalpy;
}

double Lib95::get_param_heat_c_v(double mdensity_in,
                                 double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double c_v =
        -const_R_spec_ * tau * tau *
        (get_param_d2phi_ide_dtau_dtau(mdensity_in, temperature_in) +
         get_param_d2phi_res_dtau_dtau(mdensity_in, temperature_in));

    return c_v;
}

double Lib95::get_param_heat_c_p(double mdensity_in,
                                 double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double dphi_res_ddelta =
        get_param_dphi_res_ddelta(mdensity_in, temperature_in);

    double fac_nom =
        1. + delta * dphi_res_ddelta -
        delta * tau * get_param_d2phi_res_ddelta_dtau(mdensity_in,
                                                      temperature_in);
    double fac_den =
        1. + 2. * delta * dphi_res_ddelta +
        delta * delta * get_param_d2phi_res_ddelta_ddelta(mdensity_in,
                                                          temperature_in);

    double c_p =
        const_R_spec_ *
        (fac_nom * fac_nom / fac_den -
         tau * tau * (get_param_d2phi_ide_dtau_dtau(mdensity_in,
                                                    temperature_in) +
                      get_param_d2phi_res_dtau_dtau(mdensity_in,
                                                    temperature_in)));

    return c_p;
}

double Lib95::get_param_speed_sound(double mdensity_in,
                                    double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double dphi_res_ddelta =
        get_param_dphi_res_ddelta(mdensity_in, temperature_in);

    double fac_nom =
        1. + delta * dphi_res_ddelta -
        delta * tau * get_param_d2phi_res_ddelta_dtau(mdensity_in,
                                                      temperature_in);
    double fac_den =
        tau * tau * (get_param_d2phi_ide_dtau_dtau(mdensity_in,
                                                   temperature_in) +
                     get_param_d2phi_res_dtau_dtau(mdensity_in,
                                                   temperature_in));

    double v2_s =
        const_R_spec_ * temperature_in *
        (1. + 2. * delta * dphi_res_ddelta +
         delta * delta * get_param_d2phi_res_ddelta_ddelta(mdensity_in,
                                                           temperature_in) -
         fac_nom * fac_nom / fac_den);

    return sqrt(v2_s);
}

bool Lib95::find_root_mdensity(double temperature_in,
                               double pressure_in,
                               double &mdensity_out) {
    double mden_now = mdensity_out;
    double press_now =
        get_param_pressure(mden_now, temperature_in);

    int i_iter = 0;
    bool found_root = false;
    while (!found_root) {
        i_iter += 1;

        double mden_prev = mden_now;
        double press_prev = press_now;

        double dpress_drho_prev =
            get_param_dpress_drho(mden_prev, temperature_in);

        mden_now = mden_prev +
            0.5 * (pressure_in - press_prev) / dpress_drho_prev;
        press_now =
            get_param_pressure(mden_now, temperature_in);

        if (fabs(press_now - pressure_in) <
            0.5 * eps_precision_ * fabs(press_now + pressure_in)) {
            mdensity_out = mden_now;

            found_root = true;
        }

        if (i_iter > n_iter_max_) {
            break;
        }
    }

    /*
    fprintf(stderr, "        temperature = %.8e degK, ", temperature_in);
    fprintf(stderr, "        pressure = %.8e Pa, ", pressure_in);
    fprintf(stderr, "        mdensity = %.8e kg / m^3\n", mdensity_out);
    */

    return found_root;
}

bool Lib95::find_state_coex(double temperature_in,
                            double &pressure_out,
                            double &mden_vap_out,
                            double &mden_liq_out) {
    double press_now = pressure_out;
    double mden_vap = mden_vap_out;
    double mden_liq = mden_liq_out;

    int i_iter = 0;
    bool found_state = false;
    while (!found_state) {
        i_iter += 1;

        double press_prev = press_now;

        bool found_root_vap =
            find_root_mdensity(temperature_in,
                               press_prev,
                               mden_vap);
        bool found_root_liq =
            find_root_mdensity(temperature_in,
                               press_prev,
                               mden_liq);

        if (!found_root_vap || !found_root_liq) {
            break;
        }

        double vol_vap = 1. / mden_vap;
        double vol_liq = 1. / mden_liq;

        double diff_erg_free =
            const_R_spec_ * temperature_in *
            (log(mden_liq / mden_vap) +
             get_param_phi_res(mden_liq, temperature_in) -
             get_param_phi_res(mden_vap, temperature_in));

        double press_esti =
            diff_erg_free / (vol_vap - vol_liq);

        if (fabs(press_prev - press_esti) <
            0.5 * eps_precision_ * fabs(press_prev + press_esti)) {
            pressure_out = press_prev;
            mden_vap_out = mden_vap;
            mden_liq_out = mden_liq;

            found_state = true;
        }

        if (i_iter > n_iter_max_) {
            break;
        }

        press_now = press_prev +
            0.1 * (press_esti - press_prev);
    }

    /*
    fprintf(stderr, "      temperature = %.8e degK, ", temperature_in);
    fprintf(stderr, "      pressure = %.8e Pa, ", pressure_out);
    fprintf(stderr, "      mden_vap = %.8e kg / m^3, ", mden_vap_out);
    fprintf(stderr, "      mden_liq = %.8e kg / m^3\n", mden_liq_out);
    */

    return found_state;
}

void Lib95::make_tab_coex(int nbin_in,
                          double temperature_max) {
    reset_tab_coex();

    if (temperature_max < temperature_trip_) {
        return;
    }

    nbin_coex_ = nbin_in;
    temperature_coex_min_ = temperature_trip_;
    temperature_coex_max_ = temperature_max;

    double delta_temp =
        (temperature_coex_max_ - temperature_coex_min_) /
        static_cast<double>(nbin_coex_);

    have_tab_coex_ = true;
    tab_coex_temperature_ = new double[nbin_coex_ + 1];
    tab_coex_pressure_ = new double[nbin_coex_ + 1];
    tab_coex_mden_vap_ = new double[nbin_coex_ + 1];
    tab_coex_mden_liq_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_liq_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_liq_ = new double[nbin_coex_ + 1];

    bool made_tab = true;
    if (flag_verbose_) {
        fprintf(stderr, "    temperature (degK), ");
        fprintf(stderr, "    pressure (Pa), ");
        fprintf(stderr, "    mden_vap (kg / m^3), ");
        fprintf(stderr, "    mden_liq (kg / m^3)\n");
    }
    for (int it = 0; it <= nbin_coex_; it++) {
        tab_coex_temperature_[it] =
            temperature_coex_min_ +
            delta_temp * static_cast<double>(it);

        if (it == 0) {
            tab_coex_pressure_[it] = 600.;
            tab_coex_mden_vap_[it] = 0.005;
            tab_coex_mden_liq_[it] = 1100.;
        } else {
            tab_coex_pressure_[it] = tab_coex_pressure_[it - 1];
            tab_coex_mden_vap_[it] = tab_coex_mden_vap_[it - 1];
            tab_coex_mden_liq_[it] = tab_coex_mden_liq_[it - 1];
        }

        bool found_state =
            find_state_coex(tab_coex_temperature_[it],
                            tab_coex_pressure_[it],
                            tab_coex_mden_vap_[it],
                            tab_coex_mden_liq_[it]);

        if (!found_state) {
            made_tab = false;
            break;
        }

        tab_coex_enthalpy_vap_[it] =
            get_param_enthalpy(tab_coex_mden_vap_[it],
                               tab_coex_temperature_[it]);
        tab_coex_enthalpy_liq_[it] =
            get_param_enthalpy(tab_coex_mden_liq_[it],
                               tab_coex_temperature_[it]);

        tab_coex_entropy_vap_[it] =
            get_param_entropy(tab_coex_mden_vap_[it],
                              tab_coex_temperature_[it]);
        tab_coex_entropy_liq_[it] =
            get_param_entropy(tab_coex_mden_liq_[it],
                              tab_coex_temperature_[it]);

        if (flag_verbose_) {
            fprintf(stderr, "      %.8e",
                    tab_coex_temperature_[it]);
            fprintf(stderr, "      %.8e",
                    tab_coex_pressure_[it]);
            fprintf(stderr, "      %.8e",
                    tab_coex_mden_vap_[it]);
            fprintf(stderr, "      %.8e\n",
                    tab_coex_mden_liq_[it]);
        }
    }

    if (!made_tab) {
        reset_tab_coex();
    }

    set_cspline_coex();

    return;
}

void Lib95::export_tab_coex(char *filename) {
    if (!have_tab_coex_) {
        return;
    }

    FILE *ptr_fout;
    ptr_fout = fopen(filename, "w");
    if (ptr_fout == NULL) {
        return;
    }

    fprintf(ptr_fout, "# table for coexisting phases in IAPWS95\n");
    fprintf(ptr_fout, "nbin_coex    %d\n", nbin_coex_);
    fprintf(ptr_fout, "temperature_coex_min    %e\n",
            temperature_coex_min_);
    fprintf(ptr_fout, "temperature_coex_max    %e\n",
            temperature_coex_max_);
    fprintf(ptr_fout, "# temperature (degK)");
    fprintf(ptr_fout, "    pressure (Pa)");
    fprintf(ptr_fout, "    mden_vap (kg / m^3)");
    fprintf(ptr_fout, "    mden_liq (kg / m^3)");
    fprintf(ptr_fout, "    enthalpy_vap (J / kg)");
    fprintf(ptr_fout, "    enthalpy_liq (J / kg)");
    fprintf(ptr_fout, "    entropy_vap (J / kg / degK)");
    fprintf(ptr_fout, "    entropy_liq (J / kg / degK)");
    fprintf(ptr_fout, "\n");


    for (int it = 0; it <= nbin_coex_; it++) {
        fprintf(ptr_fout, "    %.8e",
                tab_coex_temperature_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_pressure_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_mden_vap_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_mden_liq_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_enthalpy_vap_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_enthalpy_liq_[it]);
        fprintf(ptr_fout, "    %.8e",
                tab_coex_entropy_vap_[it]);
        fprintf(ptr_fout, "    %.8e\n",
                tab_coex_entropy_liq_[it]);
    }

    fclose(ptr_fout);

    return;
}

void Lib95::import_tab_coex(char *filename) {
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
    tab_coex_mden_liq_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_enthalpy_liq_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_vap_ = new double[nbin_coex_ + 1];
    tab_coex_entropy_liq_ = new double[nbin_coex_ + 1];

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
               &tab_coex_mden_liq_[it],
               &tab_coex_enthalpy_vap_[it],
               &tab_coex_enthalpy_liq_[it],
               &tab_coex_entropy_vap_[it],
               &tab_coex_entropy_liq_[it]);

        if (flag_verbose_) {
            fprintf(stderr, "    %.8e",
                    tab_coex_temperature_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_pressure_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_mden_vap_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_mden_liq_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_enthalpy_vap_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_enthalpy_liq_[it]);
            fprintf(stderr, "    %.8e",
                    tab_coex_entropy_vap_[it]);
            fprintf(stderr, "    %.8e\n",
                    tab_coex_entropy_liq_[it]);
        }
    }

    fclose(ptr_fin);

    if (!made_tab) {
        reset_tab_coex();
    }

    set_cspline_coex();

    return;
}

void Lib95::set_cspline_coex() {
    if (!have_tab_coex_) {
        return;
    }

    csp_coex_pressure_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_pressure_);
    csp_coex_mden_vap_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_mden_vap_);
    csp_coex_mden_liq_.init(nbin_coex_,
                            tab_coex_temperature_,
                            tab_coex_mden_liq_);
    csp_coex_enthalpy_vap_.init(nbin_coex_,
                                tab_coex_temperature_,
                                tab_coex_enthalpy_vap_);
    csp_coex_enthalpy_liq_.init(nbin_coex_,
                                tab_coex_temperature_,
                                tab_coex_enthalpy_liq_);
    csp_coex_entropy_vap_.init(nbin_coex_,
                               tab_coex_temperature_,
                               tab_coex_entropy_vap_);
    csp_coex_entropy_liq_.init(nbin_coex_,
                               tab_coex_temperature_,
                               tab_coex_entropy_liq_);

    return;
}

} // end namespace IAPWS
