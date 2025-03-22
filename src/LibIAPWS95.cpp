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
            coeff_res_A_[i - 55] * fn_Theta *
                (2. / coeff_res_beta_[i - 52]) *
                pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.) +
            2. * coeff_res_t_[i] * coeff_res_c_[i] *
                pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 2.) +
            4. * coeff_res_t_[i] * coeff_res_c_[i] *
                (coeff_res_c_[i] - 1.) *
                pow(fabs(delta - 1.), 2. * coeff_res_c_[i] - 2.) +
            2. * (coeff_res_A_[i - 55] / coeff_res_beta_[i - 52]) *
                  (coeff_res_A_[i - 55] / coeff_res_beta_[i - 52]) *
                pow(fabs(delta - 1.), 2. / coeff_res_beta_[i - 52] - 2.) +
            coeff_res_A_[i - 55] * fn_Theta * (4. / coeff_res_beta_[i - 52]) *
                (1. / (2. * coeff_res_beta_[i - 52]) - 1.) *
                pow(fabs(delta - 1.), 1. / coeff_res_beta_[i - 52] - 2.);
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

    double *frac_temp = new double[2];

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
            if (tab_coex_temperature_[it] > 640.) {
                frac_temp[0] =
                    (temperature_crit_ -
                     tab_coex_temperature_[it]) /
                    (temperature_crit_ -
                     tab_coex_temperature_[it - 1]);
                frac_temp[1] =
                    (tab_coex_temperature_[it] -
                     tab_coex_temperature_[it - 1]) /
                    (temperature_crit_ -
                     tab_coex_temperature_[it - 1]);

                tab_coex_pressure_[it] =
                    tab_coex_pressure_[it - 1] * frac_temp[0] +
                    pressure_crit_ * frac_temp[1];
                tab_coex_mden_vap_[it] =
                    tab_coex_mden_vap_[it - 1] * frac_temp[0] +
                    mdensity_crit_ * frac_temp[1];
                tab_coex_mden_liq_[it] =
                    tab_coex_mden_liq_[it - 1] * frac_temp[0] +
                    mdensity_crit_ * frac_temp[1];
            } else {
                tab_coex_pressure_[it] = tab_coex_pressure_[it - 1];
                tab_coex_mden_vap_[it] = tab_coex_mden_vap_[it - 1];
                tab_coex_mden_liq_[it] = tab_coex_mden_liq_[it - 1];
            }
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

    delete [] frac_temp;

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

void Lib95::set_coefficients() {
    temperature_crit_ = 647.096;
    temperature_trip_ = 273.16;
    mdensity_crit_ = 322.;
    pressure_crit_ = 22.064 * 1.0e+6;
    const_R_spec_ = 461.51805;

    coeff_ide_n_[0] = 0.;
    // coefficient ide n_i
    coeff_ide_n_[1] = -8.3204464837497;
    coeff_ide_n_[2] = 6.6832105275932;
    coeff_ide_n_[3] = 3.00632;
    coeff_ide_n_[4] = 0.012436;
    coeff_ide_n_[5] = 0.97315;
    coeff_ide_n_[6] = 1.27950;
    coeff_ide_n_[7] = 0.96956;
    coeff_ide_n_[8] = 0.24873;

    coeff_ide_gamma_[0] = 0.;
    // coefficient ide gamma_i
    coeff_ide_gamma_[1] = 0.;
    coeff_ide_gamma_[2] = 0.;
    coeff_ide_gamma_[3] = 0.;
    coeff_ide_gamma_[4] = 1.28728967;
    coeff_ide_gamma_[5] = 3.53734222;
    coeff_ide_gamma_[6] = 7.74073708;
    coeff_ide_gamma_[7] = 9.24437796;
    coeff_ide_gamma_[8] = 27.5075105;

    coeff_res_c_[0] = 0.;
    // coefficient res c_i
    coeff_res_c_[1] = 0.;
    coeff_res_c_[2] = 0.;
    coeff_res_c_[3] = 0.;
    coeff_res_c_[4] = 0.;
    coeff_res_c_[5] = 0.;
    coeff_res_c_[6] = 0.;
    coeff_res_c_[7] = 0.;
    coeff_res_c_[8] = 1.;
    coeff_res_c_[9] = 1.;
    coeff_res_c_[10] = 1.;
    coeff_res_c_[11] = 1.;
    coeff_res_c_[12] = 1.;
    coeff_res_c_[13] = 1.;
    coeff_res_c_[14] = 1.;
    coeff_res_c_[15] = 1.;
    coeff_res_c_[16] = 1.;
    coeff_res_c_[17] = 1.;
    coeff_res_c_[18] = 1.;
    coeff_res_c_[19] = 1.;
    coeff_res_c_[20] = 1.;
    coeff_res_c_[21] = 1.;
    coeff_res_c_[22] = 1.;
    coeff_res_c_[23] = 2.;
    coeff_res_c_[24] = 2.;
    coeff_res_c_[25] = 2.;
    coeff_res_c_[26] = 2.;
    coeff_res_c_[27] = 2.;
    coeff_res_c_[28] = 2.;
    coeff_res_c_[29] = 2.;
    coeff_res_c_[30] = 2.;
    coeff_res_c_[31] = 2.;
    coeff_res_c_[32] = 2.;
    coeff_res_c_[33] = 2.;
    coeff_res_c_[34] = 2.;
    coeff_res_c_[35] = 2.;
    coeff_res_c_[36] = 2.;
    coeff_res_c_[37] = 2.;
    coeff_res_c_[38] = 2.;
    coeff_res_c_[39] = 2.;
    coeff_res_c_[40] = 2.;
    coeff_res_c_[41] = 2.;
    coeff_res_c_[42] = 2.;
    coeff_res_c_[43] = 3.;
    coeff_res_c_[44] = 3.;
    coeff_res_c_[45] = 3.;
    coeff_res_c_[46] = 3.;
    coeff_res_c_[47] = 4.;
    coeff_res_c_[48] = 6.;
    coeff_res_c_[49] = 6.;
    coeff_res_c_[50] = 6.;
    coeff_res_c_[51] = 6.;
    coeff_res_c_[52] = 0.;
    coeff_res_c_[53] = 0.;
    coeff_res_c_[54] = 0.;
    // coefficient res a_i
    coeff_res_c_[55] = 3.5;
    coeff_res_c_[56] = 3.5;

    coeff_res_d_[0] = 0.;
    // coefficient res d_i
    coeff_res_d_[1] = 1.;
    coeff_res_d_[2] = 1.;
    coeff_res_d_[3] = 1.;
    coeff_res_d_[4] = 2.;
    coeff_res_d_[5] = 2.;
    coeff_res_d_[6] = 3.;
    coeff_res_d_[7] = 4.;
    coeff_res_d_[8] = 1.;
    coeff_res_d_[9] = 1.;
    coeff_res_d_[10] = 1.;
    coeff_res_d_[11] = 2.;
    coeff_res_d_[12] = 2.;
    coeff_res_d_[13] = 3.;
    coeff_res_d_[14] = 4.;
    coeff_res_d_[15] = 4.;
    coeff_res_d_[16] = 5.;
    coeff_res_d_[17] = 7.;
    coeff_res_d_[18] = 9.;
    coeff_res_d_[19] = 10.;
    coeff_res_d_[20] = 11.;
    coeff_res_d_[21] = 13.;
    coeff_res_d_[22] = 15.;
    coeff_res_d_[23] = 1.;
    coeff_res_d_[24] = 2.;
    coeff_res_d_[25] = 2.;
    coeff_res_d_[26] = 2.;
    coeff_res_d_[27] = 3.;
    coeff_res_d_[28] = 4.;
    coeff_res_d_[29] = 4.;
    coeff_res_d_[30] = 4.;
    coeff_res_d_[31] = 5.;
    coeff_res_d_[32] = 6.;
    coeff_res_d_[33] = 6.;
    coeff_res_d_[34] = 7.;
    coeff_res_d_[35] = 9.;
    coeff_res_d_[36] = 9.;
    coeff_res_d_[37] = 9.;
    coeff_res_d_[38] = 9.;
    coeff_res_d_[39] = 9.;
    coeff_res_d_[40] = 10.;
    coeff_res_d_[41] = 10.;
    coeff_res_d_[42] = 12.;
    coeff_res_d_[43] = 3.;
    coeff_res_d_[44] = 4.;
    coeff_res_d_[45] = 4.;
    coeff_res_d_[46] = 5.;
    coeff_res_d_[47] = 14.;
    coeff_res_d_[48] = 3.;
    coeff_res_d_[49] = 6.;
    coeff_res_d_[50] = 6.;
    coeff_res_d_[51] = 6.;
    coeff_res_d_[52] = 3.;
    coeff_res_d_[53] = 3.;
    coeff_res_d_[54] = 3.;
    // coefficient res b_i
    coeff_res_d_[55] = 0.85;
    coeff_res_d_[56] = 0.95;

    coeff_res_t_[0] = 0.;
    // coefficient res t_i
    coeff_res_t_[1] = -0.5;
    coeff_res_t_[2] = 0.875;
    coeff_res_t_[3] = 1.;
    coeff_res_t_[4] = 0.5;
    coeff_res_t_[5] = 0.75;
    coeff_res_t_[6] = 0.375;
    coeff_res_t_[7] = 1.;
    coeff_res_t_[8] = 4.;
    coeff_res_t_[9] = 6.;
    coeff_res_t_[10] = 12.;
    coeff_res_t_[11] = 1.;
    coeff_res_t_[12] = 5.;
    coeff_res_t_[13] = 4.;
    coeff_res_t_[14] = 2.;
    coeff_res_t_[15] = 13.;
    coeff_res_t_[16] = 9.;
    coeff_res_t_[17] = 3.;
    coeff_res_t_[18] = 4.;
    coeff_res_t_[19] = 11.;
    coeff_res_t_[20] = 4.;
    coeff_res_t_[21] = 13.;
    coeff_res_t_[22] = 1.;
    coeff_res_t_[23] = 7.;
    coeff_res_t_[24] = 1.;
    coeff_res_t_[25] = 9.;
    coeff_res_t_[26] = 10.;
    coeff_res_t_[27] = 10.;
    coeff_res_t_[28] = 3.;
    coeff_res_t_[29] = 7.;
    coeff_res_t_[30] = 10.;
    coeff_res_t_[31] = 10.;
    coeff_res_t_[32] = 6.;
    coeff_res_t_[33] = 10.;
    coeff_res_t_[34] = 10.;
    coeff_res_t_[35] = 1.;
    coeff_res_t_[36] = 2.;
    coeff_res_t_[37] = 3.;
    coeff_res_t_[38] = 4.;
    coeff_res_t_[39] = 8.;
    coeff_res_t_[40] = 6.;
    coeff_res_t_[41] = 9.;
    coeff_res_t_[42] = 8.;
    coeff_res_t_[43] = 16.;
    coeff_res_t_[44] = 22.;
    coeff_res_t_[45] = 23.;
    coeff_res_t_[46] = 23.;
    coeff_res_t_[47] = 10.;
    coeff_res_t_[48] = 50.;
    coeff_res_t_[49] = 44.;
    coeff_res_t_[50] = 46.;
    coeff_res_t_[51] = 50.;
    coeff_res_t_[52] = 0.;
    coeff_res_t_[53] = 1.;
    coeff_res_t_[54] = 4.;
    // coefficient res B_i
    coeff_res_t_[55] = 0.2;
    coeff_res_t_[56] = 0.2;

    coeff_res_n_[0] = 0.;
    // coefficient res n_i
    coeff_res_n_[1] = 0.12533547935523 * 1.0e-1;
    coeff_res_n_[2] = 0.78957634722828 * 1.0e+1;
    coeff_res_n_[3] = -0.87803203303561 * 1.0e+1;
    coeff_res_n_[4] = 0.31802509345418;
    coeff_res_n_[5] = -0.26145533859358;
    coeff_res_n_[6] = -0.78199751687981 * 1.0e-2;
    coeff_res_n_[7] = 0.88089493102134 * 1.0e-2;
    coeff_res_n_[8] = -0.66856572307965;
    coeff_res_n_[9] = 0.20433810950965;
    coeff_res_n_[10] = -0.66212605039687 * 1.0e-4;
    coeff_res_n_[11] = -0.19232721156002;
    coeff_res_n_[12] = -0.25709043003438;
    coeff_res_n_[13] = 0.16074868486251;
    coeff_res_n_[14] = -0.40092828925807 * 1.0e-1;
    coeff_res_n_[15] = 0.39343422603254 * 1.0e-6;
    coeff_res_n_[16] = -0.75941377088144 * 1.0e-5;
    coeff_res_n_[17] = 0.56250979351888 * 1.0e-3;
    coeff_res_n_[18] = -0.15608652257135 * 1.0e-4;
    coeff_res_n_[19] = 0.11537996422951 * 1.0e-8;
    coeff_res_n_[20] = 0.36582165144204 * 1.0e-6;
    coeff_res_n_[21] = -0.13251180074668 * 1.0e-11;
    coeff_res_n_[22] = -0.62639586912454 * 1.0e-9;
    coeff_res_n_[23] = -0.10793600908932;
    coeff_res_n_[24] = 0.17611491008752 * 1.0e-1;
    coeff_res_n_[25] = 0.22132295167546;
    coeff_res_n_[26] = -0.40247669763528;
    coeff_res_n_[27] = 0.58083399985759;
    coeff_res_n_[28] = 0.49969146990806 * 1.0e-2;
    coeff_res_n_[29] = -0.31358700712549 * 1.0e-1;
    coeff_res_n_[30] = -0.74315929710341;
    coeff_res_n_[31] = 0.47807329915480;
    coeff_res_n_[32] = 0.20527940895948 * 1.0e-1;
    coeff_res_n_[33] = -0.13636435110343;
    coeff_res_n_[34] = 0.14180634400617 * 1.0e-1;
    coeff_res_n_[35] = 0.83326504880713 * 1.0e-2;
    coeff_res_n_[36] = -0.29052336009585 * 1.0e-1;
    coeff_res_n_[37] = 0.38615085574206 * 1.0e-1;
    coeff_res_n_[38] = -0.20393486513704 * 1.0e-1;
    coeff_res_n_[39] = -0.16554050063734 * 1.0e-2;
    coeff_res_n_[40] = 0.19955571979541 * 1.0e-2;
    coeff_res_n_[41] = 0.15870308324157 * 1.0e-3;
    coeff_res_n_[42] = -0.16388568342530 * 1.0e-4;
    coeff_res_n_[43] = 0.43613615723811 * 1.0e-1;
    coeff_res_n_[44] = 0.34994005463765 * 1.0e-1;
    coeff_res_n_[45] = -0.76788197844621 * 1.0e-1;
    coeff_res_n_[46] = 0.22446277332006 * 1.0e-1;
    coeff_res_n_[47] = -0.62689710414685 * 1.0e-4;
    coeff_res_n_[48] = -0.55711118565645 * 1.0e-9;
    coeff_res_n_[49] = -0.19905718354408;
    coeff_res_n_[50] = 0.31777497330738;
    coeff_res_n_[51] = -0.11841182425981;
    coeff_res_n_[52] = -0.31306260323435 * 1.0e+2;
    coeff_res_n_[53] = 0.31546140237781 * 1.0e+2;
    coeff_res_n_[54] = -0.25213154341695 * 1.0e+4;
    coeff_res_n_[55] = -0.14874640856724;
    coeff_res_n_[56] = 0.31806110878444;

    // coefficient res alpha_{i + 52}
    coeff_res_alpha_[0] = 20.;
    coeff_res_alpha_[1] = 20.;
    coeff_res_alpha_[2] = 20.;

    // coefficient res gamma_{i + 52}
    coeff_res_gamma_[0] = 1.21;
    coeff_res_gamma_[1] = 1.21;
    coeff_res_gamma_[2] = 1.25;

    // coefficient res epsilon_{i + 52}
    coeff_res_epsilon_[0] = 1.;
    coeff_res_epsilon_[1] = 1.;
    coeff_res_epsilon_[2] = 1.;

    // coefficient res beta_{i + 52}
    coeff_res_beta_[0] = 150.;
    coeff_res_beta_[1] = 150.;
    coeff_res_beta_[2] = 250.;
    coeff_res_beta_[3] = 0.3;
    coeff_res_beta_[4] = 0.3;

    // coefficient res C_{i + 55}
    coeff_res_C_[0] = 28.;
    coeff_res_C_[1] = 32.;

    // coefficient res D_{i + 55}
    coeff_res_D_[0] = 700.;
    coeff_res_D_[1] = 800.;

    // coefficient res A_{i + 55}
    coeff_res_A_[0] = 0.32;
    coeff_res_A_[1] = 0.32;

    return;
}

} // end namespace IAPWS
