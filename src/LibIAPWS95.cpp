#include<stdio.h>
#include<math.h>
#include"LibIAPWS95.h"

double LibIAPWS95::get_param_phi_ide(double mdensity_in,
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

double LibIAPWS95::get_param_dphi_ide_ddelta(double mdensity_in,
                                             double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double dphi_ddelta = 1. / delta;

    return dphi_ddelta;
}

double LibIAPWS95::get_param_dphi_ide_dtau(double mdensity_in,
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

double LibIAPWS95::get_param_d2phi_ide_ddelta_ddelta(double mdensity_in,
                                                     double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_ddelta = -1. / (delta * delta);

    return d2phi_ddelta_ddelta;
}

double LibIAPWS95::get_param_d2phi_ide_ddelta_dtau(double mdensity_in,
                                                   double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double d2phi_ddelta_dtau = 0.;

    return d2phi_ddelta_dtau;
}

double LibIAPWS95::get_param_d2phi_ide_dtau_dtau(double mdensity_in,
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

double LibIAPWS95::get_param_phi_res(double mdensity_in,
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

double LibIAPWS95::get_param_dphi_res_ddelta(double mdensity_in,
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

double LibIAPWS95::get_param_dphi_res_dtau(double mdensity_in,
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

double LibIAPWS95::get_param_d2phi_res_ddelta_ddelta(double mdensity_in,
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

double LibIAPWS95::get_param_d2phi_res_ddelta_dtau(double mdensity_in,
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

double LibIAPWS95::get_param_d2phi_res_dtau_dtau(double mdensity_in,
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

double LibIAPWS95::get_param_pressure(double mdensity_in,
                                      double temperature_in) {
    double delta = mdensity_in / mdensity_crit_;
    //double tau = temperature_crit_ / temperature_in;

    double press =
        mdensity_in * const_R_spec_ * temperature_in *
        (1. + delta * get_param_dphi_res_ddelta(mdensity_in,
                                                temperature_in));

    return press;
}

double LibIAPWS95::get_param_dpress_drho(double mdensity_in,
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

double LibIAPWS95::get_param_erg_int(double mdensity_in,
                                     double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double erg_int =
        const_R_spec_ * temperature_in * tau *
        (get_param_dphi_ide_dtau(mdensity_in, temperature_in) +
         get_param_dphi_res_dtau(mdensity_in, temperature_in));

    return erg_int;
}

double LibIAPWS95::get_param_entropy(double mdensity_in,
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

double LibIAPWS95::get_param_enthalpy(double mdensity_in,
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

double LibIAPWS95::get_param_heat_c_v(double mdensity_in,
                                      double temperature_in) {
    //double delta = mdensity_in / mdensity_crit_;
    double tau = temperature_crit_ / temperature_in;

    double c_v =
        -const_R_spec_ * tau * tau *
        (get_param_d2phi_ide_dtau_dtau(mdensity_in, temperature_in) +
         get_param_d2phi_res_dtau_dtau(mdensity_in, temperature_in));

    return c_v;
}

double LibIAPWS95::get_param_heat_c_p(double mdensity_in,
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

double LibIAPWS95::get_param_speed_sound(double mdensity_in,
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
