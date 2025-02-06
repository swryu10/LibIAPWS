#include"InterCSpline.h"

void InterCSpline::init(int nbin_in,
                        double *x_in,
                        double *y_in) {
    reset();

    if (nbin_in < 1) {
        return;
    }

    nbin_ = nbin_in;
    tab_x_ = new double[nbin_ + 1];
    tab_y_ = new double[nbin_ + 1];
    tab_d2y_dx_ = new double[nbin_ + 1];

    for (int ix = 0; ix <= nbin_; ix++) {
        tab_x_[ix] = x_in[ix];
        tab_y_[ix] = y_in[ix];
    }

    bool sorted = false;
    while (!sorted) {
        bool swapped = false;
        for (int ix = 0; ix < nbin_; ix++) {
            if (tab_x_[ix] > tab_x_[ix + 1]) {
                double x_l = tab_x_[ix + 1];
                double y_l = tab_y_[ix + 1];

                double x_u = tab_x_[ix];
                double y_u = tab_y_[ix];

                tab_x_[ix] = x_l;
                tab_y_[ix] = y_l;

                tab_x_[ix + 1] = x_u;
                tab_y_[ix + 1] = y_u;

                swapped = true;
                break;
            }
        }

        if (!swapped) {
            sorted = true;
        }
    }

    xmin_ = tab_x_[0];
    xmax_ = tab_x_[nbin_];

    double *coeff_a = new double[nbin_ - 1];
    double *coeff_b = new double[nbin_ - 1];
    double *coeff_c = new double[nbin_ - 1];
    double *coeff_d = new double[nbin_ - 1];

    for (int ix = 1; ix < nbin_; ix++) {
        coeff_a[ix - 1] = (tab_x_[ix] - tab_x_[ix - 1]) / 6.;
        coeff_b[ix - 1] = (tab_x_[ix + 1] - tab_x_[ix - 1]) / 3.;
        coeff_c[ix - 1] = (tab_x_[ix + 1] - tab_x_[ix]) / 6.;
        coeff_d[ix - 1] =
            (tab_y_[ix + 1] - tab_y_[ix]) / (tab_x_[ix + 1] - tab_x_[ix]) -
            (tab_y_[ix] - tab_y_[ix - 1]) / (tab_x_[ix] - tab_x_[ix - 1]);
    }

    double *coeff_ap = new double[nbin_ - 1];
    double *coeff_bp = new double[nbin_ - 1];
    double *coeff_cp = new double[nbin_ - 1];
    double *coeff_dp = new double[nbin_ - 1];

    coeff_ap[0] = 0.;
    coeff_bp[0] = 1.;
    coeff_cp[0] = coeff_c[0] / coeff_b[0];
    coeff_dp[0] = coeff_d[0] / coeff_b[0];
    for (int ix = 1; ix < nbin_ - 1; ix++) {
        coeff_ap[ix] = 0.;
        coeff_bp[ix] = 1.;
        coeff_cp[ix] =
            coeff_c[ix] /
            (coeff_b[ix] - coeff_cp[ix - 1] * coeff_a[ix]);
        coeff_dp[ix] =
            (coeff_d[ix] - coeff_dp[ix - 1] * coeff_a[ix]) /
            (coeff_b[ix] - coeff_cp[ix - 1] * coeff_a[ix]);
    }

    tab_d2y_dx_[nbin_] = 0.;
    tab_d2y_dx_[nbin_ - 1] = coeff_dp[nbin_ - 2];
    for (int ix = nbin_ - 2; ix > 0; ix--) {
        tab_d2y_dx_[ix] =
            coeff_dp[ix - 1] - coeff_cp[ix - 1] * tab_d2y_dx_[ix + 1];
    }
    tab_d2y_dx_[0] = 0.;

    initialized_ = true;

    return;
}

double InterCSpline::get_func(double x_in,
                              double *ptr_dy_dx_out,
                              double *ptr_d2y_dx_dx_out) {
    if (!initialized_ ||
        x_in < xmin_ || x_in > xmax_) {
        if (ptr_dy_dx_out != NULL) {
            *ptr_dy_dx_out = 0.;
        }

        if (ptr_d2y_dx_dx_out != NULL) {
            *ptr_d2y_dx_dx_out = 0.;
        }

        return 0.;
    }

    int ix = 0;
    for (int jx = 0; jx < nbin_; jx++) {
        if (x_in >= tab_x_[jx] && x_in < tab_x_[jx + 1]) {
            ix = jx;
            break;
        }
    }

    double delta_x = tab_x_[ix + 1] - tab_x_[ix];
    double delta_y = tab_y_[ix + 1] - tab_y_[ix];

    double *frac_d0y = new double[2];
    frac_d0y[0] =
        (tab_x_[ix + 1] - x_in) / delta_x;
    frac_d0y[1] = 1. - frac_d0y[0];

    double *frac_d2y = new double[2];
    for (int k = 0; k < 2; k++) {
        frac_d2y[k] =
            frac_d0y[k] *
            (frac_d0y[k] * frac_d0y[k] - 1.) *
            delta_x * delta_x / 6.;
    }

    double y_out = 0.;
    for (int k = 0; k < 2; k++) {
        y_out +=
            frac_d0y[k] * tab_y_[ix + k] +
            frac_d2y[k] * tab_d2y_dx_[ix + k];
    }

    if (ptr_dy_dx_out != NULL) {
        *ptr_dy_dx_out =
            delta_y / delta_x -
            (3. * frac_d0y[0] * frac_d0y[0] - 1.) *
                delta_x * tab_d2y_dx_[ix] / 6. +
            (3. * frac_d0y[1] * frac_d0y[1] - 1.) *
                delta_x * tab_d2y_dx_[ix + 1] / 6.;
    }

    if (ptr_d2y_dx_dx_out != NULL) {
        *ptr_d2y_dx_dx_out =
            frac_d0y[0] * tab_d2y_dx_[ix] +
            frac_d0y[1] * tab_d2y_dx_[ix + 1];
    }

    return y_out;
}
