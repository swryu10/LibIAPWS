#include"InterCSpline.h"

void InterCSpline::init(int nbin_in,
                        double *x_in,
                        double *y_in) {
    if (nbin_in < 1) {
        initialized_ = false;
        return;
    }

    nbin_ = nbin_in;
    tab_x_ = new double[nbin_ + 1];
    tab_y_ = new double[nbin_ + 1];

    for (int ix = 0; ix <= nbin_; ix++) {
        tab_x_[ix] = x_in[ix];
        tab_y_[ix] = y_in[ix];
    }

    initialized_ = true;
}
