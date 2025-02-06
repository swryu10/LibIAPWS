#ifndef _INTERCSPLINE_H_
#define _INTERCSPLINE_H_

#include<stdlib.h>

/* implementation of
 * natural cubic spline interpolation */
class InterCSpline {
  private :

    int nbin_;
    double xmin_;
    double xmax_;
    double *tab_x_;
    double *tab_y_;
    double *tab_d2y_dx_;

    bool initialized_;

  public :

    InterCSpline() {
        initialized_ = false;

        return;
    }

    ~InterCSpline() {
        reset();

        return;
    }

    void reset() {
        if (!initialized_) {
            return;
        }

        delete [] tab_x_;
        delete [] tab_y_;
        delete [] tab_d2y_dx_;

        initialized_ = false;

        return;
    }

    void init(int nbin_in,
              double *x_in,
              double *y_in);

    double get_func(double x_in,
                    double *ptr_dy_dx_out = NULL,
                    double *ptr_d2y_dx_dx_out = NULL);
};

#endif
