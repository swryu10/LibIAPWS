#ifndef _INTERCSPLINE_H_
#define _INTERCSPLINE_H_

/* implementation of
 * cubic spline interpolation */
class InterCSpline {
  private :

    int nbin_;
    double *tab_x_;
    double *tab_y_;

    bool initialized_;

  public :

    InterCSpline() {
        initialized_ = false;

        return;
    }

    ~InterCSpline() {
        if (!initialized_) {
            return;
        }

        delete [] tab_x_;
        delete [] tab_y_;

        return;
    }

    void init(int nbin_in,
              double *x_in,
              double *y_in);
};

#endif
