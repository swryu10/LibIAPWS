#ifndef _BASEIAPWS_H_
#define _BASEIAPWS_H_

namespace IAPWS {

extern bool flag_verbose_;

extern int n_iter_max_;
extern double eps_precision_;

void get_complex_product(double *z_1in,
                         double *z_2in,
                         double *z_out);

void get_complex_inverse(double *z_in,
                         double *z_out);

void get_complex_log(double *z_in,
                     double *z_out);

} // end namespace IAPWS

#endif
