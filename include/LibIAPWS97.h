#ifndef _LIBIAPWS97_H_
#define _LIBIAPWS97_H_

#include<stdio.h>
#include"BaseIAPWS.h"

namespace IAPWS {

/* implementation of IAPWS R7-97(2012)
 * Revised Release on the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam */
class Lib97 {
  private :

    /* the critical temperature
     * in degK */
    double temperature_crit_;
    /* the critical mass density
     * in kg / m^3 */
    double mdensity_crit_;
    /* the critical pressure
     * in Pa */
    double pressure_crit_;
    /* specific gas constant
     * in J / kg / degK */
    double const_R_spec_;

    double temperature_ref1_;
    double pressure_ref1_;
    double tau_ref1_;
    double ppi_ref1_;

    int *coeff1_I_;
    int *coeff1_J_;
    double *coeff1_n_;

    void set_coefficients();

  public :

    Lib97() {
        temperature_crit_ = 647.096;
        mdensity_crit_ = 322.;
        pressure_crit_ = 22.064 * 1.0e+6;
        const_R_spec_ = 461.526;

        coeff1_I_ = new int[35];
        coeff1_J_ = new int[35];
        coeff1_n_ = new double[35];

        set_coefficients();

        return;
    }

    ~Lib97() {
        delete [] coeff1_I_;
        delete [] coeff1_J_;
        delete [] coeff1_n_;

        return;
    }

    void print_header(FILE *ptr_fout = stdout);

    /* parametrized specific Gibbs free energy
     * and its derivatives */
    double get_param1_gamma(double temperature_in,
                            double pressure_in);
    double get_param1_dgamma_dppi(double temperature_in,
                                  double pressure_in);
    double get_param1_dgamma_dtau(double temperature_in,
                                  double pressure_in);
    double get_param1_d2gamma_dppi_dppi(double temperature_in,
                                        double pressure_in);
    double get_param1_d2gamma_dppi_dtau(double temperature_in,
                                        double pressure_in);
    double get_param1_d2gamma_dtau_dtau(double temperature_in,
                                        double pressure_in);
};

} // end namespace IAPWS

#endif
