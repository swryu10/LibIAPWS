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

    /* specific Gibbs free energy in J / kg
     * parametrization for single-phase state */
    double get_param_g(double temperature_in,
                       double pressure_in);
    /* specific volume in m^3 / kg
     * parametrization for single-phase state */
    double get_param_vol_spec(double temperature_in,
                              double pressure_in);
    /* mass density in kg / m^3
     * parametrization for single-phase state */
    double get_param_mdensity(double temperature_in,
                              double pressure_in);
    /* specific internal energy in J / kg
     * parametrization for single-phase state */
    double get_param_erg_int(double temperature_in,
                             double pressure_in);
    /* specific entropy in J / kg / degK
     * parametrization for single-phase state */
    double get_param_entropy(double temperature_in,
                             double pressure_in);
    /* specific enthalpy in J / kg
     * parametrization for single-phase state */
    double get_param_enthalpy(double temperature_in,
                              double pressure_in);
    /* specific isobaric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_p(double temperature_in,
                              double pressure_in);
    /* specific isochoric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_v(double temperature_in,
                              double pressure_in);
    /* speed of sound in m / sec
     * parametrization for single-phase state */
    double get_param_speed_sound(double temperature_in,
                                 double pressure_in);

    int get_region(double temperature_in,
                   double pressure_in);

    /* parametrized specific Gibbs free energy
     * and its derivatives in region 1 */
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
