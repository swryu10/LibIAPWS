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

    double temperature_refB23_;
    double pressure_refB23_;

    double *coeffB23_n_;

    double temperature_ref1_;
    double pressure_ref1_;
    double tau_ref1_;
    double ppi_ref1_;

    int *coeff1_I_;
    int *coeff1_J_;
    double *coeff1_n_;

    double temperature_ref1Tph_;
    double pressure_ref1Tph_;
    double enthalpy_ref1Tph_;

    int *coeff1Tph_I_;
    int *coeff1Tph_J_;
    double *coeff1Tph_n_;

    double temperature_ref1Tps_;
    double pressure_ref1Tps_;
    double entropy_ref1Tps_;

    int *coeff1Tps_I_;
    int *coeff1Tps_J_;
    double *coeff1Tps_n_;

    double temperature_ref2_;
    double pressure_ref2_;
    double tau_ref2_res_;

    int *coeff2_ide_J_;
    double *coeff2_ide_n_;
    double coeff2mst_ide_n1_;
    double coeff2mst_ide_n2_;

    int *coeff2_res_I_;
    int *coeff2_res_J_;
    double *coeff2_res_n_;

    double temperature_ref2mst_;
    double pressure_ref2mst_;
    double tau_ref2mst_res_;

    int *coeff2mst_res_I_;
    int *coeff2mst_res_J_;
    double *coeff2mst_res_n_;

    double pressure_refB2bc_;
    double enthalpy_refB2bc_;

    double *coeffB2bc_n_;

    double temperature_ref4_;
    double pressure_ref4_;

    double *coeff4_n_;

    void set_coefficients();

  public :

    Lib97() {
        coeffB23_n_ = new double[6];

        coeff1_I_ = new int[35];
        coeff1_J_ = new int[35];
        coeff1_n_ = new double[35];

        coeff1Tph_I_ = new int[21];
        coeff1Tph_J_ = new int[21];
        coeff1Tph_n_ = new double[21];

        coeff1Tps_I_ = new int[21];
        coeff1Tps_J_ = new int[21];
        coeff1Tps_n_ = new double[21];

        coeff2_ide_J_ = new int[10];
        coeff2_ide_n_ = new double[10];

        coeff2_res_I_ = new int[44];
        coeff2_res_J_ = new int[44];
        coeff2_res_n_ = new double[44];

        coeffB2bc_n_ = new double[6];

        coeff2mst_res_I_ = new int[14];
        coeff2mst_res_J_ = new int[14];
        coeff2mst_res_n_ = new double[14];

        coeff4_n_ = new double[11];

        set_coefficients();

        return;
    }

    ~Lib97() {
        delete [] coeffB23_n_;

        delete [] coeff1_I_;
        delete [] coeff1_J_;
        delete [] coeff1_n_;

        delete [] coeff1Tph_I_;
        delete [] coeff1Tph_J_;
        delete [] coeff1Tph_n_;

        delete [] coeff1Tps_I_;
        delete [] coeff1Tps_J_;
        delete [] coeff1Tps_n_;

        delete [] coeff2_ide_J_;
        delete [] coeff2_ide_n_;

        delete [] coeff2_res_I_;
        delete [] coeff2_res_J_;
        delete [] coeff2_res_n_;

        delete [] coeff2mst_res_I_;
        delete [] coeff2mst_res_J_;
        delete [] coeff2mst_res_n_;

        delete [] coeffB2bc_n_;

        delete [] coeff4_n_;

        return;
    }

    void print_header(FILE *ptr_fout = stdout);

    /* specific Gibbs free energy in J / kg
     * parametrization for single-phase state */
    double get_param_g(double temperature_in,
                       double pressure_in,
                       bool flag_metastable = false);
    /* specific volume in m^3 / kg
     * parametrization for single-phase state */
    double get_param_vol_spec(double temperature_in,
                              double pressure_in,
                              bool flag_metastable = false);
    /* mass density in kg / m^3
     * parametrization for single-phase state */
    double get_param_mdensity(double temperature_in,
                              double pressure_in,
                              bool flag_metastable = false);
    /* specific internal energy in J / kg
     * parametrization for single-phase state */
    double get_param_erg_int(double temperature_in,
                             double pressure_in,
                             bool flag_metastable = false);
    /* specific entropy in J / kg / degK
     * parametrization for single-phase state */
    double get_param_entropy(double temperature_in,
                             double pressure_in,
                             bool flag_metastable = false);
    /* specific enthalpy in J / kg
     * parametrization for single-phase state */
    double get_param_enthalpy(double temperature_in,
                              double pressure_in,
                              bool flag_metastable = false);
    /* specific isobaric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_p(double temperature_in,
                              double pressure_in,
                              bool flag_metastable = false);
    /* specific isochoric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_v(double temperature_in,
                              double pressure_in,
                              bool flag_metastable = false);
    /* speed of sound in m / sec
     * parametrization for single-phase state */
    double get_param_speed_sound(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable = false);

    int get_region(double temperature_in,
                   double pressure_in,
                   bool flag_metastable = false);

    /* auxiliary equations for the boundary
     * between Regions 2 and 3 */
    double get_paramB23_pressure(double temperature_in);
    double get_paramB23_temperature(double pressure_in);

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

    /* backward equation for temperature
     * as function of pressure and specific enthalpy
     * in region 1 */
    double get_param1_temperature_ph(double pressure_in,
                                     double enthalpy_in);
    /* backward equation for temperature
     * as function of pressure and specific entropy
     * in region 1 */
    double get_param1_temperature_ps(double pressure_in,
                                     double entropy_in);

    /* parametrized ideal-gas part
     * of the specific Gibbs free energy
     * and its derivatives in region 2 */
    double get_param2_gamma_ide(double temperature_in,
                                double pressure_in,
                                bool flag_metastable = false);
    double get_param2_dgamma_ide_dppi(double temperature_in,
                                      double pressure_in,
                                      bool flag_metastable = false);
    double get_param2_dgamma_ide_dtau(double temperature_in,
                                      double pressure_in,
                                      bool flag_metastable = false);
    double get_param2_d2gamma_ide_dppi_dppi(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);
    double get_param2_d2gamma_ide_dppi_dtau(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);
    double get_param2_d2gamma_ide_dtau_dtau(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);

    /* parametrized residual part
     * of the specific Gibbs free energy
     * and its derivatives in region 2 */
    double get_param2_gamma_res(double temperature_in,
                                double pressure_in,
                                bool flag_metastable = false);
    double get_param2_dgamma_res_dppi(double temperature_in,
                                      double pressure_in,
                                      bool flag_metastable = false);
    double get_param2_dgamma_res_dtau(double temperature_in,
                                      double pressure_in,
                                      bool flag_metastable = false);
    double get_param2_d2gamma_res_dppi_dppi(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);
    double get_param2_d2gamma_res_dppi_dtau(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);
    double get_param2_d2gamma_res_dtau_dtau(double temperature_in,
                                            double pressure_in,
                                            bool flag_metastable = false);

    /* auxiliary equations for the boundary
     * between Regions 2b and 2c */
    double get_paramB2bc_pressure(double enthalpy_in);
    double get_paramB2bc_enthalpy(double pressure_in);

    /* parametrized vapor-liquid saturation pressure (in Pa)
     * as a function of saturation temperature (in degK) */
    double get_param4_sat_pressure(double temperature_in);
    /* parametrized vapor-liquid saturation temperature (in degK)
     * as a function of saturation pressure (in Pa) */
    double get_param4_sat_temperature(double pressure_in);
};

} // end namespace IAPWS

#endif
