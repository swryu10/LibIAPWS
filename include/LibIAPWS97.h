#ifndef _LIBIAPWS97_H_
#define _LIBIAPWS97_H_

#include<stdio.h>
#include"BaseIAPWS.h"
#include"InterCSpline.h"

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

    double temperature_low3_;
    double pressure_low3_;

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

    double temperature_ref2aTph_;
    double pressure_ref2aTph_;
    double enthalpy_ref2aTph_;

    int *coeff2aTph_I_;
    int *coeff2aTph_J_;
    double *coeff2aTph_n_;

    double temperature_ref2bTph_;
    double pressure_ref2bTph_;
    double enthalpy_ref2bTph_;

    int *coeff2bTph_I_;
    int *coeff2bTph_J_;
    double *coeff2bTph_n_;

    double temperature_ref2cTph_;
    double pressure_ref2cTph_;
    double enthalpy_ref2cTph_;

    int *coeff2cTph_I_;
    int *coeff2cTph_J_;
    double *coeff2cTph_n_;

    double temperature_ref2aTps_;
    double pressure_ref2aTps_;
    double entropy_ref2aTps_;

    double *coeff2aTps_I_;
    int *coeff2aTps_J_;
    double *coeff2aTps_n_;

    double temperature_ref2bTps_;
    double pressure_ref2bTps_;
    double entropy_ref2bTps_;

    int *coeff2bTps_I_;
    int *coeff2bTps_J_;
    double *coeff2bTps_n_;

    double temperature_ref2cTps_;
    double pressure_ref2cTps_;
    double entropy_ref2cTps_;

    int *coeff2cTps_I_;
    int *coeff2cTps_J_;
    double *coeff2cTps_n_;

    double mdensity_ref3_;
    double temperature_ref3_;

    int *coeff3_I_;
    int *coeff3_J_;
    double *coeff3_n_;

    double temperature_ref4_;
    double pressure_ref4_;

    double *coeff4_n_;

    double temperature_ref5_;
    double pressure_ref5_;

    int *coeff5_ide_J_;
    double *coeff5_ide_n_;

    int *coeff5_res_I_;
    int *coeff5_res_J_;
    double *coeff5_res_n_;

    bool have_tab_coex_;
    int nbin_coex_;
    double temperature_coex_min_;
    double temperature_coex_max_;
    double *tab_coex_temperature_;
    double *tab_coex_pressure_;
    double *tab_coex_mden_vap_;
    double *tab_coex_mden_liq_;
    double *tab_coex_enthalpy_vap_;
    double *tab_coex_enthalpy_liq_;
    double *tab_coex_entropy_vap_;
    double *tab_coex_entropy_liq_;

    InterCSpline csp_coex_pressure_;
    InterCSpline csp_coex_mden_vap_;
    InterCSpline csp_coex_mden_liq_;
    InterCSpline csp_coex_enthalpy_vap_;
    InterCSpline csp_coex_enthalpy_liq_;
    InterCSpline csp_coex_entropy_vap_;
    InterCSpline csp_coex_entropy_liq_;

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

        coeff2mst_res_I_ = new int[14];
        coeff2mst_res_J_ = new int[14];
        coeff2mst_res_n_ = new double[14];

        coeffB2bc_n_ = new double[6];

        coeff2aTph_I_ = new int[35];
        coeff2aTph_J_ = new int[35];
        coeff2aTph_n_ = new double[35];

        coeff2bTph_I_ = new int[39];
        coeff2bTph_J_ = new int[39];
        coeff2bTph_n_ = new double[39];

        coeff2cTph_I_ = new int[24];
        coeff2cTph_J_ = new int[24];
        coeff2cTph_n_ = new double[24];

        coeff2aTps_I_ = new double[47];
        coeff2aTps_J_ = new int[47];
        coeff2aTps_n_ = new double[47];

        coeff2bTps_I_ = new int[45];
        coeff2bTps_J_ = new int[45];
        coeff2bTps_n_ = new double[45];

        coeff2cTps_I_ = new int[31];
        coeff2cTps_J_ = new int[31];
        coeff2cTps_n_ = new double[31];

        coeff3_I_ = new int[41];
        coeff3_J_ = new int[41];
        coeff3_n_ = new double[41];

        coeff4_n_ = new double[11];

        coeff5_ide_J_ = new int[7];
        coeff5_ide_n_ = new double[7];

        coeff5_res_I_ = new int[7];
        coeff5_res_J_ = new int[7];
        coeff5_res_n_ = new double[7];

        set_coefficients();

        have_tab_coex_ = false;

        return;
    }

    ~Lib97() {
        reset_tab_coex();

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

        delete [] coeff2aTph_I_;
        delete [] coeff2aTph_J_;
        delete [] coeff2aTph_n_;

        delete [] coeff2bTph_I_;
        delete [] coeff2bTph_J_;
        delete [] coeff2bTph_n_;

        delete [] coeff2cTph_I_;
        delete [] coeff2cTph_J_;
        delete [] coeff2cTph_n_;

        delete [] coeff2aTps_I_;
        delete [] coeff2aTps_J_;
        delete [] coeff2aTps_n_;

        delete [] coeff2bTps_I_;
        delete [] coeff2bTps_J_;
        delete [] coeff2bTps_n_;

        delete [] coeff2cTps_I_;
        delete [] coeff2cTps_J_;
        delete [] coeff2cTps_n_;

        delete [] coeff3_I_;
        delete [] coeff3_J_;
        delete [] coeff3_n_;

        delete [] coeff4_n_;

        delete [] coeff5_ide_J_;
        delete [] coeff5_ide_n_;

        delete [] coeff5_res_I_;
        delete [] coeff5_res_J_;
        delete [] coeff5_res_n_;

        return;
    }

    void reset_tab_coex() {
        if (!have_tab_coex_) {
            return;
        }

        nbin_coex_ = 0;

        delete [] tab_coex_temperature_;
        delete [] tab_coex_pressure_;
        delete [] tab_coex_mden_vap_;
        delete [] tab_coex_mden_liq_;
        delete [] tab_coex_enthalpy_vap_;
        delete [] tab_coex_enthalpy_liq_;
        delete [] tab_coex_entropy_vap_;
        delete [] tab_coex_entropy_liq_;

        have_tab_coex_ = false;

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

    /* determine which region the system belongs
     * with given temperature and pressure */
    int get_region(double temperature_in,
                   double pressure_in,
                   bool flag_metastable = false);

    /* pressure at coexisting phase
     * in Pa (N / m^2) */
    double get_coex_pressure(double temperature_in);
    /* mass density of water vapor at coexisting phase
     * in kg / m^3 */
    double get_coex_mden_vap(double temperature_in);
    /* mass density of water liquid at coexisting phase
     * in kg / m^3 */
    double get_coex_mden_liq(double temperature_in);
    /* specific enthalpy of water vapor at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_vap(double temperature_in);
    /* specific enthalpy of water liquid at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_liq(double temperature_in);
    /* specific entropy of water vapor at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_vap(double temperature_in);
    /* specific entropy of water liquid at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_liq(double temperature_in);
    /* specific latent heat
     * in J / kg */
    double get_coex_heat_latent(double temperature_in);

    /* auxiliary equations for the boundary
     * between Regions 2 and 3 */
    double get_paramB23_pressure(double temperature_in);
    double get_paramB23_temperature(double pressure_in);


    /* mass density in kg / m^3
     * parametrization for region 1 */
    double get_param1_mdensity(double temperature_in,
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

    /* mass density in kg / m^3
     * parametrization for region 2 */
    double get_param2_mdensity(double temperature_in,
                               double pressure_in);

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

    /* backward equation for temperature
     * as function of pressure and specific enthalpy
     * in region 2 */
    double get_param2_temperature_ph(double pressure_in,
                                     double enthalpy_in);
    /* backward equation for temperature
     * as function of pressure and specific enthalpy
     * in region 2a */
    double get_param2a_temperature_ph(double pressure_in,
                                      double enthalpy_in);
    /* backward equation for temperature
     * as function of pressure and specific enthalpy
     * in region 2b */
    double get_param2b_temperature_ph(double pressure_in,
                                      double enthalpy_in);
    /* backward equation for temperature
     * as function of pressure and specific enthalpy
     * in region 2c */
    double get_param2c_temperature_ph(double pressure_in,
                                      double enthalpy_in);

    /* backward equation for temperature
     * as function of pressure and specific entropy
     * in region 2 */
    double get_param2_temperature_ps(double pressure_in,
                                     double entropy_in);
    /* backward equation for temperature
     * as function of pressure and specific entropy
     * in region 2a */
    double get_param2a_temperature_ps(double pressure_in,
                                      double entropy_in);
    /* backward equation for temperature
     * as function of pressure and specific entropy
     * in region 2b */
    double get_param2b_temperature_ps(double pressure_in,
                                      double entropy_in);
    /* backward equation for temperature
     * as function of pressure and specific entropy
     * in region 2c */
    double get_param2c_temperature_ps(double pressure_in,
                                      double entropy_in);

    /* specific Helmholtz free energy in J / kg
     * parametrization for region 3 */
    double get_param3_f(double mdensity_in,
                        double temperature_in);
    /* specific Gibbs free energy in J / kg
     * parametrization for region 3 */
    double get_param3_g(double mdensity_in,
                        double temperature_in);
    /* pressure in Pa (N / m^2)
     * parametrization for region 3 */
    double get_param3_pressure(double mdensity_in,
                               double temperature_in);
    /* derivative of pressure with respect to mass density
     * in m^2 / sec^2
     * parametrization for region 3 */
    double get_param3_dpress_drho(double mdensity_in,
                                  double temperature_in);
    /* specific internal energy in J / kg
     * parametrization for region 3 */
    double get_param3_erg_int(double mdensity_in,
                              double temperature_in);
    /* specific entropy in J / kg / degK
     * parametrization for region 3 */
    double get_param3_entropy(double mdensity_in,
                              double temperature_in);
    /* specific enthalpy in J / kg
     * parametrization for region 3 */
    double get_param3_enthalpy(double mdensity_in,
                               double temperature_in);
    /* specific isochoric heat capacity in J / kg / degK
     * parametrization for region 3 */
    double get_param3_heat_c_v(double mdensity_in,
                               double temperature_in);
    /* specific isobaric heat capacity in J / kg / degK
     * parametrization for region 3 */
    double get_param3_heat_c_p(double mdensity_in,
                               double temperature_in);
    /* speed of sound in m / sec
     * parametrization for region 3 */
    double get_param3_speed_sound(double mdensity_in,
                                  double temperature_in);

    /* parametrized specific Helmholtz free energy
     * and its derivatives in region 3 */
    double get_param3_phi(double mdensity_in,
                          double temperature_in);
    double get_param3_dphi_ddelta(double mdensity_in,
                                  double temperature_in);
    double get_param3_dphi_dtau(double mdensity_in,
                                double temperature_in);
    double get_param3_d2phi_ddelta_ddelta(double mdensity_in,
                                          double temperature_in);
    double get_param3_d2phi_ddelta_dtau(double mdensity_in,
                                        double temperature_in);
    double get_param3_d2phi_dtau_dtau(double mdensity_in,
                                      double temperature_in);

    /* parametrized vapor-liquid saturation pressure (in Pa)
     * as a function of saturation temperature (in degK) */
    double get_param4_sat_pressure(double temperature_in);
    /* parametrized vapor-liquid saturation temperature (in degK)
     * as a function of saturation pressure (in Pa) */
    double get_param4_sat_temperature(double pressure_in);

    /* parametrized ideal-gas part
     * of the specific Gibbs free energy
     * and its derivatives in region 5 */
    double get_param5_gamma_ide(double temperature_in,
                                double pressure_in);
    double get_param5_dgamma_ide_dppi(double temperature_in,
                                      double pressure_in);
    double get_param5_dgamma_ide_dtau(double temperature_in,
                                      double pressure_in);
    double get_param5_d2gamma_ide_dppi_dppi(double temperature_in,
                                            double pressure_in);
    double get_param5_d2gamma_ide_dppi_dtau(double temperature_in,
                                            double pressure_in);
    double get_param5_d2gamma_ide_dtau_dtau(double temperature_in,
                                            double pressure_in);

    /* parametrized residual part
     * of the specific Gibbs free energy
     * and its derivatives in region 5 */
    double get_param5_gamma_res(double temperature_in,
                                double pressure_in);
    double get_param5_dgamma_res_dppi(double temperature_in,
                                      double pressure_in);
    double get_param5_dgamma_res_dtau(double temperature_in,
                                      double pressure_in);
    double get_param5_d2gamma_res_dppi_dppi(double temperature_in,
                                            double pressure_in);
    double get_param5_d2gamma_res_dppi_dtau(double temperature_in,
                                            double pressure_in);
    double get_param5_d2gamma_res_dtau_dtau(double temperature_in,
                                            double pressure_in);

    /* find mass density which yields certain pressure
     * with given temperature,
     * using Newton's method for root finding */
    bool find_root3_mdensity(double temperature_in,
                             double pressure_in,
                             double &mdensity_out,
                             int sign_ini = 0);

    /* find a coexisting phase
     * based on Maxwell criterion */
    bool find_state3_coex(double temperature_in,
                          double &pressure_out,
                          double &mden_vap_out,
                          double &mden_liq_out);

    /* populate table
     * for coexisting phases */
    void make_tab_coex(int nbin_in);

    /* print out the table
     * for coexisting phases */
    void export_tab_coex(char *filename);

    /* import table for coexisting phases
     * from an external data file */
    void import_tab_coex(char *filename);

    /* initialize cubic spline interpolation
     * for thermodynamic quantities
     * at coexisting phases */
    void set_cspline_coex();
};

} // end namespace IAPWS

#endif
