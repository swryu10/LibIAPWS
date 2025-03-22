#ifndef _LIBIAPWS95_H_
#define _LIBIAPWS95_H_

#include<stdio.h>
#include"BaseIAPWS.h"
#include"InterCSpline.h"

namespace IAPWS {

/* implementation of IAPWS R6-95 (2018)
 * Revised Release on the IAPWS Formulation 1995
 * for the Thermodynamic Properties of Ordinary Water Substance
 * for General and Scientific Use */
class Lib95 {
  private :

    /* critical temperature
     * in degK */
    double temperature_crit_;
    /* triple-point temperature
     * in degK */
    double temperature_trip_;
    /* critical (mass) density
     * in kg / m^3 */
    double mdensity_crit_;
    /* critical pressure
     * in Pa (N / m^2) */
    double pressure_crit_;
    /* specific gas constant
     * in J / kg / degK */
    double const_R_spec_;

    double *coeff_ide_n_;
    double *coeff_ide_gamma_;

    double *coeff_res_c_;
    double *coeff_res_d_;
    double *coeff_res_t_;

    double *coeff_res_n_;

    double *coeff_res_alpha_;
    double *coeff_res_gamma_;
    double *coeff_res_epsilon_;

    double *coeff_res_beta_;

    double *coeff_res_C_;
    double *coeff_res_D_;
    double *coeff_res_A_;

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

    Lib95() {
        coeff_ide_n_ = new double[9];
        coeff_ide_gamma_ = new double[9];

        coeff_res_c_ = new double[57];
        coeff_res_d_ = new double[57];
        coeff_res_t_ = new double[57];
        coeff_res_n_ = new double[57];

        coeff_res_alpha_ = new double[3];
        coeff_res_gamma_ = new double[3];
        coeff_res_epsilon_ = new double[3];

        coeff_res_beta_ = new double[5];

        coeff_res_C_ = new double[2];
        coeff_res_D_ = new double[2];
        coeff_res_A_ = new double[2];

        set_coefficients();

        have_tab_coex_ = false;

        return;
    }

    ~Lib95() {
        reset_tab_coex();

        delete [] coeff_ide_n_;
        delete [] coeff_ide_gamma_;

        delete [] coeff_res_c_;
        delete [] coeff_res_d_;
        delete [] coeff_res_t_;

        delete [] coeff_res_n_;

        delete [] coeff_res_alpha_;
        delete [] coeff_res_gamma_;
        delete [] coeff_res_epsilon_;

        delete [] coeff_res_beta_;

        delete [] coeff_res_C_;
        delete [] coeff_res_D_;
        delete [] coeff_res_A_;

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

    /* the ideal-gas part of the free energy
     * and its derivatives */
    double get_param_phi_ide(double mdensity_in,
                             double temperature_in);
    double get_param_dphi_ide_ddelta(double mdensity_in,
                                     double temperature_in);
    double get_param_dphi_ide_dtau(double mdensity_in,
                                   double temperature_in);
    double get_param_d2phi_ide_ddelta_ddelta(double mdensity_in,
                                             double temperature_in);
    double get_param_d2phi_ide_ddelta_dtau(double mdensity_in,
                                           double temperature_in);
    double get_param_d2phi_ide_dtau_dtau(double mdensity_in,
                                         double temperature_in);

    /* the residual part of the free energy
     * and its derivatives */
    double get_param_phi_res(double mdensity_in,
                             double temperature_in);
    double get_param_dphi_res_ddelta(double mdensity_in,
                                     double temperature_in);
    double get_param_dphi_res_dtau(double mdensity_in,
                                   double temperature_in);
    double get_param_d2phi_res_ddelta_ddelta(double mdensity_in,
                                             double temperature_in);
    double get_param_d2phi_res_ddelta_dtau(double mdensity_in,
                                           double temperature_in);
    double get_param_d2phi_res_dtau_dtau(double mdensity_in,
                                         double temperature_in);

    /* specific Helmholtz free energy in J / kg
     * parametrization for single-phase state */
    double get_param_f(double mdensity_in,
                       double temperature_in);
    /* specific Gibbs free energy in J / kg
     * parametrization for single-phase state */
    double get_param_g(double mdensity_in,
                       double temperature_in);
    /* pressure in Pa (N / m^2)
     * parametrization for single-phase state */
    double get_param_pressure(double mdensity_in,
                              double temperature_in);
    /* derivative of pressure with respect to mass density
     * in m^2 / sec^2
     * parametrization for single-phase state */
    double get_param_dpress_drho(double mdensity_in,
                                 double temperature_in);
    /* specific internal energy in J / kg
     * parametrization for single-phase state */
    double get_param_erg_int(double mdensity_in,
                             double temperature_in);
    /* specific entropy in J / kg / degK
     * parametrization for single-phase state */
    double get_param_entropy(double mdensity_in,
                             double temperature_in);
    /* specific enthalpy in J / kg
     * parametrization for single-phase state */
    double get_param_enthalpy(double mdensity_in,
                              double temperature_in);
    /* specific isochoric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_v(double mdensity_in,
                              double temperature_in);
    /* specific isobaric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_p(double mdensity_in,
                              double temperature_in);
    /* speed of sound in m / sec
     * parametrization for single-phase state */
    double get_param_speed_sound(double mdensity_in,
                                 double temperature_in);

    /* find mass density which yields certain pressure
     * with given temperature,
     * using Newton's method for root finding */
    bool find_root_mdensity(double temperature_in,
                            double pressure_in,
                            double &mdensity_out);

    /* find a coexisting phase
     * based on Maxwell criterion */
    bool find_state_coex(double temperature_in,
                         double &pressure_out,
                         double &mden_vap_out,
                         double &mden_liq_out);

    /* populate table
     * for coexisting phases */
    void make_tab_coex(int nbin_in,
                       double temperature_max);

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

    /* pressure at coexisting phase
     * in Pa (N / m^2) */
    double get_coex_pressure(double temperature_in) {
        return csp_coex_pressure_.get_func(temperature_in);
    }
    /* mass density of water vapor at coexisting phase
     * in kg / m^3 */
    double get_coex_mden_vap(double temperature_in) {
        return csp_coex_mden_vap_.get_func(temperature_in);
    }
    /* mass density of water liquid at coexisting phase
     * in kg / m^3 */
    double get_coex_mden_liq(double temperature_in) {
        return csp_coex_mden_liq_.get_func(temperature_in);
    }
    /* specific enthalpy of water vapor at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_vap(double temperature_in) {
        return csp_coex_enthalpy_vap_.get_func(temperature_in);
    }
    /* specific enthalpy of water liquid at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_liq(double temperature_in) {
        return csp_coex_enthalpy_liq_.get_func(temperature_in);
    }
    /* specific entropy of water vapor at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_vap(double temperature_in) {
        return csp_coex_entropy_vap_.get_func(temperature_in);
    }
    /* specific entropy of water liquid at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_liq(double temperature_in) {
        return csp_coex_entropy_liq_.get_func(temperature_in);
    }
    /* specific latent heat
     * in J / kg */
    double get_coex_heat_latent(double temperature_in) {
        double h_latent =
            temperature_in *
            (csp_coex_entropy_vap_.get_func(temperature_in) -
             csp_coex_entropy_liq_.get_func(temperature_in));

        return h_latent;
    }
};

} // end namespace IAPWS

#endif
