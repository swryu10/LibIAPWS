#ifndef _LIBIAPWS06_H_
#define _LIBIAPWS06_H_

#include<stdio.h>
#include"LibIAPWS95.h"
#include"InterCSpline.h"

namespace IAPWS {

/* implementation of IAPWS R10-06 (2009)
 * Revised Release on the Equation of State 2006
 * for H2O Ice Ih */
class Lib06 {
  private :

    double temperature_trip_;
    double temperature_melt_;

    double pressure_trip_;
    double pressure_trip_num_;
    double pressure_norm_;

    double entropy_res_absol_;
    double entropy_res_Lib95_;

    double *coeff_g0_;

    double *coeff_r1_;

    double **coeff_r2_;

    double **coeff_t_;

    bool flag_s_absolute_;

    bool have_tab_coex_;
    int nbin_coex_;
    double temperature_coex_min_;
    double temperature_coex_max_;
    double *tab_coex_temperature_;
    double *tab_coex_pressure_;
    double *tab_coex_mden_vap_;
    double *tab_coex_mden_ice_;
    double *tab_coex_enthalpy_vap_;
    double *tab_coex_enthalpy_ice_;
    double *tab_coex_entropy_vap_;
    double *tab_coex_entropy_ice_;

    InterCSpline csp_coex_pressure_;
    InterCSpline csp_coex_mden_vap_;
    InterCSpline csp_coex_mden_ice_;
    InterCSpline csp_coex_enthalpy_vap_;
    InterCSpline csp_coex_enthalpy_ice_;
    InterCSpline csp_coex_entropy_vap_;
    InterCSpline csp_coex_entropy_ice_;

  public :

    Lib06() {
        temperature_trip_ = 273.16;
        temperature_melt_ = 273.152519;

        pressure_trip_ = 611.657;
        pressure_trip_num_ = 611.654771007894;
        pressure_norm_ = 101325.;

        entropy_res_absol_ = 0.18913 * 1.0e+3;
        entropy_res_Lib95_ = -0.332733756492168 * 1.0e+4;

        coeff_g0_ = new double[5];
        coeff_g0_[0] = -0.632020233335886 * 1.0e+6;
        coeff_g0_[1] = 0.655022213658955;
        coeff_g0_[2] = -0.189369929326131 * 1.0e-7;
        coeff_g0_[3] = 0.339746123271053 * 1.0e-14;
        coeff_g0_[4] = -0.556464869058991 * 1.0e-21;

        coeff_r1_ = new double[2];
        coeff_r1_[0] = 0.447050716285388 * 1.0e+2;
        coeff_r1_[1] = 0.656876847463481 * 1.0e+2;

        coeff_r2_ = new double *[3];
        for (int k = 0; k <= 2; k++) {
            coeff_r2_[k] = new double[2];
        }

        coeff_r2_[0][0] = -0.725974574329220 * 1.0e+2;
        coeff_r2_[0][1] = -0.781008427112870 * 1.0e+2;

        coeff_r2_[1][0] = -0.557107698030123 * 1.0e-4;
        coeff_r2_[1][1] = 0.464578634580806 * 1.0e-4;

        coeff_r2_[2][0] = 0.234801409215913 * 1.0e-10;
        coeff_r2_[2][1] = -0.285651142904972 * 1.0e-10;

        coeff_t_ = new double *[3];
        for (int k = 0; k <= 2; k++) {
            coeff_t_[k] = new double[2];
        }

        coeff_t_[0][0] = 0.;
        coeff_t_[0][1] = 0.;

        coeff_t_[1][0] = 0.368017112855051 * 1.0e-1;
        coeff_t_[1][1] = 0.510878114959572 * 1.0e-1;

        coeff_t_[2][0] = 0.337315741065416;
        coeff_t_[2][1] = 0.335449415919309;

        flag_s_absolute_ = false;

        have_tab_coex_ = false;

        return;
    }

    ~Lib06() {
        reset_tab_coex();

        delete [] coeff_g0_;

        delete [] coeff_r1_;

        for (int k = 0; k <= 2; k++) {
            delete [] coeff_r2_[k];
        }
        delete [] coeff_r2_;

        for (int k = 0; k <= 2; k++) {
            delete [] coeff_t_[k];
        }
        delete [] coeff_t_;

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
        delete [] tab_coex_mden_ice_;
        delete [] tab_coex_enthalpy_vap_;
        delete [] tab_coex_enthalpy_ice_;
        delete [] tab_coex_entropy_vap_;
        delete [] tab_coex_entropy_ice_;

        have_tab_coex_ = false;

        return;
    }

    void print_header(FILE *ptr_fout = stdout);

    void set_flag_s_absolute(bool flag_in) {
        flag_s_absolute_ = flag_in;
    }

    double get_param_g0(double pressure);
    double get_param_dg0_dp(double pressure);
    double get_param_d2g0_dp_dp(double pressure);

    void get_param_r(double pressure, int k,
                     double *r_out);
    void get_param_dr_dp(double pressure, int k,
                         double *dr_dp_out);
    void get_param_d2r_dp_dp(double pressure, int k,
                             double *d2r_dp_dp_out);

    void get_param_q(double temperature, int k,
                     double *q_out);
    void get_param_dq_dT(double temperature, int k,
                         double *dq_dT_out);
    void get_param_d2q_dT_dT(double temperature, int k,
                             double *d2q_dT_dT_out);

    /* parametrized specific Gibbs free energy
     * and its derivatives */
    double get_param_g(double temperature,
                       double pressure);
    double get_param_dg_dT(double temperature,
                           double pressure);
    double get_param_dg_dp(double temperature,
                           double pressure);
    double get_param_d2g_dT_dT(double temperature,
                               double pressure);
    double get_param_d2g_dT_dp(double temperature,
                               double pressure);
    double get_param_d2g_dp_dp(double temperature,
                               double pressure);

    /* mass density in kg / m^3
     * parametrization for single-phase state */
    double get_param_mdensity(double temperature_in,
                              double pressure_in);
    /* specific entropy in J / kg / degK
     * parametrization for single-phase state */
    double get_param_entropy(double temperature_in,
                             double pressure_in);
    /* specific isobaric heat capacity in J / kg / degK
     * parametrization for single-phase state */
    double get_param_heat_c_p(double temperature_in,
                              double pressure_in);
    /* specific enthalpy in J / kg
     * parametrization for single-phase state */
    double get_param_enthalpy(double temperature_in,
                              double pressure_in);
    /* specific internal energy in J / kg
     * parametrization for single-phase state */
    double get_param_erg_int(double temperature_in,
                             double pressure_in);
    /* specific Helmholtz free energy in J / kg
     * parametrization for single-phase state */
    double get_param_f(double temperature_in,
                       double pressure_in);
    /* cubic expansion coefficient in 1 / degK
     * parametrization for single-phase state */
    double get_param_coeff_alpha(double temperature_in,
                                 double pressure_in);
    /* pressure coefficient in Pa / degK
     * parametrization for single-phase state */
    double get_param_coeff_beta(double temperature_in,
                                double pressure_in);
    /* isothermal compressibility in 1 / Pa
     * parametrization for single-phase state */
    double get_param_comp_kappa_T(double temperature_in,
                                  double pressure_in);
    /* isentropic compressibility in 1 / Pa
     * parametrization for single-phase state */
    double get_param_comp_kappa_s(double temperature_in,
                                  double pressure_in);

    /* find a coexisting phase
     * based on continuity of Gibbs free energy */
    bool find_state_coex(Lib95 *ptr_lib95eos,
                         double temperature_in,
                         double &pressure_out,
                         double &mden_vap_out,
                         double &mden_ice_out);

    /* populate table
     * for coexisting phases */
    void make_tab_coex(Lib95 *ptr_lib95eos,
                       int nbin_in,
                       double temperature_min);

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
    /* mass density of H2O ice at coexisting phase
     * in kg / m^3 */
    double get_coex_mden_ice(double temperature_in) {
        return csp_coex_mden_ice_.get_func(temperature_in);
    }
    /* specific enthalpy of water vapor at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_vap(double temperature_in) {
        return csp_coex_enthalpy_vap_.get_func(temperature_in);
    }
    /* specific enthalpy of H2O ice at coexisting phase
     * in J / kg */
    double get_coex_enthalpy_ice(double temperature_in) {
        return csp_coex_enthalpy_ice_.get_func(temperature_in);
    }
    /* specific entropy of water vapor at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_vap(double temperature_in) {
        return csp_coex_entropy_vap_.get_func(temperature_in);
    }
    /* specific entropy of H2O ice at coexisting phase
     * in J / kg / degK */
    double get_coex_entropy_ice(double temperature_in) {
        return csp_coex_entropy_ice_.get_func(temperature_in);
    }
    /* specific latent heat
     * in J / kg */
    double get_coex_heat_latent(double temperature_in) {
        double h_latent =
            temperature_in *
            (csp_coex_entropy_vap_.get_func(temperature_in) -
             csp_coex_entropy_ice_.get_func(temperature_in));

        return h_latent;
    }
};

} // end namespace IAPWS

#endif
