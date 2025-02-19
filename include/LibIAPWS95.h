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
    /* temperatur of the triple point
     * in degK */
    double temperature_trip_;
    /* critical (mass) density
     * in kg / m^3 */
    double mdensity_crit_;
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

  public :

    Lib95() {
        temperature_crit_ = 647.096;
        temperature_trip_ = 273.16;
        mdensity_crit_ = 322.;
        const_R_spec_ = 461.51805;

        coeff_ide_n_ = new double[9];
        coeff_ide_n_[0] = 0.;
        // coefficient ide n_i
        coeff_ide_n_[1] = -8.3204464837497;
        coeff_ide_n_[2] = 6.6832105275932;
        coeff_ide_n_[3] = 3.00632;
        coeff_ide_n_[4] = 0.012436;
        coeff_ide_n_[5] = 0.97315;
        coeff_ide_n_[6] = 1.27950;
        coeff_ide_n_[7] = 0.96956;
        coeff_ide_n_[8] = 0.24873;

        coeff_ide_gamma_ = new double[9];
        coeff_ide_gamma_[0] = 0.;
        // coefficient ide gamma_i
        coeff_ide_gamma_[1] = 0.;
        coeff_ide_gamma_[2] = 0.;
        coeff_ide_gamma_[3] = 0.;
        coeff_ide_gamma_[4] = 1.28728967;
        coeff_ide_gamma_[5] = 3.53734222;
        coeff_ide_gamma_[6] = 7.74073708;
        coeff_ide_gamma_[7] = 9.24437796;
        coeff_ide_gamma_[8] = 27.5075105;

        coeff_res_c_ = new double[57];
        coeff_res_c_[0] = 0.;
        // coefficient res c_i
        coeff_res_c_[1] = 0.;
        coeff_res_c_[2] = 0.;
        coeff_res_c_[3] = 0.;
        coeff_res_c_[4] = 0.;
        coeff_res_c_[5] = 0.;
        coeff_res_c_[6] = 0.;
        coeff_res_c_[7] = 0.;
        coeff_res_c_[8] = 1.;
        coeff_res_c_[9] = 1.;
        coeff_res_c_[10] = 1.;
        coeff_res_c_[11] = 1.;
        coeff_res_c_[12] = 1.;
        coeff_res_c_[13] = 1.;
        coeff_res_c_[14] = 1.;
        coeff_res_c_[15] = 1.;
        coeff_res_c_[16] = 1.;
        coeff_res_c_[17] = 1.;
        coeff_res_c_[18] = 1.;
        coeff_res_c_[19] = 1.;
        coeff_res_c_[20] = 1.;
        coeff_res_c_[21] = 1.;
        coeff_res_c_[22] = 1.;
        coeff_res_c_[23] = 2.;
        coeff_res_c_[24] = 2.;
        coeff_res_c_[25] = 2.;
        coeff_res_c_[26] = 2.;
        coeff_res_c_[27] = 2.;
        coeff_res_c_[28] = 2.;
        coeff_res_c_[29] = 2.;
        coeff_res_c_[30] = 2.;
        coeff_res_c_[31] = 2.;
        coeff_res_c_[32] = 2.;
        coeff_res_c_[33] = 2.;
        coeff_res_c_[34] = 2.;
        coeff_res_c_[35] = 2.;
        coeff_res_c_[36] = 2.;
        coeff_res_c_[37] = 2.;
        coeff_res_c_[38] = 2.;
        coeff_res_c_[39] = 2.;
        coeff_res_c_[40] = 2.;
        coeff_res_c_[41] = 2.;
        coeff_res_c_[42] = 2.;
        coeff_res_c_[43] = 3.;
        coeff_res_c_[44] = 3.;
        coeff_res_c_[45] = 3.;
        coeff_res_c_[46] = 3.;
        coeff_res_c_[47] = 4.;
        coeff_res_c_[48] = 6.;
        coeff_res_c_[49] = 6.;
        coeff_res_c_[50] = 6.;
        coeff_res_c_[51] = 6.;
        coeff_res_c_[52] = 0.;
        coeff_res_c_[53] = 0.;
        coeff_res_c_[54] = 0.;
        // coefficient res a_i
        coeff_res_c_[55] = 3.5;
        coeff_res_c_[56] = 3.5;

        coeff_res_d_ = new double[57];
        coeff_res_d_[0] = 0.;
        // coefficient res d_i
        coeff_res_d_[1] = 1.;
        coeff_res_d_[2] = 1.;
        coeff_res_d_[3] = 1.;
        coeff_res_d_[4] = 2.;
        coeff_res_d_[5] = 2.;
        coeff_res_d_[6] = 3.;
        coeff_res_d_[7] = 4.;
        coeff_res_d_[8] = 1.;
        coeff_res_d_[9] = 1.;
        coeff_res_d_[10] = 1.;
        coeff_res_d_[11] = 2.;
        coeff_res_d_[12] = 2.;
        coeff_res_d_[13] = 3.;
        coeff_res_d_[14] = 4.;
        coeff_res_d_[15] = 4.;
        coeff_res_d_[16] = 5.;
        coeff_res_d_[17] = 7.;
        coeff_res_d_[18] = 9.;
        coeff_res_d_[19] = 10.;
        coeff_res_d_[20] = 11.;
        coeff_res_d_[21] = 13.;
        coeff_res_d_[22] = 15.;
        coeff_res_d_[23] = 1.;
        coeff_res_d_[24] = 2.;
        coeff_res_d_[25] = 2.;
        coeff_res_d_[26] = 2.;
        coeff_res_d_[27] = 3.;
        coeff_res_d_[28] = 4.;
        coeff_res_d_[29] = 4.;
        coeff_res_d_[30] = 4.;
        coeff_res_d_[31] = 5.;
        coeff_res_d_[32] = 6.;
        coeff_res_d_[33] = 6.;
        coeff_res_d_[34] = 7.;
        coeff_res_d_[35] = 9.;
        coeff_res_d_[36] = 9.;
        coeff_res_d_[37] = 9.;
        coeff_res_d_[38] = 9.;
        coeff_res_d_[39] = 9.;
        coeff_res_d_[40] = 10.;
        coeff_res_d_[41] = 10.;
        coeff_res_d_[42] = 12.;
        coeff_res_d_[43] = 3.;
        coeff_res_d_[44] = 4.;
        coeff_res_d_[45] = 4.;
        coeff_res_d_[46] = 5.;
        coeff_res_d_[47] = 14.;
        coeff_res_d_[48] = 3.;
        coeff_res_d_[49] = 6.;
        coeff_res_d_[50] = 6.;
        coeff_res_d_[51] = 6.;
        coeff_res_d_[52] = 3.;
        coeff_res_d_[53] = 3.;
        coeff_res_d_[54] = 3.;
        // coefficient res b_i
        coeff_res_d_[55] = 0.85;
        coeff_res_d_[56] = 0.95;

        coeff_res_t_ = new double[57];
        coeff_res_t_[0] = 0.;
        // coefficient res t_i
        coeff_res_t_[1] = -0.5;
        coeff_res_t_[2] = 0.875;
        coeff_res_t_[3] = 1.;
        coeff_res_t_[4] = 0.5;
        coeff_res_t_[5] = 0.75;
        coeff_res_t_[6] = 0.375;
        coeff_res_t_[7] = 1.;
        coeff_res_t_[8] = 4.;
        coeff_res_t_[9] = 6.;
        coeff_res_t_[10] = 12.;
        coeff_res_t_[11] = 1.;
        coeff_res_t_[12] = 5.;
        coeff_res_t_[13] = 4.;
        coeff_res_t_[14] = 2.;
        coeff_res_t_[15] = 13.;
        coeff_res_t_[16] = 9.;
        coeff_res_t_[17] = 3.;
        coeff_res_t_[18] = 4.;
        coeff_res_t_[19] = 11.;
        coeff_res_t_[20] = 4.;
        coeff_res_t_[21] = 13.;
        coeff_res_t_[22] = 1.;
        coeff_res_t_[23] = 7.;
        coeff_res_t_[24] = 1.;
        coeff_res_t_[25] = 9.;
        coeff_res_t_[26] = 10.;
        coeff_res_t_[27] = 10.;
        coeff_res_t_[28] = 3.;
        coeff_res_t_[29] = 7.;
        coeff_res_t_[30] = 10.;
        coeff_res_t_[31] = 10.;
        coeff_res_t_[32] = 6.;
        coeff_res_t_[33] = 10.;
        coeff_res_t_[34] = 10.;
        coeff_res_t_[35] = 1.;
        coeff_res_t_[36] = 2.;
        coeff_res_t_[37] = 3.;
        coeff_res_t_[38] = 4.;
        coeff_res_t_[39] = 8.;
        coeff_res_t_[40] = 6.;
        coeff_res_t_[41] = 9.;
        coeff_res_t_[42] = 8.;
        coeff_res_t_[43] = 16.;
        coeff_res_t_[44] = 22.;
        coeff_res_t_[45] = 23.;
        coeff_res_t_[46] = 23.;
        coeff_res_t_[47] = 10.;
        coeff_res_t_[48] = 50.;
        coeff_res_t_[49] = 44.;
        coeff_res_t_[50] = 46.;
        coeff_res_t_[51] = 50.;
        coeff_res_t_[52] = 0.;
        coeff_res_t_[53] = 1.;
        coeff_res_t_[54] = 4.;
        // coefficient res B_i
        coeff_res_t_[55] = 0.2;
        coeff_res_t_[56] = 0.2;

        coeff_res_n_ = new double[57];
        coeff_res_n_[0] = 0.;
        // coefficient res n_i
        coeff_res_n_[1] = 0.12533547935523 * 1.0e-1;
        coeff_res_n_[2] = 0.78957634722828 * 1.0e+1;
        coeff_res_n_[3] = -0.87803203303561 * 1.0e+1;
        coeff_res_n_[4] = 0.31802509345418;
        coeff_res_n_[5] = -0.26145533859358;
        coeff_res_n_[6] = -0.78199751687981 * 1.0e-2;
        coeff_res_n_[7] = 0.88089493102134 * 1.0e-2;
        coeff_res_n_[8] = -0.66856572307965;
        coeff_res_n_[9] = 0.20433810950965;
        coeff_res_n_[10] = -0.66212605039687 * 1.0e-4;
        coeff_res_n_[11] = -0.19232721156002;
        coeff_res_n_[12] = -0.25709043003438;
        coeff_res_n_[13] = 0.16074868486251;
        coeff_res_n_[14] = -0.40092828925807 * 1.0e-1;
        coeff_res_n_[15] = 0.39343422603254 * 1.0e-6;
        coeff_res_n_[16] = -0.75941377088144 * 1.0e-5;
        coeff_res_n_[17] = 0.56250979351888 * 1.0e-3;
        coeff_res_n_[18] = -0.15608652257135 * 1.0e-4;
        coeff_res_n_[19] = 0.11537996422951 * 1.0e-8;
        coeff_res_n_[20] = 0.36582165144204 * 1.0e-6;
        coeff_res_n_[21] = -0.13251180074668 * 1.0e-11;
        coeff_res_n_[22] = -0.62639586912454 * 1.0e-9;
        coeff_res_n_[23] = -0.10793600908932;
        coeff_res_n_[24] = 0.17611491008752 * 1.0e-1;
        coeff_res_n_[25] = 0.22132295167546;
        coeff_res_n_[26] = -0.40247669763528;
        coeff_res_n_[27] = 0.58083399985759;
        coeff_res_n_[28] = 0.49969146990806 * 1.0e-2;
        coeff_res_n_[29] = -0.31358700712549 * 1.0e-1;
        coeff_res_n_[30] = -0.74315929710341;
        coeff_res_n_[31] = 0.47807329915480;
        coeff_res_n_[32] = 0.20527940895948 * 1.0e-1;
        coeff_res_n_[33] = -0.13636435110343;
        coeff_res_n_[34] = 0.14180634400617 * 1.0e-1;
        coeff_res_n_[35] = 0.83326504880713 * 1.0e-2;
        coeff_res_n_[36] = -0.29052336009585 * 1.0e-1;
        coeff_res_n_[37] = 0.38615085574206 * 1.0e-1;
        coeff_res_n_[38] = -0.20393486513704 * 1.0e-1;
        coeff_res_n_[39] = -0.16554050063734 * 1.0e-2;
        coeff_res_n_[40] = 0.19955571979541 * 1.0e-2;
        coeff_res_n_[41] = 0.15870308324157 * 1.0e-3;
        coeff_res_n_[42] = -0.16388568342530 * 1.0e-4;
        coeff_res_n_[43] = 0.43613615723811 * 1.0e-1;
        coeff_res_n_[44] = 0.34994005463765 * 1.0e-1;
        coeff_res_n_[45] = -0.76788197844621 * 1.0e-1;
        coeff_res_n_[46] = 0.22446277332006 * 1.0e-1;
        coeff_res_n_[47] = -0.62689710414685 * 1.0e-4;
        coeff_res_n_[48] = -0.55711118565645 * 1.0e-9;
        coeff_res_n_[49] = -0.19905718354408;
        coeff_res_n_[50] = 0.31777497330738;
        coeff_res_n_[51] = -0.11841182425981;
        coeff_res_n_[52] = -0.31306260323435 * 1.0e+2;
        coeff_res_n_[53] = 0.31546140237781 * 1.0e+2;
        coeff_res_n_[54] = -0.25213154341695 * 1.0e+4;
        coeff_res_n_[55] = -0.14874640856724;
        coeff_res_n_[56] = 0.31806110878444;

        coeff_res_alpha_ = new double[3];
        // coefficient res alpha_{i + 52}
        coeff_res_alpha_[0] = 20.;
        coeff_res_alpha_[1] = 20.;
        coeff_res_alpha_[2] = 20.;

        coeff_res_gamma_ = new double[3];
        // coefficient res gamma_{i + 52}
        coeff_res_gamma_[0] = 1.21;
        coeff_res_gamma_[1] = 1.21;
        coeff_res_gamma_[2] = 1.25;

        coeff_res_epsilon_ = new double[3];
        // coefficient res epsilon_{i + 52}
        coeff_res_epsilon_[0] = 1.;
        coeff_res_epsilon_[1] = 1.;
        coeff_res_epsilon_[2] = 1.;

        coeff_res_beta_ = new double[5];
        // coefficient res beta_{i + 52}
        coeff_res_beta_[0] = 150.;
        coeff_res_beta_[1] = 150.;
        coeff_res_beta_[2] = 250.;
        coeff_res_beta_[3] = 0.3;
        coeff_res_beta_[4] = 0.3;

        coeff_res_C_ = new double[2];
        // coefficient res C_{i + 55}
        coeff_res_C_[0] = 28.;
        coeff_res_C_[1] = 32.;

        coeff_res_D_ = new double[2];
        // coefficient res D_{i + 55}
        coeff_res_D_[0] = 700.;
        coeff_res_D_[1] = 800.;

        coeff_res_A_ = new double[2];
        // coefficient res A_{i + 55}
        coeff_res_A_[0] = 0.32;
        coeff_res_A_[1] = 0.32;

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
