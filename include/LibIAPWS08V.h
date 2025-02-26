#ifndef _LIBIAPWS08V_H_
#define _LIBIAPWS08V_H_

#include<stdio.h>
#include"LibIAPWS95.h"

namespace IAPWS {

/* implementation of IAPWS R12-08 (2008)
 * Release on the IAPWS Formulation 2008
 * for the Viscosity of Ordinary Water Substance */
class Lib08V {
  private :

    /* the reference temperature
     * in degK */
    double temperature_ref_;
    /* the reference mass density
     * in kg / m^3 */
    double mdensity_ref_;
    /* the reference pressure
     * in Pa */
    double pressure_ref_;
    /* the reference temperature
     * in Pa * sec */
    double viscosity_ref_;

    double *coeff0_H_;

    double **coeff1_H_;

    double coeff2_x_mu_;
    double coeff2_q_C_;
    double coeff2_q_D_;
    double coeff2_nu_;
    double coeff2_gamma_;
    double coeff2_xi0_;
    double coeff2_Gamma0_;
    double coeff2_TR_;

    /* range in mass density (kg / m^3)
     * where the critical enhancement is applied */
    double range_visc2_mdensity_min_;
    double range_visc2_mdensity_max_;

    /* range in temperature (deg K)
     * where the critical enhancement is applied */
    double range_visc2_temperature_min_;
    double range_visc2_temperature_max_;

    Lib95 *ptr_lib95_;

  public :

    Lib08V() {
        temperature_ref_ = 647.096;
        mdensity_ref_ = 322.;
        pressure_ref_ = 22.064 * 1.0e+6;
        viscosity_ref_ = 1.0e-6;

        coeff0_H_ = new double[4];
        coeff0_H_[0] = 1.67752;
        coeff0_H_[1] = 2.20462;
        coeff0_H_[2] = 0.6366564;
        coeff0_H_[3] = -0.241605;

        coeff1_H_ = new double *[6];
        for (int i = 0; i <= 5; i++) {
            coeff1_H_[i] = new double[7];

            for (int j = 0; j <= 6; j++) {
                coeff1_H_[i][j] = 0.;
            }
        }

        coeff1_H_[0][0] = 5.20094 * 1.0e-1;
        coeff1_H_[1][0] = 8.50895 * 1.0e-2;
        coeff1_H_[2][0] = -1.08374;
        coeff1_H_[3][0] = -2.89555 * 1.0e-1;
        coeff1_H_[0][1] = 2.22531 * 1.0e-1;
        coeff1_H_[1][1] = 9.99115 * 1.0e-1;
        coeff1_H_[2][1] = 1.88797;
        coeff1_H_[3][1] = 1.26613;
        coeff1_H_[5][1] = 1.20573 * 1.0e-1;
        coeff1_H_[0][2] = -2.81378 * 1.0e-1;
        coeff1_H_[1][2] = -9.06851 * 1.0e-1;
        coeff1_H_[2][2] = -7.72479 * 1.0e-1;
        coeff1_H_[3][2] = -4.89837 * 1.0e-1;
        coeff1_H_[4][2] = -2.57040 * 1.0e-1;
        coeff1_H_[0][3] = 1.61913 * 1.0e-1;
        coeff1_H_[1][3] = 2.57399 * 1.0e-1;
        coeff1_H_[0][4] = -3.25372 * 1.0e-2;
        coeff1_H_[3][4] = 6.98452 * 1.0e-2;
        coeff1_H_[4][5] = 8.72102 * 1.0e-3;
        coeff1_H_[3][6] = -4.35673 * 1.0e-3;
        coeff1_H_[5][6] = -5.93264 * 1.0e-4;

        coeff2_x_mu_ = 0.068;
        coeff2_q_C_ = 1. / 1.9;
        coeff2_q_D_ = 1. / 1.1;
        coeff2_nu_ = 0.630;
        coeff2_gamma_ = 1.239;
        coeff2_xi0_ = 0.13;
        coeff2_Gamma0_ = 0.06;
        coeff2_TR_ = 1.5;

        range_visc2_mdensity_min_ = 245.8;
        range_visc2_mdensity_max_ = 405.3;

        range_visc2_temperature_min_ = 645.91;
        range_visc2_temperature_max_ = 650.77;

        ptr_lib95_ = NULL;

        return;
    }

    ~Lib08V() {
        delete [] coeff0_H_;

        for (int i = 0; i <= 5; i++) {
            delete [] coeff1_H_[i];
        }
        delete [] coeff1_H_;

        return;
    }

    void print_header(FILE *ptr_fout = stdout);

    void set_ptr_lib95(Lib95 *ptr_in) {
        ptr_lib95_ = ptr_in;

        return;
    }

    void set_range_visc2_mdensity(double mdensity_min,
                                  double mdensity_max) {
        if (mdensity_min < mdensity_ref_) {
            range_visc2_mdensity_min_ = mdensity_min;
        } else {
            range_visc2_mdensity_min_ = mdensity_ref_;
        }

        if (mdensity_max > mdensity_ref_) {
            range_visc2_mdensity_max_ = mdensity_max;
        } else {
            range_visc2_mdensity_max_ = mdensity_ref_;
        }

        return;
    }

    void set_range_visc2_temperature(double temp_min,
                                     double temp_max) {
        if (temp_min < temperature_ref_) {
            range_visc2_temperature_min_ = temp_min;
        } else {
            range_visc2_temperature_min_ = temperature_ref_;
        }

        if (temp_max > temperature_ref_) {
            range_visc2_temperature_max_ = temp_max;
        } else {
            range_visc2_temperature_max_ = temperature_ref_;
        }

        return;
    }

    /* returns shear viscosity
     * of (ordinary) water fluid
     * in Pa * sec */
    double get_param_viscosity(double mdensity_in,
                               double tempearture_in);

    double get_param_visc0_dimless(double tempearture_in);

    double get_param_visc1_dimless(double mdensity_in,
                                   double tempearture_in);

    double get_param_visc2_dimless(double mdensity_in,
                                   double tempearture_in);
};

} // end namespace IAPWS

#endif
