#ifndef _LIBIAPWS76_H_
#define _LIBIAPWS76_H_

/* implementation of IAPWS76-2014
 * Revised Release on Surface Tension
 * of Ordinary Water Substance */
class LibIAPWS76 {
  private :

    /* critical temperature
     * in degK */
    double temperature_crit_;

    double coeff_B_;
    double coeff_b_;
    double coeff_mu_;

  public :

    LibIAPWS76() {

        temperature_crit_ = 647.096;

        coeff_B_ = 0.2358;
        coeff_b_ = -0.625;
        coeff_mu_ = 1.256;

        return;
    }

    ~LibIAPWS76() {}

    /* returns surface tension
     * between (ordinary) liquid water and vapor
     * in J / m^2 */
    double get_tension_surf(double temperature_in);
};

#endif
