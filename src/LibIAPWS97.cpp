#include<math.h>
#include"LibIAPWS97.h"

namespace IAPWS {

void Lib97::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R7-97(2012)\n");
    fprintf(ptr_fout, "Revised Release on the IAPWS");
    fprintf(ptr_fout, " Industrial Formulation 1997\n");
    fprintf(ptr_fout, "for the Thermodynamic Properties");
    fprintf(ptr_fout, " of Water and Steam\n");
    fprintf(ptr_fout, "\n");

    return;
}

double Lib97::get_param_g(double temperature_in,
                          double pressure_in,
                          bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double fn_gamma;

    switch(n_reg) {
      case 1 :
        fn_gamma =
            get_param1_gamma(temperature_in,
                             pressure_in);
        break;
      case 2 :
        fn_gamma =
            get_param2_gamma_ide(temperature_in,
                                 pressure_in,
                                 flag_metastable) +
            get_param2_gamma_res(temperature_in,
                                 pressure_in,
                                 flag_metastable);
        break;
      case 3 :
        break;
      case 5 :
        fn_gamma =
            get_param5_gamma_ide(temperature_in,
                                 pressure_in) +
            get_param5_gamma_res(temperature_in,
                                 pressure_in);
        break;
      default :
        fn_gamma = 0.;
    }

    double g =
        const_R_spec_ * temperature_in * fn_gamma;

    return g;
}

double Lib97::get_param_vol_spec(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double ppi = 0.;
    double fn_dgamma_dppi;

    double vol_spec;

    switch(n_reg) {
      case 1 :
        ppi = pressure_in / pressure_ref1_;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        break;
      case 2 :
        if (flag_metastable) {
            ppi = pressure_in / pressure_ref2mst_;
        } else {
            ppi = pressure_in / pressure_ref2_;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        break;
      case 3 :
        break;
      case 5 :
        ppi = pressure_in / pressure_ref5_;
        fn_dgamma_dppi =
            get_param5_dgamma_ide_dppi(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dppi(temperature_in,
                                       pressure_in);
        vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        break;
      default :
        vol_spec = 0.;
    }

    return vol_spec;
}

double Lib97::get_param_mdensity(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    double mdensity =
        1. / get_param_vol_spec(temperature_in,
                                pressure_in,
                                flag_metastable);

    return mdensity;
}

double Lib97::get_param_erg_int(double temperature_in,
                                double pressure_in,
                                bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double ppi = 0.;
    double fn_dgamma_dtau;
    double fn_dgamma_dppi;

    double erg_int;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        ppi = pressure_in / pressure_ref1_;
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * fn_dgamma_dtau -
             ppi * fn_dgamma_dppi);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
            ppi = pressure_in / pressure_ref2mst_;
        } else {
            tau = temperature_ref2_ / temperature_in;
            ppi = pressure_in / pressure_ref2_;
        }
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * fn_dgamma_dtau -
             ppi * fn_dgamma_dppi);
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        ppi = pressure_in / pressure_ref5_;
        fn_dgamma_dtau =
            get_param5_dgamma_ide_dtau(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dtau(temperature_in,
                                       pressure_in);
        fn_dgamma_dppi =
            get_param5_dgamma_ide_dppi(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dppi(temperature_in,
                                       pressure_in);
        erg_int =
            const_R_spec_ * temperature_in *
            (tau * fn_dgamma_dtau -
             ppi * fn_dgamma_dppi);
        break;
      default :
        erg_int = 0.;
    }

    return erg_int;
}

double Lib97::get_param_entropy(double temperature_in,
                                double pressure_in,
                                bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_gamma;
    double fn_dgamma_dtau;

    double entropy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_gamma =
            get_param1_gamma(temperature_in,
                             pressure_in);
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        entropy =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_gamma =
            get_param2_gamma_ide(temperature_in,
                                 pressure_in,
                                 flag_metastable) +
            get_param2_gamma_res(temperature_in,
                                 pressure_in,
                                 flag_metastable);
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        entropy =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        fn_gamma =
            get_param5_gamma_ide(temperature_in,
                                 pressure_in) +
            get_param5_gamma_res(temperature_in,
                                 pressure_in);
        fn_dgamma_dtau =
            get_param5_dgamma_ide_dtau(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dtau(temperature_in,
                                       pressure_in);
        entropy =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
        break;
      default :
        entropy = 0.;
    }

    return entropy;
}

double Lib97::get_param_enthalpy(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dtau;

    double enthalpy;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        fn_dgamma_dtau =
            get_param5_dgamma_ide_dtau(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dtau(temperature_in,
                                       pressure_in);
        enthalpy =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
        break;
      default :
        enthalpy = 0.;
    }

    return enthalpy;
}

double Lib97::get_param_heat_c_p(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_d2gamma_dtau_dtau;

    double heat_c_p;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        heat_c_p =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        heat_c_p =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau;
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        fn_d2gamma_dtau_dtau =
            get_param5_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in);
        heat_c_p =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau;
        break;
      default :
        heat_c_p = 0.;
    }

    return heat_c_p;
}

double Lib97::get_param_heat_c_v(double temperature_in,
                                 double pressure_in,
                                 bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi;
    double fn_d2gamma_dppi_dtau;
    double fn_d2gamma_dtau_dtau;

    double heat_c_v;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param1_d2gamma_dppi_dppi(temperature_in,
                                         pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param1_d2gamma_dppi_dtau(temperature_in,
                                         pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        heat_c_v =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau +
            const_R_spec_ *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
            fn_d2gamma_dppi_dppi;
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_d2gamma_dppi_dppi =
            get_param2_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dppi_dtau =
            get_param2_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        heat_c_v =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau +
            const_R_spec_ *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
            fn_d2gamma_dppi_dppi;
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        fn_dgamma_dppi =
            get_param5_dgamma_ide_dppi(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dppi(temperature_in,
                                       pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param5_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param5_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param5_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in);
        heat_c_v =
            -const_R_spec_ * tau * tau * fn_d2gamma_dtau_dtau +
            const_R_spec_ *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
            (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
            fn_d2gamma_dppi_dppi;
        break;
      default :
        heat_c_v = 0.;
    }

    return heat_c_v;
}

double Lib97::get_param_speed_sound(double temperature_in,
                                    double pressure_in,
                                    bool flag_metastable) {
    int n_reg = get_region(temperature_in,
                           pressure_in,
                           flag_metastable);

    double tau = 0.;
    double fn_dgamma_dppi = 0.;
    double fn_d2gamma_dppi_dppi;
    double fn_d2gamma_dppi_dtau;
    double fn_d2gamma_dtau_dtau;

    double v2_s;

    switch(n_reg) {
      case 1 :
        tau = temperature_ref1_ / temperature_in;
        fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param1_d2gamma_dppi_dppi(temperature_in,
                                         pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param1_d2gamma_dppi_dtau(temperature_in,
                                         pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param1_d2gamma_dtau_dtau(temperature_in,
                                         pressure_in);
        v2_s = const_R_spec_ * temperature_in *
            fn_dgamma_dppi * fn_dgamma_dppi /
            ((fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
             (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
             (tau * tau * fn_d2gamma_dtau_dtau) -
             fn_d2gamma_dppi_dppi);
        break;
      case 2 :
        if (flag_metastable) {
            tau = temperature_ref2mst_ / temperature_in;
        } else {
            tau = temperature_ref2_ / temperature_in;
        }
        fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in,
                                       flag_metastable);
        fn_d2gamma_dppi_dppi =
            get_param2_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dppi_dtau =
            get_param2_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        fn_d2gamma_dtau_dtau =
            get_param2_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable) +
            get_param2_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in,
                                             flag_metastable);
        v2_s = const_R_spec_ * temperature_in *
            fn_dgamma_dppi * fn_dgamma_dppi /
            ((fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
             (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
             (tau * tau * fn_d2gamma_dtau_dtau) -
             fn_d2gamma_dppi_dppi);
        break;
      case 3 :
        break;
      case 5 :
        tau = temperature_ref5_ / temperature_in;
        fn_dgamma_dppi =
            get_param5_dgamma_ide_dppi(temperature_in,
                                       pressure_in) +
            get_param5_dgamma_res_dppi(temperature_in,
                                       pressure_in);
        fn_d2gamma_dppi_dppi =
            get_param5_d2gamma_ide_dppi_dppi(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dppi_dppi(temperature_in,
                                             pressure_in);
        fn_d2gamma_dppi_dtau =
            get_param5_d2gamma_ide_dppi_dtau(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dppi_dtau(temperature_in,
                                             pressure_in);
        fn_d2gamma_dtau_dtau =
            get_param5_d2gamma_ide_dtau_dtau(temperature_in,
                                             pressure_in) +
            get_param5_d2gamma_res_dtau_dtau(temperature_in,
                                             pressure_in);
        v2_s = const_R_spec_ * temperature_in *
            fn_dgamma_dppi * fn_dgamma_dppi /
            ((fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) *
             (fn_dgamma_dppi - tau * fn_d2gamma_dppi_dtau) /
             (tau * tau * fn_d2gamma_dtau_dtau) -
             fn_d2gamma_dppi_dppi);
        break;
      default :
        v2_s = 0.;
    }

    return sqrt(v2_s);
}

int Lib97::get_region(double temperature_in,
                      double pressure_in,
                      bool flag_metastable) {
    int n_reg = 0;

    if (flag_metastable) {
        n_reg = 2;
    } else if (temperature_in >= 273.15 &&
               temperature_in <= 623.15) {
        double press_sat =
            get_param4_sat_pressure(temperature_in);

        if (pressure_in <= press_sat) {
            n_reg = 2;
        } else if (pressure_in > press_sat &&
                   pressure_in < 1.0e+8) {
            n_reg = 1;
        }
    } else if (temperature_in > 623.15 &&
               temperature_in <= 863.15) {
        double press_B23 =
            get_paramB23_pressure(temperature_in);

        if (pressure_in <= press_B23) {
            n_reg = 2;
        } else if (pressure_in > press_B23 &&
                   pressure_in < 1.0e+8) {
            n_reg = 3;
        }
    } else if (temperature_in > 863.15 &&
               temperature_in <= 1073.15) {
        if (pressure_in > 0. &&
            pressure_in < 1.0e+8) {
            n_reg = 2;
        }
    } else if (temperature_in > 1073.15 &&
               temperature_in <= 2273.15) {
        if (pressure_in > 0. &&
            pressure_in < 50. * 1.0e+6) {
            n_reg = 5;
        }
    }

    return n_reg;
}

double Lib97::get_coex_mden_vap(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double mden_vap = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);
        double ppi = pressure_in / pressure_ref2_;

        double fn_dgamma_dppi =
            get_param2_dgamma_ide_dppi(temperature_in,
                                       pressure_in) +
            get_param2_dgamma_res_dppi(temperature_in,
                                       pressure_in);
        double vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        mden_vap = 1. / vol_spec;
    } else {
        // region 3
    }

    return mden_vap;
}

double Lib97::get_coex_mden_liq(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double mden_liq = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);
        double ppi = pressure_in / pressure_ref1_;

        double fn_dgamma_dppi =
            get_param1_dgamma_dppi(temperature_in,
                                   pressure_in);
        double vol_spec =
            (const_R_spec_ * temperature_in /
             pressure_in) *
            ppi * fn_dgamma_dppi;
        mden_liq = 1. / vol_spec;
    } else {
        // region 3
    }

    return mden_liq;
}

double Lib97::get_coex_enthalpy_vap(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double enthalpy_vap = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);
        double tau = temperature_ref2_ / temperature_in;

        double fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in);
        enthalpy_vap =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
    } else {
        // region 3
    }

    return enthalpy_vap;
}

double Lib97::get_coex_enthalpy_liq(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double enthalpy_liq = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);

        double tau = temperature_ref1_ / temperature_in;
        double fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        enthalpy_liq =
            const_R_spec_ * temperature_in *
            tau * fn_dgamma_dtau;
    } else {
        // region 3
    }

    return enthalpy_liq;
}

double Lib97::get_coex_entropy_vap(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double pressure_in =
        get_param4_sat_pressure(temperature_in);
    double entropy_vap = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);
        double tau = temperature_ref2_ / temperature_in;

        double fn_gamma =
            get_param2_gamma_ide(temperature_in,
                                 pressure_in) +
            get_param2_gamma_res(temperature_in,
                                 pressure_in);
        double fn_dgamma_dtau =
            get_param2_dgamma_ide_dtau(temperature_in,
                                       pressure_in) +
            get_param2_dgamma_res_dtau(temperature_in,
                                       pressure_in);
        entropy_vap =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
    } else {
        // region 3
    }

    return entropy_vap;
}

double Lib97::get_coex_entropy_liq(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double entropy_liq = 0.;

    if (temperature_in <= 623.15) {
        // boundary between region 1 and 2
        double pressure_in =
            get_param4_sat_pressure(temperature_in);
        double tau = temperature_ref1_ / temperature_in;

        double fn_gamma =
            get_param1_gamma(temperature_in,
                             pressure_in);
        double fn_dgamma_dtau =
            get_param1_dgamma_dtau(temperature_in,
                                   pressure_in);
        entropy_liq =
            const_R_spec_ *
            (tau * fn_dgamma_dtau - fn_gamma);
    } else {
        // region 3
    }

    return entropy_liq;
}

double Lib97::get_coex_heat_latent(double temperature_in) {
    if (temperature_in >= temperature_crit_) {
        return 0.;
    }

    double h_latent =
        temperature_in *
        (get_coex_entropy_vap(temperature_in) -
         get_coex_entropy_liq(temperature_in));

    return h_latent;
}

double Lib97::get_paramB23_pressure(double temperature_in) {
    double tau = temperature_in / temperature_refB23_;

    double ppi =
        coeffB23_n_[1] +
        coeffB23_n_[2] * tau +
        coeffB23_n_[3] * tau * tau;

    return pressure_refB23_ * ppi;
}

double Lib97::get_paramB23_temperature(double pressure_in) {
    double ppi = pressure_in / pressure_refB23_;

    double tau =
        coeffB23_n_[4] +
        sqrt((ppi - coeffB23_n_[5]) / coeffB23_n_[3]);

    return temperature_refB23_ * tau;
}

double Lib97::get_param1_gamma(double temperature_in,
                               double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double gamma = 0.;
    for (int i = 1; i <= 34; i++) {
        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        gamma += coeff1_n_[i] * dppi_pow_I * dtau_pow_J;
    }

    return gamma;
}

double Lib97::get_param1_dgamma_dppi(double temperature_in,
                                     double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double dgamma_dppi = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 1);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        dgamma_dppi -=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i]);
    }

    return dgamma_dppi;
}

double Lib97::get_param1_dgamma_dtau(double temperature_in,
                                     double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double dgamma_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_J_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 1);

        dgamma_dtau +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param1_d2gamma_dppi_dppi(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dppi_dppi = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0 ||
            coeff1_I_[i] == 1) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 2);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i]);

        d2gamma_dppi_dppi +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i] * (coeff1_I_[i] - 1));
    }

    return d2gamma_dppi_dppi;
}

double Lib97::get_param1_d2gamma_dppi_dtau(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dppi_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_I_[i] == 0 ||
            coeff1_J_[i] == 0) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i] - 1);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 1);

        d2gamma_dppi_dtau -=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_I_[i] * coeff1_J_[i]);
    }

    return d2gamma_dppi_dtau;
}

double Lib97::get_param1_d2gamma_dtau_dtau(double temperature_in,
                                           double pressure_in) {
    double tau = temperature_ref1_ / temperature_in;
    double ppi = pressure_in / pressure_ref1_;

    double dtau = tau - tau_ref1_;
    double dppi = ppi_ref1_ - ppi;

    double d2gamma_dtau_dtau = 0.;
    for (int i = 1; i <= 34; i++) {
        if (coeff1_J_[i] == 0 ||
            coeff1_J_[i] == 1) {
            continue;
        }

        double dppi_pow_I =
            get_pow_int(dppi, coeff1_I_[i]);
        double dtau_pow_J =
            get_pow_int(dtau, coeff1_J_[i] - 2);

        d2gamma_dtau_dtau +=
            coeff1_n_[i] * dppi_pow_I * dtau_pow_J *
            static_cast<double>(coeff1_J_[i] * (coeff1_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param1_temperature_ph(double pressure_in,
                                        double enthalpy_in) {
    double ppi = pressure_in / pressure_ref1Tph_;
    double eta = enthalpy_in / enthalpy_ref1Tph_;

    double fn_theta = 0.;
    for (int i = 1; i <= 20; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff1Tph_I_[i]);
        double eta_pow_J =
            get_pow_int(eta + 1., coeff1Tph_J_[i]);

        fn_theta +=
            coeff1Tph_n_[i] * ppi_pow_I * eta_pow_J;
    }

    return temperature_ref1Tph_ * fn_theta;
}

double Lib97::get_param1_temperature_ps(double pressure_in,
                                        double entropy_in) {
    double ppi = pressure_in / pressure_ref1Tps_;
    double sig = entropy_in / entropy_ref1Tps_;

    double fn_theta = 0.;
    for (int i = 1; i <= 20; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff1Tps_I_[i]);
        double sig_pow_J =
            get_pow_int(sig + 2., coeff1Tps_J_[i]);

        fn_theta +=
            coeff1Tps_n_[i] * ppi_pow_I * sig_pow_J;
    }

    return temperature_ref1Tps_ * fn_theta;
}

double Lib97::get_param2_gamma_ide(double temperature_in,
                                   double pressure_in,
                                   bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double gamma = log(ppi);

    for (int i = 1; i <= 9; i++) {
        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i]);
        gamma +=
            coeff_n * dtau_pow_J;
    }

    return gamma;
}

double Lib97::get_param2_dgamma_ide_dppi(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dgamma_dppi = 1. / ppi;

    return dgamma_dppi;
}

double Lib97::get_param2_dgamma_ide_dtau(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    double dgamma_dtau = 0.;

    for (int i = 1; i <= 9; i++) {
        if (coeff2_ide_J_[i] == 0) {
            continue;
        }

        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i] - 1);
        dgamma_dtau +=
            coeff_n * dtau_pow_J *
            static_cast<double>(coeff2_ide_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param2_d2gamma_ide_dppi_dppi(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double d2gamma_dppi_dppi = -1. / (ppi * ppi);

    return d2gamma_dppi_dppi;
}

double Lib97::get_param2_d2gamma_ide_dppi_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    //double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    return 0.;
}

double Lib97::get_param2_d2gamma_ide_dtau_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    //double ppi = pressure_in / pressure_ref2_;

    double d2gamma_dtau_dtau = 0.;

    for (int i = 1; i <= 9; i++) {
        if (coeff2_ide_J_[i] == 0 ||
            coeff2_ide_J_[i] == 1) {
            continue;
        }

        double coeff_n = coeff2_ide_n_[i];

        if (flag_metastable) {
            if (i == 1) {
                coeff_n = coeff2mst_ide_n1_;
            }

            if (i == 2) {
                coeff_n = coeff2mst_ide_n2_;
            }
        }

        double dtau_pow_J = get_pow_int(tau, coeff2_ide_J_[i] - 2);
        d2gamma_dtau_dtau +=
            coeff_n * dtau_pow_J *
            static_cast<double>(coeff2_ide_J_[i] *
                                (coeff2_ide_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param2_gamma_res(double temperature_in,
                                   double pressure_in,
                                   bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double gamma = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            gamma +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J;
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            gamma +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J;
        }
    }

    return gamma;
}

double Lib97::get_param2_dgamma_res_dppi(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double dgamma_dppi = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            dgamma_dppi +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            dgamma_dppi +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i]);
        }
    }

    return dgamma_dppi;
}

double Lib97::get_param2_dgamma_res_dtau(double temperature_in,
                                         double pressure_in,
                                         bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double dgamma_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 1);

            dgamma_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_J_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 1);

            dgamma_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_J_[i]);
        }
    }

    return dgamma_dtau;
}

double Lib97::get_param2_d2gamma_res_dppi_dppi(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dppi_dppi = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0 ||
                coeff2mst_res_I_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 2);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i]);

            d2gamma_dppi_dppi +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i] *
                                    (coeff2mst_res_I_[i] - 1));
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0 ||
                coeff2_res_I_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 2);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i]);

            d2gamma_dppi_dppi +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i] *
                                    (coeff2_res_I_[i] - 1));
        }
    }

    return d2gamma_dppi_dppi;
}

double Lib97::get_param2_d2gamma_res_dppi_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dppi_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_I_[i] == 0 ||
                coeff2mst_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 1);

            d2gamma_dppi_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_I_[i] *
                                    coeff2mst_res_J_[i]);
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_I_[i] == 0 ||
                coeff2_res_J_[i] == 0) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i] - 1);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 1);

            d2gamma_dppi_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_I_[i] *
                                    coeff2_res_J_[i]);
        }
    }

    return d2gamma_dppi_dtau;
}

double Lib97::get_param2_d2gamma_res_dtau_dtau(double temperature_in,
                                               double pressure_in,
                                               bool flag_metastable) {
    double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref2_;

    double dtau = tau - tau_ref2_res_;

    double d2gamma_dtau_dtau = 0.;
    if (flag_metastable) {
        tau = temperature_ref2mst_ / temperature_in;
        ppi = pressure_in / pressure_ref2mst_;

        dtau = tau - tau_ref2mst_res_;

        for (int i = 1; i <= 13; i++) {
            if (coeff2mst_res_J_[i] == 0 ||
                coeff2mst_res_J_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2mst_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2mst_res_J_[i] - 2);

            d2gamma_dtau_dtau +=
                coeff2mst_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2mst_res_J_[i] *
                                    (coeff2mst_res_J_[i] - 1));
        }
    } else {
        for (int i = 1; i <= 43; i++) {
            if (coeff2_res_J_[i] == 0 ||
                coeff2_res_J_[i] == 1) {
                continue;
            }

            double ppi_pow_I =
                get_pow_int(ppi, coeff2_res_I_[i]);
            double dtau_pow_J =
                get_pow_int(dtau, coeff2_res_J_[i] - 2);

            d2gamma_dtau_dtau +=
                coeff2_res_n_[i] * ppi_pow_I * dtau_pow_J *
                static_cast<double>(coeff2_res_J_[i] *
                                    (coeff2_res_J_[i] - 1));
        }
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_paramB2bc_pressure(double enthalpy_in) {
    double eta =
        enthalpy_in / enthalpy_refB2bc_;

    double ppi =
        coeffB2bc_n_[1] +
        coeffB2bc_n_[2] * eta +
        coeffB2bc_n_[3] * eta * eta;

    return pressure_refB2bc_ * ppi;
}

double Lib97::get_paramB2bc_enthalpy(double pressure_in) {
    double ppi =
        pressure_in / pressure_refB2bc_;

    double eta =
        coeffB2bc_n_[4] +
        sqrt((ppi - coeffB2bc_n_[5]) / coeffB2bc_n_[3]);

    return enthalpy_refB2bc_ * eta;
}

double Lib97::get_param2_temperature_ph(double pressure_in,
                                        double enthalpy_in) {
    double temperature_out = 0.;

    if (pressure_in < 4. * 1.0e+6) {
        temperature_out =
            get_param2a_temperature_ph(pressure_in,
                                       enthalpy_in);
    } else if (pressure_in < 6.54670 * 1.0e+6) {
        temperature_out =
            get_param2b_temperature_ph(pressure_in,
                                       enthalpy_in);
    } else {
        double enthalpy_threshold =
            get_paramB2bc_enthalpy(pressure_in);

        if (enthalpy_in > enthalpy_threshold) {
            temperature_out =
                get_param2b_temperature_ph(pressure_in,
                                           enthalpy_in);
        } else {
            temperature_out =
                get_param2c_temperature_ph(pressure_in,
                                           enthalpy_in);
        }
    }

    return temperature_out;
}

double Lib97::get_param2a_temperature_ph(double pressure_in,
                                         double enthalpy_in) {
    double ppi = pressure_in / pressure_ref2aTph_;
    double eta = enthalpy_in / enthalpy_ref2aTph_;

    double fn_theta = 0.;
    for (int i = 1; i <= 34; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff2aTph_I_[i]);
        double eta_pow_J =
            get_pow_int(eta - 2.1, coeff2aTph_J_[i]);

        fn_theta +=
            coeff2aTph_n_[i] * ppi_pow_I * eta_pow_J;
    }

    return temperature_ref2aTph_ * fn_theta;
}

double Lib97::get_param2b_temperature_ph(double pressure_in,
                                         double enthalpy_in) {
    double ppi = pressure_in / pressure_ref2bTph_;
    double eta = enthalpy_in / enthalpy_ref2bTph_;

    double fn_theta = 0.;
    for (int i = 1; i <= 38; i++) {
        double ppi_pow_I =
            get_pow_int(ppi - 2., coeff2bTph_I_[i]);
        double eta_pow_J =
            get_pow_int(eta - 2.6, coeff2bTph_J_[i]);

        fn_theta +=
            coeff2bTph_n_[i] * ppi_pow_I * eta_pow_J;
    }

    return temperature_ref2bTph_ * fn_theta;
}

double Lib97::get_param2c_temperature_ph(double pressure_in,
                                         double enthalpy_in) {
    double ppi = pressure_in / pressure_ref2cTph_;
    double eta = enthalpy_in / enthalpy_ref2cTph_;

    double fn_theta = 0.;
    for (int i = 1; i <= 23; i++) {
        double ppi_pow_I =
            get_pow_int(ppi + 25., coeff2cTph_I_[i]);
        double eta_pow_J =
            get_pow_int(eta - 1.8, coeff2cTph_J_[i]);

        fn_theta +=
            coeff2cTph_n_[i] * ppi_pow_I * eta_pow_J;
    }

    return temperature_ref2cTph_ * fn_theta;
}

double Lib97::get_param2_temperature_ps(double pressure_in,
                                        double entropy_in) {
    double temperature_out = 0.;

    if (pressure_in < 4. * 1.0e+6) {
        temperature_out =
            get_param2a_temperature_ps(pressure_in,
                                       entropy_in);
    } else if (pressure_in < 6.54670 * 1.0e+6) {
        temperature_out =
            get_param2b_temperature_ps(pressure_in,
                                       entropy_in);
    } else {
        if (entropy_in > 5.85 * 1.0e+3) {
            temperature_out =
                get_param2b_temperature_ps(pressure_in,
                                           entropy_in);
        } else {
            temperature_out =
                get_param2c_temperature_ps(pressure_in,
                                           entropy_in);
        }
    }

    return temperature_out;
}

double Lib97::get_param2a_temperature_ps(double pressure_in,
                                         double entropy_in) {
    double ppi = pressure_in / pressure_ref2aTps_;
    double sig = entropy_in / entropy_ref2aTps_;

    double fn_theta = 0.;
    for (int i = 1; i <= 46; i++) {
        double ppi_pow_I =
            pow(ppi, coeff2aTps_I_[i]);
        double sig_pow_J =
            get_pow_int(sig - 2., coeff2aTps_J_[i]);

        fn_theta +=
            coeff2aTps_n_[i] * ppi_pow_I * sig_pow_J;
    }

    return temperature_ref2aTps_ * fn_theta;
}

double Lib97::get_param2b_temperature_ps(double pressure_in,
                                         double entropy_in) {
    double ppi = pressure_in / pressure_ref2bTps_;
    double sig = entropy_in / entropy_ref2bTps_;

    double fn_theta = 0.;
    for (int i = 1; i <= 44; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff2bTps_I_[i]);
        double sig_pow_J =
            get_pow_int(10. - sig, coeff2bTps_J_[i]);

        fn_theta +=
            coeff2bTps_n_[i] * ppi_pow_I * sig_pow_J;
    }

    return temperature_ref2bTps_ * fn_theta;
}

double Lib97::get_param2c_temperature_ps(double pressure_in,
                                         double entropy_in) {
    double ppi = pressure_in / pressure_ref2cTps_;
    double sig = entropy_in / entropy_ref2cTps_;

    double fn_theta = 0.;
    for (int i = 1; i <= 30; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff2cTps_I_[i]);
        double sig_pow_J =
            get_pow_int(2. - sig, coeff2cTps_J_[i]);

        fn_theta +=
            coeff2cTps_n_[i] * ppi_pow_I * sig_pow_J;
    }

    return temperature_ref2cTps_ * fn_theta;
}

double Lib97::get_param3_f(double mdensity_in,
                           double temperature_in) {
    double f =
        const_R_spec_ * temperature_in *
        get_param3_phi(mdensity_in,
                       temperature_in);

    return f;
}

double Lib97::get_param3_g(double mdensity_in,
                           double temperature_in) {
    double g =
        get_param3_f(mdensity_in,
                    temperature_in) +
        get_param3_pressure(mdensity_in,
                            temperature_in) / mdensity_in;

    return g;
}

double Lib97::get_param3_pressure(double mdensity_in,
                                  double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    //double tau = temperature_ref3_ / temperature_in;

    double press =
        mdensity_in * const_R_spec_ * temperature_in *
        delta * get_param3_dphi_ddelta(mdensity_in,
                                       temperature_in);

    return press;
}

double Lib97::get_param3_dpress_drho(double mdensity_in,
                                     double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    //double tau = temperature_ref3_ / temperature_in;

    double dpress_drho =
        const_R_spec_ * temperature_in *
        (2. * delta * get_param3_dphi_ddelta(mdensity_in,
                                             temperature_in) +
         delta * delta *
            get_param3_d2phi_ddelta_ddelta(mdensity_in,
                                           temperature_in));

    return dpress_drho;
}

double Lib97::get_param3_erg_int(double mdensity_in,
                                 double temperature_in) {
    //double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double erg_int =
        const_R_spec_ * temperature_in * tau *
        get_param3_dphi_dtau(mdensity_in, temperature_in);

    return erg_int;
}

double Lib97::get_param3_entropy(double mdensity_in,
                                 double temperature_in) {
    //double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double entropy =
        const_R_spec_ *
        (tau * get_param3_dphi_dtau(mdensity_in,
                                    temperature_in) -
         get_param3_phi(mdensity_in, temperature_in));

    return entropy;
}

double Lib97::get_param3_enthalpy(double mdensity_in,
                                  double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double enthalpy =
        const_R_spec_ * temperature_in *
        (tau * get_param3_dphi_dtau(mdensity_in,
                                    temperature_in) +
         delta * get_param3_dphi_ddelta(mdensity_in,
                                        temperature_in));

    return enthalpy;
}

double Lib97::get_param3_heat_c_v(double mdensity_in,
                                  double temperature_in) {
    //double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double c_v =
        -const_R_spec_ * tau * tau *
        get_param3_d2phi_dtau_dtau(mdensity_in, temperature_in);

    return c_v;
}

double Lib97::get_param3_heat_c_p(double mdensity_in,
                                  double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double dphi_ddelta =
        get_param3_dphi_ddelta(mdensity_in, temperature_in);

    double fac_nom =
        delta * dphi_ddelta -
        delta * tau * get_param3_d2phi_ddelta_dtau(mdensity_in,
                                                   temperature_in);
    double fac_den =
        2. * delta * dphi_ddelta +
        delta * delta * get_param3_d2phi_ddelta_ddelta(mdensity_in,
                                                       temperature_in);

    double c_p =
        const_R_spec_ *
        (fac_nom * fac_nom / fac_den -
         tau * tau * get_param3_d2phi_dtau_dtau(mdensity_in,
                                                temperature_in));

    return c_p;
}

double Lib97::get_param3_speed_sound(double mdensity_in,
                                     double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double dphi_ddelta =
        get_param3_dphi_ddelta(mdensity_in, temperature_in);

    double fac_nom =
        delta * dphi_ddelta -
        delta * tau * get_param3_d2phi_ddelta_dtau(mdensity_in,
                                                   temperature_in);
    double fac_den =
        tau * tau * get_param3_d2phi_dtau_dtau(mdensity_in,
                                               temperature_in);

    double v2_s =
        const_R_spec_ * temperature_in *
        (2. * delta * dphi_ddelta +
         delta * delta * get_param3_d2phi_ddelta_ddelta(mdensity_in,
                                                        temperature_in) -
         fac_nom * fac_nom / fac_den);

    return sqrt(v2_s);
}

double Lib97::get_param3_phi(double mdensity_in,
                             double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double phi = coeff3_n_[1] * log(delta);
    for (int i = 2; i <= 40; i++) {
        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i]);

        phi += coeff3_n_[i] *
            delta_pow_I * tau_pow_J;
    }

    return phi;
}

double Lib97::get_param3_dphi_ddelta(double mdensity_in,
                                     double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double dphi_ddelta = coeff3_n_[1] / delta;
    for (int i = 2; i <= 40; i++) {
        if (coeff3_I_[i] == 0) {
            continue;
        }

        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i] - 1);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i]);

        dphi_ddelta +=
            coeff3_n_[i] * delta_pow_I * tau_pow_J *
            static_cast<double>(coeff3_I_[i]);
    }

    return dphi_ddelta;
}

double Lib97::get_param3_dphi_dtau(double mdensity_in,
                                   double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double dphi_dtau = 0.;
    for (int i = 2; i <= 40; i++) {
        if (coeff3_J_[i] == 0) {
            continue;
        }

        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i] - 1);

        dphi_dtau +=
            coeff3_n_[i] * delta_pow_I * tau_pow_J *
            static_cast<double>(coeff3_J_[i]);
    }

    return dphi_dtau;
}

double Lib97::get_param3_d2phi_ddelta_ddelta(double mdensity_in,
                                             double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double d2phi_ddelta_ddelta =
        -coeff3_n_[1] / (delta * delta);
    for (int i = 2; i <= 40; i++) {
        if (coeff3_I_[i] == 0 ||
            coeff3_I_[i] == 1) {
            continue;
        }

        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i] - 2);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i]);

        d2phi_ddelta_ddelta +=
            coeff3_n_[i] * delta_pow_I * tau_pow_J *
            static_cast<double>(coeff3_I_[i] *
                                (coeff3_I_[i] - 1));
    }

    return d2phi_ddelta_ddelta;
}

double Lib97::get_param3_d2phi_ddelta_dtau(double mdensity_in,
                                           double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double d2phi_ddelta_dtau = 0.;
    for (int i = 2; i <= 40; i++) {
        if (coeff3_I_[i] == 0 ||
            coeff3_J_[i] == 0) {
            continue;
        }

        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i] - 1);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i] - 1);

        d2phi_ddelta_dtau +=
            coeff3_n_[i] * delta_pow_I * tau_pow_J *
            static_cast<double>(coeff3_I_[i] *
                                coeff3_J_[i]);
    }

    return d2phi_ddelta_dtau;
}

double Lib97::get_param3_d2phi_dtau_dtau(double mdensity_in,
                                         double temperature_in) {
    double delta = mdensity_in / mdensity_ref3_;
    double tau = temperature_ref3_ / temperature_in;

    double d2phi_dtau_dtau = 0.;
    for (int i = 2; i <= 40; i++) {
        if (coeff3_J_[i] == 0 ||
            coeff3_J_[i] == 1) {
            continue;
        }

        double delta_pow_I =
            get_pow_int(delta, coeff3_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff3_J_[i] - 2);

        d2phi_dtau_dtau +=
            coeff3_n_[i] * delta_pow_I * tau_pow_J *
            static_cast<double>(coeff3_J_[i] *
                                (coeff3_J_[i] - 1));
    }

    return d2phi_dtau_dtau;
}

double Lib97::get_param4_sat_pressure(double temperature_in) {
    double fn_tau = temperature_in / temperature_ref4_;
    double fn_theta =
        fn_tau + coeff4_n_[9] / (fn_tau - coeff4_n_[10]);

    double fn_A =
        fn_theta * fn_theta +
        fn_theta * coeff4_n_[1] + coeff4_n_[2];
    double fn_B =
        fn_theta * fn_theta * coeff4_n_[3] +
        fn_theta * coeff4_n_[4] + coeff4_n_[5];
    double fn_C =
        fn_theta * fn_theta * coeff4_n_[6] +
        fn_theta * coeff4_n_[7] + coeff4_n_[8];

    double fn_beta =
        2. * fn_C / (sqrt(fn_B * fn_B - 4. * fn_A * fn_C) - fn_B);

    double press_sat =
        pressure_ref4_ * get_pow_int(fn_beta, 4);

    return press_sat;
}

double Lib97::get_param4_sat_temperature(double pressure_in) {
    double fn_beta = pow(pressure_in / pressure_ref4_, 0.25);

    double fn_E =
        fn_beta * fn_beta +
        fn_beta * coeff4_n_[3] + coeff4_n_[6];
    double fn_F =
        fn_beta * fn_beta * coeff4_n_[1] +
        fn_beta * coeff4_n_[4] + coeff4_n_[7];
    double fn_G =
        fn_beta * fn_beta * coeff4_n_[2] +
        fn_beta * coeff4_n_[5] + coeff4_n_[8];

    double fn_D =
        -2. * fn_G / (sqrt(fn_F * fn_F - 4. * fn_E * fn_G) + fn_F);
    double fn_tau =
        (coeff4_n_[10] + fn_D -
         sqrt((coeff4_n_[10] + fn_D) * (coeff4_n_[10] + fn_D) -
              4. * (coeff4_n_[9] + coeff4_n_[10] * fn_D))) / 2.;

    double temp_sat = temperature_ref4_ * fn_tau;

    return temp_sat;
}

double Lib97::get_param5_gamma_ide(double temperature_in,
                                   double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double gamma = log(ppi);
    for (int i = 1; i <= 6; i++) {
        double tau_pow_J =
            get_pow_int(tau, coeff5_ide_J_[i]);
        gamma +=
            coeff5_ide_n_[i] * tau_pow_J;
    }

    return gamma;
}

double Lib97::get_param5_dgamma_ide_dppi(double temperature_in,
                                         double pressure_in) {
    //double tau = temperature_ref2_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double dgamma_dppi = 1. / ppi;

    return dgamma_dppi;
}

double Lib97::get_param5_dgamma_ide_dtau(double temperature_in,
                                         double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    //double ppi = pressure_in / pressure_ref5_;

    double dgamma_dtau = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_ide_J_[i] == 0) {
            continue;
        }

        double tau_pow_J =
            get_pow_int(tau, coeff5_ide_J_[i] - 1);
        dgamma_dtau +=
            coeff5_ide_n_[i] * tau_pow_J *
            static_cast<double>(coeff5_ide_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param5_d2gamma_ide_dppi_dppi(double temperature_in,
                                               double pressure_in) {
    //double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double d2gamma_dppi_dppi = -1. / (ppi * ppi);

    return d2gamma_dppi_dppi;
}

double Lib97::get_param5_d2gamma_ide_dppi_dtau(double temperature_in,
                                               double pressure_in) {
    //double tau = temperature_ref5_ / temperature_in;
    //double ppi = pressure_in / pressure_ref5_;

    return 0.;
}

double Lib97::get_param5_d2gamma_ide_dtau_dtau(double temperature_in,
                                               double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    //double ppi = pressure_in / pressure_ref5_;

    double d2gamma_dtau_dtau = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_ide_J_[i] == 0 ||
            coeff5_ide_J_[i] == 1) {
            continue;
        }

        double tau_pow_J =
            get_pow_int(tau, coeff5_ide_J_[i] - 2);
        d2gamma_dtau_dtau +=
            coeff5_ide_n_[i] * tau_pow_J *
            static_cast<double>(coeff5_ide_J_[i] *
                                (coeff5_ide_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

double Lib97::get_param5_gamma_res(double temperature_in,
                                   double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double gamma = 0.;
    for (int i = 1; i <= 6; i++) {
        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i]);

        gamma +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J;
    }

    return gamma;
}

double Lib97::get_param5_dgamma_res_dppi(double temperature_in,
                                         double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double dgamma_dppi = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_res_I_[i] == 0) {
            continue;
        }

        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i] - 1);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i]);

        dgamma_dppi +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J *
            static_cast<double>(coeff5_res_I_[i]);
    }

    return dgamma_dppi;
}

double Lib97::get_param5_dgamma_res_dtau(double temperature_in,
                                         double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double dgamma_dtau = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_res_J_[i] == 0) {
            continue;
        }

        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i] - 1);

        dgamma_dtau +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J *
            static_cast<double>(coeff5_res_J_[i]);
    }

    return dgamma_dtau;
}

double Lib97::get_param5_d2gamma_res_dppi_dppi(double temperature_in,
                                               double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double d2gamma_dppi_dppi = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_res_I_[i] == 0 ||
            coeff5_res_I_[i] == 1) {
            continue;
        }

        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i] - 2);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i]);

        d2gamma_dppi_dppi +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J *
            static_cast<double>(coeff5_res_I_[i] *
                                (coeff5_res_I_[i] - 1));
    }

    return d2gamma_dppi_dppi;
}

double Lib97::get_param5_d2gamma_res_dppi_dtau(double temperature_in,
                                               double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double d2gamma_dppi_dtau = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_res_I_[i] == 0 ||
            coeff5_res_J_[i] == 0) {
            continue;
        }

        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i] - 1);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i] - 1);

        d2gamma_dppi_dtau +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J *
            static_cast<double>(coeff5_res_I_[i] *
                                coeff5_res_J_[i]);
    }

    return d2gamma_dppi_dtau;
}

double Lib97::get_param5_d2gamma_res_dtau_dtau(double temperature_in,
                                               double pressure_in) {
    double tau = temperature_ref5_ / temperature_in;
    double ppi = pressure_in / pressure_ref5_;

    double d2gamma_dtau_dtau = 0.;
    for (int i = 1; i <= 6; i++) {
        if (coeff5_res_J_[i] == 0 ||
            coeff5_res_J_[i] == 1) {
            continue;
        }

        double ppi_pow_I =
            get_pow_int(ppi, coeff5_res_I_[i]);
        double tau_pow_J =
            get_pow_int(tau, coeff5_res_J_[i] - 2);

        d2gamma_dtau_dtau +=
            coeff5_res_n_[i] * ppi_pow_I * tau_pow_J *
            static_cast<double>(coeff5_res_J_[i] *
                                (coeff5_res_J_[i] - 1));
    }

    return d2gamma_dtau_dtau;
}

void Lib97::set_coefficients() {
    temperature_crit_ = 647.096;
    mdensity_crit_ = 322.;
    pressure_crit_ = 22.064 * 1.0e+6;
    const_R_spec_ = 461.526;

    temperature_refB23_ = 1.;
    pressure_refB23_ = 1.0e+6;

    coeffB23_n_[0] = 0.;
    // coefficient boundary between regions 2 and 3 n_i
    coeffB23_n_[1] = 0.34805185628969 * 1.0e+3;
    coeffB23_n_[2] = -0.11671859879975 * 1.0e+1;
    coeffB23_n_[3] = 0.10192970039326 * 1.0e-2;
    coeffB23_n_[4] = 0.57254459862746 * 1.0e+3;
    coeffB23_n_[5] = 0.13918839778870 * 1.0e+2;

    temperature_ref1_ = 1386.;
    pressure_ref1_ = 16.53 * 1.0e+6;
    tau_ref1_ = 1.222;
    ppi_ref1_ = 7.1;

    coeff1_I_[0] = 0;
    // coefficient region 1 I_i
    coeff1_I_[1] = 0;
    coeff1_I_[2] = 0;
    coeff1_I_[3] = 0;
    coeff1_I_[4] = 0;
    coeff1_I_[5] = 0;
    coeff1_I_[6] = 0;
    coeff1_I_[7] = 0;
    coeff1_I_[8] = 0;
    coeff1_I_[9] = 1;
    coeff1_I_[10] = 1;
    coeff1_I_[11] = 1;
    coeff1_I_[12] = 1;
    coeff1_I_[13] = 1;
    coeff1_I_[14] = 1;
    coeff1_I_[15] = 2;
    coeff1_I_[16] = 2;
    coeff1_I_[17] = 2;
    coeff1_I_[18] = 2;
    coeff1_I_[19] = 2;
    coeff1_I_[20] = 3;
    coeff1_I_[21] = 3;
    coeff1_I_[22] = 3;
    coeff1_I_[23] = 4;
    coeff1_I_[24] = 4;
    coeff1_I_[25] = 4;
    coeff1_I_[26] = 5;
    coeff1_I_[27] = 8;
    coeff1_I_[28] = 8;
    coeff1_I_[29] = 21;
    coeff1_I_[30] = 23;
    coeff1_I_[31] = 29;
    coeff1_I_[32] = 30;
    coeff1_I_[33] = 31;
    coeff1_I_[34] = 32;

    coeff1_J_[0] = 0;
    // coefficient region 1 J_i
    coeff1_J_[1] = -2;
    coeff1_J_[2] = -1;
    coeff1_J_[3] = 0;
    coeff1_J_[4] = 1;
    coeff1_J_[5] = 2;
    coeff1_J_[6] = 3;
    coeff1_J_[7] = 4;
    coeff1_J_[8] = 5;
    coeff1_J_[9] = -9;
    coeff1_J_[10] = -7;
    coeff1_J_[11] = -1;
    coeff1_J_[12] = 0;
    coeff1_J_[13] = 1;
    coeff1_J_[14] = 3;
    coeff1_J_[15] = -3;
    coeff1_J_[16] = 0;
    coeff1_J_[17] = 1;
    coeff1_J_[18] = 3;
    coeff1_J_[19] = 17;
    coeff1_J_[20] = -4;
    coeff1_J_[21] = 0;
    coeff1_J_[22] = 6;
    coeff1_J_[23] = -5;
    coeff1_J_[24] = -2;
    coeff1_J_[25] = 10;
    coeff1_J_[26] = -8;
    coeff1_J_[27] = -11;
    coeff1_J_[28] = -6;
    coeff1_J_[29] = -29;
    coeff1_J_[30] = -31;
    coeff1_J_[31] = -38;
    coeff1_J_[32] = -39;
    coeff1_J_[33] = -40;
    coeff1_J_[34] = -41;

    coeff1_n_[0] = 0.;
    // coefficient region 1 n_i
    coeff1_n_[1] = 0.14632971213167;
    coeff1_n_[2] = -0.84548187169114;
    coeff1_n_[3] = -0.37563603672040 * 1.0e+1;
    coeff1_n_[4] = 0.33855169168385 * 1.0e+1;
    coeff1_n_[5] = -0.95791963387872;
    coeff1_n_[6] = 0.15772038513228;
    coeff1_n_[7] = -0.16616417199501 * 1.0e-1;
    coeff1_n_[8] = 0.81214629983568 * 1.0e-3;
    coeff1_n_[9] = 0.28319080123804 * 1.0e-3;
    coeff1_n_[10] = -0.60706301565874 * 1.0e-3;
    coeff1_n_[11] = -0.18990068218419 * 1.0e-1;
    coeff1_n_[12] = -0.32529748770505 * 1.0e-1;
    coeff1_n_[13] = -0.21841717175414 * 1.0e-1;
    coeff1_n_[14] = -0.52838357969930 * 1.0e-4;
    coeff1_n_[15] = -0.47184321073267 * 1.0e-3;
    coeff1_n_[16] = -0.30001780793026 * 1.0e-3;
    coeff1_n_[17] = 0.47661393906987 * 1.0e-4;
    coeff1_n_[18] = -0.44141845330846 * 1.0e-5;
    coeff1_n_[19] = -0.72694996297594 * 1.0e-15;
    coeff1_n_[20] = -0.31679644845054 * 1.0e-4;
    coeff1_n_[21] = -0.28270797985312 * 1.0e-5;
    coeff1_n_[22] = -0.85205128120103 * 1.0e-9;
    coeff1_n_[23] = -0.22425281908000 * 1.0e-5;
    coeff1_n_[24] = -0.65171222895601 * 1.0e-6;
    coeff1_n_[25] = -0.14341729937924 * 1.0e-12;
    coeff1_n_[26] = -0.40516996860117 * 1.0e-6;
    coeff1_n_[27] = -0.12734301741641 * 1.0e-8;
    coeff1_n_[28] = -0.17424871230634 * 1.0e-9;
    coeff1_n_[29] = -0.68762131295531 * 1.0e-18;
    coeff1_n_[30] = 0.14478307828521 * 1.0e-19;
    coeff1_n_[31] = 0.26335781662795 * 1.0e-22;
    coeff1_n_[32] = -0.11947622640071 * 1.0e-22;
    coeff1_n_[33] = 0.18228094581404 * 1.0e-23;
    coeff1_n_[34] = -0.93537087292458 * 1.0e-25;

    temperature_ref1Tph_ = 1.;
    pressure_ref1Tph_ = 1.0e+6;
    enthalpy_ref1Tph_ = 2500. * 1.0e+3;

    coeff1Tph_I_[0] = 0;
    // backward T(p,h) coefficient region 1 I_i
    coeff1Tph_I_[1] = 0;
    coeff1Tph_I_[2] = 0;
    coeff1Tph_I_[3] = 0;
    coeff1Tph_I_[4] = 0;
    coeff1Tph_I_[5] = 0;
    coeff1Tph_I_[6] = 0;
    coeff1Tph_I_[7] = 1;
    coeff1Tph_I_[8] = 1;
    coeff1Tph_I_[9] = 1;
    coeff1Tph_I_[10] = 1;
    coeff1Tph_I_[11] = 1;
    coeff1Tph_I_[12] = 1;
    coeff1Tph_I_[13] = 1;
    coeff1Tph_I_[14] = 2;
    coeff1Tph_I_[15] = 2;
    coeff1Tph_I_[16] = 3;
    coeff1Tph_I_[17] = 3;
    coeff1Tph_I_[18] = 4;
    coeff1Tph_I_[19] = 5;
    coeff1Tph_I_[20] = 6;

    coeff1Tph_J_[0] = 0;
    // backward T(p,h) coefficient region 1 J_i
    coeff1Tph_J_[1] = 0;
    coeff1Tph_J_[2] = 1;
    coeff1Tph_J_[3] = 2;
    coeff1Tph_J_[4] = 6;
    coeff1Tph_J_[5] = 22;
    coeff1Tph_J_[6] = 32;
    coeff1Tph_J_[7] = 0;
    coeff1Tph_J_[8] = 1;
    coeff1Tph_J_[9] = 2;
    coeff1Tph_J_[10] = 3;
    coeff1Tph_J_[11] = 4;
    coeff1Tph_J_[12] = 10;
    coeff1Tph_J_[13] = 32;
    coeff1Tph_J_[14] = 10;
    coeff1Tph_J_[15] = 32;
    coeff1Tph_J_[16] = 10;
    coeff1Tph_J_[17] = 32;
    coeff1Tph_J_[18] = 32;
    coeff1Tph_J_[19] = 32;
    coeff1Tph_J_[20] = 32;

    coeff1Tph_n_[0] = 0.;
    // backward T(p,h) coefficient region 1 n_i
    coeff1Tph_n_[1] = -0.23872489924521 * 1.0e+3;
    coeff1Tph_n_[2] = 0.40421188637945 * 1.0e+3;
    coeff1Tph_n_[3] = 0.11349746881718 * 1.0e+3;
    coeff1Tph_n_[4] = -0.58457616048039 * 1.0e+1;
    coeff1Tph_n_[5] = -0.15285482413140 * 1.0e-3;
    coeff1Tph_n_[6] = -0.10866707695377 * 1.0e-5;
    coeff1Tph_n_[7] = -0.13391744872602 * 1.0e+2;
    coeff1Tph_n_[8] = 0.43211039183559 * 1.0e+2;
    coeff1Tph_n_[9] = -0.54010067170506 * 1.0e+2;
    coeff1Tph_n_[10] = 0.30535892203916 * 1.0e+2;
    coeff1Tph_n_[11] = -0.65964749423638 * 1.0e+1;
    coeff1Tph_n_[12] = 0.93965400878363 * 1.0e-2;
    coeff1Tph_n_[13] = 0.11573647505340 * 1.0e-6;
    coeff1Tph_n_[14] = -0.25858641282073 * 1.0e-4;
    coeff1Tph_n_[15] = -0.40644363084799 * 1.0e-8;
    coeff1Tph_n_[16] = 0.66456186191635 * 1.0e-7;
    coeff1Tph_n_[17] = 0.80670734103027 * 1.0e-10;
    coeff1Tph_n_[18] = -0.93477771213947 * 1.0e-12;
    coeff1Tph_n_[19] = 0.58265442020601 * 1.0e-14;
    coeff1Tph_n_[20] = -0.15020185953503 * 1.0e-16;

    temperature_ref1Tps_ = 1.;
    pressure_ref1Tps_ = 1.0e+6;
    entropy_ref1Tps_ = 1.0e+3;

    coeff1Tps_I_[0] = 0;
    // backward T(p,s) coefficient region 1 I_i
    coeff1Tps_I_[1] = 0;
    coeff1Tps_I_[2] = 0;
    coeff1Tps_I_[3] = 0;
    coeff1Tps_I_[4] = 0;
    coeff1Tps_I_[5] = 0;
    coeff1Tps_I_[6] = 0;
    coeff1Tps_I_[7] = 1;
    coeff1Tps_I_[8] = 1;
    coeff1Tps_I_[9] = 1;
    coeff1Tps_I_[10] = 1;
    coeff1Tps_I_[11] = 1;
    coeff1Tps_I_[12] = 1;
    coeff1Tps_I_[13] = 2;
    coeff1Tps_I_[14] = 2;
    coeff1Tps_I_[15] = 2;
    coeff1Tps_I_[16] = 2;
    coeff1Tps_I_[17] = 2;
    coeff1Tps_I_[18] = 3;
    coeff1Tps_I_[19] = 3;
    coeff1Tps_I_[20] = 4;

    coeff1Tps_J_[0] = 0;
    // backward T(p,s) coefficient region 1 J_i
    coeff1Tps_J_[1] = 0;
    coeff1Tps_J_[2] = 1;
    coeff1Tps_J_[3] = 2;
    coeff1Tps_J_[4] = 3;
    coeff1Tps_J_[5] = 11;
    coeff1Tps_J_[6] = 31;
    coeff1Tps_J_[7] = 0;
    coeff1Tps_J_[8] = 1;
    coeff1Tps_J_[9] = 2;
    coeff1Tps_J_[10] = 3;
    coeff1Tps_J_[11] = 12;
    coeff1Tps_J_[12] = 31;
    coeff1Tps_J_[13] = 0;
    coeff1Tps_J_[14] = 1;
    coeff1Tps_J_[15] = 2;
    coeff1Tps_J_[16] = 9;
    coeff1Tps_J_[17] = 31;
    coeff1Tps_J_[18] = 10;
    coeff1Tps_J_[19] = 32;
    coeff1Tps_J_[20] = 32;

    coeff1Tps_n_[0] = 0.;
    // backward T(p,s) coefficient region 1 n_i
    coeff1Tps_n_[1] = 0.17478268058307 * 1.0e+3;
    coeff1Tps_n_[2] = 0.34806930892873 * 1.0e+2;
    coeff1Tps_n_[3] = 0.65292584978455 * 1.0e+1;
    coeff1Tps_n_[4] = 0.33039981775489;
    coeff1Tps_n_[5] = -0.19281382923196 * 1.0e-6;
    coeff1Tps_n_[6] = -0.24909197244573 * 1.0e-22;
    coeff1Tps_n_[7] = -0.26107636489332;
    coeff1Tps_n_[8] = 0.22592965981586;
    coeff1Tps_n_[9] = -0.64256463395226 * 1.0e-1;
    coeff1Tps_n_[10] = 0.78876289270526 * 1.0e-2;
    coeff1Tps_n_[11] = 0.35672110607366 * 1.0e-9;
    coeff1Tps_n_[12] = 0.17332496994895 * 1.0e-23;
    coeff1Tps_n_[13] = 0.56608900654837 * 1.0e-3;
    coeff1Tps_n_[14] = -0.32635483139717 * 1.0e-3;
    coeff1Tps_n_[15] = 0.44778286690632 * 1.0e-4;
    coeff1Tps_n_[16] = -0.51322156908507 * 1.0e-9;
    coeff1Tps_n_[17] = -0.42522657042207 * 1.0e-25;
    coeff1Tps_n_[18] = 0.26400441360689 * 1.0e-12;
    coeff1Tps_n_[19] = 0.78124600459723 * 1.0e-28;
    coeff1Tps_n_[20] = -0.30732199903668 * 1.0e-30;

    temperature_ref2_ = 540.;
    pressure_ref2_ = 1.0e+6;
    tau_ref2_res_ = 0.5;

    coeff2_ide_J_[0] = 0;
    // coefficient region 2 ideal J_i
    coeff2_ide_J_[1] = 0;
    coeff2_ide_J_[2] = 1;
    coeff2_ide_J_[3] = -5;
    coeff2_ide_J_[4] = -4;
    coeff2_ide_J_[5] = -3;
    coeff2_ide_J_[6] = -2;
    coeff2_ide_J_[7] = -1;
    coeff2_ide_J_[8] = 2;
    coeff2_ide_J_[9] = 3;

    coeff2_ide_n_[0] = 0.;
    // coefficient region 2 ideal n_i
    coeff2_ide_n_[1] = -0.96927686500217 * 1.0e+1;
    coeff2_ide_n_[2] = 0.10086655968018 * 1.0e+2;
    coeff2_ide_n_[3] = -0.56087911283020 * 1.0e-2;
    coeff2_ide_n_[4] = 0.71452738081455 * 1.0e-1;
    coeff2_ide_n_[5] = -0.40710498223928;
    coeff2_ide_n_[6] = 0.14240819171444 * 1.0e+1;
    coeff2_ide_n_[7] = -0.43839511319450 * 1.0e+1;
    coeff2_ide_n_[8] = -0.28408632460772;
    coeff2_ide_n_[9] = 0.21268463753307 * 1.0e-1;

    // coefficient region 2 ideal n_i metastable case
    coeff2mst_ide_n1_ = -0.96937268393049 * 1.0e+1;
    coeff2mst_ide_n2_ = 0.10087275970006 * 1.0e+2;

    coeff2_res_I_[0] = 0;
    // coefficient region 2 residual I_i
    coeff2_res_I_[1] = 1;
    coeff2_res_I_[2] = 1;
    coeff2_res_I_[3] = 1;
    coeff2_res_I_[4] = 1;
    coeff2_res_I_[5] = 1;
    coeff2_res_I_[6] = 2;
    coeff2_res_I_[7] = 2;
    coeff2_res_I_[8] = 2;
    coeff2_res_I_[9] = 2;
    coeff2_res_I_[10] = 2;
    coeff2_res_I_[11] = 3;
    coeff2_res_I_[12] = 3;
    coeff2_res_I_[13] = 3;
    coeff2_res_I_[14] = 3;
    coeff2_res_I_[15] = 3;
    coeff2_res_I_[16] = 4;
    coeff2_res_I_[17] = 4;
    coeff2_res_I_[18] = 4;
    coeff2_res_I_[19] = 5;
    coeff2_res_I_[20] = 6;
    coeff2_res_I_[21] = 6;
    coeff2_res_I_[22] = 6;
    coeff2_res_I_[23] = 7;
    coeff2_res_I_[24] = 7;
    coeff2_res_I_[25] = 7;
    coeff2_res_I_[26] = 8;
    coeff2_res_I_[27] = 8;
    coeff2_res_I_[28] = 9;
    coeff2_res_I_[29] = 10;
    coeff2_res_I_[30] = 10;
    coeff2_res_I_[31] = 10;
    coeff2_res_I_[32] = 16;
    coeff2_res_I_[33] = 16;
    coeff2_res_I_[34] = 18;
    coeff2_res_I_[35] = 20;
    coeff2_res_I_[36] = 20;
    coeff2_res_I_[37] = 20;
    coeff2_res_I_[38] = 21;
    coeff2_res_I_[39] = 22;
    coeff2_res_I_[40] = 23;
    coeff2_res_I_[41] = 24;
    coeff2_res_I_[42] = 24;
    coeff2_res_I_[43] = 24;

    coeff2_res_J_[0] = 0;
    // coefficient region 2 residual J_i
    coeff2_res_J_[1] = 0;
    coeff2_res_J_[2] = 1;
    coeff2_res_J_[3] = 2;
    coeff2_res_J_[4] = 3;
    coeff2_res_J_[5] = 6;
    coeff2_res_J_[6] = 1;
    coeff2_res_J_[7] = 2;
    coeff2_res_J_[8] = 4;
    coeff2_res_J_[9] = 7;
    coeff2_res_J_[10] = 36;
    coeff2_res_J_[11] = 0;
    coeff2_res_J_[12] = 1;
    coeff2_res_J_[13] = 3;
    coeff2_res_J_[14] = 6;
    coeff2_res_J_[15] = 35;
    coeff2_res_J_[16] = 1;
    coeff2_res_J_[17] = 2;
    coeff2_res_J_[18] = 3;
    coeff2_res_J_[19] = 7;
    coeff2_res_J_[20] = 3;
    coeff2_res_J_[21] = 16;
    coeff2_res_J_[22] = 35;
    coeff2_res_J_[23] = 0;
    coeff2_res_J_[24] = 11;
    coeff2_res_J_[25] = 25;
    coeff2_res_J_[26] = 8;
    coeff2_res_J_[27] = 36;
    coeff2_res_J_[28] = 13;
    coeff2_res_J_[29] = 4;
    coeff2_res_J_[30] = 10;
    coeff2_res_J_[31] = 14;
    coeff2_res_J_[32] = 29;
    coeff2_res_J_[33] = 50;
    coeff2_res_J_[34] = 57;
    coeff2_res_J_[35] = 20;
    coeff2_res_J_[36] = 35;
    coeff2_res_J_[37] = 48;
    coeff2_res_J_[38] = 21;
    coeff2_res_J_[39] = 53;
    coeff2_res_J_[40] = 39;
    coeff2_res_J_[41] = 26;
    coeff2_res_J_[42] = 40;
    coeff2_res_J_[43] = 58;

    coeff2_res_n_[0] = 0.;
    // coefficient region 2 residual n_i
    coeff2_res_n_[1] = -0.17731742473213 * 1.0e-2;
    coeff2_res_n_[2] = -0.17834862292358 * 1.0e-1;
    coeff2_res_n_[3] = -0.45996013696365 * 1.0e-1;
    coeff2_res_n_[4] = -0.57581259083432 * 1.0e-1;
    coeff2_res_n_[5] = -0.50325278727930 * 1.0e-1;
    coeff2_res_n_[6] = -0.33032641670203 * 1.0e-4;
    coeff2_res_n_[7] = -0.18948987516315 * 1.0e-3;
    coeff2_res_n_[8] = -0.39392777243355 * 1.0e-2;
    coeff2_res_n_[9] = -0.43797295650573 * 1.0e-1;
    coeff2_res_n_[10] = -0.26674547914087 * 1.0e-4;
    coeff2_res_n_[11] = 0.20481737692309 * 1.0e-7;
    coeff2_res_n_[12] = 0.43870667284435 * 1.0e-6;
    coeff2_res_n_[13] = -0.32277677238570 * 1.0e-4;
    coeff2_res_n_[14] = -0.15033924542148 * 1.0e-2;
    coeff2_res_n_[15] = -0.40668253562649 * 1.0e-1;
    coeff2_res_n_[16] = -0.78847309559367 * 1.0e-9;
    coeff2_res_n_[17] = 0.12790717852285 * 1.0e-7;
    coeff2_res_n_[18] = 0.48225372718507 * 1.0e-6;
    coeff2_res_n_[19] = 0.22922076337661 * 1.0e-5;
    coeff2_res_n_[20] = -0.16714766451061 * 1.0e-10;
    coeff2_res_n_[21] = -0.21171472321355 * 1.0e-2;
    coeff2_res_n_[22] = -0.23895741934104 * 1.0e+2;
    coeff2_res_n_[23] = -0.59059564324270 * 1.0e-17;
    coeff2_res_n_[24] = -0.12621808899101 * 1.0e-5;
    coeff2_res_n_[25] = -0.38946842435739 * 1.0e-1;
    coeff2_res_n_[26] = 0.11256211360459 * 1.0e-10;
    coeff2_res_n_[27] = -0.82311340897998 * 1.0e+1;
    coeff2_res_n_[28] = 0.19809712802088 * 1.0e-7;
    coeff2_res_n_[29] = 0.10406965210174 * 1.0e-18;
    coeff2_res_n_[30] = -0.10234747095929 * 1.0e-12;
    coeff2_res_n_[31] = -0.10018179379511 * 1.0e-8;
    coeff2_res_n_[32] = -0.80882908646985 * 1.0e-10;
    coeff2_res_n_[33] = 0.10693031879409;
    coeff2_res_n_[34] = -0.33662250574171;
    coeff2_res_n_[35] = 0.89185845355421 * 1.0e-24;
    coeff2_res_n_[36] = 0.30629316876232 * 1.0e-12;
    coeff2_res_n_[37] = -0.42002467698208 * 1.0e-5;
    coeff2_res_n_[38] = -0.59056029685639 * 1.0e-25;
    coeff2_res_n_[39] = 0.37826947613457 * 1.0e-5;
    coeff2_res_n_[40] = -0.12768608934681 * 1.0e-14;
    coeff2_res_n_[41] = 0.73087610595061 * 1.0e-28;
    coeff2_res_n_[42] = 0.55414715350778 * 1.0e-16;
    coeff2_res_n_[43] = -0.94369707241210 * 1.0e-6;

    temperature_ref2mst_ = 540.;
    pressure_ref2mst_ = 1.0e+6;
    tau_ref2mst_res_ = 0.5;

    coeff2mst_res_I_[0] = 0;
    // coefficient region 2 metastable residual I_i
    coeff2mst_res_I_[1] = 1;
    coeff2mst_res_I_[2] = 1;
    coeff2mst_res_I_[3] = 1;
    coeff2mst_res_I_[4] = 1;
    coeff2mst_res_I_[5] = 2;
    coeff2mst_res_I_[6] = 2;
    coeff2mst_res_I_[7] = 2;
    coeff2mst_res_I_[8] = 3;
    coeff2mst_res_I_[9] = 3;
    coeff2mst_res_I_[10] = 4;
    coeff2mst_res_I_[11] = 4;
    coeff2mst_res_I_[12] = 5;
    coeff2mst_res_I_[13] = 5;

    coeff2mst_res_J_[0] = 0;
    // coefficient region 2 metastable residual J_i
    coeff2mst_res_J_[1] = 0;
    coeff2mst_res_J_[2] = 2;
    coeff2mst_res_J_[3] = 5;
    coeff2mst_res_J_[4] = 11;
    coeff2mst_res_J_[5] = 1;
    coeff2mst_res_J_[6] = 7;
    coeff2mst_res_J_[7] = 16;
    coeff2mst_res_J_[8] = 4;
    coeff2mst_res_J_[9] = 16;
    coeff2mst_res_J_[10] = 7;
    coeff2mst_res_J_[11] = 10;
    coeff2mst_res_J_[12] = 9;
    coeff2mst_res_J_[13] = 10;

    coeff2mst_res_n_[0] = 0.;
    // coefficient region 2 metastable residual n_i
    coeff2mst_res_n_[1] = -0.73362260186506 * 1.0e-2;
    coeff2mst_res_n_[2] = -0.88223831943146 * 1.0e-1;
    coeff2mst_res_n_[3] = -0.72334555213245 * 1.0e-1;
    coeff2mst_res_n_[4] = -0.40813178534455 * 1.0e-2;
    coeff2mst_res_n_[5] = 0.20097803380207 * 1.0e-2;
    coeff2mst_res_n_[6] = -0.53045921898642 * 1.0e-1;
    coeff2mst_res_n_[7] = -0.76190409086970 * 1.0e-2;
    coeff2mst_res_n_[8] = -0.63498037657313 * 1.0e-2;
    coeff2mst_res_n_[9] = -0.86043093028588 * 1.0e-1;
    coeff2mst_res_n_[10] = 0.75321581522770 * 1.0e-2;
    coeff2mst_res_n_[11] = -0.79238375446139 * 1.0e-2;
    coeff2mst_res_n_[12] = -0.22888160778447 * 1.0e-3;
    coeff2mst_res_n_[13] = -0.26456501482810 * 1.0e-2;

    pressure_refB2bc_ = 1.0e+6;
    enthalpy_refB2bc_ = 1.0e+3;

    coeffB2bc_n_[0] = 0.;
    // coefficient boundary between regions 2b and 2c n_i
    coeffB2bc_n_[1] = 0.90584278514723 * 1.0e+3;
    coeffB2bc_n_[2] = -0.67955786399241;
    coeffB2bc_n_[3] = 0.12809002730136 * 1.0e-3;
    coeffB2bc_n_[4] = 0.26526571908428 * 1.0e+4;
    coeffB2bc_n_[5] = 0.45257578905948 * 1.0e+1;

    temperature_ref2aTph_ = 1.;
    pressure_ref2aTph_ = 1.0e+6;
    enthalpy_ref2aTph_ = 2. * 1.0e+6;

    coeff2aTph_I_[0] = 0;
    // backward T(p,h) coefficient region 2a I_i
    coeff2aTph_I_[1] = 0;
    coeff2aTph_I_[2] = 0;
    coeff2aTph_I_[3] = 0;
    coeff2aTph_I_[4] = 0;
    coeff2aTph_I_[5] = 0;
    coeff2aTph_I_[6] = 0;
    coeff2aTph_I_[7] = 1;
    coeff2aTph_I_[8] = 1;
    coeff2aTph_I_[9] = 1;
    coeff2aTph_I_[10] = 1;
    coeff2aTph_I_[11] = 1;
    coeff2aTph_I_[12] = 1;
    coeff2aTph_I_[13] = 1;
    coeff2aTph_I_[14] = 1;
    coeff2aTph_I_[15] = 1;
    coeff2aTph_I_[16] = 2;
    coeff2aTph_I_[17] = 2;
    coeff2aTph_I_[18] = 2;
    coeff2aTph_I_[19] = 2;
    coeff2aTph_I_[20] = 2;
    coeff2aTph_I_[21] = 2;
    coeff2aTph_I_[22] = 2;
    coeff2aTph_I_[23] = 2;
    coeff2aTph_I_[24] = 3;
    coeff2aTph_I_[25] = 3;
    coeff2aTph_I_[26] = 4;
    coeff2aTph_I_[27] = 4;
    coeff2aTph_I_[28] = 4;
    coeff2aTph_I_[29] = 5;
    coeff2aTph_I_[30] = 5;
    coeff2aTph_I_[31] = 5;
    coeff2aTph_I_[32] = 6;
    coeff2aTph_I_[33] = 6;
    coeff2aTph_I_[34] = 7;

    coeff2aTph_J_[0] = 0;
    // backward T(p,h) coefficient region 2a J_i
    coeff2aTph_J_[1] = 0;
    coeff2aTph_J_[2] = 1;
    coeff2aTph_J_[3] = 2;
    coeff2aTph_J_[4] = 3;
    coeff2aTph_J_[5] = 7;
    coeff2aTph_J_[6] = 20;
    coeff2aTph_J_[7] = 0;
    coeff2aTph_J_[8] = 1;
    coeff2aTph_J_[9] = 2;
    coeff2aTph_J_[10] = 3;
    coeff2aTph_J_[11] = 7;
    coeff2aTph_J_[12] = 9;
    coeff2aTph_J_[13] = 11;
    coeff2aTph_J_[14] = 18;
    coeff2aTph_J_[15] = 44;
    coeff2aTph_J_[16] = 0;
    coeff2aTph_J_[17] = 2;
    coeff2aTph_J_[18] = 7;
    coeff2aTph_J_[19] = 36;
    coeff2aTph_J_[20] = 38;
    coeff2aTph_J_[21] = 40;
    coeff2aTph_J_[22] = 42;
    coeff2aTph_J_[23] = 44;
    coeff2aTph_J_[24] = 24;
    coeff2aTph_J_[25] = 44;
    coeff2aTph_J_[26] = 12;
    coeff2aTph_J_[27] = 32;
    coeff2aTph_J_[28] = 44;
    coeff2aTph_J_[29] = 32;
    coeff2aTph_J_[30] = 36;
    coeff2aTph_J_[31] = 42;
    coeff2aTph_J_[32] = 34;
    coeff2aTph_J_[33] = 44;
    coeff2aTph_J_[34] = 28;

    coeff2aTph_n_[0] = 0.;
    // backward T(p,h) coefficient region 2a n_i
    coeff2aTph_n_[1] = 0.10898952318288 * 1.0e+4;
    coeff2aTph_n_[2] = 0.84951654495535 * 1.0e+3;
    coeff2aTph_n_[3] = -0.10781748091826 * 1.0e+3;
    coeff2aTph_n_[4] = 0.33153654801263 * 1.0e+2;
    coeff2aTph_n_[5] = -0.74232016790248 * 1.0e+1;
    coeff2aTph_n_[6] = 0.11765048724356 * 1.0e+2;
    coeff2aTph_n_[7] = 0.18445749355790 * 1.0e+1;
    coeff2aTph_n_[8] = -0.41792700549624 * 1.0e+1;
    coeff2aTph_n_[9] = 0.62478196935812 * 1.0e+1;
    coeff2aTph_n_[10] = -0.17344563108114 * 1.0e+2;
    coeff2aTph_n_[11] = -0.20058176862096 * 1.0e+3;
    coeff2aTph_n_[12] = 0.27196065473796 * 1.0e+3;
    coeff2aTph_n_[13] = -0.45511318285818 * 1.0e+3;
    coeff2aTph_n_[14] = 0.30919688604755 * 1.0e+4;
    coeff2aTph_n_[15] = 0.25226640357872 * 1.0e+6;
    coeff2aTph_n_[16] = -0.61707422868339 * 1.0e-2;
    coeff2aTph_n_[17] = -0.31078046629583;
    coeff2aTph_n_[18] = 0.11670873077107 * 1.0e+2;
    coeff2aTph_n_[19] = 0.12812798404046 * 1.0e+9;
    coeff2aTph_n_[20] = -0.98554909623276 * 1.0e+9;
    coeff2aTph_n_[21] = 0.28224546973002 * 1.0e+10;
    coeff2aTph_n_[22] = -0.35948971410703 * 1.0e+10;
    coeff2aTph_n_[23] = 0.17227349913197 * 1.0e+10;
    coeff2aTph_n_[24] = -0.13551334240775 * 1.0e+5;
    coeff2aTph_n_[25] = 0.12848734664650 * 1.0e+8;
    coeff2aTph_n_[26] = 0.13865724283226 * 1.0e+1;
    coeff2aTph_n_[27] = 0.23598832556514 * 1.0e+6;
    coeff2aTph_n_[28] = -0.13105236545054 * 1.0e+8;
    coeff2aTph_n_[29] = 0.73999835474766 * 1.0e+4;
    coeff2aTph_n_[30] = -0.55196697030060 * 1.0e+6;
    coeff2aTph_n_[31] = 0.37154085996233 * 1.0e+7;
    coeff2aTph_n_[32] = 0.19127729239660 * 1.0e+5;
    coeff2aTph_n_[33] = -0.41535164835634 * 1.0e+6;
    coeff2aTph_n_[34] = -0.62459855192507 * 1.0e+2;

    temperature_ref2bTph_ = 1.;
    pressure_ref2bTph_ = 1.0e+6;
    enthalpy_ref2bTph_ = 2. * 1.0e+6;

    coeff2bTph_I_[0] = 0;
    // backward T(p,h) coefficient region 2b I_i
    coeff2bTph_I_[1] = 0;
    coeff2bTph_I_[2] = 0;
    coeff2bTph_I_[3] = 0;
    coeff2bTph_I_[4] = 0;
    coeff2bTph_I_[5] = 0;
    coeff2bTph_I_[6] = 0;
    coeff2bTph_I_[7] = 0;
    coeff2bTph_I_[8] = 0;
    coeff2bTph_I_[9] = 1;
    coeff2bTph_I_[10] = 1;
    coeff2bTph_I_[11] = 1;
    coeff2bTph_I_[12] = 1;
    coeff2bTph_I_[13] = 1;
    coeff2bTph_I_[14] = 1;
    coeff2bTph_I_[15] = 1;
    coeff2bTph_I_[16] = 1;
    coeff2bTph_I_[17] = 2;
    coeff2bTph_I_[18] = 2;
    coeff2bTph_I_[19] = 2;
    coeff2bTph_I_[20] = 2;
    coeff2bTph_I_[21] = 3;
    coeff2bTph_I_[22] = 3;
    coeff2bTph_I_[23] = 3;
    coeff2bTph_I_[24] = 3;
    coeff2bTph_I_[25] = 4;
    coeff2bTph_I_[26] = 4;
    coeff2bTph_I_[27] = 4;
    coeff2bTph_I_[28] = 4;
    coeff2bTph_I_[29] = 4;
    coeff2bTph_I_[30] = 4;
    coeff2bTph_I_[31] = 5;
    coeff2bTph_I_[32] = 5;
    coeff2bTph_I_[33] = 5;
    coeff2bTph_I_[34] = 6;
    coeff2bTph_I_[35] = 7;
    coeff2bTph_I_[36] = 7;
    coeff2bTph_I_[37] = 9;
    coeff2bTph_I_[38] = 9;

    coeff2bTph_J_[0] = 0;
    // backward T(p,h) coefficient region 2b J_i
    coeff2bTph_J_[1] = 0;
    coeff2bTph_J_[2] = 1;
    coeff2bTph_J_[3] = 2;
    coeff2bTph_J_[4] = 12;
    coeff2bTph_J_[5] = 18;
    coeff2bTph_J_[6] = 24;
    coeff2bTph_J_[7] = 28;
    coeff2bTph_J_[8] = 40;
    coeff2bTph_J_[9] = 0;
    coeff2bTph_J_[10] = 2;
    coeff2bTph_J_[11] = 6;
    coeff2bTph_J_[12] = 12;
    coeff2bTph_J_[13] = 18;
    coeff2bTph_J_[14] = 24;
    coeff2bTph_J_[15] = 28;
    coeff2bTph_J_[16] = 40;
    coeff2bTph_J_[17] = 2;
    coeff2bTph_J_[18] = 8;
    coeff2bTph_J_[19] = 18;
    coeff2bTph_J_[20] = 40;
    coeff2bTph_J_[21] = 1;
    coeff2bTph_J_[22] = 2;
    coeff2bTph_J_[23] = 12;
    coeff2bTph_J_[24] = 24;
    coeff2bTph_J_[25] = 2;
    coeff2bTph_J_[26] = 12;
    coeff2bTph_J_[27] = 18;
    coeff2bTph_J_[28] = 24;
    coeff2bTph_J_[29] = 28;
    coeff2bTph_J_[30] = 40;
    coeff2bTph_J_[31] = 18;
    coeff2bTph_J_[32] = 24;
    coeff2bTph_J_[33] = 40;
    coeff2bTph_J_[34] = 28;
    coeff2bTph_J_[35] = 2;
    coeff2bTph_J_[36] = 28;
    coeff2bTph_J_[37] = 1;
    coeff2bTph_J_[38] = 40;

    coeff2bTph_n_[0] = 0.;
    // backward T(p,h) coefficient region 2b n_i
    coeff2bTph_n_[1] = 0.14895041079516 * 1.0e+4;
    coeff2bTph_n_[2] = 0.74307798314034 * 1.0e+3;
    coeff2bTph_n_[3] = -0.97708318797837 * 1.0e+2;
    coeff2bTph_n_[4] = 0.24742464705674 * 1.0e+1;
    coeff2bTph_n_[5] = -0.63281320016026;
    coeff2bTph_n_[6] = 0.11385952129658 * 1.0e+1;
    coeff2bTph_n_[7] = -0.47811863648625;
    coeff2bTph_n_[8] = 0.85208123431544 * 1.0e-2;
    coeff2bTph_n_[9] = 0.93747147377932;
    coeff2bTph_n_[10] = 0.33593118604916 * 1.0e+1;
    coeff2bTph_n_[11] = 0.33809355601454 * 1.0e+1;
    coeff2bTph_n_[12] = 0.16844539671904;
    coeff2bTph_n_[13] = 0.73875745236695;
    coeff2bTph_n_[14] = -0.47128737436186;
    coeff2bTph_n_[15] = 0.15020273139707;
    coeff2bTph_n_[16] = -0.21764114219750 * 1.0e-2;
    coeff2bTph_n_[17] = -0.21810755324761 * 1.0e-1;
    coeff2bTph_n_[18] = -0.10829784403677;
    coeff2bTph_n_[19] = -0.46333324635812 * 1.0e-1;
    coeff2bTph_n_[20] = 0.71280351959551 * 1.0e-4;
    coeff2bTph_n_[21] = 0.11032831789999 * 1.0e-3;
    coeff2bTph_n_[22] = 0.18955248387902 * 1.0e-3;
    coeff2bTph_n_[23] = 0.30891541160537 * 1.0e-2;
    coeff2bTph_n_[24] = 0.13555504554949 * 1.0e-2;
    coeff2bTph_n_[25] = 0.28640237477456 * 1.0e-6;
    coeff2bTph_n_[26] = -0.10779857357512 * 1.0e-4;
    coeff2bTph_n_[27] = -0.76462712454814 * 1.0e-4;
    coeff2bTph_n_[28] = 0.14052392818316 * 1.0e-4;
    coeff2bTph_n_[29] = -0.31083814331434 * 1.0e-4;
    coeff2bTph_n_[30] = -0.10302738212103 * 1.0e-5;
    coeff2bTph_n_[31] = 0.28217281635040 * 1.0e-6;
    coeff2bTph_n_[32] = 0.12704902271945 * 1.0e-5;
    coeff2bTph_n_[33] = 0.73803353468292 * 1.0e-7;
    coeff2bTph_n_[34] = -0.11030139238909 * 1.0e-7;
    coeff2bTph_n_[35] = -0.81456365207833 * 1.0e-13;
    coeff2bTph_n_[36] = -0.25180545682962 * 1.0e-10;
    coeff2bTph_n_[37] = -0.17565233969407 * 1.0e-17;
    coeff2bTph_n_[38] = 0.86934156344163 * 1.0e-14;

    temperature_ref2cTph_ = 1.;
    pressure_ref2cTph_ = 1.0e+6;
    enthalpy_ref2cTph_ = 2. * 1.0e+6;

    coeff2cTph_I_[0] = 0;
    // backward T(p,h) coefficient region 2c I_i
    coeff2cTph_I_[1] = -7;
    coeff2cTph_I_[2] = -7;
    coeff2cTph_I_[3] = -6;
    coeff2cTph_I_[4] = -6;
    coeff2cTph_I_[5] = -5;
    coeff2cTph_I_[6] = -5;
    coeff2cTph_I_[7] = -2;
    coeff2cTph_I_[8] = -2;
    coeff2cTph_I_[9] = -1;
    coeff2cTph_I_[10] = -1;
    coeff2cTph_I_[11] = 0;
    coeff2cTph_I_[12] = 0;
    coeff2cTph_I_[13] = 1;
    coeff2cTph_I_[14] = 1;
    coeff2cTph_I_[15] = 2;
    coeff2cTph_I_[16] = 6;
    coeff2cTph_I_[17] = 6;
    coeff2cTph_I_[18] = 6;
    coeff2cTph_I_[19] = 6;
    coeff2cTph_I_[20] = 6;
    coeff2cTph_I_[21] = 6;
    coeff2cTph_I_[22] = 6;
    coeff2cTph_I_[23] = 6;

    coeff2cTph_J_[0] = 0;
    // backward T(p,h) coefficient region 2c J_i
    coeff2cTph_J_[1] = 0;
    coeff2cTph_J_[2] = 4;
    coeff2cTph_J_[3] = 0;
    coeff2cTph_J_[4] = 2;
    coeff2cTph_J_[5] = 0;
    coeff2cTph_J_[6] = 2;
    coeff2cTph_J_[7] = 0;
    coeff2cTph_J_[8] = 1;
    coeff2cTph_J_[9] = 0;
    coeff2cTph_J_[10] = 2;
    coeff2cTph_J_[11] = 0;
    coeff2cTph_J_[12] = 1;
    coeff2cTph_J_[13] = 4;
    coeff2cTph_J_[14] = 8;
    coeff2cTph_J_[15] = 4;
    coeff2cTph_J_[16] = 0;
    coeff2cTph_J_[17] = 1;
    coeff2cTph_J_[18] = 4;
    coeff2cTph_J_[19] = 10;
    coeff2cTph_J_[20] = 12;
    coeff2cTph_J_[21] = 16;
    coeff2cTph_J_[22] = 20;
    coeff2cTph_J_[23] = 22;

    coeff2cTph_n_[0] = 0.;
    // backward T(p,h) coefficient region 2c n_i
    coeff2cTph_n_[1] = -0.32368398555242 * 1.0e+13;
    coeff2cTph_n_[2] = 0.73263350902181 * 1.0e+13;
    coeff2cTph_n_[3] = 0.35825089945447 * 1.0e+12;
    coeff2cTph_n_[4] = -0.58340131851590 * 1.0e+12;
    coeff2cTph_n_[5] = -0.10783068217470 * 1.0e+11;
    coeff2cTph_n_[6] = 0.20825544563171 * 1.0e+11;
    coeff2cTph_n_[7] = 0.61074783564516 * 1.0e+6;
    coeff2cTph_n_[8] = 0.85977722535580 * 1.0e+6;
    coeff2cTph_n_[9] = -0.25745723604170 * 1.0e+5;
    coeff2cTph_n_[10] = 0.31081088422714 * 1.0e+5;
    coeff2cTph_n_[11] = 0.12082315865936 * 1.0e+4;
    coeff2cTph_n_[12] = 0.48219755109255 * 1.0e+3;
    coeff2cTph_n_[13] = 0.37966001272486 * 1.0e+1;
    coeff2cTph_n_[14] = -0.10842984880077 * 1.0e+2;
    coeff2cTph_n_[15] = -0.45364172676660 * 1.0e-1;
    coeff2cTph_n_[16] = 0.14559115658698 * 1.0e-12;
    coeff2cTph_n_[17] = 0.11261597407230 * 1.0e-11;
    coeff2cTph_n_[18] = -0.17804982240686 * 1.0e-10;
    coeff2cTph_n_[19] = 0.12324579690832 * 1.0e-6;
    coeff2cTph_n_[20] = -0.11606921130984 * 1.0e-5;
    coeff2cTph_n_[21] = 0.27846367088554 * 1.0e-4;
    coeff2cTph_n_[22] = -0.59270038474176 * 1.0e-3;
    coeff2cTph_n_[23] = 0.12918582991878 * 1.0e-2;

    temperature_ref2aTps_ = 1.;
    pressure_ref2aTps_ = 1.0e+6;
    entropy_ref2aTps_ = 2. * 1.0e+3;

    coeff2aTps_I_[0] = 0.;
    // backward T(p,s) coefficient region 2a I_i
    coeff2aTps_I_[1] = -1.5;
    coeff2aTps_I_[2] = -1.5;
    coeff2aTps_I_[3] = -1.5;
    coeff2aTps_I_[4] = -1.5;
    coeff2aTps_I_[5] = -1.5;
    coeff2aTps_I_[6] = -1.5;
    coeff2aTps_I_[7] = -1.25;
    coeff2aTps_I_[8] = -1.25;
    coeff2aTps_I_[9] = -1.25;
    coeff2aTps_I_[10] = -1.;
    coeff2aTps_I_[11] = -1.;
    coeff2aTps_I_[12] = -1.;
    coeff2aTps_I_[13] = -1.;
    coeff2aTps_I_[14] = -1.;
    coeff2aTps_I_[15] = -1.;
    coeff2aTps_I_[16] = -0.75;
    coeff2aTps_I_[17] = -0.75;
    coeff2aTps_I_[18] = -0.5;
    coeff2aTps_I_[19] = -0.5;
    coeff2aTps_I_[20] = -0.5;
    coeff2aTps_I_[21] = -0.5;
    coeff2aTps_I_[22] = -0.25;
    coeff2aTps_I_[23] = -0.25;
    coeff2aTps_I_[24] = -0.25;
    coeff2aTps_I_[25] = -0.25;
    coeff2aTps_I_[26] = 0.25;
    coeff2aTps_I_[27] = 0.25;
    coeff2aTps_I_[28] = 0.25;
    coeff2aTps_I_[29] = 0.25;
    coeff2aTps_I_[30] = 0.5;
    coeff2aTps_I_[31] = 0.5;
    coeff2aTps_I_[32] = 0.5;
    coeff2aTps_I_[33] = 0.5;
    coeff2aTps_I_[34] = 0.5;
    coeff2aTps_I_[35] = 0.5;
    coeff2aTps_I_[36] = 0.5;
    coeff2aTps_I_[37] = 0.75;
    coeff2aTps_I_[38] = 0.75;
    coeff2aTps_I_[39] = 0.75;
    coeff2aTps_I_[40] = 0.75;
    coeff2aTps_I_[41] = 1.;
    coeff2aTps_I_[42] = 1.;
    coeff2aTps_I_[43] = 1.25;
    coeff2aTps_I_[44] = 1.25;
    coeff2aTps_I_[45] = 1.5;
    coeff2aTps_I_[46] = 1.5;

    coeff2aTps_J_[0] = 0;
    // backward T(p,s) coefficient region 2a J_i
    coeff2aTps_J_[1] = -24;
    coeff2aTps_J_[2] = -23;
    coeff2aTps_J_[3] = -19;
    coeff2aTps_J_[4] = -13;
    coeff2aTps_J_[5] = -11;
    coeff2aTps_J_[6] = -10;
    coeff2aTps_J_[7] = -19;
    coeff2aTps_J_[8] = -15;
    coeff2aTps_J_[9] = -6;
    coeff2aTps_J_[10] = -26;
    coeff2aTps_J_[11] = -21;
    coeff2aTps_J_[12] = -17;
    coeff2aTps_J_[13] = -16;
    coeff2aTps_J_[14] = -9;
    coeff2aTps_J_[15] = -8;
    coeff2aTps_J_[16] = -15;
    coeff2aTps_J_[17] = -14;
    coeff2aTps_J_[18] = -26;
    coeff2aTps_J_[19] = -13;
    coeff2aTps_J_[20] = -9;
    coeff2aTps_J_[21] = -7;
    coeff2aTps_J_[22] = -27;
    coeff2aTps_J_[23] = -25;
    coeff2aTps_J_[24] = -11;
    coeff2aTps_J_[25] = -6;
    coeff2aTps_J_[26] = 1;
    coeff2aTps_J_[27] = 4;
    coeff2aTps_J_[28] = 8;
    coeff2aTps_J_[29] = 11;
    coeff2aTps_J_[30] = 0;
    coeff2aTps_J_[31] = 1;
    coeff2aTps_J_[32] = 5;
    coeff2aTps_J_[33] = 6;
    coeff2aTps_J_[34] = 10;
    coeff2aTps_J_[35] = 14;
    coeff2aTps_J_[36] = 16;
    coeff2aTps_J_[37] = 0;
    coeff2aTps_J_[38] = 4;
    coeff2aTps_J_[39] = 9;
    coeff2aTps_J_[40] = 17;
    coeff2aTps_J_[41] = 7;
    coeff2aTps_J_[42] = 18;
    coeff2aTps_J_[43] = 3;
    coeff2aTps_J_[44] = 15;
    coeff2aTps_J_[45] = 5;
    coeff2aTps_J_[46] = 18;

    coeff2aTps_n_[0] = 0.;
    // backward T(p,s) coefficient region 2a n_i
    coeff2aTps_n_[1] = -0.39235983861984 * 1.0e+6;
    coeff2aTps_n_[2] = 0.51526573827270 * 1.0e+6;
    coeff2aTps_n_[3] = 0.40482443161048 * 1.0e+5;
    coeff2aTps_n_[4] = -0.32193790923902 * 1.0e+3;
    coeff2aTps_n_[5] = 0.96961424218694 * 1.0e+2;
    coeff2aTps_n_[6] = -0.22867846371773 * 1.0e+2;
    coeff2aTps_n_[7] = -0.44942914124357 * 1.0e+6;
    coeff2aTps_n_[8] = -0.50118336020166 * 1.0e+4;
    coeff2aTps_n_[9] = 0.35684463560015;
    coeff2aTps_n_[10] = 0.44235335848190 * 1.0e+5;
    coeff2aTps_n_[11] = -0.13673388811708 * 1.0e+5;
    coeff2aTps_n_[12] = 0.42163260207864 * 1.0e+6;
    coeff2aTps_n_[13] = 0.22516925837475 * 1.0e+5;
    coeff2aTps_n_[14] = 0.47442144865646 * 1.0e+3;
    coeff2aTps_n_[15] = -0.14931130797647 * 1.0e+3;
    coeff2aTps_n_[16] = -0.19781126320452 * 1.0e+6;
    coeff2aTps_n_[17] = -0.23554399470760 * 1.0e+5;
    coeff2aTps_n_[18] = -0.19070616302076 * 1.0e+5;
    coeff2aTps_n_[19] = 0.55375669883164 * 1.0e+5;
    coeff2aTps_n_[20] = 0.38293691437363 * 1.0e+4;
    coeff2aTps_n_[21] = -0.60391860580567 * 1.0e+3;
    coeff2aTps_n_[22] = 0.19363102620331 * 1.0e+4;
    coeff2aTps_n_[23] = 0.42660643698610 * 1.0e+4;
    coeff2aTps_n_[24] = -0.59780638872718 * 1.0e+4;
    coeff2aTps_n_[25] = -0.70401463926862 * 1.0e+3;
    coeff2aTps_n_[26] = 0.33836784107553 * 1.0e+3;
    coeff2aTps_n_[27] = 0.20862786635187 * 1.0e+2;
    coeff2aTps_n_[28] = 0.33834172656196 * 1.0e-1;
    coeff2aTps_n_[29] = -0.43124428414893 * 1.0e-4;
    coeff2aTps_n_[30] = 0.16653791356412 * 1.0e+3;
    coeff2aTps_n_[31] = -0.13986292055898 * 1.0e+3;
    coeff2aTps_n_[32] = -0.78849547999872;
    coeff2aTps_n_[33] = 0.72132411753872 * 1.0e-1;
    coeff2aTps_n_[34] = -0.59754839398283 * 1.0e-2;
    coeff2aTps_n_[35] = -0.12141358953904 * 1.0e-4;
    coeff2aTps_n_[36] = 0.23227096733871 * 1.0e-6;
    coeff2aTps_n_[37] = -0.10538463566194 * 1.0e+2;
    coeff2aTps_n_[38] = 0.20718925496502 * 1.0e+1;
    coeff2aTps_n_[39] = -0.72193155260427 * 1.0e-1;
    coeff2aTps_n_[40] = 0.20749887081120 * 1.0e-6;
    coeff2aTps_n_[41] = -0.18340657911379 * 1.0e-1;
    coeff2aTps_n_[42] = 0.29036272348696 * 1.0e-6;
    coeff2aTps_n_[43] = 0.21037527893619;
    coeff2aTps_n_[44] = 0.25681239729999 * 1.0e-3;
    coeff2aTps_n_[45] = -0.12799002933781 * 1.0e-1;
    coeff2aTps_n_[46] = -0.82198102652018 * 1.0e-5;

    temperature_ref2bTps_ = 1.;
    pressure_ref2bTps_ = 1.0e+6;
    entropy_ref2bTps_ = 0.7853 * 1.0e+3;

    coeff2bTps_I_[0] = 0;
    // backward T(p,s) coefficient region 2b I_i
    coeff2bTps_I_[1] = -6;
    coeff2bTps_I_[2] = -6;
    coeff2bTps_I_[3] = -5;
    coeff2bTps_I_[4] = -5;
    coeff2bTps_I_[5] = -4;
    coeff2bTps_I_[6] = -4;
    coeff2bTps_I_[7] = -4;
    coeff2bTps_I_[8] = -3;
    coeff2bTps_I_[9] = -3;
    coeff2bTps_I_[10] = -3;
    coeff2bTps_I_[11] = -3;
    coeff2bTps_I_[12] = -2;
    coeff2bTps_I_[13] = -2;
    coeff2bTps_I_[14] = -2;
    coeff2bTps_I_[15] = -2;
    coeff2bTps_I_[16] = -1;
    coeff2bTps_I_[17] = -1;
    coeff2bTps_I_[18] = -1;
    coeff2bTps_I_[19] = -1;
    coeff2bTps_I_[20] = -1;
    coeff2bTps_I_[21] = 0;
    coeff2bTps_I_[22] = 0;
    coeff2bTps_I_[23] = 0;
    coeff2bTps_I_[24] = 0;
    coeff2bTps_I_[25] = 0;
    coeff2bTps_I_[26] = 0;
    coeff2bTps_I_[27] = 0;
    coeff2bTps_I_[28] = 1;
    coeff2bTps_I_[29] = 1;
    coeff2bTps_I_[30] = 1;
    coeff2bTps_I_[31] = 1;
    coeff2bTps_I_[32] = 1;
    coeff2bTps_I_[33] = 1;
    coeff2bTps_I_[34] = 2;
    coeff2bTps_I_[35] = 2;
    coeff2bTps_I_[36] = 2;
    coeff2bTps_I_[37] = 3;
    coeff2bTps_I_[38] = 3;
    coeff2bTps_I_[39] = 3;
    coeff2bTps_I_[40] = 4;
    coeff2bTps_I_[41] = 4;
    coeff2bTps_I_[42] = 5;
    coeff2bTps_I_[43] = 5;
    coeff2bTps_I_[44] = 5;

    coeff2bTps_J_[0] = 0;
    // backward T(p,s) coefficient region 2b J_i
    coeff2bTps_J_[1] = 0;
    coeff2bTps_J_[2] = 11;
    coeff2bTps_J_[3] = 0;
    coeff2bTps_J_[4] = 11;
    coeff2bTps_J_[5] = 0;
    coeff2bTps_J_[6] = 1;
    coeff2bTps_J_[7] = 11;
    coeff2bTps_J_[8] = 0;
    coeff2bTps_J_[9] = 1;
    coeff2bTps_J_[10] = 11;
    coeff2bTps_J_[11] = 12;
    coeff2bTps_J_[12] = 0;
    coeff2bTps_J_[13] = 1;
    coeff2bTps_J_[14] = 6;
    coeff2bTps_J_[15] = 10;
    coeff2bTps_J_[16] = 0;
    coeff2bTps_J_[17] = 1;
    coeff2bTps_J_[18] = 5;
    coeff2bTps_J_[19] = 8;
    coeff2bTps_J_[20] = 9;
    coeff2bTps_J_[21] = 0;
    coeff2bTps_J_[22] = 1;
    coeff2bTps_J_[23] = 2;
    coeff2bTps_J_[24] = 4;
    coeff2bTps_J_[25] = 5;
    coeff2bTps_J_[26] = 6;
    coeff2bTps_J_[27] = 9;
    coeff2bTps_J_[28] = 0;
    coeff2bTps_J_[29] = 1;
    coeff2bTps_J_[30] = 2;
    coeff2bTps_J_[31] = 3;
    coeff2bTps_J_[32] = 7;
    coeff2bTps_J_[33] = 8;
    coeff2bTps_J_[34] = 0;
    coeff2bTps_J_[35] = 1;
    coeff2bTps_J_[36] = 5;
    coeff2bTps_J_[37] = 0;
    coeff2bTps_J_[38] = 1;
    coeff2bTps_J_[39] = 3;
    coeff2bTps_J_[40] = 0;
    coeff2bTps_J_[41] = 1;
    coeff2bTps_J_[42] = 0;
    coeff2bTps_J_[43] = 1;
    coeff2bTps_J_[44] = 2;

    coeff2bTps_n_[0] = 0.;
    // backward T(p,s) coefficient region 2b n_i
    coeff2bTps_n_[1] = 0.31687665083497 * 1.0e+6;
    coeff2bTps_n_[2] = 0.20864175881858 * 1.0e+2;
    coeff2bTps_n_[3] = -0.39859399803599 * 1.0e+6;
    coeff2bTps_n_[4] = -0.21816058518877 * 1.0e+2;
    coeff2bTps_n_[5] = 0.22369785194242 * 1.0e+6;
    coeff2bTps_n_[6] = -0.27841703445817 * 1.0e+4;
    coeff2bTps_n_[7] = 0.99207436071480 * 1.0e+1;
    coeff2bTps_n_[8] = -0.75197512299157 * 1.0e+5;
    coeff2bTps_n_[9] = 0.29708605951158 * 1.0e+4;
    coeff2bTps_n_[10] = -0.34406878548526 * 1.0e+1;
    coeff2bTps_n_[11] = 0.38815564249115;
    coeff2bTps_n_[12] = 0.17511295085750 * 1.0e+5;
    coeff2bTps_n_[13] = -0.14237112854449 * 1.0e+4;
    coeff2bTps_n_[14] = 0.10943803364167 * 1.0e+1;
    coeff2bTps_n_[15] = 0.89971619308495;
    coeff2bTps_n_[16] = -0.33759740098958 * 1.0e+4;
    coeff2bTps_n_[17] = 0.47162885818355 * 1.0e+3;
    coeff2bTps_n_[18] = -0.19188241993679 * 1.0e+1;
    coeff2bTps_n_[19] = 0.41078580492196;
    coeff2bTps_n_[20] = -0.33465378172097;
    coeff2bTps_n_[21] = 0.13870034777505 * 1.0e+4;
    coeff2bTps_n_[22] = -0.40663326195838 * 1.0e+3;
    coeff2bTps_n_[23] = 0.41727347159610 * 1.0e+2;
    coeff2bTps_n_[24] = 0.21932549434532 * 1.0e+1;
    coeff2bTps_n_[25] = -0.10320050009077 * 1.0e+1;
    coeff2bTps_n_[26] = 0.35882943516703;
    coeff2bTps_n_[27] = 0.52511453726066 * 1.0e-2;
    coeff2bTps_n_[28] = 0.12838916450705 * 1.0e+2;
    coeff2bTps_n_[29] = -0.28642437219381 * 1.0e+1;
    coeff2bTps_n_[30] = 0.56912683664855;
    coeff2bTps_n_[31] = -0.99962954584931 * 1.0e-1;
    coeff2bTps_n_[32] = -0.32632037778459 * 1.0e-2;
    coeff2bTps_n_[33] = 0.23320922576723 * 1.0e-3;
    coeff2bTps_n_[34] = -0.15334809857450;
    coeff2bTps_n_[35] = 0.29072288239902 * 1.0e-1;
    coeff2bTps_n_[36] = 0.37534702741167 * 1.0e-3;
    coeff2bTps_n_[37] = 0.17296691702411 * 1.0e-2;
    coeff2bTps_n_[38] = -0.38556050844504 * 1.0e-3;
    coeff2bTps_n_[39] = -0.35017712292608 * 1.0e-4;
    coeff2bTps_n_[40] = -0.14566393631492 * 1.0e-4;
    coeff2bTps_n_[41] = 0.56420857267269 * 1.0e-5;
    coeff2bTps_n_[42] = 0.41286150074605 * 1.0e-7;
    coeff2bTps_n_[43] = -0.20684671118824 * 1.0e-7;
    coeff2bTps_n_[44] = 0.16409393674725 * 1.0e-8;

    temperature_ref2cTps_ = 1.;
    pressure_ref2cTps_ = 1.0e+6;
    entropy_ref2cTps_ = 2.9251 * 1.0e+3;

    coeff2cTps_I_[0] = 0;
    // backward T(p,s) coefficient region 2c I_i
    coeff2cTps_I_[1] = -2;
    coeff2cTps_I_[2] = -2;
    coeff2cTps_I_[3] = -1;
    coeff2cTps_I_[4] = 0;
    coeff2cTps_I_[5] = 0;
    coeff2cTps_I_[6] = 0;
    coeff2cTps_I_[7] = 0;
    coeff2cTps_I_[8] = 1;
    coeff2cTps_I_[9] = 1;
    coeff2cTps_I_[10] = 1;
    coeff2cTps_I_[11] = 1;
    coeff2cTps_I_[12] = 2;
    coeff2cTps_I_[13] = 2;
    coeff2cTps_I_[14] = 2;
    coeff2cTps_I_[15] = 3;
    coeff2cTps_I_[16] = 3;
    coeff2cTps_I_[17] = 3;
    coeff2cTps_I_[18] = 4;
    coeff2cTps_I_[19] = 4;
    coeff2cTps_I_[20] = 4;
    coeff2cTps_I_[21] = 5;
    coeff2cTps_I_[22] = 5;
    coeff2cTps_I_[23] = 5;
    coeff2cTps_I_[24] = 6;
    coeff2cTps_I_[25] = 6;
    coeff2cTps_I_[26] = 7;
    coeff2cTps_I_[27] = 7;
    coeff2cTps_I_[28] = 7;
    coeff2cTps_I_[29] = 7;
    coeff2cTps_I_[30] = 7;

    coeff2cTps_J_[0] = 0;
    // backward T(p,s) coefficient region 2c J_i
    coeff2cTps_J_[1] = 0;
    coeff2cTps_J_[2] = 1;
    coeff2cTps_J_[3] = 0;
    coeff2cTps_J_[4] = 0;
    coeff2cTps_J_[5] = 1;
    coeff2cTps_J_[6] = 2;
    coeff2cTps_J_[7] = 3;
    coeff2cTps_J_[8] = 0;
    coeff2cTps_J_[9] = 1;
    coeff2cTps_J_[10] = 3;
    coeff2cTps_J_[11] = 4;
    coeff2cTps_J_[12] = 0;
    coeff2cTps_J_[13] = 1;
    coeff2cTps_J_[14] = 2;
    coeff2cTps_J_[15] = 0;
    coeff2cTps_J_[16] = 1;
    coeff2cTps_J_[17] = 5;
    coeff2cTps_J_[18] = 0;
    coeff2cTps_J_[19] = 1;
    coeff2cTps_J_[20] = 4;
    coeff2cTps_J_[21] = 0;
    coeff2cTps_J_[22] = 1;
    coeff2cTps_J_[23] = 2;
    coeff2cTps_J_[24] = 0;
    coeff2cTps_J_[25] = 1;
    coeff2cTps_J_[26] = 0;
    coeff2cTps_J_[27] = 1;
    coeff2cTps_J_[28] = 3;
    coeff2cTps_J_[29] = 4;
    coeff2cTps_J_[30] = 5;

    coeff2cTps_n_[0] = 0.;
    // backward T(p,s) coefficient region 2c n_i
    coeff2cTps_n_[1] = 0.90968501005365 * 1.0e+3;
    coeff2cTps_n_[2] = 0.24045667088420 * 1.0e+4;
    coeff2cTps_n_[3] = -0.59162326387130 * 1.0e+3;
    coeff2cTps_n_[4] = 0.54145404128074 * 1.0e+3;
    coeff2cTps_n_[5] = -0.27098308411192 * 1.0e+3;
    coeff2cTps_n_[6] = 0.97976525097926 * 1.0e+3;
    coeff2cTps_n_[7] = -0.46966772959435 * 1.0e+3;
    coeff2cTps_n_[8] = 0.14399274604723 * 1.0e+2;
    coeff2cTps_n_[9] = -0.19104204230429 * 1.0e+2;
    coeff2cTps_n_[10] = 0.53299167111971 * 1.0e+1;
    coeff2cTps_n_[11] = -0.21252975375934 * 1.0e+2;
    coeff2cTps_n_[12] = -0.31147334413760;
    coeff2cTps_n_[13] = 0.60334840894623;
    coeff2cTps_n_[14] = -0.42764839702509 * 1.0e-1;
    coeff2cTps_n_[15] = 0.58185597255259 * 1.0e-2;
    coeff2cTps_n_[16] = -0.14597008284753 * 1.0e-1;
    coeff2cTps_n_[17] = 0.56631175631027 * 1.0e-2;
    coeff2cTps_n_[18] = -0.76155864584577 * 1.0e-4;
    coeff2cTps_n_[19] = 0.22440342919332 * 1.0e-3;
    coeff2cTps_n_[20] = -0.12561095013413 * 1.0e-4;
    coeff2cTps_n_[21] = 0.63323132660934 * 1.0e-6;
    coeff2cTps_n_[22] = -0.20541989675375 * 1.0e-5;
    coeff2cTps_n_[23] = 0.36405370390082 * 1.0e-7;
    coeff2cTps_n_[24] = -0.29759897789215 * 1.0e-8;
    coeff2cTps_n_[25] = 0.10136618529763 * 1.0e-7;
    coeff2cTps_n_[26] = 0.59925719692351 * 1.0e-11;
    coeff2cTps_n_[27] = -0.20677870105164 * 1.0e-10;
    coeff2cTps_n_[28] = -0.20874278181886 * 1.0e-10;
    coeff2cTps_n_[29] = 0.10162166825089 * 1.0e-9;
    coeff2cTps_n_[30] = -0.16429828281347 * 1.0e-9;

    mdensity_ref3_ = mdensity_crit_;
    temperature_ref3_ = temperature_crit_;

    coeff3_I_[0] = 0;
    // coefficient region 3 I_i
    coeff3_I_[1] = 0;
    coeff3_I_[2] = 0;
    coeff3_I_[3] = 0;
    coeff3_I_[4] = 0;
    coeff3_I_[5] = 0;
    coeff3_I_[6] = 0;
    coeff3_I_[7] = 0;
    coeff3_I_[8] = 0;
    coeff3_I_[9] = 1;
    coeff3_I_[10] = 1;
    coeff3_I_[11] = 1;
    coeff3_I_[12] = 1;
    coeff3_I_[13] = 2;
    coeff3_I_[14] = 2;
    coeff3_I_[15] = 2;
    coeff3_I_[16] = 2;
    coeff3_I_[17] = 2;
    coeff3_I_[18] = 2;
    coeff3_I_[19] = 3;
    coeff3_I_[20] = 3;
    coeff3_I_[21] = 3;
    coeff3_I_[22] = 3;
    coeff3_I_[23] = 3;
    coeff3_I_[24] = 4;
    coeff3_I_[25] = 4;
    coeff3_I_[26] = 4;
    coeff3_I_[27] = 4;
    coeff3_I_[28] = 5;
    coeff3_I_[29] = 5;
    coeff3_I_[30] = 5;
    coeff3_I_[31] = 6;
    coeff3_I_[32] = 6;
    coeff3_I_[33] = 6;
    coeff3_I_[34] = 7;
    coeff3_I_[35] = 8;
    coeff3_I_[36] = 9;
    coeff3_I_[37] = 9;
    coeff3_I_[38] = 10;
    coeff3_I_[39] = 10;
    coeff3_I_[40] = 11;

    coeff3_J_[0] = 0;
    // coefficient region 3 J_i
    coeff3_J_[1] = 0;
    coeff3_J_[2] = 0;
    coeff3_J_[3] = 1;
    coeff3_J_[4] = 2;
    coeff3_J_[5] = 7;
    coeff3_J_[6] = 10;
    coeff3_J_[7] = 12;
    coeff3_J_[8] = 23;
    coeff3_J_[9] = 2;
    coeff3_J_[10] = 6;
    coeff3_J_[11] = 15;
    coeff3_J_[12] = 17;
    coeff3_J_[13] = 0;
    coeff3_J_[14] = 2;
    coeff3_J_[15] = 6;
    coeff3_J_[16] = 7;
    coeff3_J_[17] = 22;
    coeff3_J_[18] = 26;
    coeff3_J_[19] = 0;
    coeff3_J_[20] = 2;
    coeff3_J_[21] = 4;
    coeff3_J_[22] = 16;
    coeff3_J_[23] = 26;
    coeff3_J_[24] = 0;
    coeff3_J_[25] = 2;
    coeff3_J_[26] = 4;
    coeff3_J_[27] = 26;
    coeff3_J_[28] = 1;
    coeff3_J_[29] = 3;
    coeff3_J_[30] = 26;
    coeff3_J_[31] = 0;
    coeff3_J_[32] = 2;
    coeff3_J_[33] = 26;
    coeff3_J_[34] = 2;
    coeff3_J_[35] = 26;
    coeff3_J_[36] = 2;
    coeff3_J_[37] = 26;
    coeff3_J_[38] = 0;
    coeff3_J_[39] = 1;
    coeff3_J_[40] = 26;

    coeff3_n_[0] = 0.;
    // coefficient region 3 n_i
    coeff3_n_[1] = 0.10658070028513 * 1.0e+1;
    coeff3_n_[2] = -0.15732845290239 * 1.0e+2;
    coeff3_n_[3] = 0.20944396974307 * 1.0e+2;
    coeff3_n_[4] = -0.76867707878716 * 1.0e+1;
    coeff3_n_[5] = 0.26185947787954 * 1.0e+1;
    coeff3_n_[6] = -0.28080781148620 * 1.0e+1;
    coeff3_n_[7] = 0.12053369696517 * 1.0e+1;
    coeff3_n_[8] = -0.84566812812502 * 1.0e-2;
    coeff3_n_[9] = -0.12654315477714 * 1.0e+1;
    coeff3_n_[10] = -0.11524407806681 * 1.0e+1;
    coeff3_n_[11] = 0.88521043984318;
    coeff3_n_[12] = -0.64207765181607;
    coeff3_n_[13] = 0.38493460186671;
    coeff3_n_[14] = -0.85214708824206;
    coeff3_n_[15] = 0.48972281541877 * 1.0e+1;
    coeff3_n_[16] = -0.30502617256965 * 1.0e+1;
    coeff3_n_[17] = 0.39420536879154 * 1.0e-1;
    coeff3_n_[18] = 0.12558408424308;
    coeff3_n_[19] = -0.27999329698710;
    coeff3_n_[20] = 0.13899799569460 * 1.0e+1;
    coeff3_n_[21] = -0.20189915023570 * 1.0e+1;
    coeff3_n_[22] = -0.82147637173963 * 1.0e-2;
    coeff3_n_[23] = -0.47596035734923;
    coeff3_n_[24] = 0.43984074473500 * 1.0e-1;
    coeff3_n_[25] = -0.44476435428739;
    coeff3_n_[26] = 0.90572070719733;
    coeff3_n_[27] = 0.70522450087967;
    coeff3_n_[28] = 0.10770512626332;
    coeff3_n_[29] = -0.32913623258954;
    coeff3_n_[30] = -0.50871062041158;
    coeff3_n_[31] = -0.22175400873096 * 1.0e-1;
    coeff3_n_[32] = 0.94260751665092 * 1.0e-1;
    coeff3_n_[33] = 0.16436278447961;
    coeff3_n_[34] = -0.13503372241348 * 1.0e-1;
    coeff3_n_[35] = -0.14834345352472 * 1.0e-1;
    coeff3_n_[36] = 0.57922953628084 * 1.0e-3;
    coeff3_n_[37] = 0.32308904703711 * 1.0e-2;
    coeff3_n_[38] = 0.80964802996215 * 1.0e-4;
    coeff3_n_[39] = -0.16557679795037 * 1.0e-3;
    coeff3_n_[40] = -0.44923899061815 * 1.0e-4;

    temperature_ref4_ = 1.;
    pressure_ref4_ = 1.0e+6;

    coeff4_n_[0] = 0.;
    // coefficient region 4 n_i
    coeff4_n_[1] = 0.11670521452767 * 1.0e+4;
    coeff4_n_[2] = -0.72421316703206 * 1.0e+6;
    coeff4_n_[3] = -0.17073846940092 * 1.0e+2;
    coeff4_n_[4] = 0.12020824702470 * 1.0e+5;
    coeff4_n_[5] = -0.32325550322333 * 1.0e+7;
    coeff4_n_[6] = 0.14915108613530 * 1.0e+2;
    coeff4_n_[7] = -0.48232657361591 * 1.0e+4;
    coeff4_n_[8] = 0.40511340542057 * 1.0e+6;
    coeff4_n_[9] = -0.23855557567849;
    coeff4_n_[10] = 0.65017534844798 * 1.0e+3;

    temperature_ref5_ = 1000.;
    pressure_ref5_ = 1.0e+6;

    coeff5_ide_J_[0] = 0;
    // coefficient region 5 ideal J_i
    coeff5_ide_J_[1] = 0;
    coeff5_ide_J_[2] = 1;
    coeff5_ide_J_[3] = -3;
    coeff5_ide_J_[4] = -2;
    coeff5_ide_J_[5] = -1;
    coeff5_ide_J_[6] = 2;

    coeff5_ide_n_[0] = 0.;
    // coefficient region 5 ideal n_i
    coeff5_ide_n_[1] = -0.13179983674201 * 1.0e+2;
    coeff5_ide_n_[2] = 0.68540841634434 * 1.0e+1;
    coeff5_ide_n_[3] = -0.24805148933466 * 1.0e-1;
    coeff5_ide_n_[4] = 0.36901534980333;
    coeff5_ide_n_[5] = -0.31161318213925 * 1.0e+1;
    coeff5_ide_n_[6] = -0.32961626538917;

    coeff5_res_I_[0] = 0;
    // coefficient region 5 residual I_i
    coeff5_res_I_[1] = 1;
    coeff5_res_I_[2] = 1;
    coeff5_res_I_[3] = 1;
    coeff5_res_I_[4] = 2;
    coeff5_res_I_[5] = 2;
    coeff5_res_I_[6] = 3;

    coeff5_res_J_[0] = 0;
    // coefficient region 5 residual J_i
    coeff5_res_J_[1] = 1;
    coeff5_res_J_[2] = 2;
    coeff5_res_J_[3] = 3;
    coeff5_res_J_[4] = 3;
    coeff5_res_J_[5] = 9;
    coeff5_res_J_[6] = 7;

    coeff5_res_n_[0] = 0.;
    // coefficient region 5 residual n_i
    coeff5_res_n_[1] = 0.15736404855259 * 1.0e-2;
    coeff5_res_n_[2] = 0.90153761673944 * 1.0e-3;
    coeff5_res_n_[3] = -0.50270077677648 * 1.0e-2;
    coeff5_res_n_[4] = 0.22440037409485 * 1.0e-5;
    coeff5_res_n_[5] = -0.41163275453471 * 1.0e-5;
    coeff5_res_n_[6] = 0.37919454822955 * 1.0e-7;

    return;
}

} // end namespace IAPWS
