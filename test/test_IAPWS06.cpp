#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS06.h"

void print_functions(IAPWS::Lib06 *ptr_eos,
                     double temperature_in,
                     double pressure_in) {
    fprintf(stdout, "temperature = %f degK\n", temperature_in);
    fprintf(stdout, "pressure = %f Pa\n", pressure_in);

    fprintf(stdout, "    g = %.11e J / kg\n",
        ptr_eos->get_param_g(temperature_in, pressure_in));
    fprintf(stdout, "    dg_dp = %.11e m^3 / kg\n",
        ptr_eos->get_param_dg_dp(temperature_in, pressure_in));
    fprintf(stdout, "    dg_dT = %.11e J / kg / degK\n",
        ptr_eos->get_param_dg_dT(temperature_in, pressure_in));
    fprintf(stdout, "    d2g_dp_dp = %.11e m^3 / kg / Pa\n",
        ptr_eos->get_param_d2g_dp_dp(temperature_in, pressure_in));
    fprintf(stdout, "    d2g_dT_dp = %.11e m^3 / kg / degK\n",
        ptr_eos->get_param_d2g_dT_dp(temperature_in, pressure_in));
    fprintf(stdout, "    d2g_dT_dT = %.11e J / kg / degK^2\n",
        ptr_eos->get_param_d2g_dT_dT(temperature_in, pressure_in));
    fprintf(stdout, "    enthalpy = %.11e J / kg\n",
        ptr_eos->get_param_enthalpy(temperature_in, pressure_in));
    fprintf(stdout, "    f = %.11e J / kg\n",
        ptr_eos->get_param_f(temperature_in, pressure_in));
    fprintf(stdout, "    erg_int = %.11e J / kg\n",
        ptr_eos->get_param_erg_int(temperature_in, pressure_in));
    fprintf(stdout, "    entropy = %.11e J / kg / degK\n",
        ptr_eos->get_param_entropy(temperature_in, pressure_in));
    fprintf(stdout, "    heat_c_p = %.11e J / kg / degK\n",
        ptr_eos->get_param_heat_c_p(temperature_in, pressure_in));
    fprintf(stdout, "    mdensity = %.11e kg / m^3\n",
        ptr_eos->get_param_mdensity(temperature_in, pressure_in));
    fprintf(stdout, "    coeff_alpha = %.11e 1 / degK\n",
        ptr_eos->get_param_coeff_alpha(temperature_in, pressure_in));
    fprintf(stdout, "    coeff_beta = %.11e Pa / degK\n",
        ptr_eos->get_param_coeff_beta(temperature_in, pressure_in));
    fprintf(stdout, "    comp_kappa_T = %.11e 1 / Pa\n",
        ptr_eos->get_param_comp_kappa_T(temperature_in, pressure_in));
    fprintf(stdout, "    comp_kappa_s = %.11e 1 / Pa\n",
        ptr_eos->get_param_comp_kappa_s(temperature_in, pressure_in));

    return;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib06 iapws06eos;
    iapws06eos.print_header();

    double temp0_in = 273.16;
    double press0_in = 611.657;
    print_functions(&iapws06eos, temp0_in, press0_in);

    fprintf(stdout, "\n");

    double temp1_in = 273.152519;
    double press1_in = 101325.;
    print_functions(&iapws06eos, temp1_in, press1_in);

    fprintf(stdout, "\n");

    double temp2_in = 100.;
    double press2_in = 100. * 1.0e+6;
    print_functions(&iapws06eos, temp2_in, press2_in);

    return 0;
}
