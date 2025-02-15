#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS06.h"

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib06 iapws06eos;
    iapws06eos.print_header();

    double temp0_in = 273.16;
    double press0_in = 611.657;

    fprintf(stdout, "temperature = %f degK\n", temp0_in);
    fprintf(stdout, "pressure = %f Pa\n", press0_in);

    fprintf(stdout, "    g = %.8e J / kg\n",
        iapws06eos.get_param_g(temp0_in, press0_in));
    fprintf(stdout, "    dg_dp = %.11e m^3 / kg\n",
        iapws06eos.get_param_dg_dp(temp0_in, press0_in));
    fprintf(stdout, "    dg_dT = %.11e J / kg / degK\n",
        iapws06eos.get_param_dg_dT(temp0_in, press0_in));
    fprintf(stdout, "    d2g_dp_dp = %.11e m^3 / kg / Pa\n",
        iapws06eos.get_param_d2g_dp_dp(temp0_in, press0_in));
    fprintf(stdout, "    d2g_dT_dp = %.11e m^3 / kg / degK\n",
        iapws06eos.get_param_d2g_dT_dp(temp0_in, press0_in));
    fprintf(stdout, "    d2g_dT_dT = %.11e J / kg / degK^2\n",
        iapws06eos.get_param_d2g_dT_dT(temp0_in, press0_in));

    fprintf(stdout, "\n");

    double temp1_in = 273.152519;
    double press1_in = 101325.;

    fprintf(stdout, "temperature = %f degK\n", temp1_in);
    fprintf(stdout, "pressure = %f Pa\n", press1_in);

    fprintf(stdout, "    g = %.10e J / kg\n",
        iapws06eos.get_param_g(temp1_in, press1_in));
    fprintf(stdout, "    dg_dp = %.11e m^3 / kg\n",
        iapws06eos.get_param_dg_dp(temp1_in, press1_in));
    fprintf(stdout, "    dg_dT = %.11e J / kg / degK\n",
        iapws06eos.get_param_dg_dT(temp1_in, press1_in));
    fprintf(stdout, "    d2g_dp_dp = %.11e m^3 / kg / Pa\n",
        iapws06eos.get_param_d2g_dp_dp(temp1_in, press1_in));
    fprintf(stdout, "    d2g_dT_dp = %.11e m^3 / kg / degK\n",
        iapws06eos.get_param_d2g_dT_dp(temp1_in, press1_in));
    fprintf(stdout, "    d2g_dT_dT = %.11e J / kg / degK^2\n",
        iapws06eos.get_param_d2g_dT_dT(temp1_in, press1_in));

    fprintf(stdout, "\n");

    double temp2_in = 100.;
    double press2_in = 100. * 1.0e+6;

    fprintf(stdout, "temperature = %f degK\n", temp2_in);
    fprintf(stdout, "pressure = %f Pa\n", press2_in);

    fprintf(stdout, "    g = %.11e J / kg\n",
        iapws06eos.get_param_g(temp2_in, press2_in));
    fprintf(stdout, "    dg_dp = %.11e m^3 / kg\n",
        iapws06eos.get_param_dg_dp(temp2_in, press2_in));
    fprintf(stdout, "    dg_dT = %.11e J / kg / degK\n",
        iapws06eos.get_param_dg_dT(temp2_in, press2_in));
    fprintf(stdout, "    d2g_dp_dp = %.11e m^3 / kg / Pa\n",
        iapws06eos.get_param_d2g_dp_dp(temp2_in, press2_in));
    fprintf(stdout, "    d2g_dT_dp = %.11e m^3 / kg / degK\n",
        iapws06eos.get_param_d2g_dT_dp(temp2_in, press2_in));
    fprintf(stdout, "    d2g_dT_dT = %.11e J / kg / degK^2\n",
        iapws06eos.get_param_d2g_dT_dT(temp2_in, press2_in));

    return 0;
}
