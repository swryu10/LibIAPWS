#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS97.h"

void print_functions(IAPWS::Lib97 *ptr_eos,
                     double temperature_in,
                     double pressure_in) {
    fprintf(stdout, "temperature = %f degK\n", temperature_in);
    fprintf(stdout, "pressure = %f Pa\n", pressure_in);

    fprintf(stdout, "    v = %.9e m^3 / kg\n",
        ptr_eos->get_param_vol_spec(temperature_in, pressure_in));
    fprintf(stdout, "    erg_int = %.9e J / kg\n",
        ptr_eos->get_param_erg_int(temperature_in, pressure_in));
    fprintf(stdout, "    entropy = %.9e J / kg / degK\n",
        ptr_eos->get_param_entropy(temperature_in, pressure_in));
    fprintf(stdout, "    enthalpy = %.9e J / kg\n",
        ptr_eos->get_param_enthalpy(temperature_in, pressure_in));
    fprintf(stdout, "    heat_c_p = %.9e J / kg / degK\n",
        ptr_eos->get_param_heat_c_p(temperature_in, pressure_in));
    fprintf(stdout, "    heat_c_v = %.9e J / kg / degK\n",
        ptr_eos->get_param_heat_c_v(temperature_in, pressure_in));
    fprintf(stdout, "    speed_sound = %.9e m / sec\n",
        ptr_eos->get_param_speed_sound(temperature_in, pressure_in));

    return;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib97 iapws97eos;
    iapws97eos.print_header();

    print_functions(&iapws97eos, 300., 3. * 1.0e+6);

    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 300., 80. * 1.0e+6);

    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 500., 3. * 1.0e+6);

    return 0;
}
