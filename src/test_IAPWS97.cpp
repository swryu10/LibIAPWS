#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS97.h"

void print_functions(IAPWS::Lib97 *ptr_eos,
                     double temperature_in,
                     double pressure_in,
                     bool flag_metastable = false) {
    if (flag_metastable) {
        fprintf(stdout, "  meta-stable\n");
    }
    fprintf(stdout, "  temperature = %f degK\n", temperature_in);
    fprintf(stdout, "  pressure = %f Pa\n", pressure_in);

    fprintf(stdout, "    spec_vol = %.9e m^3 / kg\n",
        ptr_eos->get_param_vol_spec(temperature_in,
                                    pressure_in,
                                    flag_metastable));
    fprintf(stdout, "    erg_int = %.9e J / kg\n",
        ptr_eos->get_param_erg_int(temperature_in,
                                   pressure_in,
                                   flag_metastable));
    fprintf(stdout, "    entropy = %.9e J / kg / degK\n",
        ptr_eos->get_param_entropy(temperature_in,
                                   pressure_in,
                                   flag_metastable));
    fprintf(stdout, "    enthalpy = %.9e J / kg\n",
        ptr_eos->get_param_enthalpy(temperature_in,
                                    pressure_in,
                                    flag_metastable));
    fprintf(stdout, "    heat_c_p = %.9e J / kg / degK\n",
        ptr_eos->get_param_heat_c_p(temperature_in,
                                    pressure_in,
                                    flag_metastable));
    fprintf(stdout, "    heat_c_v = %.9e J / kg / degK\n",
        ptr_eos->get_param_heat_c_v(temperature_in,
                                    pressure_in,
                                    flag_metastable));
    fprintf(stdout, "    speed_sound = %.9e m / sec\n",
        ptr_eos->get_param_speed_sound(temperature_in,
                                       pressure_in,
                                       flag_metastable));

    return;
}

void print1_temperature_ph(IAPWS::Lib97 *ptr_eos,
                           double pressure_in,
                           double enthalpy_in) {
    fprintf(stdout, "  pressure = %f Pa\n", pressure_in);
    fprintf(stdout, "  enthalpy = %f J / kg\n", enthalpy_in);

    fprintf(stdout, "    T = %.9e degK\n",
        ptr_eos->get_param1_temperature_ph(pressure_in,
                                           enthalpy_in));

    return;
}

void print1_temperature_ps(IAPWS::Lib97 *ptr_eos,
                           double pressure_in,
                           double entropy_in) {
    fprintf(stdout, "  pressure = %f Pa\n", pressure_in);
    fprintf(stdout, "  entropy = %f J / kg / degK\n", entropy_in);

    fprintf(stdout, "    T = %.9e degK\n",
        ptr_eos->get_param1_temperature_ps(pressure_in,
                                           entropy_in));

    return;
}

void print2_temperature_ph(IAPWS::Lib97 *ptr_eos,
                           double pressure_in,
                           double enthalpy_in) {
    fprintf(stdout, "  pressure = %f Pa\n", pressure_in);
    fprintf(stdout, "  enthalpy = %f J / kg\n", enthalpy_in);

    fprintf(stdout, "    T = %.9e degK\n",
        ptr_eos->get_param2_temperature_ph(pressure_in,
                                           enthalpy_in));

    return;
}

void print2_temperature_ps(IAPWS::Lib97 *ptr_eos,
                           double pressure_in,
                           double entropy_in) {
    fprintf(stdout, "  pressure = %f Pa\n", pressure_in);
    fprintf(stdout, "  entropy = %f J / kg / degK\n", entropy_in);

    fprintf(stdout, "    T = %.9e degK\n",
        ptr_eos->get_param2_temperature_ps(pressure_in,
                                           entropy_in));

    return;
}

void print3_functions(IAPWS::Lib97 *ptr_eos,
                      double mdensity_in,
                      double temperature_in) {
    fprintf(stdout, "  mdensity = %f kg / m^3\n", mdensity_in);
    fprintf(stdout, "  temperature = %f degK\n", temperature_in);

    fprintf(stdout, "    pressure = %.9e Pa\n",
        ptr_eos->get_param3_pressure(mdensity_in, temperature_in));
    fprintf(stdout, "    erg_int = %.9e J / kg\n",
        ptr_eos->get_param3_erg_int(mdensity_in, temperature_in));
    fprintf(stdout, "    entropy = %.9e J / kg / degK\n",
        ptr_eos->get_param3_entropy(mdensity_in, temperature_in));
    fprintf(stdout, "    enthalpy = %.9e J / kg\n",
        ptr_eos->get_param3_enthalpy(mdensity_in, temperature_in));
    fprintf(stdout, "    heat_c_p = %.9e J / kg / degK\n",
        ptr_eos->get_param3_heat_c_p(mdensity_in, temperature_in));
    fprintf(stdout, "    heat_c_v = %.9e J / kg / degK\n",
        ptr_eos->get_param3_heat_c_v(mdensity_in, temperature_in));
    fprintf(stdout, "    speed_sound = %.9e m / sec\n",
        ptr_eos->get_param3_speed_sound(mdensity_in, temperature_in));

    return;
}

void print4_sat_pressure(IAPWS::Lib97 *ptr_eos,
                         double temperature_in) {
    fprintf(stdout, "  temperature = %f degK\n", temperature_in);

    fprintf(stdout, "    sat pressure = %.9e Pa\n",
        ptr_eos->get_param4_sat_pressure(temperature_in));

    return;
}

void print4_sat_temperature(IAPWS::Lib97 *ptr_eos,
                            double pressure_in) {
    fprintf(stdout, "  temperature = %f Pa\n", pressure_in);

    fprintf(stdout, "    sat temperature = %.9e degK\n",
        ptr_eos->get_param4_sat_temperature(pressure_in));

    return;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib97 iapws97eos;
    iapws97eos.print_header();

    fprintf(stdout, "###  BOUNDARY 23  ###\n");
    fprintf(stdout, "\n");

    double temperature_B23 = 623.15;
    double pressure_B23 = 0.165291643 * 1.0e+8;
    fprintf(stdout, "    temperature = %e degK    >",
            temperature_B23);
    fprintf(stdout, "    pressure = %.9e Pa\n",
            iapws97eos.get_paramB23_pressure(temperature_B23));
    fprintf(stdout, "    pressure = %.9e Pa    >",
            pressure_B23);
    fprintf(stdout, "    temperature = %e degK\n",
            iapws97eos.get_paramB23_temperature(pressure_B23));
    fprintf(stdout, "\n");

    fprintf(stdout, "###  BOUNDARY 2bc  ###\n");
    fprintf(stdout, "\n");

    double pressure_B2bc = 1.0e+8;
    double enthalpy_B2bc = 0.3516004323 * 1.0e+7;
    fprintf(stdout, "    enthalpy = %.9e J / kg    >",
            enthalpy_B2bc);
    fprintf(stdout, "    pressure = %.9e Pa\n",
            iapws97eos.get_paramB2bc_pressure(enthalpy_B2bc));
    fprintf(stdout, "    pressure = %.9e Pa    >",
            pressure_B2bc);
    fprintf(stdout, "    enthalpy = %.9e J / kg\n",
            iapws97eos.get_paramB2bc_enthalpy(pressure_B2bc));
    fprintf(stdout, "\n");

    fprintf(stdout, "###  REGION 1  ###\n");
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 300., 3. * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 300., 80. * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 500., 3. * 1.0e+6);
    fprintf(stdout, "\n");

    print1_temperature_ph(&iapws97eos,
                          3. * 1.0e+6, 500. * 1.0e+3);
    fprintf(stdout, "\n");

    print1_temperature_ph(&iapws97eos,
                          80. * 1.0e+6, 500. * 1.0e+3);
    fprintf(stdout, "\n");

    print1_temperature_ph(&iapws97eos,
                          80. * 1.0e+6, 1500. * 1.0e+3);
    fprintf(stdout, "\n");

    print1_temperature_ps(&iapws97eos,
                          3. * 1.0e+6, 0.5 * 1.0e+3);
    fprintf(stdout, "\n");

    print1_temperature_ps(&iapws97eos,
                          80. * 1.0e+6, 0.5 * 1.0e+3);
    fprintf(stdout, "\n");

    print1_temperature_ps(&iapws97eos,
                          80. * 1.0e+6, 3. * 1.0e+3);
    fprintf(stdout, "\n");

    fprintf(stdout, "###  REGION 2  ###\n");
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 300., 0.0035 * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 700., 0.0035 * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 700., 30. * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 450., 1. * 1.0e+6, true);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 440., 1. * 1.0e+6, true);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 450., 1.5 * 1.0e+6, true);
    fprintf(stdout, "\n");

    print2_temperature_ph(&iapws97eos,
                          0.001 * 1.0e+6, 3000. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          3. * 1.0e+6, 3000. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          3. * 1.0e+6, 4000. * 1.0e+3);
    fprintf(stdout, "\n");

    print2_temperature_ph(&iapws97eos,
                          5. * 1.0e+6, 3500. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          5. * 1.0e+6, 4000. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          25. * 1.0e+6, 3500. * 1.0e+3);
    fprintf(stdout, "\n");

    print2_temperature_ph(&iapws97eos,
                          40. * 1.0e+6, 2700. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          60. * 1.0e+6, 2700. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ph(&iapws97eos,
                          60. * 1.0e+6, 3200. * 1.0e+3);
    fprintf(stdout, "\n");

    print2_temperature_ps(&iapws97eos,
                          0.1 * 1.0e+6, 7.5 * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          0.1 * 1.0e+6, 8. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          2.5 * 1.0e+6, 8. * 1.0e+3);
    fprintf(stdout, "\n");

    print2_temperature_ps(&iapws97eos,
                          8. * 1.0e+6, 6. * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          8. * 1.0e+6, 7.5 * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          90. * 1.0e+6, 6. * 1.0e+3);
    fprintf(stdout, "\n");

    print2_temperature_ps(&iapws97eos,
                          20. * 1.0e+6, 5.75 * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          80. * 1.0e+6, 5.25 * 1.0e+3);
    fprintf(stdout, "\n");
    print2_temperature_ps(&iapws97eos,
                          80. * 1.0e+6, 5.75 * 1.0e+3);
    fprintf(stdout, "\n");

    fprintf(stdout, "###  REGION 3  ###\n");
    fprintf(stdout, "\n");

    print3_functions(&iapws97eos, 500., 650.);
    fprintf(stdout, "\n");

    print3_functions(&iapws97eos, 200., 650.);
    fprintf(stdout, "\n");

    print3_functions(&iapws97eos, 500., 750.);
    fprintf(stdout, "\n");

    fprintf(stdout, "###  REGION 4  ###\n");
    fprintf(stdout, "\n");

    print4_sat_pressure(&iapws97eos, 300.);
    print4_sat_pressure(&iapws97eos, 500.);
    print4_sat_pressure(&iapws97eos, 600.);
    fprintf(stdout, "\n");

    print4_sat_temperature(&iapws97eos, 0.1 * 1.0e+6);
    print4_sat_temperature(&iapws97eos, 1. * 1.0e+6);
    print4_sat_temperature(&iapws97eos, 10. * 1.0e+6);
    fprintf(stdout, "\n");

    fprintf(stdout, "###  REGION 5  ###\n");
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 1500., 0.5 * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 1500., 30. * 1.0e+6);
    fprintf(stdout, "\n");

    print_functions(&iapws97eos, 2000., 30. * 1.0e+6);
    fprintf(stdout, "\n");

    return 0;
}
