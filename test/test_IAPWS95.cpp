#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"LibIAPWS95.h"

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib95 iapws95eos;
    iapws95eos.print_header();

    double temp0_in = 500.;
    double mden0_in = 838.025;
    fprintf(stdout, "temperature = %f degK\n", temp0_in);
    fprintf(stdout, "mdensity = %f kg / m^3\n", mden0_in);

    fprintf(stdout, "    phi_ide = %.8e\n",
        iapws95eos.get_param_phi_ide(mden0_in, temp0_in));
    fprintf(stdout, "    dphi_ide_ddelta = %.8e\n",
        iapws95eos.get_param_dphi_ide_ddelta(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_ide_ddelta_ddelta = %.8e\n",
        iapws95eos.get_param_d2phi_ide_ddelta_ddelta(mden0_in, temp0_in));
    fprintf(stdout, "    dphi_ide_dtau = %.8e\n",
        iapws95eos.get_param_dphi_ide_dtau(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_ide_dtau_dtau = %.8e\n",
        iapws95eos.get_param_d2phi_ide_dtau_dtau(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_ide_ddelta_dtau = %.8e\n",
        iapws95eos.get_param_d2phi_ide_ddelta_dtau(mden0_in, temp0_in));

    fprintf(stdout, "    phi_res = %.8e\n",
        iapws95eos.get_param_phi_res(mden0_in, temp0_in));
    fprintf(stdout, "    dphi_res_ddelta = %.8e\n",
        iapws95eos.get_param_dphi_res_ddelta(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_res_ddelta_ddelta = %.8e\n",
        iapws95eos.get_param_d2phi_res_ddelta_ddelta(mden0_in, temp0_in));
    fprintf(stdout, "    dphi_res_dtau = %.8e\n",
        iapws95eos.get_param_dphi_res_dtau(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_res_dtau_dtau = %.8e\n",
        iapws95eos.get_param_d2phi_res_dtau_dtau(mden0_in, temp0_in));
    fprintf(stdout, "    d2phi_res_ddelta_dtau = %.8e\n",
        iapws95eos.get_param_d2phi_res_ddelta_dtau(mden0_in, temp0_in));

    fprintf(stdout, "\n");

    double temp1_in = 300.;
    fprintf(stdout, "temperature = %f degK\n", temp1_in);

    double *mden1_in = new double[3];
    mden1_in[0] = 0.9965560 * 1.0e+3;
    mden1_in[1] = 0.1005308 * 1.0e+4;
    mden1_in[2] = 0.1188202 * 1.0e+4;
    for (int i = 0; i < 3; i++) {
        fprintf(stdout, "  mdensity = %f kg / m^3\n", mden1_in[i]);
        fprintf(stdout, "    pressure = %.8e\n",
            iapws95eos.get_param_pressure(mden1_in[i], temp1_in));
        fprintf(stdout, "    heat_c_v = %.8e\n",
            iapws95eos.get_param_heat_c_v(mden1_in[i], temp1_in));
        fprintf(stdout, "    speed_sound = %.8e\n",
            iapws95eos.get_param_speed_sound(mden1_in[i], temp1_in));
        fprintf(stdout, "    entropy = %.8e\n",
            iapws95eos.get_param_entropy(mden1_in[i], temp1_in));
    }
    delete [] mden1_in;

    fprintf(stdout, "\n");

    double temp2_in = 500.;
    fprintf(stdout, "temperature = %f degK\n", temp2_in);

    double *mden2_in = new double[4];
    mden2_in[0] = 0.435;
    mden2_in[1] = 0.4532 * 1.0e+1;
    mden2_in[2] = 0.838025 * 1.0e+3;
    mden2_in[3] = 0.1084564 * 1.0e+4;
    for (int i = 0; i < 4; i++) {
        fprintf(stdout, "  mdensity = %f kg / m^3\n", mden2_in[i]);
        fprintf(stdout, "    pressure = %.8e\n",
            iapws95eos.get_param_pressure(mden2_in[i], temp2_in));
        fprintf(stdout, "    heat_c_v = %.8e\n",
            iapws95eos.get_param_heat_c_v(mden2_in[i], temp2_in));
        fprintf(stdout, "    speed_sound = %.8e\n",
            iapws95eos.get_param_speed_sound(mden2_in[i], temp2_in));
        fprintf(stdout, "    entropy = %.8e\n",
            iapws95eos.get_param_entropy(mden2_in[i], temp2_in));
    }
    delete [] mden2_in;

    fprintf(stdout, "\n");

    double temp3_in = 647.;
    fprintf(stdout, "temperature = %f degK\n", temp3_in);

    double *mden3_in = new double[1];
    mden3_in[0] = 0.358 * 1.0e+3;
    for (int i = 0; i < 1; i++) {
        fprintf(stdout, "  mdensity = %f kg / m^3\n", mden3_in[i]);
        fprintf(stdout, "    pressure = %.8e\n",
            iapws95eos.get_param_pressure(mden3_in[i], temp3_in));
        fprintf(stdout, "    heat_c_v = %.8e\n",
            iapws95eos.get_param_heat_c_v(mden3_in[i], temp3_in));
        fprintf(stdout, "    speed_sound = %.8e\n",
            iapws95eos.get_param_speed_sound(mden3_in[i], temp3_in));
        fprintf(stdout, "    entropy = %.8e\n",
            iapws95eos.get_param_entropy(mden3_in[i], temp3_in));
    }
    delete [] mden3_in;

    fprintf(stdout, "\n");

    double temp4_in = 900.;
    fprintf(stdout, "temperature = %f degK\n", temp4_in);

    double *mden4_in = new double[3];
    mden4_in[0] = 0.241;
    mden4_in[1] = 0.52615 * 1.0e+2;
    mden4_in[2] = 0.870769 * 1.0e+3;
    for (int i = 0; i < 3; i++) {
        fprintf(stdout, "  mdensity = %f kg / m^3\n", mden4_in[i]);
        fprintf(stdout, "    pressure = %.8e\n",
            iapws95eos.get_param_pressure(mden4_in[i], temp4_in));
        fprintf(stdout, "    heat_c_v = %.8e\n",
            iapws95eos.get_param_heat_c_v(mden4_in[i], temp4_in));
        fprintf(stdout, "    speed_sound = %.8e\n",
            iapws95eos.get_param_speed_sound(mden4_in[i], temp4_in));
        fprintf(stdout, "    entropy = %.8e\n",
            iapws95eos.get_param_entropy(mden4_in[i], temp4_in));
    }
    delete [] mden4_in;

    fprintf(stdout, "\n");

    double temp5_in = 275.;
    fprintf(stdout, "temperature = %f degK\n", temp5_in);

    double press5_sat = 0.698451167 * 1.0e-3 * 1.0e+6;
    fprintf(stdout, "  pressure_sat = %f Pa\n", press5_sat);
    double mden5_vap = 0.005;
    bool found_mden5_vap =
        iapws95eos.find_root_mdensity(temp5_in, press5_sat, mden5_vap);
    if (found_mden5_vap) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_vap = %.8e kg / m^3, ", mden5_vap);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden5_vap, temp5_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden5_vap, temp5_in));
    }
    double mden5_liq = 1100.;
    bool found_mden5_liq =
        iapws95eos.find_root_mdensity(temp5_in, press5_sat, mden5_liq);
    if (found_mden5_liq) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_liq = %.8e kg / m^3, ", mden5_liq);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden5_liq, temp5_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden5_liq, temp5_in));
    }

    fprintf(stdout, "\n");

    double temp6_in = 450.;
    fprintf(stdout, "temperature = %f degK\n", temp6_in);

    double press6_sat = 0.932203564 * 1.0e+6;
    fprintf(stdout, "  pressure_sat = %f Pa\n", press6_sat);
    double mden6_vap = 0.005;
    bool found_mden6_vap =
        iapws95eos.find_root_mdensity(temp6_in, press6_sat, mden6_vap);
    if (found_mden6_vap) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_vap = %.8e kg / m^3, ", mden6_vap);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden6_vap, temp6_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden6_vap, temp6_in));
    }
    double mden6_liq = 1100.;
    bool found_mden6_liq =
        iapws95eos.find_root_mdensity(temp6_in, press6_sat, mden6_liq);
    if (found_mden6_liq) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_liq = %.8e kg / m^3, ", mden6_liq);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden6_liq, temp6_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden6_liq, temp6_in));
    }

    fprintf(stdout, "\n");

    double temp7_in = 625.;
    fprintf(stdout, "temperature = %f degK\n", temp7_in);

    double press7_sat = 0.169082693 * 1.0e+2 * 1.0e+6;
    fprintf(stdout, "  pressure_sat = %f Pa\n", press7_sat);
    double mden7_vap = 0.005;
    bool found_mden7_vap =
        iapws95eos.find_root_mdensity(temp7_in, press7_sat, mden7_vap);
    if (found_mden7_vap) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_vap = %.8e kg / m^3, ", mden7_vap);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden7_vap, temp7_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden7_vap, temp7_in));
    }
    double mden7_liq = 1100.;
    bool found_mden7_liq =
        iapws95eos.find_root_mdensity(temp7_in, press7_sat, mden7_liq);
    if (found_mden7_liq) {
        fprintf(stdout, "    ");
        fprintf(stdout, "mdensity_liq = %.8e kg / m^3, ", mden7_liq);
        fprintf(stdout, "enthalpy = %.8e J / kg, ",
                iapws95eos.get_param_enthalpy(mden7_liq, temp7_in));
        fprintf(stdout, "entropy = %.8e J / kg / degK\n",
                iapws95eos.get_param_entropy(mden7_liq, temp7_in));
    }

    return 0;
}