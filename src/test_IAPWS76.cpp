#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS76.h"

int main(int argc, char *argv[]) {
    IAPWS::Lib76 iapws76sig;
    iapws76sig.print_header();

    fprintf(stdout, "  temperatue (degC)    surface tension (J / m^2)\n");

    double nbin = 74;
    double temperature_base = 273.15;
    for (int it = 0; it <= nbin; it++) {
        double temperature_degC;
        if (it == 0) {
            temperature_degC = 0.01;
        } else {
            temperature_degC = 5. * static_cast<double>(it);
        }

        double temperature_degK =
            temperature_degC + temperature_base;
        fprintf(stdout, "    %f    %.3e\n",
                temperature_degC,
                iapws76sig.get_tension_surf(temperature_degK));
    }

    return 0;
}
