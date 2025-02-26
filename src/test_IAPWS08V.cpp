#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS08V.h"

void print_viscosity(IAPWS::Lib08V *ptr_visc,
                     double mdensity_in,
                     double temperature_in) {
    fprintf(stdout, "mdensity = %f kg / m^3\n",
            mdensity_in);
    fprintf(stdout, "temperature = %f degK\n",
            temperature_in);

    fprintf(stdout, "    viscosity = %.10e Pa sec\n",
            ptr_visc->get_param_viscosity(mdensity_in,
                                          temperature_in));

    return;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::Lib95 iapws95eos;
    IAPWS::Lib08V iapws08visc;
    iapws08visc.print_header();
    iapws08visc.set_ptr_lib95(&iapws95eos);
    iapws08visc.set_range_visc2_mdensity(100., 450.);

    print_viscosity(&iapws08visc,
                    998., 298.15);
    print_viscosity(&iapws08visc,
                    1200., 298.15);

    fprintf(stdout, "\n");

    print_viscosity(&iapws08visc,
                    1000., 373.15);

    fprintf(stdout, "\n");

    print_viscosity(&iapws08visc,
                    1., 433.15);
    print_viscosity(&iapws08visc,
                    1000., 433.15);

    fprintf(stdout, "\n");

    print_viscosity(&iapws08visc,
                    1., 873.15);
    print_viscosity(&iapws08visc,
                    100., 873.15);
    print_viscosity(&iapws08visc,
                    600., 873.15);

    fprintf(stdout, "\n");

    print_viscosity(&iapws08visc,
                    1., 1173.15);
    print_viscosity(&iapws08visc,
                    100., 1173.15);
    print_viscosity(&iapws08visc,
                    400., 1173.15);

    fprintf(stdout, "\n");

    print_viscosity(&iapws08visc,
                    122., 647.35);
    print_viscosity(&iapws08visc,
                    222., 647.35);
    print_viscosity(&iapws08visc,
                    272., 647.35);
    /*
    print_viscosity(&iapws08visc,
                    322., 647.35);
    */
    print_viscosity(&iapws08visc,
                    372., 647.35);
    print_viscosity(&iapws08visc,
                    422., 647.35);

    return 0;
}
