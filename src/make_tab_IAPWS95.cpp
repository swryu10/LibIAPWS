#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS95.h"
#include"LibIAPWS97.h"

IAPWS::Lib97 *ptr_iapws97eos;

double get_pressure_sat(double temperature) {
    double press_sat =
        ptr_iapws97eos->get_param4_sat_pressure(temperature);

    return press_sat;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::n_iter_max_ = 1000000;
    IAPWS::eps_precision_ = 5. * 1.0e-8;
    IAPWS::Lib95 iapws95eos;
    iapws95eos.print_header();

    ptr_iapws97eos = new IAPWS::Lib97();
    double (*ptr_press_sat)(double) = &get_pressure_sat;

    int nbin_coex = 800;
    double temp_coex_max = 647.09;

    char filename_coex[100];
    strcpy(filename_coex, "./tab_coex_IAPWS95.txt");

    iapws95eos.make_tab_coex(nbin_coex,
                             temp_coex_max,
                             ptr_press_sat);
    fprintf(stdout, "\n");
    iapws95eos.export_tab_coex(filename_coex);
    iapws95eos.import_tab_coex(filename_coex);

    return 0;
}
