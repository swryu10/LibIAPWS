#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS06.h"

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::n_iter_max_ = 1000000;
    IAPWS::eps_precision_ = 5. * 1.0e-8;
    IAPWS::Lib95 iapws95eos;

    char filename_coex_liq[100];
    strcpy(filename_coex_liq, "./tab_coex_IAPWS95.txt");
    iapws95eos.import_tab_coex(filename_coex_liq);

    IAPWS::Lib06 iapws06eos;
    iapws06eos.print_header();

    int nbin_coex = 300;
    double temp_coex_min = 123.16;

    char filename_coex_ice[100];
    strcpy(filename_coex_ice, "./tab_coex_IAPWS06.txt");

    iapws06eos.make_tab_coex(&iapws95eos,
                             nbin_coex,
                             temp_coex_min);
    fprintf(stdout, "\n");
    iapws06eos.export_tab_coex(filename_coex_ice);
    iapws06eos.import_tab_coex(filename_coex_ice);

    int nbin_melt = 300;
    double press_melt_max = 2.1e+8;

    char filename_melt_ice[100];
    strcpy(filename_melt_ice, "./tab_melt_IAPWS06.txt");

    iapws06eos.make_tab_melt(&iapws95eos,
                             nbin_melt,
                             press_melt_max);
    fprintf(stdout, "\n");
    iapws06eos.export_tab_melt(filename_melt_ice);
    iapws06eos.import_tab_melt(filename_melt_ice);

    return 0;
}
