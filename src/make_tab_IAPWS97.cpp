#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS97.h"

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;
    IAPWS::n_iter_max_ = 1000000;
    IAPWS::eps_precision_ = 5. * 1.0e-8;
    IAPWS::Lib97 iapws97eos;
    iapws97eos.print_header();

    int nbin_coex = 100;

    char filename_coex[100];
    strcpy(filename_coex, "./tab_coex_IAPWS97.txt");

    iapws97eos.make_tab_coex(nbin_coex);
    fprintf(stdout, "\n");
    iapws97eos.export_tab_coex(filename_coex);
    iapws97eos.import_tab_coex(filename_coex);

    return 0;
}
