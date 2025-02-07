#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS95.h"

int main(int argc, char *argv[]) {
    IAPWS::Lib95 iapws95eos;

    int nbin_coex = 700;
    double temp_coex_min = 275.;
    double temp_coex_max = 625.;

    char filename_coex[100];
    strcpy(filename_coex, "./tab_coex_IAPWS95.txt");

    iapws95eos.make_tab_coex(nbin_coex,
                             temp_coex_min,
                             temp_coex_max);
    iapws95eos.export_tab_coex(filename_coex);
    iapws95eos.import_tab_coex(filename_coex);

    return 0;
}
