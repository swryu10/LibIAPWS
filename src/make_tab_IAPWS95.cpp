#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"LibIAPWS95.h"

int main(int argc, char *argv[]) {
    LibIAPWS95 iapws95eos;

    iapws95eos.make_tab_coex(700, 275., 625.);

    return 0;
}
