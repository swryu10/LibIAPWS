#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"LibIAPWS97.h"

void print_diff(int nbin_t, double *t_bin,
                int nbin_p, double *p_bin,
                double **tab_diff,
                char *filename_diff) {
    FILE *ptr_fout = fopen(filename_diff, "w");
    if (ptr_fout == NULL) {
        return;
    }

    for (int it = 0; it <= nbin_t; it++) {
        for (int ip = 0; ip <= nbin_p; ip++) {
            fprintf(ptr_fout, "    %e    %e    ",
                    t_bin[it], p_bin[ip]);
            fprintf(ptr_fout, "%e\n",
                    tab_diff[it][ip]);
        }

        if (it == nbin_t) {
            continue;
        }

        fprintf(ptr_fout, "\n");
    }

    fclose(ptr_fout);

    return;
}

int main(int argc, char *argv[]) {
    IAPWS::flag_verbose_ = true;

    IAPWS::Lib97 iapws97eos;
    iapws97eos.print_header();

    char filename_coex97[100];
    strcpy(filename_coex97, "./tab_coex_IAPWS97.txt");
    iapws97eos.import_tab_coex(filename_coex97);

    int nbin_temperature = 96;
    double temperature_min = 273.15;
    double temperature_max = 1073.15;
    double d_temperature =
        (temperature_max - temperature_min) /
        static_cast<double>(nbin_temperature);
    double *temperature_bin =
        new double[nbin_temperature + 1];
    for (int it = 0; it <= nbin_temperature; it++) {
        temperature_bin[it] = temperature_min +
            d_temperature * static_cast<double>(it);
    }

    int nbin_pressure = 128;
    double pressure_min = 611.657;
    double pressure_max = 99. * 1.0e+6;
    double *pressure_bin =
        new double[nbin_pressure + 1];
    double d_pressure =
        log(pressure_max / pressure_min);
    for (int ip = 0; ip <= nbin_pressure; ip++) {
        pressure_bin[ip] = pressure_min *
            exp(d_pressure * static_cast<double>(ip) /
                             static_cast<double>(nbin_pressure));
    }

    double **tab_bkwd_Tph =
        new double *[nbin_temperature + 1];
    double **tab_bkwd_Tps =
        new double *[nbin_temperature + 1];
    for (int it = 0; it <= nbin_temperature; it++) {
        tab_bkwd_Tph[it] =
            new double [nbin_pressure + 1];
        tab_bkwd_Tps[it] =
            new double [nbin_pressure + 1];

        for (int ip = 0; ip <= nbin_pressure; ip++) {
            tab_bkwd_Tph[it][ip] = 0.;
            tab_bkwd_Tps[it][ip] = 0.;
        }
    }

    fprintf(stderr, "temperature (K)\n");
    fprintf(stdout, "  pressure (Pa)");
    fprintf(stdout, "  n_reg");
    fprintf(stdout, "  bkwd_Tph");
    fprintf(stdout, "  bkwd_Tps\n");
    for (int it = 0; it <= nbin_temperature; it++) {
        fprintf(stdout, "  %e\n", temperature_bin[it]);
        for (int ip = 0; ip <= nbin_pressure; ip++) {
            int n_reg =
                iapws97eos.get_region(temperature_bin[it],
                                      pressure_bin[ip]);
            if (n_reg == 0) {
                continue;
            }

            double enthalpy_in =
                iapws97eos.get_param_enthalpy(temperature_bin[it],
                                              pressure_bin[ip]);
            double Tph =
                iapws97eos.get_param_temperature_ph(pressure_bin[ip],
                                                    enthalpy_in);
            tab_bkwd_Tph[it][ip] =
                (temperature_bin[it] - Tph) /
                (temperature_bin[it] + Tph);

            double entropy_in =
                iapws97eos.get_param_entropy(temperature_bin[it],
                                             pressure_bin[ip]);
            double Tps =
                iapws97eos.get_param_temperature_ps(pressure_bin[ip],
                                                    entropy_in);
            tab_bkwd_Tps[it][ip] =
                (temperature_bin[it] - Tps) /
                (temperature_bin[it] + Tps);

            fprintf(stdout, "    %e    %d",
                    pressure_bin[ip], n_reg);
            fprintf(stdout, "    %e    %e",
                    tab_bkwd_Tph[it][ip],
                    tab_bkwd_Tps[it][ip]);
            fprintf(stdout, "\n");
        }
    }

    char fname_bkwd_Tph[100];
    strcpy(fname_bkwd_Tph,
           "./bkwd_IAPWS97_Tph.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_bkwd_Tph,
               fname_bkwd_Tph);

    char fname_bkwd_Tps[100];
    strcpy(fname_bkwd_Tps,
           "./bkwd_IAPWS97_Tps.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_bkwd_Tps,
               fname_bkwd_Tps);

    for (int it = 0; it <= nbin_temperature; it++) {
        delete [] tab_bkwd_Tph[it];
        delete [] tab_bkwd_Tps[it];
    }
    delete [] tab_bkwd_Tph;
    delete [] tab_bkwd_Tps;

    return 0;
}
