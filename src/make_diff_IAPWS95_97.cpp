#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"LibIAPWS95.h"
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

    IAPWS::Lib95 iapws95eos;
    iapws95eos.print_header();

    IAPWS::Lib97 iapws97eos;
    iapws97eos.print_header();

    char filename_coex95[100];
    strcpy(filename_coex95, "./tab_coex_IAPWS95.txt");
    iapws95eos.import_tab_coex(filename_coex95);

    char filename_coex97[100];
    strcpy(filename_coex97, "./tab_coex_IAPWS97.txt");
    iapws97eos.import_tab_coex(filename_coex97);

    int nbin_temperature = 128;
    double temperature_min = 273.15;
    double temperature_max = 2273.15;
    double d_temperature =
        (temperature_max - temperature_min) /
        static_cast<double>(nbin_temperature);
    double *temperature_bin =
        new double[nbin_temperature + 1];
    for (int it = 0; it <= nbin_temperature; it++) {
        temperature_bin[it] = temperature_min +
            d_temperature * static_cast<double>(it);
    }

    int nbin_pressure = 256;
    double pressure_min = 611.657;
    double pressure_max = 99. * 1.0e+6;
    double *pressure_bin =
        new double[nbin_pressure + 1];
    double d_pressure =
        log(pressure_max / pressure_min) /
        static_cast<double>(nbin_temperature);
    for (int ip = 0; ip <= nbin_pressure; ip++) {
        pressure_bin[ip] = pressure_min *
            exp(d_pressure * static_cast<double>(ip));
    }

    double **tab_diff_mdensity =
        new double *[nbin_temperature + 1];
    double **tab_diff_erg_int =
        new double *[nbin_temperature + 1];
    double **tab_diff_entropy =
        new double *[nbin_temperature + 1];
    double **tab_diff_enthalpy =
        new double *[nbin_temperature + 1];
    double **tab_diff_heat_c_p =
        new double *[nbin_temperature + 1];
    double **tab_diff_heat_c_v =
        new double *[nbin_temperature + 1];
    for (int it = 0; it <= nbin_temperature; it++) {
        tab_diff_mdensity[it] =
            new double [nbin_pressure + 1];
        tab_diff_erg_int[it] =
            new double [nbin_pressure + 1];
        tab_diff_entropy[it] =
            new double [nbin_pressure + 1];
        tab_diff_enthalpy[it] =
            new double [nbin_pressure + 1];
        tab_diff_heat_c_p[it] =
            new double [nbin_pressure + 1];
        tab_diff_heat_c_v[it] =
            new double [nbin_pressure + 1];

        for (int ip = 0; ip <= nbin_pressure; ip++) {
            tab_diff_mdensity[it][ip] = 0.;
            tab_diff_erg_int[it][ip] = 0.;
            tab_diff_entropy[it][ip] = 0.;
            tab_diff_enthalpy[it][ip] = 0.;
            tab_diff_heat_c_p[it][ip] = 0.;
            tab_diff_heat_c_v[it][ip] = 0.;
        }
    }

    fprintf(stderr, "temperature (K)\n");
    fprintf(stdout, "  pressure (Pa)  ");
    fprintf(stdout, "n_reg  ");
    fprintf(stdout, "diff_mdensity  ");
    fprintf(stdout, "diff_erg_int  ");
    fprintf(stdout, "diff_entropy  ");
    fprintf(stdout, "diff_enthalpy  ");
    fprintf(stdout, "diff_heat_c_p  ");
    fprintf(stdout, "diff_heat_c_v");
    fprintf(stdout, "\n");
    for (int it = 0; it <= nbin_temperature; it++) {
        fprintf(stdout, "  %e\n", temperature_bin[it]);
        for (int ip = 0; ip <= nbin_pressure; ip++) {
            int n_reg =
                iapws97eos.get_region(temperature_bin[it],
                                      pressure_bin[ip]);
            if (n_reg == 0) {
                continue;
            }

            double mden97 =
                iapws97eos.get_param_mdensity(temperature_bin[it],
                                              pressure_bin[ip]);
            double mden95 = mden97;
            bool found_mden =
                iapws95eos.find_root_mdensity(temperature_bin[it],
                                              pressure_bin[ip],
                                              mden95);
            tab_diff_mdensity[it][ip] =
                (mden97 - mden95) / fabs(mden97 + mden95);

            double u97 =
                iapws97eos.get_param_erg_int(temperature_bin[it],
                                             pressure_bin[ip]);
            double u95 =
                iapws95eos.get_param_erg_int(mden95,
                                             temperature_bin[it]);
            tab_diff_erg_int[it][ip] =
                (u97 - u95) / fabs(u97 + u95);

            double s97 =
                iapws97eos.get_param_entropy(temperature_bin[it],
                                             pressure_bin[ip]);
            double s95 =
                iapws95eos.get_param_entropy(mden95,
                                             temperature_bin[it]);
            tab_diff_entropy[it][ip] =
                (s97 - s95) / fabs(s97 + s95);

            double h97 =
                iapws97eos.get_param_enthalpy(temperature_bin[it],
                                              pressure_bin[ip]);
            double h95 =
                iapws95eos.get_param_enthalpy(mden95,
                                              temperature_bin[it]);
            tab_diff_enthalpy[it][ip] =
                (h97 - h95) / fabs(h97 + h95);

            double cp97 =
                iapws97eos.get_param_heat_c_p(temperature_bin[it],
                                              pressure_bin[ip]);
            double cp95 =
                iapws95eos.get_param_heat_c_p(mden95,
                                              temperature_bin[it]);
            tab_diff_heat_c_p[it][ip] =
                (cp97 - cp95) / fabs(cp97 + cp95);

            double cv97 =
                iapws97eos.get_param_heat_c_v(temperature_bin[it],
                                              pressure_bin[ip]);
            double cv95 =
                iapws95eos.get_param_heat_c_v(mden95,
                                              temperature_bin[it]);
            tab_diff_heat_c_v[it][ip] =
                (cv97 - cv95) / fabs(cv97 + cv95);

            fprintf(stdout, "    %e    %d    ",
                    pressure_bin[ip], n_reg);
            fprintf(stdout, "%e    %e    %e    %e    %e    %e",
                    tab_diff_mdensity[it][ip],
                    tab_diff_erg_int[it][ip],
                    tab_diff_entropy[it][ip],
                    tab_diff_enthalpy[it][ip],
                    tab_diff_heat_c_p[it][ip],
                    tab_diff_heat_c_v[it][ip]);
            fprintf(stdout, "\n");
        }
    }

    char fname_diff_mdensity[100];
    strcpy(fname_diff_mdensity,
           "./diff_IAPWS95_97_mdensity.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_mdensity,
               fname_diff_mdensity);

    char fname_diff_entropy[100];
    strcpy(fname_diff_entropy,
           "./diff_IAPWS95_97_entropy.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_entropy,
               fname_diff_entropy);

    char fname_diff_enthalpy[100];
    strcpy(fname_diff_enthalpy,
           "./diff_IAPWS95_97_enthalpy.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_enthalpy,
               fname_diff_enthalpy);

    char fname_diff_erg_int[100];
    strcpy(fname_diff_erg_int,
           "./diff_IAPWS95_97_erg_int.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_erg_int,
               fname_diff_erg_int);

    char fname_diff_heat_c_p[100];
    strcpy(fname_diff_heat_c_p,
           "./diff_IAPWS95_97_heat_c_p.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_heat_c_p,
               fname_diff_heat_c_p);

    char fname_diff_heat_c_v[100];
    strcpy(fname_diff_heat_c_v,
           "./diff_IAPWS95_97_heat_c_v.txt");
    print_diff(nbin_temperature, temperature_bin,
               nbin_pressure, pressure_bin,
               tab_diff_heat_c_v,
               fname_diff_heat_c_v);

    for (int it = 0; it <= nbin_temperature; it++) {
        delete [] tab_diff_mdensity[it];
        delete [] tab_diff_erg_int[it];
        delete [] tab_diff_entropy[it];
        delete [] tab_diff_enthalpy[it];
        delete [] tab_diff_heat_c_p[it];
        delete [] tab_diff_heat_c_v[it];
    }
    delete [] tab_diff_mdensity;
    delete [] tab_diff_erg_int;
    delete [] tab_diff_entropy;
    delete [] tab_diff_enthalpy;
    delete [] tab_diff_heat_c_p;
    delete [] tab_diff_heat_c_v;

    return 0;
}
