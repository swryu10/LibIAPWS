#include<math.h>
#include"LibIAPWS97.h"

namespace IAPWS {

void Lib97::print_header(FILE *ptr_fout) {
    fprintf(ptr_fout, "\n");
    fprintf(ptr_fout, "IAPWS R7-97(2012)\n");
    fprintf(ptr_fout, "Revised Release on the IAPWS Industrial Formulation 1997\n");
    fprintf(ptr_fout, "for the Thermodynamic Properties of Water and Steam\n");
    fprintf(ptr_fout, "\n");

    return;
}

} // end namespace IAPWS
