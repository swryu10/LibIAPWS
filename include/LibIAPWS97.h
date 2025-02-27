#ifndef _LIBIAPWS97_H_
#define _LIBIAPWS97_H_

#include<stdio.h>
#include"BaseIAPWS.h"

namespace IAPWS {

/* implementation of IAPWS R7-97(2012)
 * Revised Release on the IAPWS Industrial Formulation 1997
 * for the Thermodynamic Properties of Water and Steam */
class Lib97 {
  private :


  public :

    Lib97() {

        return;
    }

    ~Lib97() {}

    void print_header(FILE *ptr_fout = stdout);
};

} // end namespace IAPWS

#endif
