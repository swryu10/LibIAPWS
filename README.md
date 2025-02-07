LibIAPWS

A C++ library to implement IAPWS (International Association for the Properties of Water and Steam) formulation.
Detailed documentations can be found in IAPWS webpage : https://iapws.org/

Currently, the following releases are implemented.
  IAPWS R1-76(2014)
    Revised Release on Surface Tension of Ordinary Water Substance
  IAPWS R6-95(2018)
    Revised Release on the IAPWS Formulation 1995
    for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use

This library can be built with cmake.
One can build at a subdirectory with the following commands.
  $ mkdir [subdirectory name]
  $ cd [subdirectory name]
  $ cmake [directory of LibIAPWS local repository]
  $ cmake --build .
