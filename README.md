# About
A C++ library to implement IAPWS (International Association for the Properties of Water and Steam) formulation.
Detailed documentations can be found in IAPWS webpage : https://iapws.org/

# Ingredients
The following releases are currently implemented.
* **IAPWS R1-76 (2014)** \
  **IAPWS::Lib76** class \
  Revised Release on Surface Tension of Ordinary Water Substance
* **IAPWS R6-95 (2018)** \
  **IAPWS::Lib95** class \
  Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use
* **IAPWS R10-06 (2009)** \
  **IAPWS::Lib06** class \
  Revised Release on the Equation of State 2006 for H2O Ice Ih
* **IAPWS R12-08 (2008)** \
  **IAPWS::Lib08V** class \
  Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance

# Build
This library can be built with cmake. \
One can build at a subdirectory with the following commands. \
&ensp;$ mkdir [subdirectory name] \
&ensp;$ cd [subdirectory name] \
&ensp;$ cmake [directory for the LibIAPWS local repository] \
&ensp;$ cmake --build .
