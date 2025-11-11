Here is the user guide for 1D Fu-liou model.

1. How to download the model codes?
Use the following command: git clone -b master --single-branch https://github.com/aerlab/courses_scripts.git

2. How to run the model?
At first, compile the fortran codes using command: make -f Makefile. An executable file "FU-wang" will be generated. Then, users can directly run this file (command: ./FU-wang) to get the radiative forcing calculation result, named "forcing.dat". Notably, once any codes (or head files) are modified, the whole model needs to be re-compiled after make clean -f Makefile.

3. Introduction of the model
FU-wang.f90 is the main code of this model. Other .f90 files contain different modules to be used. 
"zprofile" is a text file containing atmospheric profiles need to be read in the main code. The columns from left to right represent the profiles of atmospheric pressure (hPa), temperature (K), water vapor mixing ratio (kg/kg), ozone mixing ratio (kg/kg), AOD and aerosol type (1 for smoke, 2 for dust).
"forcing.dat" is the output text file listing total column AOD, downward shortwave (solar) flux at the surface, upward shortwave (solar) flux at the surface, downward shortwave (solar) flux at TOA, upward shortwave (solar) flux at TOA.
*.h files are header files included in *.f90 files. "Rad4S.h" defines some common parameters widely used in all Fortran codes. "cloud_rain.h" includes cloud properties, and "CKD.h" includes gas absorption data. "Aerosol.h" contains aerosol single-scattering optical properties for each aerosol type including asymmetry factor (a_asy), single-scattering albedo (a_ssa) and extinction cross section (a_ext). Usually, these header files contain some static data to be used.

4. Advanced usage
The subroutine "rad" in rad-FL.f90 is the main function used in the main code to calculate the radiative flux and heating rate profiles in both shortwave and longwave infrared. Users can get detailed descriptions about the input and output parameters for this function from the comments in rad-FL.f90 file.
The AOD profiles are adjusted in the current main code (Line 108-118) according to the user defined total column AOD (tau) and the ratio of dust AOD to smoke AOD. The vertical distribution of smoke/dust AOD is fixed and follows the design in "zprofile".
Another "Aerosol_3type.h" file is provided including another aerosol type, water soluble sulfuric acide relying on relative humidity. Users can apply it in the advanced application about aerosols with codes in "aero_cal-FL.f90" to be modified.   
