%D. G. Partridge, CEMPS, University of Exeter, Feb 2020
close all
clear all
clc

%! Program to solve growth of monodisperse cloud droplet population in an ascending air parcel

%!--------------------------------------------------------------------------

%!Initialise constants (from Table A1. Devenish et al., 2016)

Pi = 3.141592653589793238462643  %! In Matlab, can also simply use pi
g = 9.81                         %! Acceleration due to gravity (m s^-2)
c_pa = 1005.0                    %! Specific heat capacity of dry air (J kg^-1 K^-1)
Rho_w = 1000.0                   %! Density of liquid water (Kg m^-3)
Rho_a = 1.225 			         %! Density of air (Kg m^-3)
Eps = 0.622                      %! Ratio of molecular masses of water vapour and dry air
Lv = 2.5e6                       %! Latent heat of vapourisation (J Kg^-1)
Ra = 287.0                       %! Gas constant of dry air (J kg^-1 K^-1)
Rv = 462.0                       %! Gas constant of water vapour (J kg^-1 K^-1)
k = 0.024                        %! Thermal Conductivity of Air (J m^-1 s^-1 K^-1)
Kv = 2.21e-5                     %! Diffusivity of Water Vapour (m^2 s^-1)
T = 280                          %temperature in kelvin
qv = Rho_w/Rho_a                 %mixing ratio
Rho_v = 0.598                    %density of water vapour in kg/m^3
e = Rho_v*Rv*T 
   	 		         %! Saturation vapor pressure (provided by function svp.m)
r = 8e-7 

%%Advanced: In reality the Thermal Conductivty & Diffusivity are not constant, but depend on the temperature and pressure. 
%They also need to include kinetic effects (i.e., the effects of condensation and accommodation coefficients; e.g., Fukuta and Walter 1970).
%Assuming the constant values above is adequate for your project. If you would like to account for the temperature dependance
%please do, you can find further details in the course textbook or in "Lohmann, an introduction to clouds", eq. 7.27, page 192.
 
%--------------------------------------------------------------------------
%Define initial vertical velocity and droplet number concentration as parameters

W = 5 	    %! Vertical velocity (m/s)
N = 10      %droplet density      %! Droplet Number Conc. (#/m^3)

%!--------------------------------------------------------------------------

%Start main program here.

s = e/svp(X) - 1

A3 = (((Lv^2 * Rho_w)/(k * Rv * T^2))+((Rho_w * Rv * T)/(Kv * svp(T))))^-1;  %constant for constant temperature

drdt = A3 * (s/r) %not defined

dqldt = ( (4*Pi*Rho_w*N)/Rho_a ) * r^2 * drdt


dsdt = ( g / ( Ra * T ) ) * ( ( (Lv * Ra) / ( c_pa * Rv * T)  ) - 1 ) * W - ...
    ( ( (Lv^2) / ( c_pa * Rv * T^2) ) + 1/qv ) *  dqldt 

h(1) = 0 ; %distance from cloud base in metres
hstop = 1000 ; %metres
dh = 0.1 ; %time step
nstop = round(hstop/dh) ; 

for n = 1:nstop

    h(n+1) = h(n) + dh ;
    r = r + dh * drdt
    s = s + dh * dsdt

    r(n+1) = r(n) 
    s(n+1) = s(n)

end

%plotting








%End of main program here
%---------------------------------------------------------------------------
