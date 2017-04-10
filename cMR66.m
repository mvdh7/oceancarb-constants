function [TSO4,TSO4_err,TBr,TBr_err] = cMR66(sal)
%cMR66 Total Sulphate <TSO4> and Bromide [TBr] in seawater / mol/kg
% Input <sal> = salinity
% Source: Morris & Riley, 1966, DSR; doi:10.1016/0011-7471(66)90601-2
% "CO2 guide" refers to "Dickson, 2007, Guide to best practices for CO2
%  measurements"
% MATLAB script written by Matthew P. Humphreys, last updated 2015-01-20

% SULPHATE
SO4_Cl = 0.14000;
SO4_Cl_err = 0.00023;

SO4_mass = 32.065 + 4*15.999; % g/mol [CO2 guide]
SO4_mass_err = sqrt(0.005^2 + 4*0.001^2); % g/mol [CO2 guide]

TSO4 = (SO4_Cl/SO4_mass) * sal/1.80655; % mol/kg; factor 1.80655
% to convert sal to chlorinity is from CO2 guide, and possibly 
% from Wooster et al., 1969, L&O 14(3).
TSO4_err = TSO4 * sqrt((SO4_Cl_err/SO4_Cl).^2 ...
    + (SO4_mass_err/SO4_mass).^2);

% BROMIDE
Br_Cl = 0.003473;
Br_Cl_err = 0.000012;

Br_mass = 79.904; % g/mol [CO2 guide]
Br_mass_err = 0.001; % g/mol [CO2 guide]

TBr = (Br_Cl/Br_mass) * sal/1.80655; % mol/kg; factor 1.80655
% to convert sal to chlorinity is from CO2 guide, and possibly 
% from Wooster et al., 1969, L&O 14(3).
TBr_err = TBr * sqrt((Br_Cl_err/Br_Cl).^2 ...
    + (Br_mass_err/Br_mass).^2);

end %function cMR66