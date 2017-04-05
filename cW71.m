function [TF,TF_err] = cW71(sal)
%cW71 Total Fluoride <TF> in seawater / mol/kg
% Input <sal> = salinity
% Source: Warner, 1971, DSR; doi:10.1016/0011-7471(71)90030-1
%  NB: this is slightly different from CO2 guide recommendation
% "CO2 guide" refers to "Dickson, 2007, Guide to best practices for CO2
%  measurements"
% MATLAB script written by Matthew P. Humphreys, last updated 2015-01-20

% FLUORIDE
F_Cl = 6.75e-5;
F_Cl_err = 0.03e-5;

F_mass = 18.998; % g/mol [CO2 guide]
F_mass_err = 0.001; % g/mol [CO2 guide]

TF = (F_Cl/F_mass) * sal/1.80655; % mol/kg; factor 1.80655
% to convert sal to chlorinity is from CO2 guide, and possibly 
% from Wooster et al., 1969, L&O 14(3).
TF_err = TF * sqrt((F_Cl_err/F_Cl).^2 ...
    + (F_mass_err/F_mass).^2);

end %function cW71