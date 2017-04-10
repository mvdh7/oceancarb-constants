function [TB,TB_err] = cU74(S)
%cU74 Total Boron <TB> in seawater / micromol/kg
% Input <S> = salinity.
% Source: Uppström, 1974, Deep-Sea Res 21(2).
%          doi:10.1016/0011-7471(74)90074-6
% "CO2 guide" refers to "Dickson, 2007, Guide to best practices for CO2
%  measurements"
% MATLAB script written by Matthew P. Humphreys, last updated 2015-01-20

B_Cl = 0.232; % mg/kg / ‰
B_Cl_err = 0.005; % mg/kg / ‰

B_mass = 10.811; % g/mol [CO2 guide]
B_mass_err = 0.007; % g/mol [CO2 guide]

TB = (1e-3*B_Cl/B_mass) * S/1.80655; % mol/kg; factor 1.80655
% to convert sal to chlorinity is from CO2 guide, and possibly 
% from Wooster et al., 1969, L&O 14(3).
TB_err = TB * sqrt((B_Cl_err/B_Cl).^2 ...
    + (B_mass_err/B_mass).^2);

end %function cU74