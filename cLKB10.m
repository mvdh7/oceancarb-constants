function [TB,TB_err] = cLKB10(S)
%cLKB10 Total Boron <TB> in seawater / micromol/kg
% Input <S> = salinity.
% Source: Lee et al., 2010, GCA; doi:10.1016/j.gca.2009.12.027
% "CO2 guide" refers to "Dickson, 2007, Guide to best practices for CO2
%  measurements"
% MATLAB script written by Matthew P. Humphreys, last updated 2015-01-20

B_sal = 0.1336;     % = boron/salinity / (mg/kg)/‰
B = B_sal * S;      % = boron / mg/kg
B_RAM = 10.811e3;   % = relative atomic mass of boron / mg/mol [CO2 guide]
TB = 1e6 * B/B_RAM; % = total boron / micromol/kg

% Same calculation for uncertainty in TB, assuming perfect S
B_sal_err = 0.0005; % / (mg/kg)/‰ [Lee et al.]
B_err = B_sal_err * S; % / mg/kg
B_RAM_err = 0.007e3; % / mg/mol [CO2 guide]
TB_err = TB .* sqrt((B_err./B).^2 + (B_RAM_err/B_RAM).^2);

end %function cLKB10