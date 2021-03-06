function [gDIC, bDIC, oDIC, gALK, bALK, oALK] = ESM10_CO2SYS( ...
    PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, ...
    PRESOUT, SI, PO4, pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS)
% ESM10_CO2SYS written by Matthew.Humphreys@nioz.nl [2020-04-10]
% Calculates marine carbonate system buffer factors defined by Egleston et
%  al. (2010) Global Biogeochem. Cy. 24, GB1002, doi:10.1029/2008GB003407
%  (henceforth ESM10).
% Note that ESM10's equations have fatal typographical errors.
% The equations have been recalculated to resolve the errors, as published
%  in Richier et al. (2018) Global Change Biol., doi:10.1111/gcb.14324
%  (henceforth RAH18) and Orr et al. (2018) Mar. Chem.,
%  doi:10.1016/j.marchem.2018.10.006. The correct equations are used here.
% This version of the code rather lazily uses CO2SYS to solve the marine
%  carbonate system for the calculations. It was built with CO2SYS v1.1:
%  van Heuven et al. (2011), doi:10.3334/CDIAC/otg.CO2SYS_MATLAB_v1.1
%  The inputs to ESM10_CO2SYS are therefore identical to those for CO2SYS
%  itself. It may also work with CO2SYS other versions, but this has not
%  been tested.

% Solve the MCS using CO2SYS
co2s = CO2SYS(PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, ...
    PRESIN, PRESOUT, SI, PO4, pHSCALEIN, K1K2CONSTANTS, KSO4CONSTANTS);

% Extract relevant results
DIC = co2s(:, 2);
H = 10.^-co2s(:, 3);
bicarb = co2s(:, 6);
carb = co2s(:, 7);
CO2aq = co2s(:, 8);
B4 = co2s(:, 9);
OH = co2s(:, 10);
KB = co2s(:, 59);

% Evaluate ESM10 subfunctions (ESM10 Table 1)
S = bicarb + 4*carb + H.*B4./(KB + H) + H + OH; % see OEDG18
P = 2*CO2aq + bicarb;
Q = bicarb - H.*B4./(KB + H) - H - OH; % see RAH18
AC = bicarb + 2*carb;

% Calculate buffer factors (ESM10 Table 1)
gDIC = DIC - AC.^2./S;
bDIC = (DIC.*S - AC.^2)./AC;
oDIC = DIC - AC.*P./Q; % corrected, see RAH18
gALK = (AC.^2 - DIC.*S)./AC;
bALK = AC.^2./DIC - S;
oALK = AC - DIC.*Q./P; % corrected as for oDIC (RAH18), bicarb => Q
