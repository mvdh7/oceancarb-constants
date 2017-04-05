function [pKstarHF,valid] = cPF87(T,S)
%cPF87 Water dissociation constant pKstarT_w, Total pH scale
% Inputs: <T> = temperature / K; <S> = salinity
% Source: Perez & Fraga, 1987, Mar Chem 21(2).
%          doi:10.1016/0304-4203(87)90036-3
% Check value: T = 298.15 & S = 35 -> lnKHF = -6.09
% MATLAB script written by Matthew P. Humphreys; last updated 2015-01-20

lnbHF = -874./T - 0.111.*S.^0.5 + 9.68;
% lnbHFw_err = 0.05;

lnkHF = -lnbHF;
% lnkHFw = lnbHFw_err;

pKstarHF = -log10(exp(lnkHF));

valid = T >= 273.15+9 & T <= 273.15+33 & S >= 9 & S <= 33;

end %function cPF87