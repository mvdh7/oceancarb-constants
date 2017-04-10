function [Ksp_calc,Ksp_arag] = cM83(T,S)
%cM83 The solubility products of calcite and aragonite at one atmosphere 
% pressure in mol^2 kg^-2 [Mucci, 1983, Am J Sci]
% <T> = temperature / K; <S> = salinity
% MATLAB script written by Matthew P. Humphreys [2016-05-12].

TK = T + 273.15;

logKsp_calc = -171.9065 - 0.077993*TK + 2839.319./TK + 71.595*log10(TK) ...
    + (-0.77712 + 0.0028426*TK + 178.34./TK).*S.^0.5 ...
    - 0.07711*S + 0.0041249*S.^1.5;

logKsp_arag = -171.945 - 0.077993*TK + 2903.293./TK + 71.595*log10(TK) ...
    + (-0.068393 + 0.0017276*TK + 88.135./TK).*S.^0.5 ...
    - 0.10018*S + 0.0059415*S.^1.5;

Ksp_calc = 10.^logKsp_calc;
Ksp_arag = 10.^logKsp_arag;

end %function cM83