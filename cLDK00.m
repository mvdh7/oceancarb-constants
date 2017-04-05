function [pKstarT_C1,pKstarT_C2] = cLDK00(T,S)
%cLDK00 Carbonic acid stoichiometric equilibrium constants
%  <pKstarT_C1> and <pKstarT_C1>, Total pH scale
% Inputs: <T> = temperature / K; <S> = salinity
% Source: Lueker, Dickson & Keeling, 2000, Mar Chem 70(1-3).
%          doi:10.1016/S0304-4203(00)00022-0
% MATLAB script written by Matthew P. Humphreys; last updated 2015-01-20

pKstarT_C1 = 3633.86./T - 61.2172 + 9.6777.*log(T) ...
    - 0.011555*S + 0.0001152*S.^2;

pKstarT_C2 = 471.78./T + 25.929 - 3.16967.*log(T) ...
    - 0.01781*S + 0.0001122*S.^2;

end %function cLDK00