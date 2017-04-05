function [pKstarSWS_C1,pKstarSWS_C2] = cGP89(T,S)
%cGP89 The carbonic acid stoichiometric dissociation constants
%  <pKstarSWS_C1> and <pKstarSWS_C1> on the Seawater pH scale
%  [Goyet & Poisson, 1989, DSR].
% Empirical, synthetic seawater.
% <T> = temperature / K; <S> = salinity.
% MATLAB script written by Matthew P. Humphreys [2014-09-19].

pKstarSWS_C1 =  812.27./T + 3.356 - 0.00171*S.*log(T) + 0.000091*S.^2;
pKstarSWS_C2 = 1450.87./T + 4.604 - 0.00385*S.*log(T) + 0.000182*S.^2;

end %function cGP89