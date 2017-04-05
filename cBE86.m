function [pKstarSWS_AMP] = cBE86(T,S)
%cBE86 The 2-aminopyridine stoichiometric dissociation constant
%  <pKstarSWS_AMP> on the Seawater pH scale [Bates & Erickson, 1989, JSC].
% Empirical, synthetic seawater.
% <T> = temperature / K; <S> = salinity.
% MATLAB script written by Matthew P. Humphreys [2014-09-19].

pKstarSWS_AMP = 2498.31./T - 15.3274 + 2.4050*log(T) ...
    + S.*(0.012929 - 2.9417e-5*T);

end %function cBE86