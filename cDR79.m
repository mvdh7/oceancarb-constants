function pKstarF_HF = cDR79(T,S)
%DR79 The HF stoichiometric dissociation constant <pKstarF_HF> 
%  on the Free pH scale [Dickson & Riley, 1979, Mar Chem].
% Ionic strength <I> calculated from DOE94.
% Empirical, synthetic seawater.
% <T> = temperature / K; <S> = salinity
% MATLAB script written by Matthew P. Humphreys [2014-09-19].

I = 19.924 * S ./ (1000 - 1.005 * S);

lnKF = 1590.2./T - 12.641 + 1.525*sqrt(I) + log(1 - 0.001005*S);

pKstarF_HF = -log10(exp(lnKF));

end %function cDR79