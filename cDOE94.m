function [pKstarT_W,I,pKstarF_HSO4] = cDOE94(T,S)
%cDOE94 The ion product of water <pKstarT_W> on the Total pH scale,
%  the HSO4- stoichiometric dissociation constant <pKstarF_HSO4> 
%  on the Free pH scale, and ionic strength from salinity <I> [DOE, 1994],
%  as reprinted by Zeebe & Wolf Gladrow [2001].
% <T> = temperature / K; <S> = salinity.
% 
% MATLAB script written by Matthew P. Humphreys [2014-09-19].

%% Ionic strength
I = 19.924 * S ./ (1000 - 1.005 * S);

%% HSO4
lnKS =  - 4276.1./T + 141.328 -  23.093*log(T)            ...
     + (-13856  ./T + 324.57  -  47.986*log(T)) .* I.^0.5 ...
     + ( 35474  ./T - 771.54  + 114.723*log(T)) .* I      ...
     - (  2698  ./T)                            .* I.^1.5 ...
     + (  1776  ./T)                            .* I.^2   ...
     + log(1 - 0.001005*S);

pKstarF_HSO4 = -log10(exp(lnKS));

%% H2O
lnKW = 148.96502 - 13847.26./T - 23.6521*log(T) ...
    + (118.67./T - 5.977 + 1.0495*log(T)) .* S.^0.5 - 0.01615*S;

pKstarT_W = -log10(exp(lnKW));

end %function cDOE94