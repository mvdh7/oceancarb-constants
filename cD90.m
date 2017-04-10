function [pKstarT_B,pKstarF_HSO4,valid] = cD90(T,S)
%D90 Boric acid & bisulphate eq'm constants
% Sources:
%  Boric acid: Dickson, 1990, Deep-Sea Res Pt A 37(5).
%               doi:10.1016/0198-0149(90)90004-F [D90a]
%  Bisulphate: Dickson, 1990, J Chem Thermodyn 22(2).
%               doi:10.1016/0021-9614(90)90074-Z [D90b]
%         via: Dickson et al., 2007. Guide to best practices for CO2
%               measurements, PICES Special Publication 3. [DSC07]
% Inputs: <T> = temperature / K; <S> = salinity.
% MATLAB script written by Matthew P. Humphreys; last updated 2015-01-20

%% BORIC ACID - Total pH scale
% D90a eq'n 23
lnKB = (-8966.90   -  2890.53  *S.^0.5 - 77.942   *S ...
        + 1.728*S.^1.5 - 0.0996*S.^2) ./ T ...
     +    148.0248 +   137.1942*S.^0.5 +   1.62142*S ...
     - (   24.4344 +    25.085 *S.^0.5 +   0.2474 *S) .* log(T) ...
     +      0.053105*S.^0.5 .* T;
 
pKstarT_B = -log10(exp(lnKB));

%% BISULPHATE - Free pH scale
% D90b is paywalled, so equations are from DSC07:
I = 19.942*S ./ (1000 - 1.005*S); % Ionic strength
lnKS = -4276.1./T + 141.328 - 23.093*log(T) ...
    + ((-13856./T) + 324.57 - 47.986*log(T)) .* I.^0.5 ...
    + ((35474./T) - 771.54 + 114.723*log(T)) .* I ...
    - (2698./T) .* I.^1.5 ...
    + (1776./T) .* I.^2 ...
    + log(1 - 0.001005*S);
pKstarF_HSO4 = -log10(exp(lnKS));

%% VALIDITY
% D90a and D90b both use the same T & S ranges:
valid = T >= 273.15 & T <= 318.15 & S >= 5  & S <= 45;

end %function cD90