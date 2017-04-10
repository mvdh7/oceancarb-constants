function pKstarF_HSO4 = cWM13(T,S)
%cWM13 The HSO4- stoichiometric dissociation constant <pKstarF_HSO4> 
%  on the Free pH scale [Waters & Millero, 2013, Mar Chem].
% Theoretical, Pitzer seawater model.
% <T> = temperature / K; <S> = salinity.
% MATLAB script written by Matthew P. Humphreys [2014-09-19].

%% Equation 29
% Constants (from Corrigendum, Table 6)
c1  =  4.24666         ;
c2  = -0.152671        ;
c3  =  2.67059   * 1e-2;
c4  = -4.2128    * 1e-5;
c5  =  0.2542181       ;
c6  = -5.09534   * 1e-3;
c7  =  7.1589    * 1e-4;
c8  = -2.91179   * 1e-3;
c9  =  2.09968   * 1e-5;
c10 = -4.03724   * 1e-5;
% Equation
log10KK = (c1 + c2*T + c3*T.*log(T) + c4*T.^2) .* S.^0.5 ...
        + (c5 + c6*T + c7*T.*log(T))           .* S      ...
        + (c8 + c9*T)                          .* S.^1.5 ...
        +  c10                                  * S.^2   ;

%% Equation (30) 
% Constants (taken directly from Clegg et al. [1994, J Chem Soc F Trans])
a1  =    562.69486        ;
a2  = -  102.5154         ;
a3  = -    1.117033 * 1e-4;
a4  =      0.2477538      ;
a5  = -13273.75           ;
% Equation
log10Ko = a1 + a2*log(T) + a3*T.^2 + a4*T + a5./T;

%% Calculate pKstar
pKstarF_HSO4 = -(log10KK + log10Ko);

end %function cWM13