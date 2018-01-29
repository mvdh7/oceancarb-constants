function pw = cWP02(tk)
% Calculate vapour pressure of pure water for 273 < tk/K < 647
%  following Wagner & Pruss (2002) via Dickson et al. (2007) CO2 guide
% Check value: pw = 3169.8 Pa @ 298.15 K
% Written by Matthew P. Humphreys [2018-01-19]

Tc = 647.096; % K
pc = 22.064; % MPa

theta = 1 - tk/Tc;

a1 = - 7.85951783;
a2 =   1.84408259;
a3 = -11.7866497 ;
a4 =  22.6807411 ;
a5 = -15.9618719 ;
a6 =   1.80122502;

ln_pp = Tc./tk .* (a1*theta + a2*theta.^1.5 + a3*theta.^3 ...
    + a4*theta.^3.5 + a5*theta.^4 + a6*theta.^7.5);

pw = pc*exp(ln_pp)*1e6; % Pa

end %function cWP02
