function [aou,aoupct,o2s] = cGG92(T,S,O2)
% calculates saturated oxygen from T and S in umol/kg
% based on Garcia & Gordon (1992)

% 'Combined fit' constants (umol/kg)
A0 = 5.80818;
A1 = 3.20684;
A2 = 4.11890;
A3 = 4.93845;
A4 = 1.01567;
A5 = 1.41575;
B0 = -0.00701211;
B1 = -0.00725958;
B2 = -0.00793334;
B3 = -0.00554491;
C0 = -0.000000132412;

% Scaled temperature <Ts>
Ts = log((298.15-T)./(273.15+T));

% Calculate saturated oxygen <o2s>
o2s = exp(A0 + A1*Ts + A2*Ts.^2 + A3*Ts.^3 + A4*Ts.^4 + A5*Ts.^5 ...
    + S.*(B0 + B1*Ts + B2*Ts.^2 + B3*Ts.^3) + C0*S.^2);

% Calculate AOU
aou = o2s - O2;
aoupct = 100 * aou./o2s;

end