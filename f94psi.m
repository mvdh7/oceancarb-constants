function PSI = f94psi(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMP,PRES,SI,PHOS)
% Frankignoulle et al. (1994)

co2s = CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMP,TEMP,PRES,PRES, ...
    SI,PHOS,3,10,3);

AC = (co2s(:,6) + 2*co2s(:,7)) * 1e-6;
H  = co2s(:,28) * 1e-6;
K1 = co2s(:,54);
K2 = co2s(:,55);
KH = co2s(:,58);
KB = co2s(:,59);
TB = co2s(:,79) * 1e-6;

P = KB.*TB./(H + KB).^2 + KH./H.^2 + 1;
Q = H + 2*K2;
R = AC.*Q + 2*K2.*AC + P.*Q.*H;
S = (K2.*AC + P.*Q.*(H+K2) + R.*H./K1) ...
    ./ -(Q.*(1 + 2*H./K1).*H);

PSI = -(R + 2*Q.*S.*H) ...
    .* (K2.*AC + P.*Q.*(H+K2) + R.*H./K1) ...
    ./ (R.*S.*Q.^2);

end %function f94psi