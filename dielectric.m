function D = dielectric(tk,pres)
% Dielectric constant of water
% tk = temperature / K
% pres = pressure / bar

%% Pitzer book chapter - Appendix E
u1 =  3.4279e+2;
u2 = -5.0866e-3;
u3 =  9.4690e-7;
u4 = -2.0525e00;
u5 =  3.1159e+3;
u6 = -1.8289e+2;
u7 = -8.0325e+3;
u8 =  4.2142e+6;
u9 =  2.1417e00;

e1000 = u1 * exp(u2*tk + u3*tk.^2);

uC = u4 + u5./(u6 + tk);
uB = u7 + u8./tk + u9*tk;

D = e1000 + uC .* log((uB + pres) ./ (uB + 1000));

%% Hi-T approach of Pitzer
% % rho = water density / g /cm3
% % equation follows DIPPR105
% Aw = 0.14395;
% Bw = 0.0112;
% Cw = 649.727;
% Dw = 0.05107;
% 
% rho = 1e-3 * Aw ./ Bw.^(1 + (1 - tk/Cw).^Dw);
% 
% g = 1 + 2.68*rho + 6.69*rho.^5 .* ((565./tk).^0.3 - 1);
% 
% RMM_H2O = 18.01528;
% ALPHA = 1.444e-24;
% MU = 1.84e-18;
% 
% F = (4*pi*avogadro*rho / (3*RMM_H2O)) ...
%     .* (ALPHA + MU^2*g ./ (3*boltzmann*1e7*tk));
% % 1e7 factor is to convert Boltzmann constant into erg/K (ie CGS units)
% 
% D = (1 + 9*F + sqrt((1 + 9*F).^2 + 8)) / 4;

end %function mcs_dielectric
