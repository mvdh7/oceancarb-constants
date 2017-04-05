function [pKstarT_w,pKstarT_P1,pKstarT_P2,pKstarT_P3] = cDSC07(T,S)
%cDSC07 Water & phosphoric acid dissociation constants, Total pH scale
% Inputs: <T> = temperature / K; <S> = salinity
% Source: Dickson et al., 2007. Guide to best practices for CO2
%          measurements, PICES Special Publication 3.
% Check values: T = 298.15 & S = 35 -> lnkw = -30.434, lnkP1 = -3.71, 
%                lnkP2 = -13.727, lnkP3 = -20.24
% MATLAB script written by Matthew P. Humphreys; last updated 2015-01-20

% These equations are all taken from DSC07, but they are originally from
%  M95; DSC07 have subtracted 0.015 from the constant term in each case,
%  to approximately convert pH SWS to Total.

%% WATER
% kw = [H+][OH-]
lnkw = (-13847.26./T) + 148.9652 - 23.6521*log(T) ...
    + ((118.67./T) - 5.977 + 1.0495*log(T)) .* S.^0.5 ...
    - 0.01615 * S;
pKstarT_w = -log10(exp(lnkw));

%% PHOSPHORIC ACID
% These equations, first presented by M95, are based on a composite of data
%  from KP67, DR79 and JW79. NB: the equations below are from DSC07.
% kP1 = [H+][H2PO4-]/[H3PO4]
lnkP1 = -4576.752./T + 115.525 - 18.453*log(T) ...
    + (-106.736./T + 0.69171) .* S.^0.5 ...
    + (-0.65643./T - 0.01844) .* S;
pKstarT_P1 = -log10(exp(lnkP1));
% kP2 = [H+][HPO42-]/[H2PO4-]
lnkP2 = -8814.715./T + 172.0883 - 27.927*log(T) ...
    + (-160.34./T + 1.3566) .* S.^0.5 ...
    + (0.37335./T - 0.05778) .* S;
pKstarT_P2 = -log10(exp(lnkP2));
% kP3 = [H+][PO43-]/[HPO42-]
lnkP3 = -3070.75./T - 18.141 ...
    + (17.27039./T + 2.81197) .* S.^0.5 ...
    + (-44.99486./T - 0.09984) .* S;
pKstarT_P3 = -log10(exp(lnkP3));

%% REFERENCES
%  KP67: Kester & Pytkowicz, 1967, Limnol Oceanogr 12(2).
%  DR79: Dickson & Riley, 1979, Mar Chem 7(2).
%         doi:10.1016/0304-4203(79)90002-1
%  JW79: Johansson & Wedborg, 1979, Mar Chem 8(1).
%         doi:10.1016/0304-4203(79)90032-X
%   M95: Millero, 1995, Geochim Cosmochim Acta 59(4).
%         doi:10.1016/0016-7037(94)00354-O
% DSC07: Dickson et al., 2007. Guide to best practices for CO2
%         measurements, PICES Special Publication 3.

end %function cDSC07