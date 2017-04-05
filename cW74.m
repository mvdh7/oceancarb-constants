function [k0,ln_k0] = cW74(t,s)
% 2016-04-27
% Citation: Weiss (1974) Mar Chem 2
% k0 units: mol /kg /atm

tk = t + 273.15;

ln_k0 = 9345.17./tk - 60.2409 + 23.3585*log(tk/100) ...
    + s.*(0.023517 - 0.00023656*tk + 0.0047036*(tk/100).^2);
k0 = exp(ln_k0);

end %function cW74