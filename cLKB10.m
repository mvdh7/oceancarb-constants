function TB = cLKB10(sal)
%calk_LKB10 Total boron <TB> in seawater / mol/kg
% Input <sal> = salinity; output <TB> = total boron / mol/kg
% Total boron = [boric acid] + [borate(-)]

B_sal = 0.1336;     % = boron/salinity / (mg/kg)/‰ (LKB10)
B = B_sal * sal;    % = boron / mg/kg
B_RAM = 10.811e3;   % = relative atomic mass of boron / mg/mol (DSC07)
TB = B/B_RAM;       % = total boron / mol/kg

end %function cLKB10
