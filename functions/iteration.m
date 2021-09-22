function y = iteration(p,Power_supplied)
% from graph shown
Power = 2.45e2*(p/10).^1.48e-2; % for 1 mole/s kW*(s/mol)
moles = Power_supplied./Power;

M = 18; % g/mol, molar mass H2O
mp = M*moles*1.e-3; % mass flow rate, kg/s

%% Injection
R = 8.31; % J/molK

% Gases mass flow rate
Mmol_H2 = 2.01588;
Mmol_O2 = 31.9988;
mO2 = Mmol_O2/(2*(Mmol_H2+0.5*Mmol_O2))*mp;
% Specific heat ratio 

% Inlet
Temp = 298;
vO2 = 35; % m/s
roO2 = p*1.e5/(R*1.e3/Mmol_O2*Temp); % kg/m^3
AinjO2 = mO2./(roO2.*vO2);   % m^2

% injector design
Cd_O2 = 0.85; % from book (Vince)
DP_injO2 = 8./roO2.*(mO2./(Cd_O2.*AinjO2)).^2; % Injection plate pressure loss

y = p-DP_injO2*1.e-5-4.65;

end







