clear all
close all

LOCALPATH=pwd;
PARENTPATH=fileparts(LOCALPATH);
COMMONPATH='\CODE';
PATH=[PARENTPATH,COMMONPATH];
addpath(genpath(PATH));

%% data
run('CEA_interpolation')
clear p Mmol x Tcc

pt = 5.6382;
Power_supplied = 1.2; 

g0 = 9.81;
eps = 250;

de_max = 3e-2; % m
m = 0.5; % kg of water

% from graph shown
Power = 2.45e2*(pt/10).^1.48e-2; % for 1 mole/s kW*(s/mol)
moles = Power_supplied./Power;

% with less pressure, # moles increases

M = 18; % g/mol, molar mass H2O
mp = M*moles*1.e-3; % mass flow rate, kg/s

%% Tanks
tanks = tank(pt); % giving water storage pressure

%% Valves
R = 8.31; % J/molK

% Gases mass flow rate
Mmol_H2 = 2.01588;
Mmol_O2 = 31.9988;

mH2 = Mmol_H2/(Mmol_H2+0.5*Mmol_O2)*mp;
mO2 = Mmol_O2/(2*(Mmol_H2+0.5*Mmol_O2))*mp;

OF = Mmol_O2/(2*Mmol_H2);
cp_O2 = 0.918;
cp_H2 = 14.31;
% Specific heat ratio 
k_O2 = cp_O2/(cp_O2-R/Mmol_O2); 
k_H2 = cp_H2/(cp_H2-R/Mmol_H2);

Temp = 298;
roH2 = pt*1.e5/(R*1.e3/Mmol_H2*Temp); % pressure of the tank
roO2 = pt*1.e5/(R*1.e3/Mmol_O2*Temp); % kg/m^3

% %%ON_OFF VALVE
% Deltap_H2O = 0.01; %bar   COming from interpolation from datasheet
% %%CHECK VALVE
% SG_H2 = roH2/997;  %Ratio btw Density H2 and density Ration [-]  
% SG_O2 = roO2/997;  %Ratio btw Density O2 and density Ration [-]
% lm2gpm = 0.26417287472922;  % [l/min]->[gallons/min] Conversion
% Q_H2 = mH2./roH2*1000*60*lm2gpm; % H2 mass flow rate [gallons/min]
% Q_O2 = mO2./roO2*1000*60*lm2gpm; % O2 mass flow rate [gallons/min]
% CV = 0.19;% Fluid Coefficient [gallons/min]
% Deltap_H2 = ((SG_H2.*Q_H2.^2)./CV.^2); %Pressure drop [psi]
% Deltap_O2 = ((SG_O2.*Q_O2.^2)./CV.^2);
% Deltap_H2 = Deltap_H2*0.0689;  %Pressure conversion [psi]->[bar]
% Deltap_O2 = Deltap_O2*0.0689;

%% Injection
AinjH2 = 5.1760e-07; % m^2
AinjO2 = 2.5880e-07; % m^2
vH2 = mH2/(roH2*AinjH2); % m/s
vO2 = mO2/(roO2*AinjO2); % m/s
aH2 = sqrt(k_H2*R*1e3./Mmol_H2.*Temp);
aO2 = sqrt(k_O2*R*1e3./Mmol_O2.*Temp);
MH2 = vH2./aH2;
MO2 = vO2./aO2;

% injection losses
Cd_O2 = 0.85;
Cd_H2 = 0.2133;
DP_injO2 = 8./roO2.*(mO2./(Cd_O2.*AinjO2)).^2; % Injection plate pressure loss
DP_injH2 = 8./roH2.*(mH2./(Cd_H2.*AinjH2)).^2; % Both flows require the same losses
N_oxfu = 1; % configuration decided with swirling

% conversion
DP_injO2 = DP_injO2*1.e-5;
DP_injH2 = DP_injH2*1.e-5;
pcc = pt-(DP_injH2+DP_injO2)/2;

%% CEA 
% H2, O2, pcc

Tcc = ppval(f_T, pcc);
Mmol = ppval(f_M, pcc);
k = ppval(f_g, pcc);
vt = ppval(f_Sv, pcc); % [m/s]
ro_t = ppval(f_rho,pcc);

mu = [0.97057  1.04  1.0613  1.0738  1.0827  1.0895  1.0951  1.0998  1.1039  1.1075  1.1107]*1.e-4; % dependent only with T --> constant
f_mu = spline([1 10 20 30 40 50 60 70 80 90 100],mu);
mu = ppval(f_mu, pcc);

fi = k.*(2./(k+1)).^((k+1)./(k-1)); 
cstar = sqrt(R*1e3*Tcc./Mmol./fi);
At = 3.2969e-07;
Ae = eps.*At; % m^2
de = sqrt(Ae*4/pi)*1e3; % mm
dt = sqrt(At*4/pi)*1e3; % mm
Ae_max = de_max^2*pi/4;

%% CC design
L_star = 0.71; % overestimated, m
Mx = 0.006;
CC = CC_design(At,L_star,k,Mx);

%% Nozzle
% Characteristics
alfa = 15; % degrees
beta = 45;  % degrees
rcrt = 1.2; % curvature throat/throat radius

dcc = 6.4333;

mi = mu./ro_t;
Re_t = vt.*dt.*1.e-3./mi;
Re = sqrt(1/rcrt).*Re_t;

lambda = 0.5*(1+cosd(alfa)); % divergence loss factor
C_d = 1-((k+1)./2).^(0.75).*(3.266-(2.128./(k+1))).*Re.^(-0.5)+0.9428.*(k-1).*(k+2).*Re.^(-1)./(k+1).^0.5;
mp_real = mp.*C_d;

%% Performances
% pressure fed system chosen instead of turbopump due to the low dimensions
PeP0_guess = 1.e-6;
options = optimset('Display','off');
PeP0 = fsolve(@(x)sqrt(k*(2/(k+1))^((k+1)/(k-1)))/((x)^(1/k)*sqrt((2*k)/(k-1)*(1-x^((k-1)/k))))-eps,PeP0_guess,options);
Pe = PeP0.*pcc;

% ct = lambda*sqrt(2.*k.^2./(k-1).*(2./(k+1)).^((k+1)./(k-1)).*(1-(PeP0).^((k-1)./k)))+(Pe/pcc).*eps;
ct = lambda*sqrt(fi.*(2.*k./(k-1)).*(1-(PeP0).^((k-1)./k)))+(Pe/pcc).*eps;
Isp = cstar.*ct/g0;
T_real = Isp.*mp_real*g0; % difference due to C_d of nozzle in mp
t = m./mp;
vexit = sqrt(2*k./(k-1)*1.e3*R./Mmol.*Tcc.*(1-(PeP0.^((k-1)./k))));
Mexit = sqrt(2*k./(k-1).*(1./PeP0).^((k-1)./k));
T2 = vexit.*mp_real+Ae.*Pe*1.e5; % to confirm
Itot = T_real.*t;