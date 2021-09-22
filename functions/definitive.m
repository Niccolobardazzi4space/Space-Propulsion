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

pcc = 4.65;
% p=linspace(1,100,100);
Power_supplied = 1;

options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-16);
pt = fsolve(@(x) iteration(x,Power_supplied),6,options);

g0 = 9.81;
eps = 250;

de_max = 3e-2; % m
m = 0.5; % kg of water

% from graph shown
Power = 2.45e2*(pcc/10).^1.48e-2; % for 1 mole/s kW*(s/mol)
moles = Power_supplied./Power;

% with less pressure, # moles increases

M = 18; % g/mol, molar mass H2O
mp = M*moles*1.e-3; % mass flow rate, kg/s

%% Tanks
tanks = tank(pt); % giving water storage pressure

%% Injection
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

% Inlet
Temp = 298;
vH2 = 35; % m/s
vO2 = 35; % m/s
roH2 = pt*1.e5/(R*1.e3/Mmol_H2*Temp); % pressure of the tank
roO2 = pt*1.e5/(R*1.e3/Mmol_O2*Temp); % kg/m^3
aH2 = sqrt(k_H2*R*1e3./Mmol_H2.*Temp);
aO2 = sqrt(k_O2*R*1e3./Mmol_O2.*Temp);
MH2 = vH2./aH2;
MO2 = vO2./aO2;
AinjH2 = mH2./(roH2.*vH2);   % m^2
dinjH2 = (4*AinjH2/pi).^0.5; % m
AinjO2 = mO2./(roO2.*vO2);   % m^2
dinjO2 = (4*AinjO2/pi).^0.5; % m

% injector design
Cd_O2 = 0.85; % from book (Vince)
DP_injO2 = 8./roO2.*(mO2./(Cd_O2.*AinjO2)).^2; % Injection plate pressure loss
DP_injH2 = DP_injO2; % Both flows require the same losses
Cd_H2 = mH2./(AinjH2).*sqrt(8./(roH2.*DP_injH2));
N_oxfu = 1; % configuration decided with swirling

%% VALVES
Afeed_H2=(((pi/4)).*(4.*dinjH2).^2); %Area of H2 feed lines m^2
Dfeed_H2=sqrt(4.*Afeed_H2./pi); %Diameter of H2 fedd lines   m
Afeed_O2=(((pi./4)).*(4.*dinjO2).^2); %Area of O2 feed lines m^2
Dfeed_O2=sqrt(4.*Afeed_O2./pi); %Diameter of H2 fedd lines   m
D_H2O = 6.35/1000; %m Fixed to the arianegroup valve
A_H2O = pi.*D_H2O.^2/4; %Area of water feed line (before PEM) m^2
% v_H2O = (roH2.*Afeed_H2.*v_H2+roO2.*Afeed_O2.*v_O2)./(997.*A_H2O); %m/s ...
% v_H2O_check = mp./(997.*A_H2O);

%%ON_OFF VALVE
bar_to_psi = 14.5038; %conversion factor bar to psi
lm2gpm = 0.26417287472922;  % [l/min]->[gallons/min] Conversion
%Data
deltaP_data_onoff = 0.15*bar_to_psi; %psi, pressure drop taken from datasheet
mdot_data_onoff = 4.5*1e-3; %kg/s, mass flow rate corresponding to the pressure drop taken from datasheet H2O
Q_data_onoff = mdot_data_onoff/997*1000*60*lm2gpm; %flow rate of the pressure drop
SG_data_onoff = 1; %Ratio btw Density H2O and density of fluid [-]  
Cv_onoff = sqrt(SG_data_onoff*Q_data_onoff^2/deltaP_data_onoff); %flow coefficient of ON/OFF valve

%Our engine
Q_H2O = mp./997*1000*60*lm2gpm; % H2O volume flow rate [gallons/min]
SG_H2O = 1;
deltaP_onoff = ((SG_H2O.*Q_H2O.^2)/Cv_onoff.^2) ./ bar_to_psi; %actual P loss ONOFF bar

%%CHECK VALVE
%Data
SG_air = 1.225/997; %Air used in graph
scfm_to_gpm = 74.8052; %conversion factor scfm to gpm
Q_data_check = 10*scfm_to_gpm; %flow rate of the pressure drop graph, gpm
deltaP_data_check = 4; %pressure drop indentified on the graph
Cv_check = sqrt(SG_air*Q_data_check^2/deltaP_data_check); %Cv of the check valve

%Our engine
SG_H2 = roH2/997;  %Ratio btw Density H2 and density of hydrogen [-]  
SG_O2 = roO2/997;  %Ratio btw Density O2 and density of oxygen [-]
Q_H2 = mH2./roH2.*1000.*60.*lm2gpm; % H2 mass flow rate [gallons/min]
Q_O2 = mO2./roO2.*1000.*60.*lm2gpm; % O2 mass flow rate [gallons/min]
deltaP_H2_check = ((SG_H2.*Q_H2.^2)./Cv_check^2) ./ bar_to_psi; %actual P loss check bar
deltaP_O2_check = ((SG_O2.*Q_O2.^2)./Cv_check^2) ./ bar_to_psi; %actual P loss check bar

% figure()
% plot(p,deltaP_onoff)
% title('ON/OFF Valve Pressure Loss','Interpreter','latex','FontSize',14)
% xlabel('Pressure of storage tank [bar]','Interpreter','latex','FontSize',14)
% ylabel('Pressure Loss water line[bar]','Interpreter','latex','FontSize',14)
% figure()
% plot(p,deltaP_H2_check)
% title('Check Valve H2 Pressure Loss','Interpreter','latex','FontSize',14)
% xlabel('Pressure of storage tank [bar]','Interpreter','latex','FontSize',14)
% ylabel('Pressure Loss H2 line [bar]','Interpreter','latex','FontSize',14)
% figure()
% plot(p,deltaP_O2_check)
% title('Check Valve O2 Pressure Loss','Interpreter','latex','FontSize',14)
% xlabel('Pressure of storage tank [bar]','Interpreter','latex','FontSize',14)
% ylabel('Pressure Loss O2 line [bar]','Interpreter','latex','FontSize',14)

%% CEA 
% H2, O2, pcc

% conversion
% DP_injO2 = DP_injO2*1.e-5;
% Dp = Deltap_O2+DP_injO2;
% p = p-Dp;

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
At = cstar.*mp./(pcc*1e5); % m^2
Ae = eps.*At; % m^2
de = sqrt(Ae*4/pi)*1e3; % mm
dt = sqrt(At*4/pi)*1e3; % mm
Ae_max = de_max(1)^2*pi/4;

%% CC design
L_star = 0.71; % overestimated, m
Mx = 0.006;
CC = CC_design(At,L_star,k,Mx);

%% Nozzle
% Characteristics
alfa = 15; % degrees
beta = 45;  % degrees
rcrt = 1.2; % curvature throat/throat radius

dcc = CC.dc;

mi = mu./ro_t;
Re_t = vt.*dt.*1.e-3./mi;
Re = sqrt(1/rcrt).*Re_t;

L_c = (de/2-dt/2)/tand(alfa); % Nozzle length [mm]
L_d = (dcc/2-dt/2)/tand(beta); % Nozzle length [mm]
L = L_c+L_d;
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
T = ct.*At.*pcc.*1e5;
T_real = Isp.*mp_real*g0; % difference due to C_d of nozzle in mp
t = m./mp;
vexit = sqrt(2*k./(k-1)*1.e3*R./Mmol.*Tcc.*(1-(PeP0.^((k-1)./k))));
Mexit = sqrt(2*k./(k-1).*(1./PeP0).^((k-1)./k));
T2 = vexit.*mp_real+Ae.*Pe*1.e5; % to confirm
Itot = T_real.*t;

