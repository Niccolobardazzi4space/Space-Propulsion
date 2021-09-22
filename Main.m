%% Cubesat primary propulsion with water propellant
% 
%  CREATED BY:             Person Code:           Matriculation Number:
%    Niccolò Bardazzi          10800456                         963039
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  SPACE PROPULSION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

LOCALPATH=pwd;
PARENTPATH=fileparts(LOCALPATH);
COMMONPATH='\CODE';
PATH=[PARENTPATH,COMMONPATH];
addpath(genpath(PATH));

%% data
run('CEA_interpolation')

p = linspace(1,100,1000);
g0 = 9.81;
eps = ones(1,length(p))*250;

de_max = ones(1,length(p))*3e-2; % m
m = 0.5; % kg of water

% from graph shown
Power = 2.45e2*(p/10).^1.48e-2; % for 1 mole/s kW*(s/mol)
Power_supplied = 1; % 1kW
moles = Power_supplied./Power;

% with less pressure, # moles increases

M = 18; % g/mol, molar mass H2O
mp = M*moles*1.e-3; % mass flow rate, kg/s

pcc = p*1; % losses due to injection plate and check valves

%% Tanks
tanks = tank(p); % giving water storage pressure

%% Injection
R = 8.31; % J/molK

% Gases mass flow rate
Mmol_H2 = 2.01588;
Mmol_O2 = 31.9988;
mH2 = Mmol_H2/(Mmol_H2+Mmol_O2)*mp;
mO2 = Mmol_O2/(Mmol_H2+Mmol_O2)*mp;
OF = Mmol_O2/(2*Mmol_H2);
cp_O2 = 0.918;
cp_H2 = 14.31;
% Specific heat ratio 
k_O2 = cp_O2/(cp_O2-R/Mmol_O2); 
k_H2 = cp_H2/(cp_H2-R/Mmol_H2);

% Inlet
Temp = 298;
vH2 = ones(1,length(p))*15;
vO2 = ones(1,length(p))*15;
roH2 = p*1.e5/(R*1.e3/Mmol_H2*Temp);
roO2 = p*1.e5/(R*1.e3/Mmol_O2*Temp); % kg/m^3
aH2 = sqrt(k_H2*R*1e3./Mmol_H2.*Temp);
aO2 = sqrt(k_O2*R*1e3./Mmol_O2.*Temp);
MH2 = vH2./aH2;
MO2 = vO2./aO2;
AinjH2 = mH2./(roH2.*vH2);   % m^2
dinjH2 = (4*AinjH2/pi).^0.5; % m
AinjO2 = mO2./(roO2.*vO2);   % m^2
dinjO2 = (4*AinjO2/pi).^0.5; % m

% injector design
Cd_O2 = ones(1,length(p))*0.85; % from book (Vince)
DP_injO2 = 8./roO2.*(mO2./(Cd_O2.*AinjO2)).^2; % Injection plate pressure loss
DP_injH2 = DP_injO2; % Both flows require the same losses
Cd_H2 = mH2./(AinjH2).*sqrt(8./(roH2.*DP_injH2));
N_oxfu = 1; % configuration decided with swirling

%% VALVES
%%ON_OFF VALVE
Deltap_H2O = 0.01; %bar   COming from interpolation from datasheet
%%CHECK VALVE
SG_H2 = roH2/997;  %Ratio btw Density H2 and density Ration [-]  
SG_O2 = roO2/997;  %Ratio btw Density O2 and density Ration [-]
lm2gpm = 0.26417287472922;  % [l/min]->[gallons/min] Conversion
Q_H2 = mH2./roH2*1000*60*lm2gpm; % H2 mass flow rate [gallons/min]
Q_O2 = mO2./roO2*1000*60*lm2gpm; % O2 mass flow rate [gallons/min]
CV = 0.19;% Fluid Coefficient [gallons/min]
Deltap_H2 = ((SG_H2.*Q_H2.^2)./CV.^2); %Pressure drop [psi]
Deltap_O2 = ((SG_O2.*Q_O2.^2)./CV.^2);
Deltap_H2 = Deltap_H2*0.0689;  %Pressure conversion [psi]->[bar]
Deltap_O2 = Deltap_O2*0.0689;

%% CEA 
% H2, O2, pcc

% conversion
% DP_injO2 = DP_injO2*1.e-5;
% Dp = Deltap_O2+DP_injO2;
% p = p-Dp;

Tcc = ppval(f_T, p);
Mmol = ppval(f_M, p);
k = ppval(f_g, p);
vt = ppval(f_Sv, p); % [m/s]
ro_t = ppval(f_rho,x);

mu = [0.97057  1.04  1.0613  1.0738  1.0827  1.0895  1.0951  1.0998  1.1039  1.1075  1.1107]*1.e-4; % dependent only with T --> constant
f_mu = spline([1 10 20 30 40 50 60 70 80 90 100],mu);
mu = ppval(f_mu, p);

fi = k.*(2./(k+1)).^((k+1)./(k-1)); %Is this correct? The first k.^2
cstar = sqrt(R*1e3*Tcc./Mmol./fi);
At = cstar.*mp./(pcc*1e5); % m^2
Ae = eps.*At; % m^2
de = sqrt(Ae*4/pi)*1e3; % mm
dt = sqrt(At*4/pi)*1e3; % mm
Ae_max = de_max(1)^2*pi/4;

%% CC design
L_star = 0.635; % overestimated, m
Mx = 0.01;
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
PeP0_guess = ones(1,length(p))*1.e-6;
PeP0 = zeros(1,length(p));
options = optimset('Display','off');
for i = 1:length(p)
    PeP0(i) = fsolve(@(x)sqrt(k(i)*(2/(k(i)+1))^((k(i)+1)/(k(i)-1)))/((x)^(1/k(i))*sqrt((2*k(i))/(k(i)-1)*(1-x^((k(i)-1)/k(i)))))-eps(i),PeP0_guess(i),options);
end
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

%% Plot
figure(1)
subplot(1,3,1)
plot(p,T_real,'k')
title('Thrust')
xlabel('Pressure [bar]')
ylabel('T [N]')
hold on 
plot(p,T,'--g')
hold off
axis tight
legend('real','ideal')

subplot(1,3,2)
plot(p,Itot,'k')
axis tight
title('Total impulse')
xlabel('Pressure [bar]')
ylabel('I_{tot} [N\cdot s]')

subplot(1,3,3)
plot(p,Isp,'k')
axis tight
title('Specific impulse')
xlabel('Pressure [bar]')
ylabel('I_{sp} [s]')

figure(2)
subplot(2,1,1)
plot(p,C_d)
axis tight
title('Discharge coefficient')
xlabel('Pressure [bar]')
ylabel('C_d')
subplot(2,1,2)
Lambda = ones(1,length(p))*lambda;
plot(p,Lambda)
axis tight
title('Divergence losses')
ylim([0.97 1])
xlabel('Pressure [bar]')
ylabel('\lambda')

figure(3)
subplot(2,3,1)
plot(p,Tcc)
axis tight
title('Combustion chamber temperature')
xlabel('Pressure [bar]')
ylabel('Tcc [K]')

subplot(2,3,2)
plot(p,Mmol,'k')
axis tight
title('Molar mass')
xlabel('Pressure [bar]')
ylabel('M_{mol} [g/mol]')

subplot(2,3,3)
plot(p,k,'g')
axis tight
title('Specific heat ratio')
xlabel('Pressure [bar]')
ylabel('k')

subplot(2,3,4)
plot(p,vt,'r')
axis tight
title('Throat velocity')
xlabel('Pressure [bar]')
ylabel('v_t [m/s]')

subplot(2,3,5)
plot(p,ro_t,'c')
axis tight
title('Throat density')
xlabel('Pressure [bar]')
ylabel('\rho_t [kg/m^3]')

subplot(2,3,6)
plot(p,mu,'m')
axis tight
title('Dynamic viscosity')
xlabel('Pressure [bar]')
ylabel('\mu [Pa\cdot s]')


figure(4)
subplot(1,2,1)
plot(p,ct)
axis tight
title('Thrust coefficient')
xlabel('Pressure [bar]')
ylabel('c_T')
subplot(1,2,2)
plot(p,cstar)
axis tight
title('Characteristic velocity')
xlabel('Pressure [bar]')
ylabel('c^* [m/s]')


% %% Hypothesis
% Isp = 325;
% T = 0.3;
% mp = T/Isp/g0;
% At = 3e-2/150;


%% Pareto
% Parameters for Pareto
param.eps = eps(1);
f_mpreal = spline(p,mp_real); param.f_mpreal = f_mpreal;
f_mp = spline(p,mp);          param.f_mp = f_mp;
param.f_T = f_T;
param.f_M = f_M;
param.f_g = f_g;
param.m = m;

% Optimization parameters
A = [];  b = []; Aeq = []; beq = [];
lb = 1;
ub = 100;
nonlcon = [];
initialpoints = 10; % Initial point for pareto front

ParetoObjective = @(x)pareto_water(x,param);

options = optimoptions(@paretosearch,'UseParallel',true,'Display','iter','InitialPoints',initialpoints,...
                        'PlotFcn',{@psplotparetof},...
                        'ParetoSetSize',1000);
                   
[x,Fval,exitFlag,Output] = paretosearch(ParetoObjective,length(lb),A,b,Aeq,beq,lb,ub,nonlcon,options);       

xlabel('$D_t [mm]$','Interpreter','latex','FontSize',14)
ylabel('$I_{tot} [N\cdot s]$ ','Interpreter','latex','FontSize',14)
zlabel('$T_{cc} [K]$','Interpreter','latex','FontSize',14)
grid on
title('Pareto Front')

p = polyfit(Fval(:,1),Fval(:,2),2);
x1 = linspace(Fval(1,1),Fval(end,1),100);
y1 = polyval(p,x1);
hold on
plot(x1,y1,'--b')



