function [tanks] = tank(p)
% This function sizes spherical tanks of fuel and of pressurizant gas
% input:  - pressure of fuel tank [bar]
%         
%
% output: - thickness [mm]
%         - diameter [mm]
%         - radius [mm]
%         - Surface [mm^2]
%         - Mass [g]


% Water
m = 0.5; % kg
roH2O = 997; % water density, kg/m^3
T = 273; % freezing point, K 

%Helium
k=1.67; %specific heat ratio of helium
T1=310; %K, initial He temperature
T2=298; %K, final He temperature
pi_h=p*(T2/T1)^(-k/(k-1)); %initial pressure into helium tank, bar

% Tank material
%Carbon fiber
%F = 0.895e9; % tensile strength, Pa
%roCF = 1550; % carbon fiber density, kg/m^3
%Steel
F= 0.862*1.e9;
roST=7830;
%roT=roCF;
roT=roST;


% Thickness,diameter and volume of tanks
V = m/roH2O; % m^3
V = V*1.02; %volume increasing for safety
V_helium_f=V/(1-(p/pi_h)*(T1/T2)); %final volume of helium in the tank
V_helium_i=V_helium_f-V; %Volume of helium tanks
tanks.water.r = ((3/4)*V/pi)^(1/3); %radius of water tank m
tanks.helium.r= ((3/4)*V_helium_i/pi)^(1/3); % radius of helium tank m
tanks.water.d = 2*tanks.water.r; %diameter of water tank m
tanks.helium.d = 2*tanks.helium.r; %diameter of helium  m
pb_w = 2*p*1.e5; %pressure at which water tank explodes (burst pressure) Pa
pb_h= 2*pi_h*1.e5; %pressure at which helium tank explodes (burst pressure) Pa
tanks.water.t = pb_w*tanks.water.r/(2*F); %thickness of water tank m
tanks.helium.t= pb_h*tanks.helium.r/(2*F); %thickness of helium tank m
tanks.helium.p_initial = pi_h;
tanks.helium.p_final = p;
tanks.water.p = p;


% Mass and surface of tanks
tanks.water.S = 4*pi*tanks.water.r^2; % surface of water tank
tanks.helium.S = 4*pi*tanks.helium.r^2; % surface of helium tank
tanks.water.M = tanks.water.S*tanks.water.t*roT; %mass of water tank
tanks.helium.M= tanks.helium.S*tanks.water.t*roT; %mass of helium tank

%Results
tanks.water.r=tanks.water.r*1e3; %mm
tanks.helium.r=tanks.helium.r*1e3; %mm
tanks.water.d=tanks.water.d*1e3; %mm
tanks.helium.d=tanks.helium.d*1e3; %mm
tanks.water.S=tanks.water.S*1e6; %mm^2
tanks.helium.S=tanks.helium.S*1e6; %mm^2
tanks.water.t=tanks.water.t*1e3; %mm
tanks.helium.t=tanks.helium.t*1e3; %mm
tanks.water.M=tanks.water.M*1e3; %g
tanks.helium.M=tanks.helium.M*1e3; %g

end
    