function CC = CC_design(At,L_star,gamma,Mx)

%function CC = CC_design(At,L_star,eps_C,mp,T,MM)
%CC_DESIGN Preliminary design of the combustion chamber
% Creator: Pietro Gazzi (pietrogazzi01@gmail.com)
%   input : 
%           - At       Throat area               [m^2]
%           - L_star   CC characteristic length  [m]
%           - eps_C    Chamber construction      [Ac/At]
%           - gamma    Heat ratio
%           - Mx       Interface mach number
%           [- mp       Mass flow rate            [kg/s]]
%           [- T        Flame Temperature         [K]]
%           [- MM       Molar mass                [kg/kmol]]



My = 1;     %Throat Mach Number 
eps_C = Mx/My*sqrt(((1+((gamma-1)/2)*My^2)/(1+((gamma-1)/2)*Mx^2))^((gamma+1)/(gamma-1)));
%R = 8.32*10^3;
Ac = (1/eps_C)*At;         %Chamber Area [m^2]
Vc = L_star*At;        %CC volume [m^3]
%rho = P/(R/MM*T);     %fluid density
L_c = Vc./Ac;  %CC length [m]



%%Output Struct
CC.Vc = Vc;
CC.L_c = L_c*10^3;
CC.Ac = Ac;
CC.eps_C=1/eps_C;
CC.dc = sqrt(4*Ac/pi)*10^3;

end

