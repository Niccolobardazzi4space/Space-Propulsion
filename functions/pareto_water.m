function y = pareto_water(x,param)

p = x; % t_flyby = x(2); t_arr = x(3);

% constants
R = 8.31;
g0 = 9.81;

eps = param.eps;
m = param.m;
f_T = param.f_T;
f_M = param.f_M;
f_g = param.f_g;
f_mpreal = param.f_mpreal;
f_mp = param.f_mp;

pcc = p; % losses due to injection plate and check valves

Tcc = ppval(f_T, pcc);
Mmol = ppval(f_M, pcc);
k = ppval(f_g, pcc);
mp = ppval(f_mp, pcc);
mp_real = ppval(f_mpreal, pcc);

fi = k*(2/(k+1))^((k+1)/(k-1));
cstar = sqrt(R*1e3*Tcc/Mmol/fi);
At = cstar*mp/(pcc*1e5);

PeP0_guess = 1.e-7;
options = optimset('Display','off');
PeP0 = fsolve(@(x)sqrt(k*(2/(k+1))^((k+1)/(k-1)))/((x)^(1/k)*sqrt((2*k)/(k-1)*(1-x^((k-1)/k))))-eps,PeP0_guess,options);

Pe = PeP0*pcc;
ct = (0.983*sqrt(fi.*(2.*k./(k-1)).*(1-(PeP0).^((k-1)./k)))+(Pe/pcc).*eps);
Isp = cstar.*ct/g0;
T_real = Isp.*mp_real*g0;
t = m/mp;
Itot = T_real*t;

dt = sqrt(At*4/pi)*1.e3;

y(1) = dt; y(2) = Itot; y(3) = Tcc;
end