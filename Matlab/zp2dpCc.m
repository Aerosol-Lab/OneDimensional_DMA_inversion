%% Mobility Diameter from Electrical Mobility Calculation using Cunningham Slip Correction Factor (iteration)
function [dpList, dp_FM] = zp2dpCc(ZList)
dpList = zeros(length(ZList),1);
for dpi = 1:length(ZList)
    Zp = ZList(dpi);
n = 1; %assume singly charged
e = 1.6022E-19; %electron charge
mu = 1.82E-5; %dynamic viscosity
mair = 28.97*1.66053892E-27;
kb = 1.3806488E-23;
rho = 1.24;
eps = 1.36;
T = 298;

if (1/Zp)<5E6
    coeff = (mair*pi/(8*kb*T))^0.5 * 3*e/(pi*rho*eps);
    dp0 = (coeff/Zp)^0.5 - 0.3E-9; % initial guess for FM regime
    dp_FM = dp0;
else
    dp0 = n*e/(3*pi*mu*Zp) - 0.3E-9; %initial guess fo continuum regime at mobility diameter, dp.
end

% Cunningham Slip Correction Factor (currently guess at formula)
a1 = 1.257;
a2 = 0.4;
a3 = 0.55;
i = 1;
Kn = 2*66.5E-9/(dp0 + 0.3E-9); %Knudsen number
cFactor = 1 + Kn*(a1 + a2*exp(-2*a3/Kn)); %Slip correction factor
dp = n*e*cFactor/(3*pi*mu*Zp) - 0.3E-9;
dpMAT(i) = dp0;
dpMAT(i+1) = dp;

while abs(dp0-dp)/dp0 > 0.0001 && i < 1000
    i = i + 1;
    dp0 = (dp + dp0)/2;
    Kn = 2*66.5E-9/(dp0+0.3E-9); %Knudsen number
    cFactor = 1 + Kn*(a1 + a2*exp(-2*a3/Kn));
    dp = n*e*cFactor/(3*pi*mu*Zp) - 0.3E-9;
    dpMAT(i) = dp;
end
dpList(dpi) = dp;
end
