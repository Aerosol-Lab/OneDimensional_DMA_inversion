function [DMATfer] = KernelCalc_InvMobility_dp(invZpArray, Qsh, Qa,V,T,L,R1,R2)

%% NDMA Diffusive Transfer Function and Initial Parameters ****
kb = 1.38064852e-23;
e = 1.6022E-19;
n = 1; % assuming singly charged.

ZpStar = (Qsh/(2*pi*L*V))*log(R2/R1); % Mobility calculation based on DMA V

%%%%% Transfer Function Calculation:

y = (R1/R2)^2;
kap = L*R2/(R2^2 - R1^2);
Iy = (0.25*(1-y^2)*(1-y)^2 + (5/18)*(1 - y^3)*(1 - y)*log(y) + (1/12)*(1-y^4)*(log(y))^2)/((1-y)*(-0.5*(1+y)*log(y) - (1-y))^2);
B = (Qa/Qsh);
delta = 0; % delta is 0 if the aerosol and sampling flowrates are equal, which is always the case in this experiment.
Z_p = (invZpArray.*ZpStar).^(-1); % array of non-dimensional mobilities
GDMA = 4*(1 + B)^2*(Iy + (2*(1+B)*kap)^(-2))/(1-y);
sig = (GDMA*log(R2/R1)*kb*T/(n*e*V))^0.5*(Z_p).^0.5;
sigStar = sig(5)/(Z_p(5)^0.5);

DMATfer = (2^0.5*B*(1-delta))^(-1)*sig.*(eps1((Z_p - (1 + B))./(2^0.5*sig)) + eps1((Z_p - (1 - B))./(2^0.5*sig)) - 2.*eps1((Z_p - 1)./(2^0.5*sig)));
% The line above computes the diffusive transfer function for the DMA
% based on Stolzenberg's 1988 analysis for previously defined Zp range

% ZpPlot = transpose(Zp);
% dpPlot = transpose(dp);
% sigPlot = transpose(sig);
% Z_pPlot = transpose(Z_p);
% DMATferPlot = transpose(DMATfer);
% sigStar = sig(5)/(Z_p(5)^0.5);

%% Kernel Calculation

G = DMATfer;
G = transpose(G);

negTransform = G > 0;
G = negTransform.*G; % these two lines eliminate negative values
smallTransform = G > 1E-14;
G = smallTransform.*G;