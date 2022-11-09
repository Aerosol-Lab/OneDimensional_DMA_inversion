function l = getLambda(dp,V)

Qaerosol = 0.28; %flow in LPM
Qshealth = 6.0; %flow in LPM

P = 101325; % Operating Pressure
T = 300; % Operating Temperature
%Assuming air is the gas used.
e = 1.6022E-19;

L = 0.44369; % TSI spec
R1 = 0.00937; % from TSI spec sheet: http://cires1.colorado.edu/jimenez-group/Manuals/SMPS_manual
R2 = 0.01961; % 

%% unit conversion to standard SI units
Qaerosol = Qaerosol * 1.6666e-5;
Qshealth = Qshealth * 1.6666e-5;

%% Derived quantities
S = 110.4;
muref = 1.716e-5;
Tref = 273.15;
mu = muref*((T/Tref)^(3.0/2.0))*(Tref + S)/(T + S);
mfp0 = 67.3e-9;
T0 = 296.15;
p0 = 101325.0;
mfp = mfp0* (T / T0)* (p0 / P)* (1.0 + 110.4 / T0) / (1.0 + 110.4 / T);

invZpKernelInt = 3*pi*mu.*(dp + 0.3E-9)./(e*Cc(dp,mfp));  %Inverse mobility for singly charged particles
ZpStarDebug = (Qshealth./(2*pi*L*V)).*log(R2/R1); % Mobility calculation based on DMA V


%% Calculate Charge fraction
% Charged Fraction: calc charged fraction of particles for bipolar diffusion charging based on wiedensohler (for subr f_charge_W). Gopalakrishnan not used because fits don't go to 1 nm (only 10 nm)
numCharges = 3; % can be used to generalize charge distribution. 
[f0,f1,f2] = f_charge_W(dp);

fList = [f0;f1;f2];
%% Calculate tranfer function DMA
thetaDMA =0;

for(i = 2:numCharges)
    thetaDMAfi = fList(i,:).*KernelCalc_InvMobility_dp((1/(i-1)).*invZpKernelInt, Qshealth, Qaerosol,V,T,L,R1,R2); %calc Kernels KernelCalc_InvMobility_dp(invZpArray, Qsh, Qa,V,T,L,R1,R2) 1/(i-1) accounts for doubly charged mobility reduction
    thetaDMA = thetaDMA + thetaDMAfi;
end

l = thetaDMA;
end