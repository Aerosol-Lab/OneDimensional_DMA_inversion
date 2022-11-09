%% Cunningham Slip Correction Factor

function ans = Cc(dp,mfp)

A1 = 1.257;
A2 = 0.4;
A3 = 0.55;

Kn = 2*mfp./dp;
ans = 1 + Kn.*(A1 + A2*exp(-2*A3./Kn));