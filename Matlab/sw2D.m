DMA1V = y;
DMA2V = x;

Qaerosol = 0.28; %flow in LPM
Qshealth = 3.0; %flow in LPM

P = 101325; % Operating Pressure
T = 300; % Operating Temperature
%Assuming air is the gas used.


L = 0.44369; % TSI spec
R1 = 0.00937; % from TSI spec sheet: http://cires1.colorado.edu/jimenez-group/Manuals/SMPS_manual
R2 = 0.01961; % 
Area = 3.1415*(R2^2-R1^2);

%% unit conversion to standard SI units
Qaerosol = Qaerosol * 1.6666e-5;
Qshealth = Qshealth * 1.6666e-5;

zpDMA1 = Qshealth*log(R2/R1)./(2.0*3.1415*L*DMA1V);
zpDMA2 = Qshealth*log(R2/R1)./(2.0*3.1415*L*DMA2V);

dp2 = []
for i = 1:length(zpDMA2)
    dp2(i) = zp2dpCc(zpDMA2(i));
end

mfp = 67e-9;
c = [];
for i = 1:length(dp2)
c(i) = Cc(dp2(i),mfp);
end
e = 1.609e-19;
z = [];
for i = 1:length(zpDMA1)
    z(i,:) = zpDMA1(i).*3.0*3.1415*1.8e-5.*dp2./c./e;
end

[X,Y] = meshgrid(dp2*1e9,log(zpDMA1*10000))

contourf(X,Y,Z,100,'edgecolor','none')
xlabel('DMA 2 (Size) [nm]')
ylabel('DMA 1 (Mobility) [cm V^-1 s^-1]')
colorbar

