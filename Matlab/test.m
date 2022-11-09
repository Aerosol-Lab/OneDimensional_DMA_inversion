%% Inversion DEMO
close all
clear all
%Settings

%Experiment space
numVoltages = 31; % number of experimental conditions
voltMin = 10;
voltMax = 10e3;
voltages = exp(linspace(log(voltMin),log(voltMax),numVoltages));

Qaerosol = 0.28; %flow in LPM
Qshealth = 6.0; %flow in LPM
L = 0.44369; % TSI spec
R1 = 0.00937; % from TSI spec sheet: http://cires1.colorado.edu/jimenez-group/Manuals/SMPS_manual
R2 = 0.01961; % 
Qaerosol = Qaerosol * 1.6666e-5;
Qshealth = Qshealth * 1.6666e-5;

ZpStarDebug = (Qshealth./(2*pi*L*voltages)).*log(R2/R1); % Mobility calculation based on DMA V
dpExperiment = zp2dpCc(ZpStarDebug);
dpMin2 = min(dpExperiment);%1e-9;
dpMax2 = max(dpExperiment);%1e-6;

%phantom settings
dpMin = 1e-9;
dpMax = 1e-6;
numParticles =1000;
Ntotal = 1;%1e2*(1e6); %concentration
sigma = 1.3;
dpMedian = 200e-9;

%% create phantom function 
dpBins = exp(linspace(log(dpMin),log(dpMax),numParticles+1));
dpPhantom = exp(0.5*log(dpBins(2:end).*dpBins(1:end-1)));
n = (Ntotal./(log(sigma)*sqrt(2*pi))).*exp(-(log(dpPhantom/dpMedian).^2)/(2.0*log(sigma)^2)); %distribution function
n2 = (Ntotal./(log(1.1)*sqrt(2*pi))).*exp(-(log(dpPhantom/(100e-9)).^2)/(2.0*log(1.1)^2)); %distribution function
n3 = (Ntotal./(log(1.4)*sqrt(2*pi))).*exp(-(log(dpPhantom/(50e-9)).^2)/(2.0*log(1.4)^2)); %distribution function
n = n+ 0.2*n2 + 0.1*n3;
figure
hold on
plot(dpPhantom*1e9,n*1e-6)
hold off
set(gca,'XScale','log')
xlabel('Particle Diameter [nm]')
ylabel('Distriubtion Function "dn/dlogdp" [cm^-3]')

%% N_i = int (lambda_i*n) -inf to inf
%lambda is a continous function for some discreate condition i

lambdaList =zeros(numVoltages,numParticles);
for i = 1:numVoltages
lambdaList(i,:) = getLambda(dpPhantom,voltages(i));
end

%calculate "experimental" result
N = zeros(1,numVoltages);

for i = 1:numVoltages
N(i)  = trapz(dpPhantom,lambdaList(i,:).*n./dpPhantom); %Add 1/dpPhantom to change integration from log dp to dp. Could also change dpPhantom to log(dpPhantom)
end

figure
hold on
plot(dpPhantom*1e9,n*1e-6)
hold off
ylabel('Distriubtion Function "dn/dlogdp" [cm^-3]')
yyaxis right



plot(dpExperiment*1e9,N,'sq','MarkerEdgeColor','k','MarkerFaceColor','r')
ylabel('Concentration [cm^-3]')
yRange = ylim;
yRange = yRange -min(yRange);
ylim(yRange)
hold off
set(gca,'XScale','log')
xlabel('Particle Diameter [nm]')

%% On to inversion
%Now we have "experimental data" N(i) we can preform inversion and see if
%we get that same result as our original distribution

%% make A
%dpResult = dpPhantom;

numParticles = numVoltages*10;
%dpBins = [dpExperiment; dpExperiment(end)*1.5]; %exp(linspace(log(dpMin),log(dpMax*10),numParticles+1));
dpBins = exp(linspace(log(dpMin2/1.2),log(dpMax2*1.2),numParticles+1));
dpResult =exp(0.5*log(dpBins(2:end).*dpBins(1:end-1))); %dpPhantom;

A = zeros(numVoltages,numParticles);
b = N';

L2d = eye(numParticles,numParticles)*2;
L2d(1,1) = 1;
for i = 2:numParticles-1
L2d(i,i-1) = -1;
L2d(i,i+1) = -1;
end
L2d(numParticles,numParticles) = 1;

LtL2d = transpose(L2d)*L2d;

% Get A
for i = 1:numVoltages
    for j = 1:numParticles
        logdpMinij = log(dpBins(j));
        logdpMaxij = log(dpBins(j+1));
        resolution = 500;
        xIntegrate = linspace(logdpMinij,logdpMaxij,resolution);
        lambdaIntegrate =  getLambda(exp(xIntegrate),voltages(i)); 
        integrand = trapz(xIntegrate,lambdaIntegrate);
           A(i,j) = integrand;
    end
end

%LtL2d = eye(numParticles,numParticles);
%% Inversion 
%x = solveInversion(A,b,lambda,LtL2d);

%%


updateFigure(dpPhantom,dpResult,n,dpExperiment,N,A,b,LtL2d,lambdaList,numVoltages);


function updateFigure(dpPhantom,dpResult,n,dpExperiment,N,A,b,LtL2d,lambdaList,numVoltages)
lambda = 1e-2;
[x,r,r2,r3] = solveInversion(A,b,lambda,LtL2d);
Nrecalculated = reCalcN(dpResult,dpPhantom,x,lambdaList,numVoltages);

fig = figure;
s11=subplot(1,2,1);
left_color = [0 0 0];
right_color = [0 0 0];
set(s11,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis right
hold on
box on
plot(dpExperiment*1e9,N,'sq','MarkerEdgeColor','k','MarkerFaceColor','r')
p2 = plot(dpExperiment*1e9,Nrecalculated,'^','MarkerEdgeColor','k','MarkerFaceColor','m')
hold off
ylabel('Concentration [cm^-3]')
yRange = ylim;
yRange = yRange -min(yRange);
ylim(yRange)
hold off
set(gca,'XScale','log')
xlabel('Particle Diameter [nm]')
yyaxis left
hold on

plot(dpPhantom*1e9,n*1e-6)
xlocal = x;
p = plot(dpResult*1e9,xlocal*1e-6,'--k')
yRange = ylim;
yRange = max(yRange,0);
ylim(yRange);
hold off
ylabel('Distriubtion Function "dn/dlogdp" [cm^-3]')
%p.XDataSource = 'xlocal';

legend('Phantom Dist','Invetered Dist','Phantom Data','Inverted Data')
s = subplot(1,2,2);
box on
hold on
xlabel('Residual (Ax-b)')
ylabel('Solution Residual (Lx)')
%Lcurve generate once
[Residual, RegPerameter] = genLCurve(A,b,LtL2d);
p3 = plot(Residual,RegPerameter);
set(s,'XScale','log')
set(s,'YScale','log')
p4 = plot(r,r3,'sq','MarkerEdgeColor','k','MarkerFaceColor','r');
t1 = text(r,r3,["\lambda " num2str(lambda,'%.2e')])
hold off

KEY_IS_PRESSED = 0;
contLoop = 0;

set(gcf, 'KeyPressFcn', @myKeyPressFcn)
Key ='';
log10Lambda = log10(lambda);
deltalog10Lambda = 0.1;
drawnow

disp("Use left and right arrows to increase or decrease regularlization parameter \lambda. q to quit")
while or(~KEY_IS_PRESSED,contLoop)
      drawnow
      %disp(["Value ", num2str(10^log10Lambda)])

      %using strcmp for string comparison if comparison is true = 1
      switchString = num2str(Key);
      switch switchString
          case 'leftarrow'
              KEY_IS_PRESSED = 0;
              contLoop = 1;
              log10Lambda = max(log10Lambda - deltalog10Lambda,-10);
              Key ='';
              
          case 'rightarrow'
              KEY_IS_PRESSED = 0;
              contLoop = 1;
              log10Lambda = min(log10Lambda + deltalog10Lambda,1);
              Key ='';
              
          case 'q'
              KEY_IS_PRESSED = 1;
              contLoop = 0;
              Key ='';
              disp('Quit ME!')
          otherwise
              KEY_IS_PRESSED = 0;
              contLoop = 0;
              Key ='';
      end
       
      if(contLoop == 1)
          [xlocal,r,r2,r3] = solveInversion(A,b,10^(log10Lambda),LtL2d);
          set(p,'ydata',xlocal*1e-6);
          Nrec =  reCalcN(dpResult,dpPhantom,xlocal,lambdaList,numVoltages);
          set(p2,'ydata',Nrec);
          set(p4,'xdata',r,'ydata',r3);
          rangeX = log(get(s, 'XLim'));
          rangeY = log(get(s, 'YLim'));
          factorX = 0.1*abs(rangeX(2)-rangeX(1));
          factorY = 0.1*abs(rangeY(2)-rangeY(1));
          set(t1,'Position',[exp(log(r)+factorX) exp(log(r3)+factorY)],'String',["\lambda " num2str(10^(log10Lambda),'%.2e')]);
          refreshdata
          drawnow
      end
end

disp('loop ended')

 function myKeyPressFcn(hObject, event)
      Key=get(hObject,'CurrentKey');
 end

end

function [x,r,r2,r3] = solveInversion(A,b,lambda,L)
    Aprime = transpose(A)*A + lambda^2*L;
    bprime = transpose(A)*b;
    
    x = linsolve(Aprime,bprime);
    r = (b-A*x);
    r = sqrt(sum(r.^2)/length(r));
    r2 = (bprime-Aprime*x);
    r2 = sqrt(sum(r2.^2)/length(r2));
    r3 = (L*x);
    r3 = sqrt(sum(r3.^2)/length(r3));

end

function  Nrecalculated = reCalcN(dpResult,dpPhantom,x,lambdaList,numVoltages)

xinterp = interp1(dpResult,x,dpPhantom,'linear','extrap');
Nrecalculated = zeros(1,numVoltages);
for i = 1:numVoltages
Nrecalculated(i)  = trapz(dpPhantom,lambdaList(i,:).*xinterp./dpPhantom); %Add 1/dpPhantom to change integration from log dp to dp. Could also change dpPhantom to log(dpPhantom)
end
end

function  [Residual, RegPerameter] = genLCurve(A,b,L)

res = 500; 
Lmin = 1e-10; 
Lmax = 10;
Llist = exp(linspace(log(Lmin),log(Lmax),res));

Residual = zeros(length(Llist),1);
RegPerameter = zeros(1,length(Llist));
for i = 1: length(Llist)
[x,r,r2,r3] = solveInversion(A,b,Llist(i),L);
Residual(i) = r;
RegPerameter(i) = r3;
end
end