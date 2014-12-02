

PF = @PAL_Weibull;

alpha = 0.12;
beta = 6;
gamma = 0.5;
lambda = 0.2;

c50_y = (((1-lambda)-gamma)/2)+gamma;

params = [alpha,beta,gamma,lambda];

xCurve = logspace(.001,2,10000);
xCurveLOG = log10(xCurve);
xCurve = xCurve./100;

yCurve = PF(params,xCurve);

c50_x = PF(params,c50_y,'Inverse');

noiseSigStair = .02;
noiseSigConstStim = .03;


minCont = 0.015;
maxCont = 0.8;
threshCont = c50_x+normrnd(0.02,noiseSigStair);

cont2 = 10^(log10(minCont*100)+.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont3 = 10^(log10(threshCont*100)-.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont5 = 10^(log10(threshCont*100)+(1/6)*(log10(maxCont*100) - log10(threshCont*100)))/100;
cont6 = 10^(log10(maxCont*100)-.25*(log10(maxCont*100) - log10(threshCont*100)))/100;


contLevels = [minCont,cont2,cont3,threshCont,cont5,cont6,maxCont];

behavior = PF(params,contLevels) + normrnd(0,noiseSigConstStim,size(contLevels,1),size(contLevels,2));


OutOfNum = 56.*ones(size(behavior));
NumPos = round(behavior.*OutOfNum);
paramsFree = [1,1,0,1];
lapseLimits = [.1,.3];
fit_params = PAL_PFML_Fit(contLevels,NumPos,OutOfNum,params,paramsFree, PF, 'lapseLimits', lapseLimits);


% plots

axisVect = [min(xCurveLOG),max(xCurveLOG),.2,1.01];

figure(1),clf,hold on
plot(log10(contLevels.*100),behavior,'bo','LineWidth',2)

plot(xCurveLOG,PF(params,xCurve),'b--','LineWidth',2)
plot(xCurveLOG,PF(fit_params,xCurve),'r--','LineWidth',2)


plot([min(xCurveLOG),log10(100*c50_x)],[c50_y,c50_y],'k-.');
plot([log10(100*c50_x), log10(100*c50_x)], [axisVect(3),c50_y],'k-.');

title('Simulated curve & data')
ylabel('proportion correct')
xlabel('contrast')
set(gca,'XTick',log10(contLevels*100),'XTickLabel',contLevels*100,'YTick',.4:.1:1)
axis(axisVect)
hold off


