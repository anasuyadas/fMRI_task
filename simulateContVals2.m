

PF = @PAL_Weibull;

alpha = 0.2;
beta = 3;
gamma = 0.5;
lambda = 0.209;

c50_y = (((1-lambda)-gamma)/2)+gamma;

params = [alpha,beta,gamma,lambda];

xCurve = logspace(.001,2,10000);
xCurveLOG = log10(xCurve);
xCurve = xCurve./100;

yCurve = PF(params,xCurve);

c50_x = PF(params,c50_y,'Inverse');
thresh_x = PF(params,.709,'Inverse');

noiseSigStair = .02;
noiseSigConstStim = .0;


minCont = 0.015;
maxCont = 0.8;
threshCont = thresh_x+normrnd(0.02,noiseSigStair);

cont2 = 10^(log10(minCont*100)+.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont3 = 10^(log10(threshCont*100)-(1/4)*(log10(threshCont*100) - log10(minCont*100)))/100;
cont5 = 10^(log10(threshCont*100)-(1/8)*(log10(threshCont*100) - log10(minCont*100)))/100;
cont6 = 10^(log10(maxCont*100)-.25*(log10(maxCont*100) - log10(threshCont*100)))/100;


contLevels = [minCont,cont2,cont3,cont5,threshCont,cont6,maxCont];

probCorr = PF(params,contLevels);

numTrials = 56;
behavior = zeros(size(probCorr));
for trial = 1:numTrials
    behavior = behavior + (rand(size(probCorr))+normrnd(0,noiseSigConstStim,size(behavior,1),size(behavior,2)) < probCorr);
end
behavior

OutOfNum = numTrials.*ones(size(behavior));
NumPos = behavior;
paramsFree = [1,1,0,1];
lapseLimits = [0,.5];
fit_params = PAL_PFML_Fit(contLevels,NumPos,OutOfNum,params,paramsFree, PF, 'lapseLimits', lapseLimits);
c50FIT_y = (((1-fit_params(4))-fit_params(3))/2)+fit_params(3);
c50FIT_x = PF(fit_params,c50FIT_y,'Inverse');


% plots

axisVect = [min(xCurveLOG),max(xCurveLOG),.4,1];

figure(1),clf,hold on
plot(log10(contLevels.*100),NumPos./OutOfNum,'ko','LineWidth',2)

plot(xCurveLOG,PF(params,xCurve),'b--','LineWidth',2)
plot(xCurveLOG,PF(fit_params,xCurve),'r--','LineWidth',2)


plot([min(xCurveLOG),log10(100*c50_x)],[c50_y,c50_y],'k-.');
plot([log10(100*c50_x), log10(100*c50_x)], [axisVect(3),c50_y],'k-.');

plot([min(xCurveLOG),log10(100*c50FIT_x)],[c50FIT_y,c50FIT_y],'k:');
plot([log10(100*c50FIT_x), log10(100*c50FIT_x)], [axisVect(3),c50FIT_y],'k:');

legend('Data','Underlying Curve','Fitted Curve','c50 - Underlying', 'c50 - Fitted')
title('Simulated curve & data')
ylabel('proportion correct')
xlabel('contrast')
set(gca,'XTick',log10(contLevels*100),'XTickLabel',contLevels*100,'YTick',.4:.1:1)
axis(axisVect)
hold off


