threshCont = .07303;

minCont = .015;
maxCont = .80;


cont2 = 10^(log10(minCont*100)+.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont3 = 10^(log10(threshCont*100)-.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont4 = 10^(log10(threshCont*100)-(1/8)*(log10(threshCont*100) - log10(minCont*100)))/100;
cont6 = 10^(log10(maxCont*100)-.25*(log10(maxCont*100) - log10(threshCont*100)))/100;


contLevels = [minCont,cont2,cont3,cont4threshCont,cont6,maxCont];

figure,clf, hold on
plot(log10(contLevels.*100),[.5,.5,.66,.75,.84,1,1],'bo','LineWidth',2)
axis([minCont,2,.4,1.1])
hold off