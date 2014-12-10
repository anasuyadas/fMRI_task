


numBlocks = 6;
participant = 'MO_day2';

fileNums = [1:6]; %which file numbers area the constant stimulus blocks?

numContrasts = 7;

numCorr.loc = cell(4,1);
numResponses.loc = cell(4,1);
for i = 1:4
        numCorr.loc{i} = zeros(1,numContrasts);
        numResponses.loc{i} = zeros(1,numContrasts);
end

numCorr.diagonal = cell(2,1);
numResponses.diagonal = cell(2,1);
for i = 1:2
        numCorr.diagonal{i} =  zeros(1,numContrasts);
        numResponses.diagonal{i} =  zeros(1,numContrasts);
end

diagonal = repmat([2,2],[1,5]);


addpath(sprintf('data/%s',participant));
for block = 1:numBlocks
    currExpFile = dir(sprintf('data/%s/*stim0%d.mat',participant,fileNums(block)));
    load(sprintf('%s',currExpFile.name));
    
    exp = getTaskParameters(myscreen,task);
    
    excludeTrials = stimulus.fixationBreakTrialVect-1;
    
    diagThisBlock = diagonal(block);
    for trialNum = 1:exp.nTrials
        index = exp.randVars.trialIndex(trialNum);
        if ~isnan(exp.response(trialNum))
            if ~(sum(excludeTrials==trialNum))
                %by location
                numCorr.loc{stimulus.randVars.targetLocation(index)}(stimulus.randVars.contrast(index)) =...
                    numCorr.loc{stimulus.randVars.targetLocation(index)}(stimulus.randVars.contrast(index)) +...
                   (stimulus.randVars.targetOrientation(index)==exp.response(trialNum));
                
                numResponses.loc{stimulus.randVars.targetLocation(index)}(stimulus.randVars.contrast(index)) =...
                    numResponses.loc{stimulus.randVars.targetLocation(index)}(stimulus.randVars.contrast(index)) + 1;
            
                %by diagonal
                numCorr.diagonal{diagThisBlock}(stimulus.randVars.contrast(index)) = ...
                    numCorr.diagonal{diagThisBlock}(stimulus.randVars.contrast(index)) +...
                   (stimulus.randVars.targetOrientation(index)==exp.response(trialNum));
                
                numResponses.diagonal{diagThisBlock}(stimulus.randVars.contrast(index)) = ...
                    (numResponses.diagonal{diagThisBlock}(stimulus.randVars.contrast(index)) + 1);
            end
        end
    end
    
    fprintf('Completed block %d out of %d \n',block,numBlocks);
    
end

perfStruct{1} = struct('numCorr',numCorr,'numResponses',numResponses);

%% fit curves

palamedesPath = '/Users/purplab/Documents/MATLAB/palamedes1_6_0/Palamedes';

addpath(genpath(palamedesPath));

contLevels = stimulus.contrasts;

PF = @PAL_Weibull;

alphaguess = 0:.1:1;
betaguess = .5:1:30;
gammaguess = 0.5;
lambdaguess = 0:.01:.5;

paramsFree = [1 1 0 1];

lapseLimits = [0, .5];

searchGrid = struct('alpha', alphaguess,'beta', betaguess,'gamma',gammaguess, 'lambda',lambdaguess);


for locFit = 1:4
    NumPos = perfStruct{1}.numCorr.loc{locFit};
    OutOfNum = perfStruct{1}.numResponses.loc{locFit};
    
   
    [fitPercCorr.loc{locFit}.params,...
     fitPercCorr.loc{locFit}.LL,...
     fitPercCorr.loc{locFit}.exitflag]...
                              = PAL_PFML_Fit(contLevels,...
                                             NumPos,...
                                             OutOfNum,...
                                             searchGrid, paramsFree, PF, 'lapseLimits', lapseLimits);
                                         
    fitPercCorr.loc{locFit}.c50PERFORMANCE = (((1-fitPercCorr.loc{locFit}.params(4))-.5)/2)+.5;                                   
          
    fitPercCorr.loc{locFit}.c50 =...
        PF(fitPercCorr.loc{locFit}.params,...
           fitPercCorr.loc{locFit}.c50PERFORMANCE,'Inverse');
       
    fprintf('Completed FIT for location %d out of %d \n',locFit,4);
    fprintf('Success = %d \n \n',fitPercCorr.loc{locFit}.exitflag)
    
end


for diagonalFit = 1:2
    NumPos = perfStruct{1}.numCorr.diagonal{diagonalFit};
    OutOfNum = perfStruct{1}.numResponses.diagonal{diagonalFit};
    
    [fitPercCorr.diagonal{diagonalFit}.params,...
     fitPercCorr.diagonal{diagonalFit}.LL,...
     fitPercCorr.diagonal{diagonalFit}.exitflag]...
                              = PAL_PFML_Fit(contLevels,...
                                             NumPos,...
                                             OutOfNum,...
                                             searchGrid, paramsFree, PF, 'lapseLimits', lapseLimits);
                        
                                         
    fitPercCorr.diagonal{diagonalFit}.c50PERFORMANCE = (((1-fitPercCorr.diagonal{diagonalFit}.params(4))-.5)/2)+.5;                                   
                                         
    fitPercCorr.diagonal{diagonalFit}.c50 =...
        PF(fitPercCorr.diagonal{diagonalFit}.params,...
           fitPercCorr.diagonal{diagonalFit}.c50PERFORMANCE,'Inverse');
                                         
    fprintf('Completed FIT for diagonal %d out of %d \n',diagonalFit,2);
    fprintf('Success = %d \n \n',fitPercCorr.diagonal{diagonalFit}.exitflag)
end

perfStruct{1}.fitPercCorr = fitPercCorr;


%% plot
xCurve = logspace(.001,2,10000);
xCurveLOG = log10(xCurve);
xCurve = xCurve./100;

axisVect = [min(xCurveLOG),max(xCurveLOG),.2,1.01];

colorVect = [0,0,1;1,0,0;0,0,1;1,0,0];

for locPlot = 1:4

c50 = log10(perfStruct{1}.fitPercCorr.loc{locPlot}.c50*100);
c50PERF = perfStruct{1}.fitPercCorr.loc{locPlot}.c50PERFORMANCE;
    
figure(locPlot),clf,hold on
plot(log10(contLevels.*100),...
     perfStruct{1}.numCorr.loc{locPlot}./perfStruct{1}.numResponses.loc{locPlot},...
     'o','LineWidth',2,'Color',colorVect(locPlot,:))

plot(xCurveLOG,PF(perfStruct{1}.fitPercCorr.loc{locPlot}.params,xCurve),...
     '--','LineWidth',2,'Color',colorVect(locPlot,:))


plot([min(xCurveLOG),c50],[c50PERF,c50PERF],'k-.');
plot([c50, c50], [axisVect(3),c50PERF],'k-.');


title(sprintf('Loc %d',locPlot))
ylabel('proportion correct')
xlabel('contrast')
set(gca,'XTick',log10(contLevels*100),'XTickLabel',contLevels*100,'YTick',.4:.1:1)
axis(axisVect)
hold off

end





for diagonalPlot = 1:2

c50 = log10(perfStruct{1}.fitPercCorr.diagonal{diagonalPlot}.c50*100);
c50PERF = perfStruct{1}.fitPercCorr.diagonal{diagonalPlot}.c50PERFORMANCE;
    
    
    
figure(diagonalPlot+4),clf,hold on
plot(log10(contLevels.*100),...
     perfStruct{1}.numCorr.diagonal{diagonalPlot}./perfStruct{1}.numResponses.diagonal{diagonalPlot},...
     'o','LineWidth',2,'Color',colorVect(diagonalPlot,:))

plot(xCurveLOG,PF(perfStruct{1}.fitPercCorr.diagonal{diagonalPlot}.params,xCurve),...
     '--','LineWidth',2,'Color',colorVect(diagonalPlot,:))

plot([min(xCurveLOG),c50],[c50PERF,c50PERF],'k-.');
plot([c50, c50], [axisVect(3),c50PERF],'k-.');

title(sprintf('Diagonal %d',diagonalPlot))
ylabel('proportion correct')
xlabel('contrast')
set(gca,'XTick',log10(contLevels*100),'XTickLabel',contLevels*100,'YTick',.4:.1:1)
axis(axisVect)
hold off



end