numBlocks = 2;
participant = 'idTestALL_7';

fileNums = ['05';'06']; %which file numbers area the constant stimulus blocks?


addpath(sprintf('data/%s',participant));
for block = 1:numBlocks
    currExpFile = dir(sprintf('data/%s/*stim%s.mat',participant,fileNums(block,:)));
    load(sprintf('%s',currExpFile.name));
    
    
    numReversals = length(stimulus.stair.reversals);
    numTrials = length(stimulus.stair.strength);
    threshEst = mean(stimulus.stair.strength(stimulus.stair.reversals((numReversals-6):end)));
    
    axisVect = [0,numTrials,0,max(stimulus.stair.strength)+.01];
    
    figure, clf, hold on
    plot(stimulus.stair.strength)
    plot(stimulus.stair.reversals,stimulus.stair.strength(stimulus.stair.reversals),'bo')
    
    plot([0,length(stimulus.stair.strength)], [threshEst,threshEst],'k-.')
     axis(axisVect)
     hold off


end