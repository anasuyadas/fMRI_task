
particpant = 'idTestALL_5';

oriThresh = 1.34;%oriThreshEst;
contThresh = .065;%contThreshEst;
numBlocks = 2;


diagonals = [1,2];
diagonals = repmat(diagonals,[1,numBlocks/length(diagonals)]);

for block = 1:numBlocks
    pl_MRI_constStim(particpant,'Eye',1,'indTilt',oriThresh,'indContrast',contThresh,'diagonal',diagonals(block))
    
end
