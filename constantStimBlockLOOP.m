
particpant = 'AD';
  
trainDiagonal = 1;

oriThresh = oriThreshEst;
contThresh = contThreshEst;
numBlocks = 6;


diagonals = [1,2];
diagonals = repmat(diagonals,[1,numBlocks/length(diagonals)]);

for block = 1:numBlocks
    pl_MRI_constStim(particpant,'Eye',1,'diagonal',diagonals(block),...
                     'indTilt',oriThresh,'indContrast',contThresh,'trainDiagonal',trainDiagonal)
    
end
