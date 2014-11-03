
particpant = 'HL';
  
 
 
oriThresh = oriThreshEst;
contThresh = contThreshEst;
numBlocks = 6;


diagonals = [2,1];
diagonals = repmat(diagonals,[1,numBlocks/length(diagonals)]);

for block = 1:numBlocks
    pl_MRI_constStim(particpant,'Eye',1,'indTilt',oriThresh,'indContrast',contThresh,'diagonal',diagonals(block))
    
end
