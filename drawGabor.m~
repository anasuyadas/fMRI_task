function drawGabor(desiredContrast,position,orientation,sf,trialIndex)
% drawGaborPedCdeltaC
%
%        $Id: drawGabor.m, v 1 2014/10/17 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,orientation,sf)
%    purpose: draw a gabor stimulus on the screen with a specified contrast

global stimulus;

if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
    disp(sprintf('[drawGabor] Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
    keyboard
   
end

mglBltTexture(stimulus.tex{trialIndex},position,0,0,orientation);

% disp(sprintf('orientation for this trial is %f',desiredContrast));
end