function drawGabor(desiredContrast,position,orientation,sf)
% drawGaborPedCdeltaC
%
%        $Id: drawGabor.m, v 1 2007/01/18 19:40:56 ?? ?? $
%      usage: drawGabor(desiredContrast,position,orientation,sf)
%    purpose: draw a gabor stimulus on the screen with a specified contrast
%             within a given clut range (it finds the closest contrast to
%             the requested one within the available clut range)

global stimulus;
% now find closest matching contrast we can display with this gamma table
displayContrastNum = min(round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast),stimulus.nDisplayContrasts);
% disp(sprintf('Desired contrast: %0.4f Actual contrast: %0.4f Difference: %0.4f',desiredContrast,stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts),desiredContrast-stimulus.currentMaxContrast*(displayContrastNum/stimulus.nDisplayContrasts)));
if round(stimulus.nDisplayContrasts*desiredContrast/stimulus.currentMaxContrast)>stimulus.nDisplayContrasts
    disp(sprintf('[drawGabor] Desired contrast out of range %0.2f > %0.2f',desiredContrast,stimulus.currentMaxContrast));
    keyboard
end

mglBltTexture(stimulus.tex{sf}(displayContrastNum+1),position,0,0,orientation); %mglBltTexture(texture,position,hAlignment,vAlignment,rotation)

disp(sprintf('orientation for this trial is %d',orientation));
end