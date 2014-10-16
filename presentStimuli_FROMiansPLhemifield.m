function presentStimuli(wPtr, trialAngles, locations, cont, phase_offset)
global params; 
% Inpute:
% Angles - is an array of angles to present each gabor
% locations - is a matrix of n*2 location of each corresponding stimulus
% This function will draw all the stimuli rotated by angles in the specific
% location
if (nargin<4), contrast = 255*ones(length(angles),1); end
if (nargin<5), phase_offset = zeros(length(angles),1); end
if (size(trialAngles,2) ~= params.stim.num), error('In order to present all stimuli correctly, the amount of angles should match the amount of stimuli pressented.'); end
if (size(locations,3) ~= params.stim.num), error('In order to present all stimuli correctly, the amount of locations should match the amount of stimuli pressented.'); end

%Compute grating
grating = genGrating(params.stim.sizePix, params.stim.cyclesPerImage, phase_offset);
mask = My2DGauss_nonSym(params.stim.sizePix ,0,2);
gabor = grating.*mask;

for i =1:params.stim.num
    texture = params.stim.bkColor + 128*cont(i).*gabor;
    
    %display grating
    rect =  CenterRectOnPointd(params.stim.rectPix, locations(i,1), locations(i,2));%rect of stimulus, to be centered at location
    
    
    if params.stim.colorTest == 1
        Screen('FrameOval', wPtr, params.stim.color(trialAngles,:), rect, 3, 3);
    else
    textureIndex=Screen('MakeTexture', wPtr, texture); 
    Screen('DrawTexture', wPtr, textureIndex, [], rect, trialAngles(i));
    end
end

