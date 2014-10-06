function myscreen = localizer_2stims(observer,Eye)

%%% Pilot with only valid and invalid conditions


global stimulus;
global MGL;

% check arguments
if ~any(nargin == 2)
    help transientAttention
    return
end

thisdir = pwd;
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'data'))
    disp('Making data directory');
    mkdir('data');
end

% make an observer directory if necessary
datadirname = fullfile(thisdir,'data',observer);
if ~isdir(datadirname);
    disp(sprintf('Making observer directory %s',datadirname));
    mkdir(datadirname);
end

disp(sprintf('[ DATA ] saving data in: %s',datadirname));

stimulus = [];
stimulus.EyeTrack=Eye;

% clearing old variables:
clear task myscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen

myscreen = initScreen('CMU_CRT');
myscreen.datadir = datadirname;
myscreen.allowpause = 0;
myscreen.saveData = -2;
myscreen.background=.5;

if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.segmin =     [8.5 8.5];  % segments: 1:right, 2:down, 3:left 4:up
task{1}.segmax =     [8.5 8.5];
task{1}.segquant =   [0 0];
task{1}.getResponse = [0 0];
task{1}.synchToVol = ones(size(task{1}.segquant));
% task{1}.synchToVol(end) = 1;
task{1}.fudgeLastVolume = 1;

n_repeats = 14; % ((9/1.5)*2*14)*1.5/60 = 168TR = 4.2mins

[contrast, location, repeat] = ndgrid(1,1,1:n_repeats);
task{1}.numTrials = length(location(:));
random_order = randperm(task{1}.numTrials);
task{1}.randVars.contrast = contrast(random_order);
task{1}.randVars.targetLocation = location(random_order); %one of the 4 positions
task{1}.randVars.uniform.targetOrientation = 1:2;
task{1}.randVars.uniform.distractorOrientation1 = 1:2;
task{1}.randVars.uniform.distractorOrientation2 = 1:2;
task{1}.randVars.uniform.distractorOrientation3 = 1:2;
task{1}.randVars.uniform.random_cueonly_cond = 1:4;
task{1}.randVars.len_ = task{1}.numTrials;
stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.LocationIndices=unique(location);

task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task);

myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
clear stimulus.tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = StartSegmentCallback(task, myscreen)
global stimulus;

task.thistrial.count = 0;

%down
task.thistrial.targetLocation1 = 1;
task.thistrial.targetLocation2 = 2;

task.thistrial.phase1 = randsample(1:20:360,1);
task.thistrial.phase2 = randsample(1:20:360,1);

stimulus.Tilt = randsample(-6:2:6,1);
stimulus.rotation = [stimulus.Tilt -stimulus.Tilt];

stimulus.tmp.targetLocation1(1)  = stimulus.eccentricity_h*[stimulus.locations{task.thistrial.targetLocation1}(1)];
stimulus.tmp.targetLocation1(2)  = stimulus.eccentricity_v*[stimulus.locations{task.thistrial.targetLocation1}(2)];

stimulus.tmp.targetLocation2(1)  = stimulus.eccentricity_h*[stimulus.locations{task.thistrial.targetLocation2}(1)];
stimulus.tmp.targetLocation2(2)  = stimulus.eccentricity_v*[stimulus.locations{task.thistrial.targetLocation2}(2)];

stimulus.FixationStarted=0;

mglClearScreen(stimulus.grayColor);
setGammaTable(1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;

task.thistrial.count = task.thistrial.count + 1;

if task.thistrial.count == 25
    task.thistrial.phase1 = randsample(1:20:360,1);
    task.thistrial.phase2 = randsample(1:20:360,1);
    stimulus.Tilt = randsample(-6:2:6,1);
    stimulus.rotation = [stimulus.Tilt -stimulus.Tilt];
    task.thistrial.count = 0;
end

mglClearScreen(stimulus.grayColor);

for cross=1:4
    mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
end

if (task.thistrial.thisseg == 1)
    drawGabor(stimulus.contrasts(task.thistrial.contrast),myscreen,stimulus.tmp.targetLocation1, stimulus.rotation(task.thistrial.targetOrientation),1,task.thistrial.phase1);
    drawGabor(stimulus.contrasts(task.thistrial.contrast),myscreen,stimulus.tmp.targetLocation2, stimulus.rotation(task.thistrial.targetOrientation),1,task.thistrial.phase2);
elseif (task.thistrial.thisseg == 2)
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the observer's response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor);
if ~task.thistrial.gotResponse
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function drawGabor(desiredContrast,myscreen,position,orientation,sf,p)
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

mglBltTexture(stimulus.phaseGrating{p}.tex{sf}(displayContrastNum+1),position,0,0,orientation); %mglBltTexture(texture,position,hAlignment,vAlignment,rotation)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to create a gamma table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setGammaTable(maxContrast)
global stimulus;

% set the reserved colors
gammaTable(1:size(stimulus.reservedColors,1),1:size(stimulus.reservedColors,2))=stimulus.reservedColors;
% create the gamma table
cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
luminanceVals = cmin:((cmax-cmin)/(stimulus.nGratingColors-1)):cmax;
% now get the linearized range
redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');

% set the gamma table
gammaTable((stimulus.minGratingColors+1):256,:)=[redLinearized;greenLinearized;blueLinearized]';
% set the gamma table
mglSetGammaTable(gammaTable);
% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)
global MGL;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gabors
stimulus.width = 4;                                  % in deg
stimulus.height = 4;
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg

stimulus.DistractorRotation = [0];
stimulus.init = 1;

stimulus.spf = 4;                % in cpd
stimulus.theorientation = 90;      % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity_h = 5;  % in deg
stimulus.eccentricity_v = 2.65;  % in deg

% fixation
stimulus.FCwidth = 0.3;
stimulus.FClinewidth = 3;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate befire trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=1;
stimulus.edgeDist=0;%presents stimuli after this duration when fixation detected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames and locations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% locations 1:10  starts at upper right, moves around clockwise
stimulus.locations = {[cosd(45), -sind(45)],[-cosd(45), -sind(45)],[-cosd(45), sind(45)],[cosd(45), sind(45)]}; %4 locations: x,y coordinates are specified here. cosd is in degrees.. This will be multiplied by eccentricity
for i=1:length(stimulus.locations)
    stimulus.placeholders{i}= [stimulus.eccentricity_h*stimulus.locations{i} 0 0 0 0];
    stimulus.placeholders{i}(2,:)=[stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[(stimulus.cornerDist-stimulus.edgeDist) stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist (stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(3,:)=[stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[(stimulus.cornerDist-stimulus.edgeDist) -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist -(stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(4,:)=[stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-(stimulus.cornerDist-stimulus.edgeDist) -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist -(stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(5,:)=[stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-(stimulus.cornerDist-stimulus.edgeDist) stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist (stimulus.cornerDist-stimulus.edgeDist)]];
end

stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];

row=[4 5 2 3];
signs{1}=[-1 -1];signs{2}=[-1 1];signs{3}=[1 1];signs{4}=[1 -1];

for i=1:4
    stimulus.centralx{i}(3:4)=.3*[stimulus.locations{i}];
end

stimulus.contrasts =1;

stimulus.nReservedColors=size(stimulus.reservedColors,1);
stimulus.nGratingColors = 256-(2*floor(stimulus.nReservedColors/2)+1);
stimulus.minGratingColors = 2*floor(stimulus.nReservedColors/2)+1;
stimulus.midGratingColors = stimulus.minGratingColors+floor(stimulus.nGratingColors/2);
stimulus.maxGratingColors = 255;
stimulus.deltaGratingColors = floor(stimulus.nGratingColors/2);

% to set up color values
stimulus.black = [0 0 0];
stimulus.white = [1/255 1/255 1/255];
stimulus.green = [0 160 0];
stimulus.blue = [0 0 160];
stimulus.greencorrect = [0 255 0];
stimulus.redincorrect = [255 0 0];
stimulus.orangenoanswer = [255 215 0];
stimulus.grey = [.025 .025 .025];
stimulus.background = [stimulus.midGratingColors/255 stimulus.midGratingColors/255 stimulus.midGratingColors/255];
stimulus.fixation.color = [0; .6; 0]'; % green


% calculate a grating, a gaussian envelope (gaussian is in the alpha
% channel), and a mask (for now, just high-contrast random noise)
for thisSF = 1:length(stimulus.spf)      %only one spatial frequency
    for p = 1:20:360
        stimulus.gratingMatrix{thisSF}.phaseGrating{p} = mglMakeGrating(stimulus.width,stimulus.height,stimulus.spf(thisSF),stimulus.theorientation,p);
    end
end

%%% Makes a grating windowed with raised cosine
xpxpcm = myscreen.screenWidth/myscreen.displaySize(1);
ypxpcm = myscreen.screenHeight/myscreen.displaySize(2);

xpxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*xpxpcm);
ypxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*ypxpcm);

for p = 1:20:360
    res = mkR([size(stimulus.gratingMatrix{1}.phaseGrating{p},1) size(stimulus.gratingMatrix{1}.phaseGrating{p},2)]);
    stimulus.sizedg = 3;
    [Xtbl,Ytbl] = rcosFn(30, (stimulus.sizedg*xpxpdeg)/2, [1 0]); %2 = sharp transition (edge effect?) / 50 = radius of the circle in pixel
    stimulus.phaseGrating{p}.grating(:,:,4) = 255*pointOp(res, Ytbl, Xtbl(1), Xtbl(2)-Xtbl(1), 0);
end

% making the texture for all the Gabor stimuli:
disppercent(-inf,'Calculating gabors');
for thisSF = 1:length(stimulus.spf)
    for thisContrast = 0:stimulus.deltaGratingColors
        for p = 1:20:360
            %% stimulus.texture
            stimulus.phaseGrating{p}.grating(:,:,1) = stimulus.midGratingColors+stimulus.gratingMatrix{thisSF}.phaseGrating{p}*thisContrast;
            stimulus.phaseGrating{p}.grating(:,:,2) = stimulus.phaseGrating{p}.grating(:,:,1);
            stimulus.phaseGrating{p}.grating(:,:,3) = stimulus.phaseGrating{p}.grating(:,:,1);
            stimulus.phaseGrating{p}.tex{thisSF}(thisContrast+1) = mglCreateTexture(stimulus.phaseGrating{p}.grating);
            disppercent(thisContrast/stimulus.deltaGratingColors);
        end
    end
end
disppercent(inf);
stimulus.nDisplayContrasts = stimulus.deltaGratingColors;
disppercent(inf);

% calculate gray color
stimulus.grayColor = stimulus.background; %stimulus.midGratingColors/255;

stimulus.namecond = 'exo';
