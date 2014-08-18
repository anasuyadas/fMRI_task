function myscreen = exo_mri_pilot_2stim(observer,contrast,IndTilt,Eye)

%%% Pilot with only valid and invalid conditions
% The raised cosine reqires matlabPyrTools
% http://www.cns.nyu.edu/~lcv/software.php
% psuedorandomize across all combination of stims+ITI


global stimulus;
global MGL;

% check arguments
if ~any(nargin == 3)
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
stimulus.Tilt=IndTilt;

% clearing old variables:
clear task myscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen

stimulus.EyeTrack=Eye;
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
task{1}.segmin =     [0.1 .06 .04 .1 .3 .3 .64 .03 1];  % segments: 1:fixation, 2:cue 3: ISI 4:stimulus,5: post-stim ISI, 6:response cue 7:response, 8:feedback dur, 9:ITI code
task{1}.segmax =     [0.1 .06 .04 .1 .3 .3 .64 .03 1];  % ITI is variable (btw. 1 and 1.5), rest is constant
task{1}.segquant =   [0 0 0 0 0 0 0 0 0]; % I guess, ITI varies in steps of 0.25
task{1}.getResponse = [0 0 0 0 0 0 1 0 0]; % responses are allowed during response intervals
task{1}.synchToVol = zeros(size(task{1}.segquant));
task{1}.synchToVol(end) = 0;
task{1}.fudgeLastVolume = 1;

n_repeats = 2; %  trials per block = 48 trials. Number of volumes = (48)+(16*2)+(16*3)+(16*4) = 192 = 4.8mins = 4mins and 48secs

% CHANGE: here as our conditions are very different
[contrast, iti, location, diagonal, repeat] = ndgrid(1:3,1:3,1:2,1:2,1:n_repeats);

% CHANGE: num trials different
task{1}.numTrials = length(location(:));
random_order = randperm(task{1}.numTrials);

task{1}.randVars.contrast = contrast(random_order);

% CHANGE: 
task{1}.randVars.pre_post = pre_post(random_order);

% CHANGE: we have just one cue type
task{1}.randVars.CueCondition = CueCondition(random_order); %1=valid, 2=invalidb% vector from 1:10 with 1:4:valid, 5:8 invalid, 9:10 blank
task{1}.randVars.iti = iti(random_order); %1, 2 or 3 TR


task{1}.randVars.targetLocation = location(random_order); %one of the 4 positions

% We only need to specify target orientation
task{1}.randVars.uniform.targetOrientation = 1:2;

% CHANGE: we dont need this
task{1}.randVars.uniform.distractorOrientation1 = 1:2;
task{1}.randVars.uniform.distractorOrientation2 = 1:2;
task{1}.randVars.uniform.distractorOrientation3 = 1:2;

% are these blank trials?
task{1}.randVars.uniform.random_cueonly_cond = 1:4;


task{1}.randVars.len_ = task{1}.numTrials;
stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.LocationIndices=unique(location);

%%% Randomise the position of the cue in the cue only condition
% count_cue_only = 1;
% task{1}.count_cue_only = count_cue_only;
% random_cueonly_cond = [];
% for a = 1:(task{1}.numTrials/10*2)
%     for nn = 1:length([1 2 3 4])
%         random_cueonly_cond = [random_cueonly_cond;nn];
%     end
% end
% random_cueonly_cond = random_cueonly_cond(1:(task{1}.numTrials/10*2));
% task{1}.random_cueonly_cond = shufflerows(random_cueonly_cond,1000);
%%%

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
% segments: 1:ITI,   2:fixation,    3:stimulus, 4:response
global stimulus;

if (task.thistrial.thisseg == 9) % ITI
    stimulus.trialend = stimulus.trialend + 1;
elseif (task.thistrial.thisseg == 1) % fixation
    
    if task.thistrial.iti == 1
        iti = 3;
    elseif task.thistrial.iti == 2
        iti = 4.5;
    elseif task.thistrial.iti == 3
        iti = 6;
    end
    
    if task.thistrial.pre_post==1% pre-cue
        task.thistrial.order_cue = 2;
        task.thistrial.order_stim = 4;
        task.thistrial.seglen = [0.1 .067 .053 .05 .1 .6 .5 .03 iti];
    elseif task.thistrial.pre_post==2% post-cue
        task.thistrial.order_cue = 2;
        task.thistrial.order_stim = 4;
        task.thistrial.seglen = [0.1 .05 .053 .067 .1 .6 .5 .03 iti];
    end;
    
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
    stimulus.tmp.distractorIndices=stimulus.LocationIndices(stimulus.LocationIndices~=task.thistrial.targetLocation);
    for Locs=1:length(stimulus.tmp.distractorIndices)
        stimulus.tmp.distractorLocations{Locs}= stimulus.eccentricity*[stimulus.locations{stimulus.tmp.distractorIndices(Locs)}];
    end
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.respcueLocation{task.thistrial.targetLocation}; %if polygon
    stimulus.tmp.respcueLocation=task.thistrial.targetLocation; %if central x
    stimulus.tmp.WedgeStart=stimulus.CueWedges(task.thistrial.targetLocation);
    for i=1:4
        stimulus.tmp.NeutralcueLocation{i}=stimulus.NeutralcueLocation{i};
    end
    if task.thistrial.CueCondition==1%valid
        stimulus.tmp.ExocueLocation=stimulus.ExocueLocation{task.thistrial.targetLocation};
    elseif task.thistrial.CueCondition==2 %invalid % CHANGE
        possibleLoc = [1 2];
        potentialLoc = possibleLoc(possibleLoc~=task.thistrial.targetLocation);
        stimulus.tmp.ExocueLocation=stimulus.ExocueLocation{randsample(potentialLoc,1)};
    end
elseif (task.thistrial.thisseg == 8) % response
    stimulus.trialnum = stimulus.trialnum + 1;
    if ~task.thistrial.gotResponse
        mglPlaySound(stimulus.noanswer);
    end;
end

mglClearScreen(stimulus.grayColor);
setGammaTable(1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;

mglClearScreen(stimulus.grayColor);%###

if (task.thistrial.thisseg == 9) % ITI
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
elseif (task.thistrial.thisseg == 1) % FIXATION
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
    if stimulus.EyeTrack
        ep=myscreen.eyetracker.eyepos;
        if (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist && ~stimulus.FixationStarted
            stimulus.FixationStart=mglGetSecs;
            stimulus.FixationStarted=1;
        elseif (sqrt(ep(end,1)^2+ep(end,2)^2))>=stimulus.TrialStartFixDist && stimulus.FixationStarted
            stimulus.FixationStarted=0;
        elseif (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist
            stimulus.FixationDur=mglGetSecs(stimulus.FixationStart);
            if stimulus.FixationDur >=stimulus.TrialStartFixDur
                task = jumpSegment(task);
            end
        end
    else
        
        %task = jumpSegment(task);
    end
elseif (task.thistrial.thisseg == task.thistrial.order_cue) % Exo cue
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
    mglLines2(stimulus.tmp.ExocueLocation(1), stimulus.tmp.ExocueLocation(2), stimulus.tmp.ExocueLocation(3), stimulus.tmp.ExocueLocation(4),6,stimulus.white);
elseif (task.thistrial.thisseg == 3) % ISI
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
elseif (task.thistrial.thisseg == task.thistrial.order_stim) % STIMULUS
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
    if stimulus.EyeTrack
        ep=myscreen.eyetracker.eyepos;
        if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
            stimulus.FixationBreak(stimulus.trialnum)=1;
        end
    end
    drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.targetLocation, stimulus.rotation(task.thistrial.targetOrientation), 1);
    for Loc=1:length(stimulus.tmp.distractorIndices)
        eval(sprintf('drawGabor(stimulus.contrasts(task.thistrial.contrast),stimulus.tmp.distractorLocations{Loc}, stimulus.rotation(task.thistrial.distractorOrientation%g), 1);',Loc));
    end
elseif (task.thistrial.thisseg == 5) % ISI
    %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
    for cross=1:4
        mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
    end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
elseif (task.thistrial.thisseg == 6) % RESPONSE CUE but no answer
    mglLines2(stimulus.centralx{stimulus.tmp.respcueLocation}(1), stimulus.centralx{stimulus.tmp.respcueLocation}(2), stimulus.centralx{stimulus.tmp.respcueLocation}(3), stimulus.centralx{stimulus.tmp.respcueLocation}(4),8,stimulus.black);
    
    possibleLoc = [1 2 3 4];
        potentialLoc = possibleLoc(possibleLoc~=stimulus.tmp.respcueLocation);
        if ~task.thistrial.gotResponse
            for cross=potentialLoc
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
            end
        else
            if stimulus.tmp.response == 1
                for cross=potentialLoc
                    mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
                end
            elseif stimulus.tmp.response == 0
                for cross=potentialLoc
                    mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.redincorrect);
                end
            elseif stimulus.tmp.response == 3
                for cross=potentialLoc
                    mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
                end
            end
        end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
elseif (task.thistrial.thisseg == 7) % The response cue disappears and they have 500ms to answer
    if ~task.thistrial.gotResponse
        %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
        for cross=1:4
            mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
        end
    else
        if stimulus.tmp.response == 1
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
            end
        elseif stimulus.tmp.response == 0
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.redincorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.redincorrect);
            end
        elseif stimulus.tmp.response == 3
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
            end
        end
    end
    for Gabor=stimulus.LocationIndices
        for corner=2:5
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,3),stimulus.placeholders{Gabor}(corner,4),1,stimulus.white);
            mglLines2(stimulus.placeholders{Gabor}(corner,1),stimulus.placeholders{Gabor}(corner,2),stimulus.placeholders{Gabor}(corner,5),stimulus.placeholders{Gabor}(corner,6),1,stimulus.white);
        end
    end
elseif (task.thistrial.thisseg == 8)
    if ~task.thistrial.gotResponse
        %mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white);
        for cross=1:4
            mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.white);
        end
    else
        if stimulus.tmp.response == 1
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
            end
        elseif stimulus.tmp.response == 0
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.redincorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.redincorrect);
            end
        elseif stimulus.tmp.response == 3
            %             mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
            for cross=1:4
                mglLines2(stimulus.centralx{cross}(1), stimulus.centralx{cross}(2), stimulus.centralx{cross}(3), stimulus.centralx{cross}(4),8,stimulus.greencorrect);
            end
        end
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the observer's response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor); %###
if ~task.thistrial.gotResponse
    
    % check response correct or not
    if task.thistrial.CueCondition==9 || task.thistrial.CueCondition==10 %cue-only
        stimulus.tmp.response = task.thistrial.whichButton == 3; %press 3 to have the same motor response as in the main conditions
    else
        stimulus.tmp.response = task.thistrial.whichButton == (task.thistrial.targetOrientation); %1 for left and 2 for right
    end;
    
    % give feeback:
    if stimulus.tmp.response == 1
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
        %mglPlaySound(stimulus.CorrectSound);
    elseif stimulus.tmp.response == 0
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.redincorrect);
        %mglPlaySound(stimulus.IncorrectSound);
    elseif stimulus.tmp.response == 3
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
        %mglPlaySound(stimulus.CorrectSound);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to draw the gabor stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% check the table
% plot(stimulus.linearizedGammaTable.redTable,'k');
% hold on
% plot(256*(0.25:0.5/250:0.75),redLinearized,'r');

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
stimulus.width = 4;%stimulus.gaussSdx*7;             % in deg
stimulus.height = 4;%stimulus.gaussSdy*7;            % in deg
%stimulus.gaussSdx = 1; %0.8;  %0.3; %0.5; %1; % stimulus.width/7;                % in deg
%stimulus.gaussSdy = 1; %0.8;  %0.3; %0.5; %1; % stimulus.height/7;               % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg

stimulus.rotation = [stimulus.Tilt -stimulus.Tilt]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.DistractorRotation = [0];
stimulus.init = 1;

stimulus.sf = 4;                % in cpd
stimulus.orientation = 90;      % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity = 4.6;    % in deg

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
    stimulus.placeholders{i}= [stimulus.eccentricity*stimulus.locations{i} 0 0 0 0];
    stimulus.placeholders{i}(2,:)=[stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[(stimulus.cornerDist-stimulus.edgeDist) stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist (stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(3,:)=[stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[(stimulus.cornerDist-stimulus.edgeDist) -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[stimulus.cornerDist -(stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(4,:)=[stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-(stimulus.cornerDist-stimulus.edgeDist) -stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist -(stimulus.cornerDist-stimulus.edgeDist)]];
    stimulus.placeholders{i}(5,:)=[stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-(stimulus.cornerDist-stimulus.edgeDist) stimulus.cornerDist] stimulus.placeholders{i}(1,1:2)+[-stimulus.cornerDist (stimulus.cornerDist-stimulus.edgeDist)]];
end
% Note: locations should specify a vector of length 1, i.e.,
%                       sum(locations.^2)==1
stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];

% CHANGE: 
stimulus.precue.size  = 1.0;% for the peripheral cues 
stimulus.precue.color = [1 1 1]; % white
stimulus.cueDist=1; %degs from outer edge of gabor
stimulus.CueWedges=[0 90 180 270];

stimulus.ExoCueSize=.16;

% Response cue
stimulus.rcue.XLocation{1} = [-.6; 0; -0.4; 0];
stimulus.rcue.YLocation{1} = [ 0; .6; 0; -.6];
stimulus.rcue.XLocation{2} = [.7; 0; 0.4; 0];
stimulus.rcue.YLocation{2} = [ 0; .6; 0; -.6];
stimulus.rcue.color        = [0; .6; 0]; % green
row=[4 5 2 3];
signs{1}=[-1 -1];signs{2}=[-1 1];signs{3}=[1 1];signs{4}=[1 -1];
for i=1:4
    stimulus.NeutralcueLocation{i}(1:2)=.8*[stimulus.locations{i}];
    stimulus.NeutralcueLocation{i}(3:4)=1.05*[stimulus.locations{i}];
    %stimulus.NeutralcueLocation{i}(1:2) = (stimulus.locations{i} * (stimulus.eccentricity+3)); 
    %stimulus.NeutralcueLocation{i}(3:4) = (stimulus.locations{i} * (stimulus.eccentricity+3)); 
end

for i=1:4
    %stimulus.centralx{i}(1:2)=.8*[stimulus.locations{i}];
    stimulus.centralx{i}(3:4)=.3*[stimulus.locations{i}];
    %stimulus.NeutralcueLocation{i}(1:2) = (stimulus.locations{i} * (stimulus.eccentricity+3)); 
    %stimulus.NeutralcueLocation{i}(3:4) = (stimulus.locations{i} * (stimulus.eccentricity+3)); 
end


% %%% If peripheral neutral cues
% stimulus.NeutralcueLocation{1}(1) = stimulus.NeutralcueLocation{1}(1)-1.99;
% stimulus.NeutralcueLocation{1}(3) = stimulus.NeutralcueLocation{1}(3)-1.82;
% stimulus.NeutralcueLocation{2}(1) = stimulus.NeutralcueLocation{2}(1)-1.99;
% stimulus.NeutralcueLocation{2}(3) = stimulus.NeutralcueLocation{2}(3)-1.82;
% stimulus.NeutralcueLocation{3}(1) = stimulus.NeutralcueLocation{3}(1)+1.99;
% stimulus.NeutralcueLocation{3}(3) = stimulus.NeutralcueLocation{3}(3)+1.82;
% stimulus.NeutralcueLocation{4}(1) = stimulus.NeutralcueLocation{4}(1)+1.99;
% stimulus.NeutralcueLocation{4}(3) = stimulus.NeutralcueLocation{4}(3)+1.82;

for i=1:4; 
%     stimulus.ExocueLocation{i}(1:2) = (stimulus.locations{i} * (stimulus.eccentricity+3.5)); 
%     stimulus.ExocueLocation{i}(3:4) = (stimulus.locations{i} * (stimulus.eccentricity+3.5)); 
    stimulus.ExocueLocation{i}(1:2) = (stimulus.locations{i}); 
    stimulus.ExocueLocation{i}(3:4) = (stimulus.locations{i}); 
end

% CHANGE: do we need placeholders
stimulus.ExocueLocation{1}(1) = stimulus.placeholders{1}(1,1)-.35;
stimulus.ExocueLocation{1}(3) = stimulus.placeholders{1}(1,1)+.35;
stimulus.ExocueLocation{2}(1) = stimulus.placeholders{2}(1,1)-.35;
stimulus.ExocueLocation{2}(3) = stimulus.placeholders{2}(1,1)+.35;
stimulus.ExocueLocation{3}(1) = stimulus.placeholders{3}(1,1)-.35;
stimulus.ExocueLocation{3}(3) = stimulus.placeholders{3}(1,1)+.35;
stimulus.ExocueLocation{4}(1) = stimulus.placeholders{4}(1,1)-.35;
stimulus.ExocueLocation{4}(3) = stimulus.placeholders{4}(1,2)+.35;

stimulus.respcueLocation{1}=[0.35 -0.55;0.55 -0.35;0.7 -0.7];
stimulus.respcueLocation{2}=[-0.35 -0.55;-0.55 -0.35;-0.7 -0.7];
stimulus.respcueLocation{3}=[-0.35 0.55;-0.55 0.35;-0.7 0.7];
stimulus.respcueLocation{4}=[0.35 0.55;0.55 0.35;0.7 0.7];

stimulus.contrasts =.07;

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
for thisSF = 1:length(stimulus.sf)      %only one spatial frequency
    gratingMatrix{thisSF} = mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf(thisSF),stimulus.orientation,stimulus.phase);
end

%%% Makes a gabor (grating windowed with a gaussian window)
% grating(:,:,4) = 255*mglMakeGaussian(stimulus.width,stimulus.height,stimulus.gaussSdx,stimulus.gaussSdy);

%%% Makes a grating windowed with raised cosine
xpxpcm = myscreen.screenWidth/myscreen.displaySize(1);
ypxpcm = myscreen.screenHeight/myscreen.displaySize(2);

xpxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*xpxpcm);
ypxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*ypxpcm);

res = mkR([size(gratingMatrix{1},1) size(gratingMatrix{1},2)]);
stimulus.sizedg = 3;
[Xtbl,Ytbl] = rcosFn(30, (stimulus.sizedg*xpxpdeg)/2, [1 0]); %2 = sharp transition (edge effect?) / 50 = radius of the circle in pixel
grating(:,:,4) = 255*pointOp(res, Ytbl, Xtbl(1), Xtbl(2)-Xtbl(1), 0);


% making the texture for all the Gabor stimuli:
disppercent(-inf,'Calculating gabors');
for thisSF = 1:length(stimulus.sf)
    for thisContrast = 0:stimulus.deltaGratingColors
        %% stimulus.texture
        grating(:,:,1) = stimulus.midGratingColors+gratingMatrix{thisSF}*thisContrast;
        grating(:,:,2) = grating(:,:,1);
        grating(:,:,3) = grating(:,:,1);
        stimulus.tex{thisSF}(thisContrast+1) = mglCreateTexture(grating);
        disppercent(thisContrast/stimulus.deltaGratingColors);
    end
end
disppercent(inf);
stimulus.nDisplayContrasts = stimulus.deltaGratingColors;
disppercent(inf);

% calculate gray color
stimulus.grayColor = stimulus.background; %stimulus.midGratingColors/255;

% sounds
stimulus.noanswer = find(strcmp(MGL.soundNames,'Ping'));
stimulus.CorrectSound = find(strcmp(MGL.soundNames,'Submarine'));
stimulus.IncorrectSound = find(strcmp(MGL.soundNames,'Basso'));

stimulus.namecond = 'exo';
