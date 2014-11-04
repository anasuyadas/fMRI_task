function myscreen = pl_MRI_stair_Ori(observer,varargin)

%%% Pilot with only valid and invalid conditions
% The raised cosine reqires matlabPyrTools
% http://www.cns.nyu.edu/~lcv/software.php
% psuedorandomize across all combination of stims+ITI


global stimulus;
global MGL;

eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye'}));

if ieNotDefined('indContrast'),indContrast = .8;end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('indTilt'),indTilt = 10;end % default tilt
if ieNotDefined('Eye'),Eye = 0;end % no eye-tracking
if ieNotDefined('cueType'),cueType = 0;end

%contLevels = makeContLevels(indContrast); %  min = 1.5% , max = 80%

 
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
stimulus.Tilt=indTilt;

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
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize)
if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
    myscreen.eyetracker.savedata = true;%%%%% TO ADD FOR ONLINE EYETRACKING
    myscreen.eyetracker.data = [1 1 1 0];%%%%% TO ADD FOR ONLINE EYETRACKING
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 0;
task{1}.segmin =     [0.5 .06 .04 .1 .4 .4 1 .03 2];  % segments: 1:fixation, 2:cue 3: ISI 4:stimulus,5: post-stim ISI, 6:response cue 7:response, 8:feedback dur, 9:ITI code
task{1}.segmax =     [0.5 .06 .04 .1 .4 .4 1 .03 2];  
task{1}.segquant =   [0 0 0 0 0 0 0 0 0]; % I guess, ITI varies in steps of 0.25
task{1}.getResponse = [0 0 0 0 0 0 1 0 0]; % responses are allowed during response intervals



n_repeats = 15; %num trials = num(conts*oris*locs*repeats)

if diagonal == 1  
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[1,3],1:n_repeats);
else 
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[2,4],1:n_repeats);

end



task{1}.numTrials = length(location(:)); % n*n_repeats
task{1}.origNumTrials = length(location(:)); % n*n_repeats
random_order = randperm(task{1}.numTrials);
 
stimulus.randVars.targetLocation = location(random_order); %one of the 2 positions
stimulus.randVars.contrast = contrast(random_order);
stimulus.randVars.targetOrientation = ori(random_order);


task{1}.randVars.len_ = task{1}.numTrials;
task{1}.randVars.trialIndex = random_order;


stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.FixationBreakCurrent = 0;
stimulus.FixationBreakRecent= 0;
stimulus.trialAttemptNum = 1;
stimulus.numFixBreaks = 0;
stimulus.fixationBreakTrialVect = 0;
stimulus.LocationIndices=unique(location);
stimulus.upDated = 1;
stimulus.fixBreakTRACKindex = 0;
stimulus.testFix1 = 0;
stimulus.firstFixBreak = 0;
stimulus.increasedAttemptNum= 0;

stimulus.indTilt=indTilt;
stimulus.preCue.type = cueType;

task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback,@recalibrateCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAIRCASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set individual tilt

stair.upRule = 1;
stair.downRule = 3;
stair.startThresh = 10;
stair.stepSize = 2;
stair.minStepSize = .1;
stair.halfRule = 'levitt';

stimulus.stair = upDownStaircase(stair.upRule,stair.downRule,stair.startThresh,[stair.stepSize, stair.minStepSize],stair.halfRule);
stimulus.stair.minThreshold = .04; stimulus.stair.maxThreshold = 10;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,indContrast);

myscreen = eyeCalibDisp(myscreen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;

while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task,myscreen,phaseNum] = updateTaskHack(task,myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end
clear stimulus.tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = endTask(myscreen,task);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TASK 1: function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = StartSegmentCallback(task, myscreen)
% segments: 1:ITI,   2:fixation,    3:stimulus, 4:response
global stimulus

if (task.thistrial.thisseg == 9) % ITI
    stimulus.trialend = stimulus.trialend + 1;
    stimulus.increasedAttemptNum = 0;
    stimulus.testFix1 = 0;
elseif (task.thistrial.thisseg == 1) % fixation
    iti = .6;%task.thistrial.iti;
    task.thistrial.seglen =[0.5 .06 .04 .1 .4 .4 1 .03 2];
    %need to make sure that there are only two locations per run
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{stimulus.randVars.targetLocation(task.thistrial.trialIndex)}];
    
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.respcueLocation{stimulus.randVars.targetLocation(task.thistrial.trialIndex)}; %if polygon
    stimulus.tmp.respcueLocation=stimulus.randVars.targetLocation(task.thistrial.trialIndex); %if central x
    %stimulus.tmp.WedgeStart=stimulus.CueWedges(task.thistrial.targetLocation);
    
    %just neutral cues - no exo cues
    for i=1:2
        stimulus.tmp.preCueNeutLocation{i}=stimulus.preCueNeutLocation{i};
    end
    if ~stimulus.testFix1 
        stimulus.FixationBreak(task.trialnum) = 0;
        stimulus.FixationBreakCurrent = 0;
        stimulus.updateCurrent = 1;
        stimulus.upDated = 0;
        stimulus.testFix1  = 1;
    end
    
    if (1 < task.trialnum) && ~stimulus.increasedAttemptNum
        stimulus.trialAttemptNum = stimulus.trialAttemptNum+1;
        stimulus.increasedAttemptNum = 1;
    end
    
elseif (task.thistrial.thisseg == 8) % response
    stimulus.trialnum = stimulus.trialnum + 1;

end

mglClearScreen(stimulus.grayColor);
setGammaTable(1);
end




%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;

mglClearScreen(stimulus.grayColor);%###
stimulus.trialend = task.numTrials;
if (task.thistrial.thisseg == 9) % ITI
    stimulus.testFix1 = 0;
    drawFixation(task);
elseif (task.thistrial.thisseg == 1) % Initial Fixation
    
    
    drawFixation(task);
    
    if ~stimulus.testFix1
        stimulus.FixationBreak(task.trialnum) = 0;
        stimulus.FixationBreakCurrent = 0;
        stimulus.updateCurrent = 1;
        stimulus.testFix1  = 1;
    end
    
    if stimulus.EyeTrack && ~stimulus.FixationBreakCurrent, fixCheck(myscreen,task); end
elseif (task.thistrial.thisseg == 2) % Pre Cue
    drawFixation(task);
    
    if stimulus.EyeTrack && ~stimulus.FixationBreakCurrent, fixCheck(myscreen,task); end
    if ~stimulus.FixationBreakCurrent  || ~stimulus.EyeTrack
    drawPreCue(stimulus.randVars.targetLocation(task.thistrial.trialIndex));
    end
    
elseif (task.thistrial.thisseg == 3) % ISI 1
    drawFixation(task);
    if stimulus.EyeTrack && ~stimulus.FixationBreakCurrent, fixCheck(myscreen,task); end
    
elseif (task.thistrial.thisseg == 4) % Stimulus
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
    % the contrast value is the threshold itself
    if ~stimulus.FixationBreakCurrent  || ~stimulus.EyeTrack
        drawGabor(stimulus.contrasts(stimulus.randVars.contrast(task.thistrial.trialIndex)),...
              stimulus.tmp.targetLocation,...
              ((stimulus.rotation(stimulus.randVars.targetOrientation(task.thistrial.trialIndex))*stimulus.stair.threshold)),1,...
              task.thistrial.trialIndex);

    end
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation(task);
    if stimulus.EyeTrack && ~stimulus.FixationBreakCurrent, fixCheck(myscreen,task); end
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation(task);
    if stimulus.EyeTrack && ~stimulus.FixationBreakCurrent, fixCheck(myscreen,task); end
    if ~stimulus.FixationBreakCurrent  || ~stimulus.EyeTrack
        drawRespCue(stimulus.randVars.targetLocation(task.thistrial.trialIndex)); % has to be a positive integer
    end

elseif (task.thistrial.thisseg == 7) % Resp Window
    drawFixation(task);
    
elseif (task.thistrial.thisseg == 8) % Feedback
    drawFixation(task);
    
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% response call back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task, myscreen] = responseCallback(task,myscreen)
global stimulus;
mglClearScreen(stimulus.grayColor); %###
if ~task.thistrial.gotResponse
    
    % check response correct or not
stimulus.tmp.response = task.thistrial.whichButton == (stimulus.randVars.targetOrientation(task.thistrial.trialIndex)); %1 for left and 2 for right
    
end
if ~stimulus.FixationBreakCurrent || ~stimulus.EyeTrack
    stimulus.stair = upDownStaircase(stimulus.stair,stimulus.tmp.response);
end
disp(sprintf('threshold for this trial is %s',stimulus.stair.threshold));
end

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % staircase call back
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [task,stimulus] = staircaseCallback(task,stimulus)
% global stimulus;
% stimulus.stair = upDownStaircase(stimulus.stair,stimulus.tmp.response);
% 
% end
function drawRespCue(loc)
    global stimulus
    
    mglLines2(stimulus.respcueLocation{loc}(1), stimulus.respcueLocation{loc}(3),...
              stimulus.respcueLocation{loc}(2), stimulus.respcueLocation{loc}(4),stimulus.respCue.width,stimulus.black);
    
end

%% recalibrateCallback
function [task, myscreen] = recalibrateCallback(task,myscreen)
global stimulus


if stimulus.FixationBreakCurrent
    
    if  stimulus.numFixBreaks < 2
        
        stimulus.FixationBreakRecent = 0;
        stimulus.recalib(stimulus.trialAttemptNum) = 0;
    elseif (stimulus.fixationBreakTrialVect(stimulus.numFixBreaks) - stimulus.fixationBreakTrialVect(stimulus.numFixBreaks-1)) < 3
        
        if stimulus.FixationBreakRecent < 3
            
            stimulus.FixationBreakRecent = stimulus.FixationBreakRecent+1;
            stimulus.recalib(stimulus.trialAttemptNum) = 0;
        else
            stimulus.FixationBreakRecent = 0;
            stimulus.recalib(stimulus.trialAttemptNum) = 1;
%             myscreen = eyeCalibDisp(myscreen);
            eyeCalibDisp(myscreen);            
            
        end
            
    else
        stimulus.FixationBreakRecent = 0;
        stimulus.recalib(stimulus.trialAttemptNum) = 0;
    end
else
    stimulus.recalib(stimulus.trialAttemptNum) = 0;
end
end
