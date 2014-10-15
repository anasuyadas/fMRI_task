function myscreen = pl_MRI_constStim(observer,varargin)

%%% Pilot with only valid and invalid conditions
% The raised cosine reqires matlabPyrTools
% http://www.cns.nyu.edu/~lcv/software.php
% psuedorandomize across all combination of stims+ITI


global stimulus;
global MGL;


% check arguments
% if ~any(nargin == 3)
%     help transientAttention
%     return
% % end
% 
eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye'}));

if ieNotDefined('indContrast'),indContrast = .4;end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('indTilt'),indTilt = 10;end % default tilt
if ieNotDefined('Eye'),Eye = 0;end % no eye-tracking
if ieNotDefined('cueType'),cueType = 0;end

contLevels = makeContLevels(indContrast); %  min = 1.5% , max = 80%

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
mglVisualAngleCoordinates(156.5,[36.5,26.5]) %distance from screen
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
task{1}.segmin =     [0.1 .06 .04 .1 .3 .3 .8 .03 2];  % segments: 1:fixation, 2:cue 3: ISI 4:stimulus,5: post-stim ISI, 6:response cue 7:response, 8:feedback dur, 9:ITI code
task{1}.segmax =     [0.1 .06 .04 .1 .3 .3 .8 .03 2];  
task{1}.segquant =   [0 0 0 0 0 0 0 0 0]; % I guess, ITI varies in steps of 0.25
task{1}.getResponse = [0 0 0 0 0 0 1 0 0]; % responses are allowed during response intervals



n_repeats = 2;%  trials per block n= 36; 3contrast*3ITIs*2location 
% Number of volumes = (n)+(n/3*2)+(n/3*3)+(n/3*4).
%n_repeats will have to be adjusted depending on our TR to keep block
%length approximately ~5minutes

if diagonal == 1  
    [contrast,ori,location,trialNum] = ndgrid(1:7,1:2,[1,3],1:n_repeats);
%     [contrast2,ori2,location2] = ndgrid(1:7,1:2,[1,3]);
else 
    [contrast,ori,location,trialNum] = ndgrid(1:7,1:2,[2,4],1:n_repeats);
%     [contrast2,ori2,location2] = ndgrid(1:7,1:2,[2,4]);
end
%contrast =3 is blank trials. We wants on ~10% of total trials to be blank
%trials. Re-assign 4 out of 6 blank trials to be non-blank stim containing
%trials
% contrast(3,:,[1:2],1)=1;
% contrast(3,:,[1:2],2)=2; 


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
stimulus.trialAttemptNum = 0;
stimulus.numFixBreak = 0;
stimulus.fixationBreakTrialVect = 0;
stimulus.LocationIndices=unique(location);

stimulus.indTilt=indTilt;
stimulus.preCue.type = cueType;

task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback,@recalibrateCallback);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,contLevels);
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
elseif (task.thistrial.thisseg == 1) % fixation
    iti = .6;%task.thistrial.iti;
    task.thistrial.seglen =[0.1 .06 .04 0.1 .3 .3 .8 .03 iti];
    %need to make sure that there are only two locations per run
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{stimulus.randVars.targetLocation(task.thistrial.trialIndex)}];%[stimulus.locations{task.thistrial.targetLocation}];
    
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.respcueLocation{stimulus.randVars.targetLocation(task.thistrial.trialIndex)}; %if polygon
    stimulus.tmp.respcueLocation=stimulus.randVars.targetLocation(task.thistrial.trialIndex); %if central x
    %stimulus.tmp.WedgeStart=stimulus.CueWedges(task.thistrial.targetLocation);
    
    %just neutral cues - no exo cues
    for i=1:2
        stimulus.tmp.preCueNeutLocation{i}=stimulus.preCueNeutLocation{i};
    end
    
    
elseif (task.thistrial.thisseg == 8) % response
    stimulus.trialnum = stimulus.trialnum + 1;
%     if ~task.thistrial.gotResponse
%         mglPlaySound(stimulus.noanswer);
%     end;
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
    drawFixation(task);
elseif (task.thistrial.thisseg == 1) % Initial Fixation
    stimulus.FixationBreak(task.trialnum) = 0;
    stimulus.FixationBreakCurrent = 0;
    
    drawFixation(task);
    %disp(sprintf('total num of trials so far %f',task.trialnum))
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
elseif (task.thistrial.thisseg == 2) % Pre Cue
    drawFixation(task);
    
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
    if ~stimulus.FixationBreakCurrent  || ~stimulus.EyeTrack
    drawPreCue(stimulus.randVars.targetLocation(task.thistrial.trialIndex));
    end
    
elseif (task.thistrial.thisseg == 3) % ISI 1
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
    
elseif (task.thistrial.thisseg == 4) % Stimulus
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
    % the contrast value is the threshold itself
    if ~stimulus.FixationBreakCurrent  || ~stimulus.EyeTrack
    drawGabor(stimulus.contrasts(stimulus.randVars.contrast(task.thistrial.trialIndex)),...
              stimulus.tmp.targetLocation,...
              ((stimulus.rotation(stimulus.randVars.targetOrientation(task.thistrial.trialIndex))*stimulus.indTilt)),1,...
              task.thistrial.trialIndex);
    end
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck(myscreen,task); end
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
%% makeContLevels
function [contLevels] = makeContLevels(indThresh,varargin)
%makes a vector of contrast values (between 0 and 1), dependent on the
%contrast threshold 'indThresh'
%max & min contrast


% threshCont = .10;
threshCont = indThresh;

if ieNotDefined('minCont'), minCont = .015; end
if ieNotDefined('maxCont'), maxCont = .80; end


cont2 = 10^(log10(minCont*100)+.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont3 = 10^(log10(threshCont*100)-.25*(log10(threshCont*100) - log10(minCont*100)))/100;
cont5 = 10^(log10(threshCont*100)+.25*(log10(maxCont*100) - log10(threshCont*100)))/100;
cont6 = 10^(log10(maxCont*100)-.25*(log10(maxCont*100) - log10(threshCont*100)))/100;


contLevels = [minCont,cont2,cont3,threshCont,cont5,cont6,maxCont];
end

%% recalibrateCallback
function [task, myscreen] = recalibrateCallback(task,myscreen)
global stimulus

stimulus.trialAttemptNum = stimulus.trialAttemptNum+1;

if stimulus.FixationBreakCurrent
    
    stimulus.numFixBreak = stimulus.numFixBreak+1;
    stimulus.fixationBreakTrialVect(stimulus.numFixBreak) = stimulus.trialAttemptNum;
    
    if  stimulus.numFixBreak < 2
        stimulus.FixationBreakRecent = 0;
    elseif (stimulus.fixationBreakTrialVect(stimulus.numFixBreak) - stimulus.fixationBreakTrialVect(stimulus.numFixBreak-1)) < 3
        
        if stimulus.FixationBreakRecent < 3
            
            stimulus.FixationBreakRecent = stimulus.FixationBreakRecent+1;
            
        else stimulus.FixationBreakRecent
            
            stimulus.FixationBreakRecent = 0;
            myscreen = eyeCalibDisp(myscreen,'Press SPACEBAR, then C to calibrate');
            
        end
    else
        stimulus.FixationBreakRecent = 0;
    end
end

end
