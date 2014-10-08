function myscreen = pl_MRI_stair_con(observer,varargin)

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

if ieNotDefined('indContrast'),indContrast = .8;end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('indTilt'),indTilt = 10;end % default tilt
if ieNotDefined('Eye'),Eye = 0;end % no eye-tracking
if ieNotDefined('cueType'),cueType = 0;end

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
mglVisualAngleCoordinates(57) %distance from screen
if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
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



n_repeats = 15;%  trials per block n= 36; 3contrast*3ITIs*2location 
% Number of volumes = (n)+(n/3*2)+(n/3*3)+(n/3*4).
%n_repeats will have to be adjusted depending on our TR to keep block
%length approximately ~5minutes

if diagonal == 1  
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[1,3],1:n_repeats);
else 
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[2,4],1:n_repeats);
end
%contrast =3 is blank trials. We wants on ~10% of total trials to be blank
%trials. Re-assign 4 out of 6 blank trials to be non-blank stim containing
%trials
% contrast(3,:,[1:2],1)=1;
% contrast(3,:,[1:2],2)=2; 


task{1}.numTrials = length(location(:)); % n*n_repeats
random_order = randperm(task{1}.numTrials);
 
task{1}.randVars.targetLocation = location(random_order); %one of the 2 positions
task{1}.randVars.len_ = task{1}.numTrials;
task{1}.randVars.contrast = contrast(random_order);
task{1}.randVars.targetOrientation = ori(random_order);

stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.LocationIndices=unique(location);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAIRCASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set individual tilt

stimulus.indTilt=indTilt;
stimulus.preCue.type = cueType;
stair.upRule = 1;
stair.downRule = 2;
stair.startThresh = 40;
stair.stepSize = 10;
stair.minStepSize = .1;
stair.halfRule = 'levitt';

stimulus.stair = upDownStaircase(stair.upRule,stair.downRule,stair.startThresh,[stair.stepSize, stair.minStepSize],stair.halfRule);
stimulus.stair.minThreshold = .5; stimulus.stair.maxThreshold = 50;


task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);%initStimulus('stimulus',myscreen,indContrast,diagonal);
stimulus = myInitStimulus(stimulus,myscreen,task,indContrast);

myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
    % update the task
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task,myscreen,phaseNum] = updateTask(task,myscreen,phaseNum);
    
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
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
    
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.respcueLocation{task.thistrial.targetLocation}; %if polygon
    stimulus.tmp.respcueLocation=task.thistrial.targetLocation; %if central x
    %stimulus.tmp.WedgeStart=stimulus.CueWedges(task.thistrial.targetLocation);
    
    %just neutral cues - no exo cues
    for i=1:2
        stimulus.tmp.preCueNeutLocation{i}=stimulus.preCueNeutLocation{i};
    end
    
    
elseif (task.thistrial.thisseg == 8) % response
    stimulus.trialnum = stimulus.trialnum + 1;
    if ~task.thistrial.gotResponse
        mglPlaySound(stimulus.noanswer);
    end;
end

mglClearScreen(stimulus.grayColor);
setGammaTable(1);
end


%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;

mglClearScreen(stimulus.grayColor);%###

if (task.thistrial.thisseg == 9) % ITI
    drawFixation(task);
    
elseif (task.thistrial.thisseg == 1) % Initial Fixation
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
elseif (task.thistrial.thisseg == 2) % Pre Cue
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    drawPreCue(task.thistrial.targetLocation);
    
elseif (task.thistrial.thisseg == 3) % ISI 1
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 4) % Stimulus
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    % the contrast value is the threshold itself
    drawGabor(stimulus.stair.threshold/100,...
              stimulus.tmp.targetLocation,...
              ((stimulus.rotation(task.thistrial.targetOrientation)*stimulus.indTilt)),1);
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    drawRespCue(task.thistrial.targetLocation); % has to be a positive integer

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
    if task.thistrial.contrast ==  3 %cue-only
        stimulus.tmp.response = task.thistrial.whichButton == 3; %press 3 to have the same motor response as in the main conditions
    else
        stimulus.tmp.response = task.thistrial.whichButton == (task.thistrial.targetOrientation); %1 for left and 2 for right
    end;
    
end
stimulus.stair = upDownStaircase(stimulus.stair,stimulus.tmp.response);
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
