function myscreen = pl_MRI_localizer(observer,varargin)

%%% main function that runs the task in the scanner using only a neutral
%%% cue
% The raised cosine reqires matlabPyrTools
% http://www.cns.nyu.edu/~lcv/software.php
% psuedorandomize across all combination of stims+ITI


global stimulus;
global MGL;

mglOpen
% check arguments
% if ~any(nargin == 3)
%     help transientAttention
%     return
% % end
% 
% eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye'}));

if ieNotDefined('indContrast'),indContrast = [0.5];end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('Eye'),Eye = 0;end % no eye-tracking

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

% clearing old variables:
clear task myscreen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initalize the screen

stimulus.EyeTrack=Eye;
myscreen = initScreen('disp2');
myscreen.datadir = datadirname;
myscreen.allowpause = 0;
myscreen.saveData = -2;
myscreen.background=.5;

if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}.waitForBacktick = 1;
task{1}.segmin =     [12 12];  % run the localizer for each diagonal in each block
task{1}.segmax =     [12 12];  % this assumes a 2s TR... will need to be changed if we change our TR. Basically 16TRs/cycle 
task{1}.segquant =   [0 0];
task{1}.getResponse = [0 0];
task{1}.synchToVol = ones(size(task{1}.segquant));
% task{1}.synchToVol(end) = 1;
task{1}.fudgeLastVolume = 1;

n_repeats = 3;


[contrast,location,repeats] = ndgrid(1,1:2,1:n_repeats);



%contrast =3 is blank trials. We wants on ~10% of total trials to be blank
%trials. Re-assign 4 out of 6 blank trials to be non-blank stim containing
%trials

task{1}.numTrials = length(location(:)); % n*n_repeats
random_order = randperm(task{1}.numTrials);
 
task{1}.randVars.targetLocation = location(random_order); %one of the 2 positions
task{1}.randVars.len_ = task{1}.numTrials;

task{1}.randVars.contrast = contrast(random_order);

%we want a block design. 
% task{1}.randVars.uniform.targetOrientation = 1:2;
% task{1}.randVars.iti= iti(random_order)
% task{1}.randVars.iti= task{1}.randVars.iti.*1.5 % replace if TR changes

stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.LocationIndices=unique(location);


task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen); %initStimulus('stimulus',myscreen,indContrast,diagonal);
stimulus = myInitStimulus(stimulus,myscreen,task,indContrast);

myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the fixation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the first task to be the fixation staircase task
global fixStimulus;
if ~easyFixTask
  % default values
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1;
  fixStimulus.fixLineWidth = 3;
  fixStimulus.stimTime = 0.4;
  fixStimulus.responseTime = 1;
else
  % make cross bigger and task slower
  fixStimulus.diskSize = 0.5;
  fixStimulus.fixWidth = 1+1*easyFixTask;
  fixStimulus.fixLineWidth = 3+2*easyFixTask;
  fixStimulus.stimTime = 0.4+0.4*easyFixTask;
  fixStimulus.responseTime = 1+1*easyFixTask;
end
global fixStimulus
fixStimulus.pos = [xOffset yOffset];
[task{1} myscreen] = fixStairInitTask(myscreen);
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


task.thistrial.targetLocation1 = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
task.thistrial.targetLocation2 = 2;

 (task.thistrial.thisseg == 1) % fixation
    iti =task.thistrial.iti;%iti = .6;
    task.thistrial.seglen =[0.1 .06 .04 .1 .3 .3 .64 .03 iti];
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
    
    
(task.thistrial.thisseg == 8) % response
    stimulus.trialnum = stimulus.trialnum + 1;
    if ~task.thistrial.gotResponse
       %mglPlaySound(stimulus.noanswer);
    end;


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
%     if stimulus.EyeTrack, fixCheck; end
elseif (task.thistrial.thisseg == 2) % Pre Cue
    drawFixation(task);
%     if stimulus.EyeTrack, fixCheck; end
    drawPreCue(task.thistrial.targetLocation);
    
elseif (task.thistrial.thisseg == 3) % ISI 1 
    drawFixation(task);
%     if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 4) % Stimulus
    drawFixation(task);
%     if stimulus.EyeTrack, fixCheck; end
    drawGabor(stimulus.contrasts(task.thistrial.contrast),...
              stimulus.tmp.targetLocation,...
              (stimulus.orientation+(stimulus.rotation(task.thistrial.targetOrientation)*stimulus.indTilt)) ,1);
          % the above line of code adds or subtracts the tilt from the base
          % orientation
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation(task);
%     if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation(task);
%     if stimulus.EyeTrack, fixCheck; end
    drawRespCue(task.thistrial.targetLocation);

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
%
if ~task.thistrial.gotResponse
    
    % check response correct or not
    if task.thistrial.contrast ==  3 %cue-only
        stimulus.tmp.response = task.thistrial.whichButton == 3; %press 3 to have the same motor response as in the main conditions
    else
        stimulus.tmp.response = task.thistrial.whichButton == (task.thistrial.targetOrientation); %1 for left and 2 for right
    end;
    
end

end

function drawRespCue(loc)
    global stimulus
    
    mglLines2(stimulus.respcueLocation{loc}(1), stimulus.respcueLocation{loc}(3),...
              stimulus.respcueLocation{loc}(2), stimulus.respcueLocation{loc}(4),stimulus.respCue.width,stimulus.black);
    
end
