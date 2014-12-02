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
eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye','easyFixTask'}));

if ieNotDefined('indContrast'),indContrast = [0.5];end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('Eye'),Eye = 0;end % no eye-tracking
if ieNotDefined('easyFixTask'),easyFixTask = 1;end

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
myscreen = initScreen();
myscreen.background = 'gray';
myscreen.datadir = datadirname;
myscreen.allowpause = 0;
myscreen.saveData = -2;
myscreen.background=.5;
mglVisualAngleCoordinates(myscreen.displayDistance,myscreen.displaySize)

if stimulus.EyeTrack
    myscreen = eyeCalibDisp(myscreen);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin =     [1.75 1.75 1.75 1.75 1.75  1.75 1.7];  % run the localizer for each diagonal in each block
task{1}{1}.segmax =     [1.75 1.75 1.75 1.75 1.75  1.75 1.71];  % this assumes a 2s TR... will need to be changed if we change our TR. Basically 16TRs/cycle 
task{1}{1}.segquant =   [0 0 0 0 0 0 0 0 0];
task{1}{1}.getResponse = [0 0 0 0 0 0 0 0 0];
task{1}{1}.synchToVol = [1 0 0 0 0 0 0];

task{1}{1}.fudgeLastVolume = 1;

n_repeats =5;

[contrast,location,repeats] = ndgrid(1,1:2,1:n_repeats);

%contrast =3 is blank trials. We wants on ~10% of total trials to be blank
%trials. Re-assign 4 out of 6 blank trials to be non-blank stim containing
%trials

task{1}{1}.numTrials = length(location(:)); % n*n_repeats
random_order = randperm(task{1}{1}.numTrials);

task{1}{1}.randVars.len_ = task{1}{1}.numTrials;
task{1}{1}.randVars.trialIndex = random_order;

stimulus.randVars.targetLocation = location(random_order); %one of the 2 positions
stimulus.randVars.contrast = contrast(random_order);

%we want a block design. 
% task{1}{1}.randVars.uniform.targetOrientation = 1:2;
% task{1}{1}.randVars.iti= iti(random_order)
% task{1}{1}.randVars.iti= task{1}{1}.randVars.iti.*1.5 % replace if TR changes

stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.LocationIndices=unique(location);


task{1}{1}.random = 1;
[task{1}{1}, myscreen] = initTask(task{1}{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen); %initStimulus('stimulus',myscreen,indContrast,diagonal);
stimulus = myInitStimulus(stimulus,myscreen,task{1},indContrast);

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
  fixStimulus.fixLineWidth = 1;
  fixStimulus.stimTime = 0.4;
  fixStimulus.responseTime = 1;
else
  % make cross bigger and task slower
  fixStimulus.diskSize = 30;
  fixStimulus.fixWidth =.4+1*easyFixTask;
  fixStimulus.fixLineWidth = 2+2*easyFixTask;
  fixStimulus.stimTime = 0.4+0.4*easyFixTask;
  fixStimulus.responseTime = 1+1*easyFixTask;
end
global fixStimulus
fixStimulus.pos = [0 0];
[task{2} myscreen] = fixStairInitTask(myscreen);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
mglSimulateRun(1.75,154,0)
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    % runs automatically the task, you only need to change: StartSegmentCallback,DrawStimulusCallback,responseCallback
    [task{1},myscreen,phaseNum] = updateTask(task{1},myscreen,phaseNum);
    [task{2} myscreen] = updateTask(task{2},myscreen,1);
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
mglClearScreen(stimulus.grayColor);
if (task.thistrial.thisseg == 1)
    stimulus.targetOrientation = randsample(0:25:180,1)
     if (stimulus.randVars.targetLocation(task.thistrial.trialIndex) == 1) 
       stimulus.tmp.targetLocation1= stimulus.eccentricity*[stimulus.locations{1}];
        stimulus.tmp.targetLocation2= stimulus.eccentricity*[stimulus.locations{3}];
     else 
         stimulus.tmp.targetLocation1 = stimulus.eccentricity*[stimulus.locations{2}];
         stimulus.tmp.targetLocation2= stimulus.eccentricity*[stimulus.locations{4}];
   
     end
else
     stimulus.targetOrientation= stimulus.targetOrientation + 180/7;
end

mglClearScreen(stimulus.grayColor);
setGammaTable(1);
end


%%
function [task, myscreen] = DrawStimulusCallback(task, myscreen)
global stimulus;

mglClearScreen(stimulus.grayColor);%###


%     if stimulus.EyeTrack, fixCheck; end
     drawGabor(stimulus.contrasts(stimulus.randVars.contrast(task.thistrial.trialIndex)),...
              stimulus.tmp.targetLocation1,...
              (stimulus.targetOrientation),1,...
              task.thistrial.trialIndex);
          
       drawGabor(stimulus.contrasts(stimulus.randVars.contrast(task.thistrial.trialIndex)),...
              stimulus.tmp.targetLocation2,...
              (stimulus.targetOrientation),1,...
              task.thistrial.trialIndex);
          % the above line of code adds or subtracts the tilt from the base
          % orientation


end


