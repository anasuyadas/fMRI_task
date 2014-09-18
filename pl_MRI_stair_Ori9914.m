function myscreen = pl_MRI_stair_Ori(observer,varargin)

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
% end

eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye'}));

if ieNotDefined('indContrast'),indContrast = [.8];end % initialize some default contrast vals
if ieNotDefined('diagonal'),diagonal = 1;end % default diagonal. Can be zero or 1. diagonal 1: upper right+ lower left; diagonal 2: lower right + upper left. THIS NEEDSS TO BE DOUBLE CHECKED
if ieNotDefined('indTilt'),indTilt = 5;end % default tilt
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
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[1,4],1:n_repeats);
else 
    [contrast,ori,location,trialNum] = ndgrid(1,1:2,[2,3],1:n_repeats);
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


task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback,@staicaseCallback);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,indContrast,diagonal);

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
global stimulus;

if (task.thistrial.thisseg == 9) % ITI
    stimulus.trialend = stimulus.trialend + 1;
elseif (task.thistrial.thisseg == 1) % fixation
    iti = task.thistrial.iti;
    task.thistrial.seglen =[0.1 .06 .04 .1 .3 .3 .64 .03 iti];
    %need to make sure that there are only two locations per run
    stimulus.tmp.targetLocation  = stimulus.eccentricity*[stimulus.locations{task.thistrial.targetLocation}];
    
    stimulus.FixationStarted=0;
    %response cue
    stimulus.tmp.respcueLocation=stimulus.respcueLocation{task.thistrial.targetLocation}; %if polygon
    stimulus.tmp.respcueLocation=task.thistrial.targetLocation; %if central x
    stimulus.tmp.WedgeStart=stimulus.CueWedges(task.thistrial.targetLocation);
    
    %just neutral cues - no exo cues
    for i=1:2
        stimulus.tmp.NeutralcueLocation{i}=stimulus.NeutralcueLocation{i};
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
    drawFixation;
    
elseif (task.thistrial.thisseg == 1) % Initial Fixation
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
elseif (task.thistrial.thisseg == 2) % Pre Cue
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
    drawPreCue;
    
elseif (task.thistrial.thisseg == 3) % ISI 1
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 4) % Stimulus
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
    drawGabor(stimulus.contrasts(task.thistrial.contrast),...
              stimulus.tmp.targetLocation,...
              stimulus.rotation(task.thistrial.targetOrientation)*stimulus.stair.threshold,1);
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation;
    if stimulus.EyeTrack, fixCheck; end
    drawRespCue(stimulus.tmp.targetLocation);

elseif (task.thistrial.thisseg == 7) % Resp Window
    drawFixation;
    
elseif (task.thistrial.thisseg == 8) % Feedback
    drawFixation;

end
    

end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,contrast,diagonal)
global MGL;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;


xpxpcm = myscreen.screenWidth/myscreen.displaySize(1);
ypxpcm = myscreen.screenHeight/myscreen.displaySize(2);

xpxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*xpxpcm);
ypxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*ypxpcm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gabors
stimulus.width = 4;%stimulus.gaussSdx*7;             % in deg
stimulus.height = 4;%stimulus.gaussSdy*7;            % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg

stimulus.rotation = [1 -1]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.init = 1;

stimulus.sf = 4;                % in cpd
stimulus.orientation = 0;      % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity = 4.6;    % in deg

stimulus.locations = {[-cosd(45),sind(45)],[cosd(45), sind(45)],[cosd(45), -sind(45)],[-cosd(45), -sind(45)]};

for loc = 1:4;
    stimulus.locationsPix{loc,1} = [stimulus.locations{loc}(1)*stimulus.eccentricity*xpxpdeg,...
                                   stimulus.locations{loc}(2)*stimulus.eccentricity*ypxpdeg];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.FCwidth = 0.3;
stimulus.FClinewidth = 3;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate befire trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=1;
stimulus.edgeDist=0;%presents stimuli after this duration when fixation detected



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set parameters -> in degrees
respCue.length = 1;
respCue.width =.1;
respCue.startEcc = .3;

stimulus.respcueLocation{1}= respCue.startEcc+...
                             [];
stimulus.respcueLocation{2}=[-0.35 -0.55;-0.55 -0.35;-0.7 -0.9];
stimulus.respcueLocation{3}=[-0.35 0.55;-0.55 0.35;-0.7 0.9];
stimulus.respcueLocation{4}=[0.35 0.55;0.55 0.35;0.7 0.9];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters -> in degrees
preCue.length = .8;
preCue.height = .2;
preCue.distToStim = .4;

neutCue.distToFixation = .3 + stimulus.FCwidth; %distance from center


%convert to pixels
preCue.x1 = (-preCue.length/2)*xpxpdeg;
preCue.x2 = (preCue.length/2)*xpxpdeg;

preCue.y1 = (stimulus.height+preCue.distToStim)*ypxpdeg;
preCue.y2 = (stimulus.height+preCue.distToStim+preCue.height)*ypxpdeg;

% make matrices in stimulus structure
stimulus.preCueExgLocation{1} = (stimulus.locationsPix{1}+...
                                 [preCue.x1, preCue.y1;...
                                  preCue.x2, preCue.y1;...
                                  preCue.x2, preCue.y2;...
                                  preCue.x1, preCue.y2]);

stimulus.preCueExgLocation{2} = (stimulus.locationsPix{2}+...
                                 [preCue.x1, preCue.y1;...
                                  preCue.x2, preCue.y1;...
                                  preCue.x2, preCue.y2;...
                                  preCue.x1, preCue.y2]);

stimulus.preCueExgLocation{3} = (stimulus.locationsPix{3}+...
                                 [preCue.x1, -preCue.y1;...
                                  preCue.x2, -preCue.y1;...
                                  preCue.x2, -preCue.y2;...
                                  preCue.x1, -preCue.y2]);

stimulus.preCueExgLocation{4} = (stimulus.locationsPix{4}+...
                                 [preCue.x1, -preCue.y1;...
                                  preCue.x2, -preCue.y1;...
                                  preCue.x2, -preCue.y2;...
                                  preCue.x1, -preCue.y2]);


stimulus.preCueNeutLocation{1} = [-preCue.length/4, neutCue.distToFixation;...
                                   preCue.length/4, neutCue.distToFixation;...
                                   preCue.length/4, neutCue.distToFixation+preCue.height;...
                                  -preCue.length/4, neutCue.distToFixation+preCue.height];

stimulus.preCueNeutLocation{1} = [-preCue.length/4, -neutCue.distToFixation;...
                                   preCue.length/4, -neutCue.distToFixation;...
                                   preCue.length/4, -(neutCue.distToFixation+preCue.height);...
                                  -preCue.length/4, -(neutCue.distToFixation+preCue.height)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAIRCASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stair.upRule = 1;
stair.downRule = 3;
stair.startThresh = log10(5);
stair.stepSize = .5;
stair.minStepSize = .05;
stair.halfRule = 'levitt';

stimulus.stair = upDownStaircase(stair.upRule,stair.downRule,stair.startThresh,[stair.stepSize, stair.minStepSize],stair.halfRule);


end

%% draw fixation
function drawFixation
    global stimulus
    global task
    
    if ~task.thistrial.gotResponse
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white); %if there is no response & it is not the response window or feedback window, just present white fixation
    
    elseif (task.thistrial.thisseg == 7 && ~task.thistrial.gotResponse); %No fixation before response
    
    elseif ((task.thistrial.thisseg == 7 && task.thistrial.gotResponse) || (task.thistrial.thisseg == 8 && task.thistrial.gotResponse)) %Once there is a response, present fixation of a different color
        if stimulus.tmp.response == 1, mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
        elseif stimulus.tmp.response == 0, mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.stimulus.redincorrect); end
    
    elseif (task.thistrial.thisseg == 8 && ~task.thistrial.gotResponse) %If there is no response, present orange fixation
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.orangenoanswer);
    end
end
%% fixCheck
function fixCheck
    global stimulus
    global task
    
    ep = myscreen.eyetracker.eyepos;
    
    if task.thistrial.thisseg == 1;
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
        if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
            stimulus.FixationBreak(stimulus.trialnum)=1;
        end
    end
end

%%
function drawPreCue
    global stimulus
    global task
    
    if stimulus.preCue == 0
        
    elseif stimulus.preCue == 1
        
    end

end
%%
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
end

%%

function drawRespCue
    global task
    global stimulus
    task.thistrial.targetLocation
    mglLines2(stimulus.centralx{stimulus.tmp.respcueLocation}(1), stimulus.centralx{stimulus.tmp.respcueLocation}(2), stimulus.centralx{stimulus.tmp.respcueLocation}(3), stimulus.centralx{stimulus.tmp.respcueLocation}(4),8,stimulus.black);
    
end
