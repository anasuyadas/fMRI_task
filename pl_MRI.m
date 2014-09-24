function myscreen = pl_MRI(observer,varargin)

%%% Pilot with only valid and invalid conditions
% The raised cosine reqires matlabPyrTools
% http://www.cns.nyu.edu/~lcv/software.php
% psuedorandomize across all combination of stims+ITI


global stimulus;
global MGL;

mglVisualAngleCoordinates(57,[42 26]); %distance from screen, height & width of monitor
% check arguments
% if ~any(nargin == 3)
%     help transientAttention
%     return
% % end
% 
% eval(evalargs(varargin,0,0,{'indContrast','diagonal','IndTilt','Eye'}));

if ieNotDefined('indContrast'),indContrast = [0 0.8 1];end % initialize some default contrast vals
if ieNotDefined('indOri'),indOri = 3;end % initialize some default contrast vals
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
myscreen = initScreen('disp2');
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



n_repeats = 3;%  trials per block n= 36; 3contrast*3ITIs*2location 
% Number of volumes = (n)+(n/3*2)+(n/3*3)+(n/3*4).
%n_repeats will have to be adjusted depending on our TR to keep block
%length approximately ~5minutes
if diagonal == 1  
    [contrast, iti, ori,location,repeats] = ndgrid(1:3,1:3,1:2,1:2,1:n_repeats);
else 
    [contrast, iti, ori,location,repeats] = ndgrid(1:3,1:3,1:2,3:4,1:n_repeats);
end
%contrast =3 is blank trials. We wants on ~10% of total trials to be blank
%trials. Re-assign 4 out of 6 blank trials to be non-blank stim containing
%trials
contrast(3,:,[1:2],1)=1;
contrast(3,:,[1:2],2)=2; 


task{1}.numTrials = length(location(:)); % n*n_repeats
random_order = randperm(task{1}.numTrials);
 
task{1}.randVars.targetLocation = location(random_order); %one of the 2 positions
task{1}.randVars.len_ = task{1}.numTrials;
task{1}.randVars.contrast = contrast(random_order);
task{1}.randVars.targetOrientation = ori(random_order);
task{1}.randVars.iti= iti(random_order)
task{1}.randVars.iti= task{1}.randVars.iti.*1.5 % replace if TR changes

stimulus.trialend = 0;
stimulus.trialnum=1;
stimulus.FixationBreak=zeros(1,length(location(:)));
stimulus.LocationIndices=unique(location);


task{1}.random = 1;
[task{1}, myscreen] = initTask(task{1},myscreen,@StartSegmentCallback,@DrawStimulusCallback,@responseCallback);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myscreen = initStimulus('stimulus',myscreen); %initStimulus('stimulus',myscreen,indContrast,diagonal);
stimulus = myInitStimulus(stimulus,myscreen,task,indContrast,indOri);

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
    drawGabor(stimulus.contrasts(task.thistrial.contrast),...
              stimulus.tmp.targetLocation,...
              stimulus.rotation(task.thistrial.targetOrientation) ,1);
    
elseif (task.thistrial.thisseg == 5) % ISI 2
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
    
elseif (task.thistrial.thisseg == 6) % Resp Cue
    drawFixation(task);
    if stimulus.EyeTrack, fixCheck; end
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
function [task,stimulus] = responseCallback(task,stimulus)
%
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
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,contrast,orientation)
global MGL;

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;


stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];

stimulus.contrasts =contrast;

stimulus.nReservedColors=size(stimulus.reservedColors,1);
stimulus.nGratingColors = 256-(2*floor(stimulus.nReservedColors/2)+1);
stimulus.minGratingColors = 2*floor(stimulus.nReservedColors/2)+1;
stimulus.midGratingColors = stimulus.minGratingColors+floor(stimulus.nGratingColors/2);
stimulus.maxGratingColors = 255;
stimulus.deltaGratingColors = floor(stimulus.nGratingColors/2);

stimulus.nDisplayContrasts = stimulus.deltaGratingColors;


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
stimulus.grayColor = stimulus.background;

% audio
stimulus.noanswer = find(strcmp(MGL.soundNames,'Ping'));
stimulus.CorrectSound = find(strcmp(MGL.soundNames,'Submarine'));
stimulus.IncorrectSound = find(strcmp(MGL.soundNames,'Basso'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gabors
stimulus.width = 4;%stimulus.gaussSdx*7;             % in deg
stimulus.height = 4;%stimulus.gaussSdy*7;            % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg


stimulus.rotation = [1 -1]*orientation; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.init = 1;

stimulus.sf = 4;                % in cpd
stimulus.orientation = 0;       % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity = 4.6;    % in deg

stimulus.locations = {[-cosd(45),sind(45)];[cosd(45), -sind(45)];[cosd(45), sind(45)];[-cosd(45), -sind(45)]};

stimulus.locationsEcc = cell(4,1);
for loc = 1:4;
    stimulus.locationsEcc{loc} = [stimulus.locations{loc}(1)*stimulus.eccentricity,...
                                   stimulus.locations{loc}(2)*stimulus.eccentricity];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make stim texture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for thisSF = 1:length(stimulus.sf)      %only one spatial frequency
    gratingMatrix{thisSF} = mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf(thisSF),stimulus.orientation,stimulus.phase);
end

res = mkR([size(gratingMatrix{1},1) size(gratingMatrix{1},2)]);
stimulus.sizedg = 3;
[Xtbl,Ytbl] = rcosFn(30, (stimulus.sizedg)/2, [1 0]); %2 = sharp transition (edge effect?) / 50 = radius of the circle in pixel
grating(:,:,4) = 255*pointOp(res, Ytbl, Xtbl(1), Xtbl(2)-Xtbl(1), 0);



disppercent(-inf,'Calculating gabors');
for thisSF = 1:length(stimulus.sf)
    for thisContrast = 0:stimulus.deltaGratingColors
        % stimulus.texture
        grating(:,:,1) = stimulus.midGratingColors+gratingMatrix{thisSF}*thisContrast;
        grating(:,:,2) = grating(:,:,1);
        grating(:,:,3) = grating(:,:,1);
        stimulus.tex{thisSF}(thisContrast+1) = mglCreateTexture(grating);
        disppercent(thisContrast/stimulus.deltaGratingColors);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.FCwidth = 0.5;
stimulus.FClinewidth = 3;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate befire trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=1;
stimulus.edgeDist=0;%presents stimuli after this duration when fixation detected



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set parameters -> in degrees
respCue.length = .5;
respCue.width =3;
respCue.startEcc = 0.5;

stimulus.respcueLocation = cell(4,1);
for loc = 1:4
        stimulus.respcueLocation{loc} = [stimulus.locations{loc}(1)*respCue.startEcc,...
                                        stimulus.locations{loc}(2)*respCue.startEcc;...
                        
                                        stimulus.locations{loc}(1)*(respCue.startEcc+respCue.length),...
                                        stimulus.locations{loc}(2)*(respCue.startEcc+respCue.length)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters -> in degrees
stimulus.preCue = 0; %always neutral for staircase

preCue.length = .5;
preCue.height = .2;
preCue.distToStim = .5 + (stimulus.height/2); %distance from center

neutCue.distToFixation = 0.5 + stimulus.FCwidth; %distance from center



% make matrices in stimulus structure
stimulus.preCueExgLocation = cell(4,1);
for loc = 1:4
    if stimulus.locationsEcc{loc}(2) > 0 
        stimulus.preCueExgLocation{loc} = [(stimulus.locationsEcc{loc}(1))-preCue.length/2,...
                                           (stimulus.locationsEcc{loc}(2))+preCue.distToStim;...
                               
                                           (stimulus.locationsEcc{loc}(1))+preCue.length/2,...
                                           (stimulus.locationsEcc{loc}(2))+preCue.distToStim];
    else
         stimulus.preCueExgLocation{loc} =[(stimulus.locationsEcc{loc}(1))-preCue.length/2,...
                                           (stimulus.locationsEcc{loc}(2))-preCue.distToStim;...
                        
                                           (stimulus.locationsEcc{loc}(1))+preCue.length/2,...
                                           (stimulus.locationsEcc{loc}(2))-preCue.distToStim];
    end
end


stimulus.preCueNeutLocation{1} = [-preCue.length/4, neutCue.distToFixation;...
                                   preCue.length/4, neutCue.distToFixation];

stimulus.preCueNeutLocation{2} = [-preCue.length/4, -neutCue.distToFixation;...
                                  -preCue.length/4, -neutCue.distToFixation];

                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAIRCASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stair.upRule = 1;
stair.downRule = 3;
stair.startThresh = log10(5);
stair.stepSize = .8;
stair.minStepSize = .02;
stair.halfRule = 'levitt';

stimulus.stair = upDownStaircase(stair.upRule,stair.downRule,stair.startThresh,[stair.stepSize, stair.minStepSize],stair.halfRule);


end

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
end

%% draw fixation
function drawFixation(task)
    global stimulus

    
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
function drawPreCue(loc)
    global stimulus
    
    if stimulus.preCue == 0
        mglLines2(stimulus.preCueNeutLocation{loc}(1),stimulus.preCueNeutLocation{loc}(3),...
                  stimulus.preCueNeutLocation{loc}(2),stimulus.preCueNeutLocation{loc}(4),stimulus.FCwidth,stimulus.white);
    elseif stimulus.preCue == 1
        mglLines2(stimulus.preCueExgLocation{loc}(1),stimulus.preCueExgLocation{loc}(3),...
                  stimulus.preCueExgLocation{loc}(2),stimulus.preCueExgLocation{loc}(4),stimulus.FCwidth,stimulus.white);
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

function drawRespCue(loc)
    global stimulus
    
    mglLines2(stimulus.respcueLocation{loc}(1), stimulus.respcueLocation{loc}(3),...
              stimulus.respcueLocation{loc}(2), stimulus.respcueLocation{loc}(4),stimulus.FCwidth,stimulus.black);
    
end
