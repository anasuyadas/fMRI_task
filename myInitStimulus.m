function stimulus = myInitStimulus(stimulus,myscreen,task,contrast)
global MGL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;


% xpxpcm = myscreen.screenWidth/myscreen.displaySize(1);
% ypxpcm = myscreen.screenHeight/myscreen.displaySize(2);
%
% xpxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*xpxpcm);
% ypxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*ypxpcm);
% 
% centerpix = [myscreen.screenWidth/2,myscreen.screenHeight/2];

stimulus.contrasts =contrast;

stimulus.frameThick = .08;
stimulus.reservedColors = [0 0 0; 1 1 1; 0 .6 0];


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
stimulus.width = 10;%4;%stimulus.gaussSdx*7;             % in deg
stimulus.height = 10;%4;%stimulus.gaussSdy*7;            % in deg
stimulus.gaussSdx = stimulus.width/7;                % in deg
stimulus.gaussSdy = stimulus.height/7;               % in deg
stimulus.sizedg = 20;%should be reset to 3degs

stimulus.rotation = [-1 1]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.init = 1;

stimulus.sf = 4;                % in cpd
stimulus.orientation = 90;       % in deg
stimulus.phase = 0;             % in deg
stimulus.eccentricity = 4.6;    % in deg

stimulus.locations = {[-cosd(45),sind(45)];[cosd(45), sind(45)];[cosd(45), -sind(45)];[-cosd(45), -sind(45)]};

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
stimulus.FCwidth = 1;
stimulus.FClinewidth = 3;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate befire trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=1;
stimulus.edgeDist=0;%presents stimuli after this duration when fixation detected



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set parameters -> in degrees
respCue.length = 0.5;
respCue.startEcc = 0.5;

stimulus.respCue.width =3;
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
stimulus.preCue.type = 0; %always neutral for staircase
respCue.width =3;

preCue.length = .8;
preCue.height = .2;
preCue.distToStim = .4 + (stimulus.height/2); %distance from center

neutCue.distToFixation = .3 + stimulus.FCwidth; %distance from center

stimulus.preCue.width =3;

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
                                  preCue.length/4, -neutCue.distToFixation];                           
end