function stimulus = myInitStimulusCONTstair(stimulus,myscreen,task,contrast)
global MGL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% let's get the linearized gamma table
stimulus.linearizedGammaTable = mglGetGammaTable;
stimulus.linearizedGammaTable.redTable(1:3) = 0; % this is just to provisionally deal with what appears to be some bug: the first value in each of these gamma tables is a NaN
stimulus.linearizedGammaTable.greenTable(1:3) = 0;
stimulus.linearizedGammaTable.blueTable(1:3) = 0;


stimulus.xpxpcm = myscreen.screenWidth/myscreen.displaySize(1);
stimulus.ypxpcm = myscreen.screenHeight/myscreen.displaySize(2);

stimulus.xpxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*stimulus.xpxpcm);
stimulus.ypxpdeg = ceil(tan(2*pi/360)*myscreen.displayDistance*stimulus.ypxpcm);

centerpix = [myscreen.screenWidth/2,myscreen.screenHeight/2];

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
stimulus.white = [255 255 255];
stimulus.green = [0 160 0];
stimulus.blue = [0 0 160];
stimulus.greencorrect = [0 200 20];
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
stimulus.sizedg = 3;%should be reset to 3degs

stimulus.rotation = [1,-1]; % this is the tilt orientation of the gabor stimulus from vertical in Degrees
stimulus.init = 1;

stimulus.sf = 4;                % in cpd
stimulus.orientation = 0;       % in deg
stimulus.phase = 180.*rand(1,task{1}.numTrials); % in deg
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
%ONLY for first trial
gratingMatrix = mglMakeGrating(stimulus.width,stimulus.height,stimulus.sf,90,stimulus.phase(1));


 res = mkR([size(gratingMatrix,1) size(gratingMatrix,2)]);
 
 [Xtbl,Ytbl] = rcosFn(size(gratingMatrix,1)/5,stimulus.sizedg*stimulus.xpxpdeg/2, [1 0]);%(stimulus.sizedg)/2, [1 0]); %1st argument is eidth pixels => MAKE INTO VARIABLE
 grating(:,:,4) = 255*pointOp(res, Ytbl, Xtbl(1), Xtbl(2)-Xtbl(1), 0);
 
 

% stimulus.texture
grating(:,:,1) = stimulus.midGratingColors+gratingMatrix*(127*stimulus.stair.threshold);
grating(:,:,2) = grating(:,:,1);
grating(:,:,3) = grating(:,:,1);
stimulus.tex{1} = mglCreateTexture(grating);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimulus.FCwidth = .6;
stimulus.FClinewidth = 3;
stimulus.TrialStartFixDist=2; %2 degree radius in which to fixate before trial starts
stimulus.TrialStartFixDur=.25;
stimulus.cornerDist=1;
stimulus.edgeDist=0;%presents stimuli after this duration when fixation detected

stimulus.FCloc = cell(2,1);

stimulus.FCloc{1} = [  cosd(135)*(stimulus.FCwidth),  sind(135)*(stimulus.FCwidth);...
                      -cosd(135)*(stimulus.FCwidth), -sind(135)*(stimulus.FCwidth)];
                         
stimulus.FCloc{2} = [ cosd(45)*(stimulus.FCwidth),  sind(45)*(stimulus.FCwidth);...
                     -cosd(45)*(stimulus.FCwidth), -sind(45)*(stimulus.FCwidth)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set parameters -> in degrees
respCue.length = stimulus.FCwidth;
respCue.startEcc = 0;

stimulus.respCue.width =stimulus.FClinewidth;
stimulus.respcueLocation = cell(4,1);

% for loc = 1:4
%         stimulus.respcueLocation{loc} = [0, 0;...
%                                          stimulus.locations{loc}(1)*(respCue.length),...
%                                          stimulus.locations{loc}(2)*(respCue.length)];
% end

stimulus.respcueLocation{1} = [stimulus.FCloc{1}(1,:);0,0];
stimulus.respcueLocation{2} = [stimulus.FCloc{2}(1,:);0,0];
stimulus.respcueLocation{3} = [0,0;stimulus.FCloc{1}(2,:)];
stimulus.respcueLocation{4} = [0,0;stimulus.FCloc{2}(2,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRE CUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set parameters -> in degrees
% stimulus.preCue.type = 1; %should be always neutral for staircase
respCue.width =3;

preCue.length = 1.2;
preCue.height = .2;
preCue.distToStim = .2 + (stimulus.height/2); %distance from center

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