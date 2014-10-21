function fixCheck(myscreen,task)
    global stimulus
  
    ep = myscreen.eyetracker.eyepos;


%     if task.thistrial.thisseg == 1;
%         if (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist && ~stimulu
% s.FixationStarted
%             stimulus.FixationStart=mglGetSecs;
%             stimulus.FixationStarted=1;
%         elseif (sqrt(ep(end,1)^2+ep(end,2)^2))>=stimulus.TrialStartFixDist && stimulus.FixationStarted
%             stimulus.FixationStarted=0;
%         elseif (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist
%             stimulus.FixationDur=mglGetSecs(stimulus.FixationStart);
%             if stimulus.FixationDur >=stimulus.TrialStartFixDur
%             end
%         end
%         
%     else
        if (sqrt(ep(end,1)^2+ep(end,2)^2)) > stimulus.TrialStartFixDist && ~stimulus.FixationBreakCurrent
            stimulus.FixationBreak(task.trialnum)=1;
            stimulus.FixationBreakCurrent = 1;
            stimulus.updateCurrent = 0;
            
            if ~stimulus.firstFixBreak
                stimulus.fixBreakTRACKindex(1) = task.thistrial.trialIndex;
                stimulus.fixBreakTRACKseg(1) = task.thistrial.thisseg;
                stimulus.firstFixBreak = 1;
                stimulus.numFixBreaks = 1;
                stimulus.fixationBreakTrialVect(stimulus.numFixBreaks) = stimulus.trialAttemptNum;
            else
                stimulus.fixBreakTRACKindex(1,length(stimulus.fixBreakTRACKindex)+1) = task.thistrial.trialIndex;
                stimulus.fixBreakTRACKseg(1,length(stimulus.fixBreakTRACKseg)+1) = task.thistrial.thisseg;
                stimulus.numFixBreaks = stimulus.numFixBreaks+1;
                stimulus.fixationBreakTrialVect(stimulus.numFixBreaks) = stimulus.trialAttemptNum;
            end
            
        end
end