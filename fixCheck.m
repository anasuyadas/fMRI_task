function fixCheck(myscreen,task)
    global stimulus
  
    ep = myscreen.eyetracker.eyepos;
    
%     if task.thistrial.thisseg == 1;
%         if (sqrt(ep(end,1)^2+ep(end,2)^2))<=stimulus.TrialStartFixDist && ~stimulus.FixationStarted
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
        if (sqrt(ep(end,1)^2+ep(end,2)^2))>stimulus.TrialStartFixDist
            stimulus.FixationBreak(task.trialnum)=1;
        end
%     end
end