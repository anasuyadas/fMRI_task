%% draw fixation
function drawFixation(task)
    global stimulus
    
    mglTextSet('Arial',24,[1 1 1],0,0,0,1,0,0,0);
    thisText = mglText('please fixate');
    
    if stimulus.FixationBreakCurrent && (1 <= task.thistrial.thisseg < 9)  %(stimulus.FixationBreak(task.trialnum))  %If there is a fixn break, present 
       mglBltTexture(thisText,[0 0],'left','top');
    
    elseif ~task.thistrial.gotResponse
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.white); %if there is no response & it is not the response window or feedback window, just present white fixation
    
    elseif (task.thistrial.thisseg == 7 && ~task.thistrial.gotResponse); %No fixation before response
    
    elseif ((task.thistrial.thisseg == 7 && task.thistrial.gotResponse) || (task.thistrial.thisseg == 8 && task.thistrial.gotResponse)) %Once there is a response, present fixation of a different color
        if stimulus.tmp.response == 1, mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.greencorrect);
        elseif stimulus.tmp.response == 0, mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.redincorrect); end
    
    elseif (task.thistrial.thisseg == 8 && ~task.thistrial.gotResponse)  %If there is no response, present orange fixation
        mglFixationCross(stimulus.FCwidth,stimulus.FClinewidth,stimulus.orangenoanswer);
    
    end
end