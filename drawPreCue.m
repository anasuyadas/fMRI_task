function drawPreCue(loc)
    global stimulus
    
    if stimulus.preCue.type == 0
        mglLines2(stimulus.preCueNeutLocation{1}(1),stimulus.preCueNeutLocation{1}(3),...
                  stimulus.preCueNeutLocation{1}(2),stimulus.preCueNeutLocation{1}(4),stimulus.preCue.width,stimulus.white);
        mglLines2(stimulus.preCueNeutLocation{2}(1),stimulus.preCueNeutLocation{2}(3),...
                  stimulus.preCueNeutLocation{2}(2),stimulus.preCueNeutLocation{2}(4),stimulus.preCue.width,stimulus.white);
    elseif stimulus.preCue.type == 1
        mglLines2(stimulus.preCueExgLocation{loc}(1),stimulus.preCueExgLocation{loc}(3),...
                  stimulus.preCueExgLocation{loc}(2),stimulus.preCueExgLocation{loc}(4),stimulus.preCue.width,stimulus.white);
    end

end