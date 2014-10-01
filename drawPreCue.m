function drawPreCue(loc)
    global stimulus
    
    if stimulus.preCue == 0
        mglLines2(stimulus.preCueNeutLocation{1}(1),stimulus.preCueNeutLocation{1}(3),...
                  stimulus.preCueNeutLocation{1}(2),stimulus.preCueNeutLocation{1}(4),1,stimulus.white);
        mglLines2(stimulus.preCueNeutLocation{2}(1),stimulus.preCueNeutLocation{2}(3),...
                  stimulus.preCueNeutLocation{2}(2),stimulus.preCueNeutLocation{2}(4),1,stimulus.white);
    elseif stimulus.preCue == 1
        mglLines2(stimulus.preCueExgLocation{loc}(1),stimulus.preCueExgLocation{loc}(3),...
                  stimulus.preCueExgLocation{loc}(2),stimulus.preCueExgLocation{loc}(4),1,stimulus.white);
    end

end