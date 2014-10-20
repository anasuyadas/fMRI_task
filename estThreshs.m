

oriThreshEst1 = mean(stimulus.stair.strength(stimulus.stair.reversals((length(stimulus.stair.reversals)-6):end)));
%%
oriThreshEst2 = mean(stimulus.stair.strength(stimulus.stair.reversals((length(stimulus.stair.reversals)-6):end)));

oriThreshEst = mean([oriThreshEst1,oriThreshEst2]);




%%

contThreshEst1 = mean(stimulus.stair.strength(stimulus.stair.reversals((length(stimulus.stair.reversals)-6):end)));
%%
contThreshEst2 = mean(stimulus.stair.strength(stimulus.stair.reversals((length(stimulus.stair.reversals)-6):end)));

contThreshEst = mean([contThreshEst1,contThreshEst2]);