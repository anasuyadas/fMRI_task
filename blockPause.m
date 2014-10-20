function blockPause(block,numBlocks)
mglTextSet('Helvetica',24,[1 1 1],0,0,0,1,0,0,0);
thisText = mglText(sprintf('You have completed block %d out of %d',block,numBlocks));

myscreen = mglBltTexture(thisText,[0 0],'left','top');

waitSecs(15)

end