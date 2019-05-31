function plotConcentrationGradient(data, t)

  figure

  uTime = unique(data.time);
  tmpDiff = abs(uTime - t);
  tClosest = uTime(find(tmpDiff == min(tmpDiff),1)); 

  index = find(data.time == tClosest);

  %% Locate growth cones
  gcID = setdiff(data.ID(index),data.parentID(index));

  for j = 1:length(gcID)
    gcIdx = index(find(data.ID(index) == gcID(j)));
    dist = data.dist(gcIdx);
    conc = data.tubulinConc(gcIdx);

    curIdx = gcIdx;

    compCtr = 0;

    while(data.parentID(curIdx) ~= -1)    
      parentIdx = index(find(data.ID(index) == data.parentID(curIdx)));

      dist = [dist; data.dist(parentIdx)];
      conc = [conc; data.tubulinConc(parentIdx)];

      curIdx = parentIdx;

      compCtr = compCtr + 1;

      if(compCtr > 1e6)
        disp('WARNING: Did we get stuck in an infinite loop?')
        beep
        return
      end

    end

      

    plot(dist*1e6,conc,'k-')
    hold on

  end

  title(sprintf('Tubulin concentration at time %d hours', tClosest/3600))
  xlabel('Distance to soma (\mum)')
  ylabel('Concentration (mM)')

end
