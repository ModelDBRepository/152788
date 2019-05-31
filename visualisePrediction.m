% This script visualises the growth of the neuron, used for the prediction
% simulations.

close all, format compact

dataPath = 'output/';

predictGC = unique([]);

dataFile = 'GC-slavery-predict';

for i = 1:length(predictGC)
  dataFile = sprintf('%s-%d', dataFile, predictGC(i));
end

dataFile = strcat(dataFile,'.txt');

morphFile = 'input/Ramaker-neuronMorph.txt';
gcFiles = { 'input/Ramaker-GC1-growthspeed.txt', ...
	    'input/Ramaker-GC2-growthspeed.txt', ...
	    'input/Ramaker-GC3-growthspeed.txt' };


if(~exist('data') ...
   | length(predictGC) ~= length(oldPredictGC) ...
   | nnz(predictGC-oldPredictGC) > 0)

  % Read new data
  data = readData([dataPath dataFile]);

  % Saved for future reference
  oldPredictGC = predictGC;

  % Read in morphology and growth cones

  tmp = load(morphFile);
  morphData.x = tmp(:,1);
  morphData.y = tmp(:,2);

  gcData = struct('time',[],'x',[],'y',[]);

  for i = 1:length(gcFiles)
    tmp = load(gcFiles{i});
    tmpData.time = tmp(:,1);
    tmpData.x = tmp(:,2);
    tmpData.y = tmp(:,3);

    gcData(i) = tmpData;
  end

end


%%%% Now we got to plot it...

T = 100000;

gcBranchX = {};
gcBranchY = {};

% Loop through the growth cones
for i = 1:length(gcData)
  % Find relevant time points
  idx = find(gcData(i).time <= T);

  if(isempty(idx))
    continue
  end

  % If a newer time point is closer to the point two steps back than the
  % old point, then remove it.

  branchMorphX = gcData(i).x(idx(1:2));
  branchMorphY = gcData(i).y(idx(1:2));

  for j = 3:length(idx)

    % Only process point if it is different from previous point
    if(branchMorphX(end) ~= gcData(i).x(idx(j)) ...
       | branchMorphY(end) ~= gcData(i).y(idx(j)))

      oldDir = [branchMorphX(end)-branchMorphX(end-1), ...
     	        branchMorphY(end)-branchMorphY(end-1)];

      newDir = [gcData(i).x(idx(j))-branchMorphX(end-1), ...
	        gcData(i).y(idx(j))-branchMorphY(end-1)];

      while(sum(oldDir.*newDir) < 0 & length(branchMorphX) > 2)
        % Remove previous point
        branchMorphX(end) = [];
        branchMorphY(end) = [];

        oldDir = [branchMorphX(end)-branchMorphX(end-1), ...
        	  branchMorphY(end)-branchMorphY(end-1)];

        newDir = [gcData(i).x(idx(j))-branchMorphX(end-1), ...
	          gcData(i).y(idx(j))-branchMorphY(end-1)];

      end

      branchMorphX(end+1) = gcData(i).x(idx(j));
      branchMorphY(end+1) = gcData(i).y(idx(j));    

    end

  end

  gcBranchX{i} = branchMorphX;
  gcBranchY{i} = branchMorphY;

end

plot(morphData.x,morphData.y,'k-')
hold on
for i=1:length(gcBranchX)
%  gcDist = sqrt((morphData.x-gcBranchX{i}(1))^2 ...
%		+(morphData.y-gcBranchY{i}(1))^2);
%
%  linkIdx = find(gcDist == min(gcDist),1);


  plot(gcBranchX{i},gcBranchY{i},'b-')
end
