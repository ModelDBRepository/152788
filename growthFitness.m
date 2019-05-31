% This function calculates how similar the growth of the simulated
% neuron is to that of the experimental neuron.

function errorValue = growthFitness(simFile,expFiles)

  function expData = readExpData(expFiles)

    for i = 1:length(expFiles)

      tmp = load(expFiles{i});
      expData(i).time = tmp(:,1);
      expData(i).v = tmp(:,2);
      expData(i).dist = tmp(:,3);
 
    end

  end

  if(~exist('expFiles'))
    expFiles = { 'input/Ramaker/Ramaker-GC1-growthspeed.txt', ...
		 'input/Ramaker/Ramaker-GC2-growthspeed.txt', ...
		 'input/Ramaker/Ramaker-GC3-growthspeed.txt' };
  end

  simData = readData(simFile);  

  experimentData = readExpData(expFiles);

  % Find the growth cones
  gcID = setdiff(simData.ID,simData.parentID);

  % We do not know which one is which though... ack.
  % For the Ramaker data we know the following:
  % GC1 starts earliest.
  % GC3 is lowest at 8e4 seconds

  for i = 1:length(gcID)
    tFirst(i) = simData.time(find(simData.ID == gcID(i),1));
  end


  [foo idx] = sort(tFirst);
  gcID = gcID(idx);


  idx2 = find(simData.ID == gcID(2) & simData.time > 8e4,1);
  idx3 = find(simData.ID == gcID(3) & simData.time > 8e4,1);

  if(simData.dist(idx2) < simData.dist(idx3))
    gcID([2 3]) = gcID([3 2]);
  end

gcID

  % How do I match times best...

  gcError = NaN*gcID;

  for iG = 1:length(gcID)
    % Find all values for respective growth cone
    idx = find(simData.ID == gcID(iG));

    tIdx = idx(find(min(experimentData(iG).time) <= simData.time(idx) ...
		    & simData.time(idx) < max(experimentData(iG).time)));
    % Simulated distances at experimental time points
    % Experimental data at simulated time points
    expD = interp1(experimentData(iG).time, ...
		   experimentData(iG).dist, ...
		   simData.time(tIdx));

    gcError(iG) = sqrt(sum((expD - simData.dist(tIdx)).^2));

  end

  errorValue = sum(gcError);

  fprintf('Fitness value %d\n', errorValue)

end
