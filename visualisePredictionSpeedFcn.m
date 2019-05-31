% This script visualises the growth of the neuron, used for the prediction
% simulations.

%close all, format compact
   function visualisePredictionSpeedFcn(dataFile, errorValue)

%dataPath = 'output/';
%predictGC = unique([]);
%dataFile = 'GC-slavery-predict';
%
%for i = 1:length(predictGC)
%  dataFile = sprintf('%s-%d', dataFile, predictGC(i));
%end
%
%dataFile = strcat(dataFile,'.txt');

morphFile = 'input/Ramaker/Ramaker-980625-neuronMorph.txt';
gcFiles = { 'input/Ramaker/Ramaker-980625-GC-1.txt', ...
	    'input/Ramaker/Ramaker-980625-GC-2.txt', ...
	    'input/Ramaker/Ramaker-980625-GC-3.txt' };



  % Read new data
  % data = readData([dataPath dataFile]);
  data = readData(dataFile);

  % Saved for future reference
  % oldPredictGC = predictGC;

  % Read in morphology and growth cones

  tmp = load(morphFile);
  morphData.x = tmp(:,1);
  morphData.y = tmp(:,2);

  gcData = struct('time',[],'speed',[],'dist',[]);

  maxTime = -inf;
  
  for i = 1:length(gcFiles)
    tmp = load(gcFiles{i});
    tmpData.time  = tmp(:,1);
    tmpData.speed = tmp(:,2);
    tmpData.dist  = tmp(:,3);

    gcData(i) = tmpData;
    
    if(max(tmpData.time) > maxTime)
      maxTime = max(tmpData.time);
    end
  end

  % Add endpoints
  for i = 1:length(gcData)
    if(gcData(i).time < maxTime)
      disp('Added end point, assuming same dist as previous point')
      gcData(i).time = [gcData(i).time; maxTime];
      gcData(i).speed = [gcData(i).speed; gcData(i).speed(end)];
      gcData(i).dist = [gcData(i).dist; gcData(i).dist(end)];      
    end
  end
  

%%%% Now we got to plot it...

gcColor = [1 0 0; 0 1 0; 0 0 1];
gcID = setdiff(data.ID,data.parentID);

for i = 1:length(gcID)
  gcIDparent(i) = data.parentID(find(gcID(i) == data.ID,1,'last'));
end

nComp = zeros(size(gcID));
for i=1:length(gcID)
  nComp(i) = nnz(data.ID == gcID(i));
  fprintf('%d: %d\n', gcID(i), nComp(i))
end

if(0)
  plotConcentrationGradient(data,2*3600);
  plotConcentrationGradient(data,7*3600);
  plotConcentrationGradient(data,12*3600);
  plotConcentrationGradient(data,15*3600);
  plotConcentrationGradient(data,27*3600);
end

% Plot GC concentration evolution

figure
for i = 1:length(gcFiles)
  idx = find(data.ID == gcID(i));
  plot(data.time(idx)/3600,data.tubulinConc(idx),'k-')
  hold on  
end
xlabel('Time (hours)')
ylabel('Concentration (mM)')
title(sprintf('Error: %d', errorValue))


% Find the branch points.

endTime = max(data.time);
idx = find(data.time == endTime);

pIDsorted = sort(data.parentID(idx));
sameParentID = pIDsorted(find(diff(pIDsorted) == 0));
branchParentID = unique(sameParentID);

branchChildID = [];
for i = 1:length(branchParentID)
  branchChildID = [branchChildID, ...
	          data.ID(idx(find(data.parentID(idx) ...
				   == branchParentID(i))))];
end

figure
for i = 1:length(branchChildID)
  idx = find(data.ID == branchChildID(i));
  plot(data.time(idx)/3600,data.flux(idx),'k-','linewidth',i)
  hold on

end

plot([min(data.time(idx)) max(data.time(idx))]/3600,[0 0],'r-')

xlabel('Time (hours)','fontsize',24)
ylabel('Flux at branch','fontsize',24)
set(gca,'fontsize',20)



% Plot flux

figure
for i = 1:length(gcFiles)

  %idx = find(data.ID == gcID(i));
  idx = find(data.ID == gcIDparent(i));

  plot(data.time(idx)/3600,data.flux(idx),'k-')
  hold on  
end
xlabel('Time (hours)')
ylabel('Flux')
title(sprintf('Error: %d', errorValue))


% Plot distances

figure
for i = 1:length(gcFiles)

  plot(gcData(i).time/3600,gcData(i).dist*1e6,'k--')

  idx = find(data.ID == gcID(i));
  hold on  
  plot(data.time(idx)/3600,data.dist(idx)*1e6,'color',gcColor(i,:))
end
xlabel('Time (hours)','fontsize',24)
ylabel('Distance (\mum)','fontsize',24)
set(gca,'fontsize',20)
title(sprintf('Error: %d', errorValue))




if(0)
  figure,

  uID = unique(data.ID); 
  for i = 1:length(uID), 
    idx = find(data.ID==uID(i)); 
    plot(data.time(idx),data.dist(idx),'k')
    hold on

    if(nnz(diff(data.dist(idx)) > 0) > 5)
      fprintf('GC? id: %d\n', uID(i))
    end

  end
  xlabel('Time (s)')
  ylabel('Distance (m)')

  allGCTime = [];

  for i=1:length(gcID)

    idx = find(data.ID == gcID(i));
    gcTime{i} = data.time(idx);
    gcDist{i} = data.dist(idx) - data.dist(idx(1));
    allGCTime = union(allGCTime, gcTime{i});

  end

  totalGCdist = zeros(size(allGCTime));

  for i = 1:length(gcID)
    for j = 1:length(gcTime{i})
      idx = find(allGCTime == gcTime{i}(j));
      totalGCdist(idx) = totalGCdist(idx) + gcDist{i}(j);
    end
  end

  figure
  plot(allGCTime/3600,totalGCdist*1e6,'k-')
  xlabel('Time (hours)','fontsize',24)
  ylabel('Total axon length (\mum)','fontsize',24)
  set(gca,'fontsize',20)
  box off
end

% Function end
end
