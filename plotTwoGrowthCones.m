% The purpose of this script is to check if the diffusion is correct
%
% It compares the solution given by iGrow with the analytical solution.
%

clear all, close all

%dataPath = '/home/hjorth/models/growth/iGrow/output/';
dataPath = 'output/';
dataFile = 'two-growth-cones-output.txt';

data = readData([dataPath dataFile]);

time = data.time;
ID = data.ID;
parentID = data.parentID;
tubulinConc = data.tubulinConc;
dist = data.dist;


allIDs = unique(ID);
allTime = unique(time);
allConc = NaN*zeros(length(allTime),length(allIDs));

for i=1:length(allIDs)
  for j = 1:numel(allTime)
    c = tubulinConc(ID == allIDs(i) & time == allTime(j));
    if(~isempty(c))
      allConc(j,i) = c;
    end
  end
end

figure
plot(allTime, allConc)
xlabel('Time (s)')
ylabel('Concentration (mM)')


%%% We also want a plot showing the concentration as a function of distance

C = colormap;

pLeg = {};

uTime = unique(time);

showT = linspace(min(time),max(time),5);
allDistTime = [];
figure

for i=1:length(showT)

  tDiff = abs(uTime - showT(i));
  tClosest = uTime(find(tDiff == min(tDiff),1));
  allDistTime(end+1,1) = tClosest;

  for iC = 1:3
    lineCol(iC) = interp1(linspace(0,1,size(C,1)),C(:,iC),tClosest/max(uTime));
  end

  idx = find(time == tClosest);

  for j=1:length(idx)
    lineStart = dist(idx(j));
    concStart = tubulinConc(idx(j));
    jdx = idx(find(parentID(idx(j)) == ID(idx)));

    if(~isempty(jdx))
      lineEnd = dist(jdx);
      concEnd = tubulinConc(jdx); 
      p = plot([lineStart lineEnd], [concStart concEnd], ...
               'color', lineCol);
      hold on
      pAll(i) = p;
    end
  end

  pLeg{i} = sprintf('%ds',showT(i));
end

% The distance should be to the center of the compartment

xlabel('Distance (\mum)')
ylabel('Concentration (mM)')
legend(pAll,pLeg,'location','best')

saveas(gcf,'output/pics/two-growth-cones-test.pdf','pdf')




% Do a plot that shows distance vs concentration for the growth cones
%

% 1. Find the IDs of the growth cones
% Here a growth cone is defined as a compartment without a child

uID = unique(ID);
gcID = [];

for i=1:length(uID)
  idx = find(uID(i) == parentID);
  if(isempty(idx))
    gcID(end+1) = uID(i);
  end
end


% 2. Plot them

figure

for i=1:length(gcID)

  for iC = 1:3
    lineCol(iC) = interp1(linspace(0,1,size(C,1)),C(:,iC),i/length(gcID));
  end

  idx = find(ID == gcID(i));
  plot(dist(idx),tubulinConc(idx),'color',lineCol)
  hold on

end
xlabel('Arc length (m)')
ylabel('Concentration (mM)')



figure

subplot(2,1,1)
for i=1:length(gcID)

  for iC = 1:3
    lineCol(iC) = interp1(linspace(0,1,size(C,1)),C(:,iC),i/length(gcID));
  end

  idx = find(ID == gcID(i));
  plot(time(idx),tubulinConc(idx),'color',lineCol)
  hold on
end

xlabel('Time (s)')
ylabel('Concentration (mM)')

subplot(2,1,2)
for i=1:length(gcID)

  for iC = 1:3
    lineCol(iC) = interp1(linspace(0,1,size(C,1)),C(:,iC),i/length(gcID));
  end

  idx = find(ID == gcID(i));
  plot(time(idx),dist(idx),'color',lineCol)
  hold on
end

xlabel('Time (s)')
ylabel('Length (m)')
