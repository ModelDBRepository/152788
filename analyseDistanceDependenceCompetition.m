% This script and analyseDistanceDependence.m looks at the influence
% of one modified growth cone has on its neighbouring growth cones.

close all

if(exist('data') ~= 1)

  switch(1)

    case 1


      saveFigBase = 'FIGS/Martine-cell1';
      
      fileNameBase = 'DATA/cell1.swc-NOGCmod-out.txt';
      
      GCnum = 0:79; %0:36 %setdiff(0:54,[37 38 42]); %0:40

      for i = 1:length(GCnum)
        fileName{i} = sprintf('DATA/cell1.swc-GCmod-%d-out.txt',GCnum(i));      
      end


    case 4

      saveFigBase = 'FIGS/Martine-cell4';      
      
      fileNameBase = 'DATA/cell4.swc-NOGCmod-out.txt';

      GCnum = 0:46; %0:24; %0:11; % 39?

      for i = 1:length(GCnum)
        fileName{i} = sprintf('DATA/cell4.swc-GCmod-%d-out.txt',GCnum(i));
      end


    end

  dataRef = readData(fileNameBase);

  for i=1:length(fileName)
    data(i) = readData(fileName{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nextFig = 1;


if(1)
  % Make dendrogram figure

  for i = 1:length(data)
    figure(nextFig)
    set(gcf,'visible','off')

    switch(6)
      case 1
        % All neurites same diameter in fig
        makeDendrogram(dataRef,'end','k-',5)
        makeDendrogram(data(i),'end','r-',2,GCnum(i)+1)
      case 2
        % Scale the neurites
        makeDendrogram(dataRef,'end','k-',-10)
        makeDendrogram(data(i),'end','r-',-2,GCnum(i)+1)

      case 3
        makeDendrogram(dataRef,'end','k-',-10)
        makeDendrogram(data(i),'end','gradient',-2,GCnum(i)+1)

      case 4
        makeDendrogram(dataRef,500,'gradient',-10)

      case 5
        makeDendrogram(dataRef,'end','-',-10,NaN,[1 1 1]*0.5)
        makeDendrogram(dataRef, 0,'k-',-10)
        makeDendrogram(data(i),'end','gradient',-2,GCnum(i)+1)

      case 6
        makeDendrogram(data(i),'end','gradient',-2,GCnum(i)+1)
        makeDendrogram(dataRef, 0,'growthcones',-1,NaN,[1 1 1]*0.7)
        makeDendrogram(dataRef,'end','growthcones',-1,NaN,[1 1 1]*0)

        
    end

    if(0)
      % Plotting reference line, do not use t=0, use the first datapoint 
      % larger than 0
      t = unique(dataRef.time);

      makeDendrogram(dataRef,t(2),'b-',1)
    end

    title(fileName{i})
    figName = strcat(fileName{i},'-dendrogram.pdf');
    saveas(gcf,figName,'pdf');

    nextFig = nextFig + 1;
    
    if(strcmpi(get(gcf,'visible'),'off'))
      disp('Closing all')
      close all
    end
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start by identifying the compartment IDs between runs for GCs. 
% Assume compartments are in the same order in reference and in modified
% simulations.
%
% !!! OLD, new code retains the compartment IDs between saves.

GCidRefTMP = setdiff(dataRef.ID,dataRef.parentID); % Incorrectly ordered
%GCidx = NaN*ones(size(GCidRef));

idx = find(dataRef.time == 0);
GCidx = idx(find(ismember(dataRef.ID(idx),GCidRefTMP)));
GCidRef = dataRef.ID(GCidx);

time = unique(dataRef.time);
deltaLref = NaN*zeros(length(GCidRef), length(time));
deltaL = NaN*zeros(length(GCidRef),length(time),length(data));
deltaLrel = NaN*zeros(length(GCidRef),length(time),length(data));

for i = 1:length(GCidRef)
  deltaLref(i,:) = dataRef.dist(dataRef.ID == GCidRef(i));
end

GCdist0 = dataRef.dist(GCidx);

GCIDlookup = zeros(length(GCidRef),length(data));

for i = 1:length(data)
  GCIDlookup(:,i) = data(i).ID(GCidx);

  % Verify that these are indeed growth cones
  id = setdiff(data(i).ID,data(i).parentID);

  if(length(intersect(GCIDlookup(:,i),id)) ~= length(id))
    disp('Mismatch between reference and modified simulations')
    disp('Are you using the wrong files?')
    keyboard
  end

end

for i = 1:length(data)
  for j = 1:length(GCidRef)
    tmp = data(i).dist(find(data(i).ID == GCIDlookup(j,i)));
    deltaL(j,1:numel(tmp),i) = tmp;
    if(numel(tmp) < size(deltaL,2))
      fprintf(['Data %d: Growth cone %d disappeared, keeping last distance ' ...
               'constant.\n'], i, GCidRef(i))
      deltaL(j,numel(tmp)+1:end,i) = tmp(end);
    end
  end

  deltaLrel(:,:,i) = (deltaL(:,:,i) - deltaLref) ./ deltaLref;

end

% Define competition as the ratio of length change in the modified growth
% cone divided by the length change in all other growth cones at T=end.


for i = 1:length(data)
  gIdx = GCnum(i)+1;

  deltaLmodGC(i) = deltaL(gIdx,end,i) - deltaLref(gIdx,end);
  restIdx = setdiff(1:length(GCidRef),gIdx);
  deltaLrestGC(i) = sum(deltaL(restIdx,end,i) - deltaLref(restIdx,end));

end

%plot(GCdist0(1:length(data)),-deltaLmodGC./deltaLrestGC,'k*')


pIdx = 30;

figure
plot(time/3600,deltaLrel(:,:,pIdx),'k-', ...
     [min(time) max(time)]/3600, [0 0], 'r--')
xlabel('Time (hours)','fontsize',20)
ylabel('Growth (relative to reference simulation','fontsize',20)
set(gca,'fontsize',16);
box off

fName = sprintf('%s-growth-gc-%d.pdf',saveFigBase, pIdx);
saveas(gcf,fName,'pdf');

figure
p = plot(time/3600, 1e6*(deltaL(:,:,pIdx)-deltaLref), 'k', ...
	 time/3600, 1e6*sum(deltaL(:,:,pIdx)-deltaLref,1), 'k--', ...
	 [min(time) max(time)]/3600, [0 0], 'r--');
legend(p([1 end-1 end]),'Growth cone growth','Summed growth','0', ...
       'location','best')
xlabel('Time (hours)','fontsize',20)
ylabel('Growth, relative to reference (micrometers)','fontsize',20)
set(gca,'fontsize',16)
box off

fName = sprintf('%s-growth-rel-to-ref.pdf',saveFigBase);
saveas(gcf,fName,'pdf');


figure
plot(GCdist0(1:length(data))*1e6,-deltaLrestGC./deltaLmodGC,'k.','markersize',16)
%xlabel('Distance from soma of modified growth cone (\mum)','fontsize',20)
%ylabel('Total retraction length / modified growth cone elongation','fontsize',20)
xlabel('Distance from soma (\mum)','fontsize',30)
ylabel('Relative retraction','fontsize',30)
set(gca,'fontsize',25);
box off

fName = sprintf('%s-growth-ret-div-by-elongation.pdf',saveFigBase);
saveas(gcf,fName,'pdf');
