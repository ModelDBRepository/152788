% This function reads in multiple simulations and show the result in a
% aggregated figure.

% clear all
close all, format compact


dataPath = 'output/';
  
dataFileMask = 'Fig2C-diffusion-and-active-transport-D-%d-A-%d-X-%d-Y-%d.output.perturbed.txt';
transientFile = 'Fig2C-diffusion-and-active-transport-D-%d-A-%d-X-%d-Y-%d.output.original.txt';

summaryFile = 'input/Fig2C-diffusion-and-active-transport-summary.txt';
summaryDataMask = '%d input/Fig2C-diffusion-and-active-transport-D-%d-A-%d-X-%d-Y-%d.input %f %f %f %f\n';

if(~exist('data'))

  baselineData = struct('time',[], 'ID', [], 'parentID', [], ...
                        'x1', [], 'y1', [], 'z1', [], ...
                        'x2', [], 'y2', [], 'z2', [], ...
                        'r', [], 'tubulinConc', [], 'flux',[],'dist', []);

  data = struct('time',[], 'ID', [], 'parentID', [], ...
                'x1', [], 'y1', [], 'z1', [], ...
                'x2', [], 'y2', [], 'z2', [], ...
                'r', [], 'tubulinConc', [], 'flux', [], 'dist', []);

  runID = [];
  GCidBase = [];
  GCidData = [];
  branchIDData = [];

  ctr = 1;

  nD = 1;
  nA = 20;
  nX = 1;
  nY = 1;

  for iD = 1:nD
    for iX = 1:nX
      for iY = 1:nY
        for iA = 1:nA

          baselineData(ctr) = readData([dataPath sprintf(transientFile,iD,iA,iX,iY)]);
          data(ctr) = readData([dataPath sprintf(dataFileMask,iD,iA,iX,iY)]);
          runID(ctr,:) = [iD, iA, iX, iY];

          % Get ID of growth cones. Only compartments that are not
          % parents to any other compartments.
        
          GCid(ctr,:) = setdiff(baselineData(ctr).ID, ...
                                baselineData(ctr).parentID); 
          
          % The growth cone with the lower ID is the one that is
          % getting perturbed.
          
          % Get ID of compartments just past branch point
          endTime = max(data(ctr).time);
          lastIdx = find(data(ctr).time == endTime);
          allParents = unique(data(ctr).parentID(lastIdx));
          n = hist(data(ctr).parentID(lastIdx),allParents);
          branchParent = allParents(find(n == 2));
          
          branchID(ctr,:) = ...
              data(ctr).ID(lastIdx(data(ctr).parentID(lastIdx) ...
                                   == branchParent)); 
          
          % Make sure that the branchIDs match the growth cone IDs
          compID = GCid(ctr,1);
          compIdx = lastIdx(data(ctr).ID(lastIdx) == compID);
          parentID = data(ctr).parentID(compIdx);
          
          while(compID ~= branchID(ctr,1) ...
                & compID ~= branchID(ctr,2))
            
            compID = parentID;
            compIdx = lastIdx(data(ctr).ID(lastIdx) == compID);
            parentID = data(ctr).parentID(compIdx);
            
            assert(parentID ~= -1); % If this happens we never found
                                    % the branch points.
          end
          
          if(compID == branchID(ctr,1))
            % IDs are in the right order
            
          elseif(compID == branchID(ctr,2))
            % They are flipped, fix it!
            branchID(ctr,:) = branchID(ctr,[2 1]);
          else
            disp('This should not happen!')
            keyboard
          end
          
          
          ctr = ctr + 1;
          
        end
      end
    end
  end

  % Read the summary file

   workerID = [];
   fileTag = [];
   runDiffusion = [];
   runActTransp = [];
   runX = [];
   runY = [];

   fid = fopen(summaryFile,'r');
   str = fgets(fid);
   ctr = 1;
   while(str ~= -1)
     [workerID fileTag(ctr,1), fileTag(ctr,2), fileTag(ctr,3), fileTag(ctr,4), ...
       runDiffusion(ctr), runActTransp(ctr), runX(ctr), runY(ctr)] = ...
       strread(str,summaryDataMask);
     str = fgets(fid);
     ctr = ctr + 1;
   end
   fclose(fid);

   % Now match the fileTag to the runID

   remapIdx = [];  

   for i = 1:size(runID,1)
     idx = find(runID(i,1) == fileTag(:,1) ...
                & runID(i,2) == fileTag(:,2) ...
                & runID(i,3) == fileTag(:,3) ...
                & runID(i,4) == fileTag(:,4));
     
     if(length(idx) ~= 1)
       disp('More than one row with the same ID in summary file')
       keyboard
     end

     remapIdx(i) = idx;
   end

   % Remap the indexes
   runActTransp = runActTransp(remapIdx);
   runDiffusion = runDiffusion(remapIdx);
   runX = runX(remapIdx);
   runY = runY(remapIdx);
   clear workerID fileTag % To avoid confusion

else
  disp('Data already loaded, skipping reading. Clear data to reread')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User is responsible for picking values matching the saved step size
tAfter = 108000; %10000;

% Here we are making the assumption that we do not vary X or Y
assert(all(all(runID(:,3:4) == 1)))

% Find the location of the growth cone after tAfter seconds
% for baseline case, and for perturbed case

tAfterClosest = zeros(numel(data),1);


for i = 1:numel(data)
  
  tU = unique(data(i).time);
  [maxError,closestTimeIdx] = min(abs(tU - (tAfter+min(tU))));
  tAfterClosest(i) = tU(closestTimeIdx);
  assert(maxError < tAfter * 0.01)

  baseDist(i,1) = baselineData(i).dist(find(baselineData(i).ID == GCid(i,1) ...
                                            & baselineData(i).time == tAfterClosest(i)));
  baseDist(i,2) = baselineData(i).dist(find(baselineData(i).ID == GCid(i,2) ...
                                            & baselineData(i).time == tAfterClosest(i)));    
  
  pertDist(i,1) = data(i).dist(find(data(i).ID == GCid(i,1) ...
                                    & data(i).time == tAfterClosest(i)));
  pertDist(i,2) = data(i).dist(find(data(i).ID == GCid(i,2)...
                                    & data(i).time == tAfterClosest(i)));  

end


[actTransp,aIdx] = sort(runActTransp);

fig = figure;

p = semilogx(actTransp,1e6*baseDist(aIdx,1),'-k', ...
             actTransp,1e6*baseDist(aIdx,2),'k-', ...
             actTransp,1e6*pertDist(aIdx,1),'r-', ...
             actTransp,1e6*pertDist(aIdx,2),'b-', ...
             'linewidth',1);

legend(p(2:4),'Control','Perturbed +50%','Neighbour','location','northwest');
xlabel('Fraction of tubulin actively transported','fontsize',30)
%ylabel('Neurite length (micrometer)', 'fontsize',30)
ylabel('Neurite length (\mum)', 'fontsize',30)

title(sprintf('%d hours after pertubation',tAfter/3600))
set(gca,'fontsize',25)
axis tight
box off

saveas(gcf,'FIGS/Fig2C-diffusion-and-active-transport.pdf','pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




plotIdx = numel(data) - 5;
plotIdxRef = 1;

figure
hold on

% Control, with weak and strong active transport

idx = find(baselineData(plotIdxRef).ID == GCid(plotIdxRef,1));
plot(baselineData(plotIdxRef).time(idx)/3600, ...
     baselineData(plotIdxRef).dist(idx)*1e6, ...
     'k--','linewidth',1);
idx = find(baselineData(plotIdxRef).ID == GCid(plotIdxRef,1));
p(1) = plot(baselineData(plotIdxRef).time(idx)/3600, ...
            baselineData(plotIdxRef).dist(idx)*1e6, ...
            'k--', 'linewidth',1);

idx = find(baselineData(plotIdx).ID == GCid(plotIdx,1));
plot(baselineData(plotIdx).time(idx)/3600, ...
     baselineData(plotIdx).dist(idx)*1e6, ...
     'k-','linewidth',1);
idx = find(baselineData(plotIdx).ID == GCid(plotIdx,1));
p(1) = plot(baselineData(plotIdx).time(idx)/3600, ...
            baselineData(plotIdx).dist(idx)*1e6, ...
            'k-', 'linewidth',1);

% Branch with polymerisation increased & neighbour for the two cases

idx = find(data(plotIdxRef).ID == GCid(plotIdxRef,1));
p(2) = plot(data(plotIdxRef).time(idx)/3600,...
            data(plotIdxRef).dist(idx)*1e6, ...
            'r--', 'linewidth',1);

idx = find(data(plotIdxRef).ID == GCid(plotIdxRef,2));
p(3) = plot(data(plotIdxRef).time(idx)/3600, ...
            data(plotIdxRef).dist(idx)*1e6, ...
            'b--','linewidth',1);


idx = find(data(plotIdx).ID == GCid(plotIdx,1));
p(2) = plot(data(plotIdx).time(idx)/3600,...
            data(plotIdx).dist(idx)*1e6, ...
            'r-', 'linewidth',1);

idx = find(data(plotIdx).ID == GCid(plotIdx,2));
p(3) = plot(data(plotIdx).time(idx)/3600, ...
            data(plotIdx).dist(idx)*1e6, ...
            'b-','linewidth',1);

legend(p,'Control','Perturbed +50%','Neighbour','location','northwest');
%legend(p,'Control','Perturbed +50%','Neighbour','location','southeast');
xlabel('Time (hours)','fontsize',30)
% ylabel('Neurite length (micrometers)','fontsize',30)
ylabel('Neurite length (\mum)','fontsize',30)
set(gca,'fontsize',25)

title(sprintf('Diffusion: %.1d m^2/s   Active transport: %.1d and %.1d',...
              runDiffusion(plotIdx), runActTransp(plotIdx), runActTransp(plotIdxRef)))

saveas(gcf,'FIGS/Fig2C-diffusion-and-active-transport-growth-example.pdf','pdf')












