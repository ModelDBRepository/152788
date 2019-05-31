% For figure 2A
% 
% This function does not handle active transport simulations, only
% those using just diffusion. See !!! notes below if it is modified
% for active transport also.
%
% This function reads in multiple simulations and show the result in a
% aggregated figure.

% clear all
close all, format compact


% simName = 'diffusion-only';
simName = 'Fig2A-diffusion-only';

numDiff = 10;
numX = 1;
numY = 1;

dataPath = 'output/';
  
dataFileMask = strcat(simName,'-dR-%d-X-%d-Y-%d.output.perturbed.txt');
transientFile = strcat(simName,'-dR-%d-X-%d-Y-%d.output.original.txt');

summaryFile = sprintf('input/%s-summary.txt',simName);
summaryDataMask = strcat('%d input/', strcat(simName,'-dR-%d-X-%d-Y-%d.input %f %f %f\n'));

if(~exist('data'))

  % Read in both baseline data, and disrupted data
  
  baselineData = struct('time',[], 'ID', [], 'parentID', [], ...
                        'x1', [], 'y1', [], 'z1', [], ...
                        'x2', [], 'y2', [], 'z2', [], ...
                        'r', [], 'tubulinConc', [], 'flux',[],'dist', []);

  data = struct('time',[], 'ID', [], 'parentID', [], ...
                'x1', [], 'y1', [], 'z1', [], ...
                'x2', [], 'y2', [], 'z2', [], ...
                'r', [], 'tubulinConc', [], 'flux', [], 'dist', []);

  runID = [];
  GCid = [];
  branchID = [];

  
  ctr = 1;

  for iD = 1:numDiff
    for iX = 1:numX
      for iY = 1:numY

        % Load files
        baselineData(ctr) = readData([dataPath sprintf(transientFile,iD,iX,iY)]);
        data(ctr) = readData([dataPath sprintf(dataFileMask,iD,iX,iY)]);
        runID(ctr,:) = [iD, iX, iY];

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

  % Read the summary file

   workerID = [];
   fileTag = [];
   runDiffusion = [];
   runX = [];
   runY = [];

   fid = fopen(summaryFile,'r');
   str = fgets(fid);
   ctr = 1;
   while(str ~= -1)
     % !!! This needs updating to handle active transport!!
     [workerID fileTag(ctr,1), fileTag(ctr,2), fileTag(ctr,3), ...
      runDiffusion(ctr), runX(ctr), runY(ctr)] = ...
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
		     & runID(i,3) == fileTag(:,3));

     if(length(idx) ~= 1)
       disp('More than one row with the same ID in summary file')
       keyboard
     end

     remapIdx(i) = idx;
   end

   % Remap the indexes
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
assert(all(all(runID(:,2:3) == 1)))

% Find the location of the growth cone after tAfter seconds
% for baseline case, and for perturbed case

tAfterClosest = zeros(numel(data),1);

for i = 1:numel(data)
  
  tU = unique(data(i).time);
  [maxError,closestTimeIdx] = min(abs(tU - (tAfter+min(tU))));
  tAfterClosest(i) = tU(closestTimeIdx);
  assert(maxError < tAfter * 0.1)

  baseDist(i,1) = baselineData(i).dist(find(baselineData(i).ID == GCid(i,1) ...
                                            & baselineData(i).time == tAfterClosest(i)));
  baseDist(i,2) = baselineData(i).dist(find(baselineData(i).ID == GCid(i,2) ...
                                            & baselineData(i).time == tAfterClosest(i)));    
  
  pertDist(i,1) = data(i).dist(find(data(i).ID == GCid(i,1) ...
                                    & data(i).time == tAfterClosest(i)));
  pertDist(i,2) = data(i).dist(find(data(i).ID == GCid(i,2)...
                                    & data(i).time == tAfterClosest(i)));  

end


% Sort them so they are in sequence
[diffusion,dIdx] = sort(runDiffusion);

fig = figure;
p = semilogx(diffusion,1e6*baseDist(dIdx,1),'-k', ...
             diffusion,1e6*baseDist(dIdx,2),'k-', ...
             diffusion,1e6*pertDist(dIdx,1),'r-', ...
             diffusion,1e6*pertDist(dIdx,2),'b-', ...
             'linewidth',1);

legend(p(2:4),'Control','Perturbed +50%','Neighbour','location','northwest');
xlabel('Diffusion (m^2/s)','fontsize',30)
%ylabel('Neurite length (micrometer)', 'fontsize',30)
ylabel('Neurite length (\mum)', 'fontsize',30)

title(sprintf('%d hours after pertubation',tAfter/3600))
set(gca,'fontsize',25)

box off

saveas(gcf,'FIGS/Fig2A-diffusion-only.pdf','pdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plotIdx = 4; %8; % 4 is closer to what we normally do

figure
hold on

idx = find(baselineData(plotIdx).ID == GCid(plotIdx,1));
plot(baselineData(plotIdx).time(idx)/3600, ...
     baselineData(plotIdx).dist(idx)*1e6, ...
     'k-','linewidth',1);
idx = find(baselineData(plotIdx).ID == GCid(plotIdx,1));
p(1) = plot(baselineData(plotIdx).time(idx)/3600, ...
            baselineData(plotIdx).dist(idx)*1e6, ...
            'k-', 'linewidth',1);

idx = find(data(plotIdx).ID == GCid(plotIdx,1));
p(2) = plot(data(plotIdx).time(idx)/3600,...
            data(plotIdx).dist(idx)*1e6, ...
            'r-', 'linewidth',1);

idx = find(data(plotIdx).ID == GCid(plotIdx,2));
p(3) = plot(data(plotIdx).time(idx)/3600, ...
            data(plotIdx).dist(idx)*1e6, ...
            'b-','linewidth',1)

legend(p,'Control','Perturbed +50%','Neighbour','location','northwest');
xlabel('Time (hours)','fontsize',30)
%ylabel('Neurite length (micrometers)','fontsize',30)
ylabel('Neurite length (\mum)','fontsize',30)
set(gca,'fontsize',25)

title(sprintf('Diffusion: %.1d m^2/s',runDiffusion(plotIdx)))

saveas(gcf,'FIGS/Fig2A-diffusion-only-growth-example.pdf','pdf')



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot all growth cone lengths as a function of time

figure
for i = 1:numel(data)
  hold on
  idx = find(baselineData(i).ID == GCid(i,1));
  plot(baselineData(i).time(idx),baselineData(i).dist(idx),'k-')
  idx = find(baselineData(i).ID == GCid(i,2));  
  plot(baselineData(i).time(idx),baselineData(i).dist(idx),'k-')  
  
  idx = find(data(i).ID == GCid(i,1));
  plot(data(i).time(idx),data(i).dist(idx),'r-');
  idx = find(data(i).ID == GCid(i,2));
  plot(data(i).time(idx),data(i).dist(idx),'b-');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












