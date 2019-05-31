% This function read in the summary files then calculates the fitness 
% for all the values.

clear all, close all, format compact

filenameMask = ...
  'output/GridSearch/predict-growth-speed-Di-%d-A-%d-SC-%d-PR-%d-DR-%d-De-%d.txt';

summaryFile = ...
  'input/GridSearch/predict-growth-speed-summary.txt';

summaryDataMask = '%d input/GridSearch/predict-growth-speed-Di-%d-A-%d-SC-%d-PR-%d-DR-%d-De-%d.txt %f %f %f %f %f %f';


  % !!!! This must match
  GCtoPredict = 2; % 0, 1 or 2
  fprintf('Assuming predicting GC %d\n', GCtoPredict)

  altGC = 1; % GC to verify with

  growthConeFiles = { 'input/Ramaker/Ramaker-980625-GC-1.txt', ...
		      'input/Ramaker/Ramaker-980625-GC-2.txt', ...
		      'input/Ramaker/Ramaker-980625-GC-3.txt' };


  % First read in summary data

  workerID = [];
  fileTag = [];
  runDiffusion = [];
  runActTransp = [];
  runSomaConc = [];
  runPolyRate = [];
  runDepolyRate = [];
  runTubulinDecay = [];

  fid = fopen(summaryFile,'r');
  str = fgets(fid);
  ctr = 1;
  while(str ~= -1)

    % Remove end of line
    if(str(end) == char(10))
      str = str(1:end-1);
    end

    try
      [workerID fileTag(ctr,1), fileTag(ctr,2), ...
	        fileTag(ctr,3), fileTag(ctr,4), ...
	        fileTag(ctr,5), fileTag(ctr,6) ...
       runDiffusion(ctr), runActTransp(ctr), ...
       runSomaConc(ctr), runPolyRate(ctr), ...
       runDepolyRate(ctr), runTubulinDecay(ctr)] = ...
         strread(str,summaryDataMask);
    catch exception
      % I messed up somewhere... oops.
      getReport(exception)
      keyboard
    end

    str = fgets(fid);
    ctr = ctr + 1;
  end
 
  fclose(fid);


% Next, calculate the fitness value for all of them

  fitError = NaN*zeros(size(workerID));

  for i = 1:size(fileTag,1)
    
    fileName = sprintf(filenameMask,fileTag(i,:));
    allFiles{i} = fileName;

  end

  matlabpool open

  % Parallelise the slow step
  parfor i = 1:size(fileTag,1)

    fileName = sprintf(filenameMask,fileTag(i,:));
    fitError(i) = growthFitness(fileName,growthConeFiles);

  end

  matlabpool close


figure
plot(runDiffusion,fitError,'k.')
xlabel('Diffusion')
ylabel('Error')

figure
plot(runActTransp,fitError,'k.')
xlabel('Active transport')
ylabel('Error')

figure
plot(runSomaConc, fitError, 'k.')
xlabel('Soma tubulin concentration')
ylabel('Error')

figure
plot(runPolyRate,fitError,'k.')
xlabel('Polymerisation rate')
ylabel('Error')

figure
plot(runDepolyRate, fitError, 'k.')
xlabel('Depolymerisation rate')
ylabel('Error')

figure
plot(runTubulinDecay, fitError, 'k.')
xlabel('Tubulin decay')
ylabel('Error')

uDiff = unique(runDiffusion);
uDecay = unique(runTubulinDecay);

[diffAll,decayAll] = meshgrid(runDiffusion,runTubulinDecay);

tagIdx = [2 3 4 5];
uTag = unique(fileTag(:,tagIdx),'rows');

for j = 1:size(uTag,1)
  errorAll = NaN*ones(size(diffAll));

  for i = 1:numel(diffAll)
    idx = find(diffAll(i) == runDiffusion(:) ...
	       & decayAll(i) == runTubulinDecay(:) ... 
	       & fileTag(:,tagIdx(1)) == uTag(j,1) ...
	       & fileTag(:,tagIdx(2)) == uTag(j,2) ...
	       & fileTag(:,tagIdx(3)) == uTag(j,3) ...
	       & fileTag(:,tagIdx(4)) == uTag(j,4));
     errorAll(i) = fitError(idx);

   end
 
   if(0)
     figure
     surf(diffAll,decayAll,errorAll)
     xlabel('Diffusion')
     ylabel('Decay')
     title(num2str(uTag(j,:)))
  end
end


% Plot the best fit...
idx = find(fitError == min(fitError));

if(length(idx) > 3)
  disp('Multiple files with same fittness. Only showing the first three first')
end

for i = 1:min(length(idx),3)
  visualisePredictionSpeedFcn(allFiles{idx(i)}, fitError(idx(i)));
end




%%%%%% Run the simulation with another growth cone to check

inFile = strrep(allFiles{idx(1)},'output','input');


%[newFile, fakeSummary] = verifyParameterPredictionOnOtherGrowthCone(inFile,1);
[newFile, fakeSummary] = verifyParameterPredictionOnOtherGrowthCone(inFile,altGC);
%runStr = sprintf('python runSimWorkerPredictSpeed.py %s 1 > regrow-output.txt', fakeSummary)
runStr = sprintf('python runSimWorkerPredictSpeed.py %s 1', fakeSummary)
system(runStr)

visualisePredictionSpeedFcn(newFile, growthFitness(newFile));
