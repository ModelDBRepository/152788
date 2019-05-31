% This code writes the files necessary to run the predict speed simulations
% with a set of different parameters.
%
%

%
% setupPredictSpeedSimulation
% python runSimWorkerPredictSpeed.py input/GridSearch/predict-growth-speed-summary.txt 1
% python runSimWorkerPredictSpeedParallel.py input/GridSearch/predict-growth-speed-summary.txt
% calculatePredictionFitness

function setupPredictSpeedSimulation()

  % 1.3e-9 probably upper limit, assuming water and molecule with
  % diameter 0.5nm

  %diffRate = linspace(1e-13,1e-9,11); %linspace(1e-13,1e-9,5);
  %actRate = [0 1 10 100 1000] * 440e-9*6e-3;

  % Tried with 0.25e-9 and higher for diffusion constant, and that was
  % waaay too high.

  %diffRate = linspace(1e-13,1e-10,25)
  diffRate = linspace(1e-13,0.5e-10,15)
  actRate = 440e-9*6e-3 * [0 1 10 100];

  %tubulinSomaConc = 10e-3 *linspace(0.51,0.55,5); %[0.55 0.75 1] %[5e-3 10e-3];
  tubulinSomaConc = linspace(5.5e-3,50e-3,5)

  % These two are only used when predicting growth
  tubulinPolymerisationRate = [1]*33e-6/(5e-3*3600); 
  tubulinDepolymerisationRate = [1]*33e-6/3600;

  %growthSpeedAt2xC = linspace(5e-6,50e-6,5);
  %zeroGrowthConc = 5e-3;
  %
  %tubulinPolymerisationRate = growthSpeedAt2xC/(zeroGrowthConc*3600);
  %tubulinDepolymerisationRate = growthSpeedAt2xC/3600;

  tubulinDegradation = 5.67e-7 * [1 10 100 1e3] %[0.1 1 10 100 1000 1e4];
  % tubulinProductionRate = [1]*10e-3*(4*pi*(10e-6)^2) % !!! Double check

  clockEndTime = 1.6e5;

  morphFile = 'input/Ramaker/Ramaker-980625-neuronMorph.txt';

  % !!!!!!!!!!!!!
  % Need to change GCtoPredict also in calculatePredictionFitness
  GCtoPredict = 2; % 0, 1 or 2


  predictGC = sprintf('[%d]', GCtoPredict); 

  growthConeFiles = ['[''input/Ramaker/Ramaker-980625-GC-1.txt'',' ...
		     '''input/Ramaker/Ramaker-980625-GC-2.txt'',' ...
		     '''input/Ramaker/Ramaker-980625-GC-3.txt'']'];

  % This one is used if we want to see if parameter set is possible
  %noNegativeConcentrations = 'False'; % This allows negative concentrations

  % This one is used when we test fitness
  noNegativeConcentrations = 'True'; % No neg allowed, stop growth instead

  nWorkers = 1;
  nJobs = length(diffRate) ...
        * length(actRate) ...
        * length(tubulinSomaConc) ...
        * length(tubulinPolymerisationRate) ...
        * length(tubulinDepolymerisationRate) ...
        * length(tubulinDegradation); % ...

        % * length(tubulinProductionRate);

  jobID = mod(0:nJobs-1,nWorkers)+1;

  inputFilenameMask = ...
    'input/GridSearch/predict-growth-speed-Di-%d-A-%d-SC-%d-PR-%d-DR-%d-De-%d.txt';
  outputFilenameMask = ...
    'output/GridSearch/predict-growth-speed-Di-%d-A-%d-SC-%d-PR-%d-DR-%d-De-%d.txt';

  fidSum = fopen('input/GridSearch/predict-growth-speed-summary.txt','w');
  ctr = 1;

  for iDi = 1:length(diffRate)
    for iA = 1:length(actRate)
      for iSC = 1:length(tubulinSomaConc)
        for iPR = 1:length(tubulinPolymerisationRate)
          for iDR = 1:length(tubulinDepolymerisationRate)
            for iDe = 1:length(tubulinDegradation)
	      % for iP = 1:length(tubulinProductionRate)


  	        inFilename = sprintf(inputFilenameMask, ...
                                 iDi, iA, iSC, iPR, iDR, iDe);
            outFilename = sprintf(outputFilenameMask, ...
                                  iDi, iA, iSC, iPR, iDR, iDe);

            fid = fopen(inFilename,'w');
            fprintf(fid,'Experiment.tubulinDiffusionConstant = %d\n', ...
                    diffRate(iDi));
            fprintf(fid,'Experiment.tubulinActiveTransportRate = %d\n', ...
                    actRate(iA));
            fprintf(fid,'Experiment.tubulinConcentrationSoma = %d\n', ...
                    tubulinSomaConc(iSC));
            fprintf(fid,'Experiment.neuriteGrowthPoly = %d\n', ...
                        tubulinPolymerisationRate(iPR));
            fprintf(fid,'Experiment.neuriteGrowthDepoly = %d\n', ...
                    tubulinDepolymerisationRate(iDR));
            fprintf(fid,'Experiment.tubulinDegradationConstant = %d\n', ...
                    tubulinDegradation(iDe));
            fprintf(fid,'self.saveFileName = "%s"\n', outFilename);
            fprintf(fid,'self.clockEnd = %d\n', clockEndTime);
            fprintf(fid,'Experiment.morphFile = "%s"\n', morphFile);
            fprintf(fid,'Experiment.growthConeFiles = %s\n', ...
                    growthConeFiles);
            fprintf(fid,'Experiment.predictGCnumber = %s\n', predictGC);
            fprintf(fid,'self.solver.noNegativeConcentrations=%s\n', ...
                    noNegativeConcentrations);
            
            fprintf(fid, 'self.clampSomaConcentration = True\n'); %!!!!
            
            
            % Not saving production rate, calculated
            % in experimentPredictGrowthSpeed.py
            

            fclose(fid);
            
            fprintf(fidSum,'%d %s %d %d %d %d %d %d\n', ...
                    jobID(ctr), ...
                    inFilename, ...
                    diffRate(iDi), ...
                    actRate(iA), ...
                    tubulinSomaConc(iSC), ...
                    tubulinPolymerisationRate(iPR), ...
                    tubulinDepolymerisationRate(iDR), ...
                    tubulinDegradation(iDe));
            
            
            ctr = ctr + 1;
            
            % end
            end
          end
        end
      end
    end
  end

  fclose(fidSum);

end
