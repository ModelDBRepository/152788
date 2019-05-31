% This matlab scripts sets up the simulation parameters for searching
% a range of parameters within the space.
%
%

% First scenario:
% We want to start with diffusion only, within range of 1e-14 to 1e-9 m^2/s
%
% The branching is as normal, one primary branch that is X micrometers long
% which then split into two branches that are Y micrometers each.
% 

function runSimFig2A()

  disp('Preparing diffusion only simulation, Y morphology')

  % For biological moleculs, range 1e-11 to 1e-10 (Wikipedia):
  % http://en.wikipedia.org/wiki/Fick's_laws_of_diffusion
  %
  diffRange = linspace(1e-12,1e-9,10); % Lower end was 1e-13 before
  % diffRange = unique([2e-10, logspace(-12,-9,10)]); % Lower end was 1e-13 before  
  xRange = 50e-6; %linspace(10e-6,100e-6,3);
  yRange = 50e-6; %linspace(10e-6,100e-6,3);
 
  nWorkers = 3;
  nJobs = length(diffRange)*length(xRange)*length(yRange);

  jobID = mod(0:nJobs-1,nWorkers)+1;

  inputFilenameMask = 'input/Fig2A-diffusion-only-dR-%d-X-%d-Y-%d.input';
  outputFilenameMask = 'output/Fig2A-diffusion-only-dR-%d-X-%d-Y-%d.output';

  fidSum = fopen('input/FIG2A-diffusion-only-summary.txt','w');
  ctr = 1;

  for i = 1:length(diffRange)
    for j = 1:length(xRange)
      for k = 1:length(yRange)

        inFilename = sprintf(inputFilenameMask,i,j,k);
        outFilename = sprintf(outputFilenameMask,i,j,k);

        fid = fopen(inFilename,'w');
        fprintf(fid,'Experiment.tubulinDiffusionConstant = %d\n', diffRange(i));
        fprintf(fid,'Experiment.tubulinActiveTransportRate = %d\n', 0);
        fprintf(fid,'self.distA = %d\n', xRange(j));
        fprintf(fid,'self.distB = %d\n', yRange(k));
        fprintf(fid,'self.saveFileName = "%s"\n', outFilename);
        
        % !!! Added soma concetration clamp
        fprintf(fid, 'self.clampSomaConcentration = True\n');
        fprintf(fid, 'self.clockEnd = 5e5\n');
        fprintf(fid,'self.polyRateModifier = 1.5\n');
        
        fprintf(fid, 'Experiment.tubulinConcentrationSoma = 5.5e-3\n');

        fclose(fid);

        % First column is the ID of the worker that is responsible for
        % running this simulation, second column is the file name of the
        % info file that has all parameters etc.
     
        fprintf(fidSum,'%d %s %d %d %d\n', ...
                jobID(ctr), inFilename, diffRange(i), xRange(j), yRange(k));

        ctr = ctr + 1;

      end
    end
  end

  fclose(fidSum);
 
  disp('python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 1')
  disp('python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 2')
  disp('python runSimWorker.py input/FIG2A-diffusion-only-summary.txt 3')

end
