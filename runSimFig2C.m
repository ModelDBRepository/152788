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

function runSimFig2C()

  disp('Preparing diffusion and active transport simulation, Y morphology')

  % actRange = logspace(-3,2,20) * 440e-9*6e-3;
  actRange = logspace(-2,2,20) * 440e-9*6e-3;  
  diffRange = 1e-11; 
  xRange = 50e-6; %linspace(10e-6,100e-6,3);
  yRange = 50e-6; %linspace(20e-6,100e-6,3);
 
  nWorkers = 3;
  nJobs = length(diffRange)*length(xRange)*length(yRange)*length(actRange);

  jobID = mod(0:nJobs-1,nWorkers)+1;

  inputFilenameMask = 'input/Fig2C-diffusion-and-active-transport-D-%d-A-%d-X-%d-Y-%d.input';
  outputFilenameMask = 'output/Fig2C-diffusion-and-active-transport-D-%d-A-%d-X-%d-Y-%d.output';

  fidSum = fopen('input/Fig2C-diffusion-and-active-transport-summary.txt','w');
  ctr = 1;

  for i = 1:length(diffRange)
    for j = 1:length(xRange)
      for k = 1:length(yRange)
        for m = 1:length(actRange)

          inFilename = sprintf(inputFilenameMask,i,m,j,k);
          outFilename = sprintf(outputFilenameMask,i,m,j,k);
          
          fid = fopen(inFilename,'w');
          fprintf(fid,'Experiment.tubulinDiffusionConstant = %d\n', diffRange(i));
          fprintf(fid,'Experiment.tubulinActiveTransportRate = %d\n', actRange(m));
          fprintf(fid,'self.distA = %d\n', xRange(j));
          fprintf(fid,'self.distB = %d\n', yRange(k));
          fprintf(fid,'self.saveFileName = "%s"\n', outFilename);

          % !!! Added soma clamp!
          fprintf(fid, 'self.clampSomaConcentration = True\n');
          fprintf(fid, 'self.clockEnd = 5e5\n');
          fprintf(fid,'self.polyRateModifier = 1.5\n');
        
          fprintf(fid, 'Experiment.tubulinConcentrationSoma = 5.5e-3\n');
                    
          fclose(fid);

          % First column is the ID of the worker that is responsible for
          % running this simulation, second column is the file name of the
          % info file that has all parameters etc.
     
          fprintf(fidSum,'%d %s %d %d %d %d\n', ...
                  jobID(ctr), inFilename, diffRange(i), actRange(m), xRange(j), yRange(k));

          ctr = ctr + 1;

        end
      end
    end
  end

  fclose(fidSum);
 
  disp('python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 1')
  disp('python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 2')
  disp('python runSimWorker.py input/Fig2C-diffusion-and-active-transport-summary.txt 3')


end
  