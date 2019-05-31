%
% This code uses the previously found parameters for growth cone prediction
% and checks if it works on other growth cones also (where non-predicted GC 
% are slaved).
%

function [outFile, fakeSummaryFile] = verifyParameterPredictionOnOtherGrowthCone(bestParameterFile,newGC,expSet)

  if(~exist('expSet'))
    expSet = 1;
  end

  fid = fopen(bestParameterFile,'r');
  newFile = strrep(bestParameterFile,'.txt', ...
		   sprintf('-newGC-%d.txt', newGC));

  fidOut = fopen(newFile,'w');

  % We want to read in the best parameter set, and change which growth
  % cone should be predicted.

  str = fgets(fid);

  while(str ~= -1)

    if(nnz(strfind(str,'Experiment.predictGCnumber')))
      % Replace the GC number info
      str = sprintf('Experiment.predictGCnumber = [%d]\n', newGC);
    end

    if(nnz(strfind(str,'Experiment.growthConeFiles')))

%      switch(newGC)    
%        case 0
%          growthConeFiles = ['[''input/Ramaker-GC1-growthspeed.txt'',' ...
%		  	     '''input/Ramaker-GC2-growthspeed.txt'',' ...
%			     '''input/Ramaker-GC3-growthspeed.txt'']'];
%        case 1
%          growthConeFiles = ['[''input/Ramaker-GC1-growthspeed.txt'',' ...
%	    	     	     '''input/Ramaker-GC2-growthspeed-predict2.txt'',' ...
%			     '''input/Ramaker-GC3-growthspeed.txt'']'];
%        case 2
%          growthConeFiles = ['[''input/Ramaker-GC1-growthspeed.txt'',' ...
%			     '''input/Ramaker-GC2-growthspeed.txt'',' ...
%			     '''input/Ramaker-GC3-growthspeed-predict3.txt'']';
%        otherwise
%          disp('Unknown GC')
%          growthConeFiles = ['[''input/Ramaker-GC1-growthspeed.txt'',' ...
%		  	     '''input/Ramaker-GC2-growthspeed.txt'',' ...
%			     '''input/Ramaker-GC3-growthspeed.txt'']'];
%
%      end

       switch(expSet)
         case 1

           growthConeFiles = ['[''input/Ramaker/Ramaker-980625-GC-1.txt'',' ...
			      '''input/Ramaker/Ramaker-980625-GC-2.txt'',' ...
			      '''input/Ramaker/Ramaker-980625-GC-3.txt'']'];

         case 2

           growthConeFiles = ['[''input/Ramaker/Ramaker-980513-GC-1.txt'',' ...
			      '''input/Ramaker/Ramaker-980513-GC-2.txt'',' ...
			      '''input/Ramaker/Ramaker-980513-GC-3.txt'']'];

       end


      str = sprintf('Experiment.growthConeFiles = %s\n', growthConeFiles);

    end


    if(nnz(strfind(str,'output/')))
      str = strrep(str,'.txt', sprintf('-newGC-%d.txt', newGC));
      outFile = strrep(str(strfind(str,'=')+1:end),'"','');

      while(outFile(1) == ' ')
        outFile = outFile(2:end);
      end

      if(outFile(end) == char(10))
        outFile = outFile(1:end-1);
      end
 
   end

    fprintf(fidOut,str);

    str = fgets(fid);

  end

  fclose(fid);
  fclose(fidOut);

  fakeSummaryFile = 'input/Swarm/fakeSummary.txt';
  fidSum = fopen(fakeSummaryFile,'w');

  fprintf(fidSum, '1 %s -1 -1 -1 -1 -1 -1\n', newFile);

end
