function data = readData(fileName)

  fprintf('Reading %s\n', fileName)

  data = struct('time',[], 'ID', [], 'parentID', [], ...
 	        'x1', [], 'y1', [], 'z1', [], ...
	        'x2', [], 'y2', [], 'z2', [], ...
	        'r', [], 'tubulinConc', [], ...
		'flux', [], 'dist', []);


  fid = fopen(fileName,'r');

  % Check how many lines there are

  nLines = 0;
  tmp = 1;

  while(tmp ~= -1)
    tmp = fgets(fid);
    nLines = nLines + 1;
  end

  % Reopen file
  fclose(fid);
  fid = fopen(fileName,'r');

  % Subtract the empty line at the end and the header line
  nLines = nLines - 2;

  header = fgets(fid);

  readFlux = 0;

  if(strcmp(header(1:end-1),'time;id;parent id;start coords;end coords;radie;dist;tubulin;flux'))
    readFlux = 1;
    disp('Data file contains net in flux from parent, reading.')
  elseif(~strcmp(header(1:end-1),'time;id;parent id;start coords;end coords;radie;dist;tubulin'))
    disp('Header format is not what was expected')
  end

  fprintf('Reading %d lines of data.\n',nLines)

  % Preallocate
  time = zeros(nLines,1);
  ID = zeros(nLines,1);
  parentID = zeros(nLines,1);
  x1 = zeros(nLines,1);
  y1 = zeros(nLines,1);
  z1 = zeros(nLines,1);
  x2 = zeros(nLines,1);
  y2 = zeros(nLines,1);
  z2 = zeros(nLines,1);
  r = zeros(nLines,1);
  dist = zeros(nLines,1);
  tubulinConc = zeros(nLines,1);

  if(readFlux)
    flux = zeros(nLines,1);
  end

  tmp = fgets(fid);
  i = 1;

  while(tmp ~= -1)
    if(mod(i,10000) == 0)
      fprintf('Reading line %d\n',i)
    end

    try
      % Remove linebreak at end of line, char(10)
      if(~isempty(tmp) & (tmp(end) == '\n' | tmp(end) == char(10)))
        tmp = tmp(1:end-1);
      end


      if(readFlux)
        [time(i,1), ID(i,1), parentID(i,1), ...
         x1(i,1),y1(i,1),z1(i,1), ...
         x2(i,1),y2(i,1),z2(i,1), ...
         r(i,1), dist(i,1), ...
	 tubulinConc(i,1), flux(i,1)] = ...
          strread(tmp,'%f;%f;%f;(%f,%f,%f);(%f,%f,%f);%f;%f;%f;%f');
      else
        [time(i,1), ID(i,1), parentID(i,1), ...
         x1(i,1),y1(i,1),z1(i,1), ...
         x2(i,1),y2(i,1),z2(i,1), ...
         r(i,1), dist(i,1), ...
         tubulinConc(i,1)] = ...
          strread(tmp,'%f;%d;%d;(%f,%f,%f);(%f,%f,%f);%f;%f;%f');
      end
    catch exception
      disp(getReport(exception))
      disp('Something went wrong.... argh!')
      disp('Check if two workers were writing to the same file.')
      keyboard
    end

    i = i + 1;
    tmp = fgets(fid);

  end

  fprintf('Read %d lines total\n',i)

  if(nnz(tubulinConc < 0))
    fprintf('%d positions with negative tubulin concentration\n', ...
	    nnz(tubulinConc < 0))
    beep
  end

  if(nnz(diff(time) < 0))
    disp('This file is corrupted, the time points are not in order.')
    disp('Did two workers write to the same file?')
  end

  % Store everything in a struct

  data.time = time;
  data.ID = ID;
  data.parentID = parentID;
  data.x1 = x1;
  data.y1 = y1;
  data.z1 = z1;
  data.x2 = x2;
  data.y2 = y2;
  data.z2 = z2;
  data.r = r;
  data.tubulinConc = tubulinConc;
  data.dist = dist;

  if(readFlux)
    data.flux = flux;
  end

  fclose(fid);

end
