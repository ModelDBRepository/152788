% Reads the SWC file, use plotSWC to show it.
%
% Send bug reports to:
% johannes.hjorth@cncr.vu.nl


function morphData = readSWC(fileName)

  fprintf('Reading %s\n', fileName)
  fid = fopen(fileName,'r');

  str = fgets(fid);  

  morphData.id = [];
  morphData.type = [];
  morphData.x = [];
  morphData.y = [];
  morphData.z = [];
  morphData.r = [];
  morphData.parent = [];

  while(str ~= -1)

    if(str(1) == '#')
      % Comment line
      str = fgets(fid);  
      continue
    end

    dataStr = strread(str,'%f');
   
    morphData.id(end+1)     = dataStr(1);
    morphData.type(end+1)   = dataStr(2);
    morphData.x(end+1)      = dataStr(3);
    morphData.y(end+1)      = dataStr(4);
    morphData.z(end+1)      = dataStr(5);
    morphData.r(end+1)      = dataStr(6);
    morphData.parent(end+1) = dataStr(7);

    str = fgets(fid);  

  end

end
