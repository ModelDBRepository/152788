function exportToPovRay(filePath, fileName, speedup, timePoint)

  if(~exist('timePoint'))
    timePoint = [];
	end

  radieScale = 5;
  cameraZoomIn = 3;

  data = [];

  if(~exist('filePath'))
    filePath = 'output/';
  end

  if(~exist('fileName'))
    fileName = 'output.txt';
  end

  if(~exist('speedup'))
    speedup = 50;
  end


  function importFile(filePath, fileName)
    disp(sprintf('Reading %s%s', filePath, fileName))

    fid = fopen([filePath fileName],'r');
    firstLine = fgetl(fid);

    tmpTokens = textscan(firstLine,'%s','delimiter',';');
    data.columnTokens = tmpTokens{1};

    %for i=1:length(data.columnTokens)
    %  disp(data.columnTokens{i})
    %end
 
    dataLine = fgetl(fid);

    % Number of conc variables stored
    nConc = length(data.columnTokens) - 8;

    formatString = '%d;%d;%d;(%n,%n,%n);(%n,%n,%n);%n';
    for i=1:nConc
      formatString = [formatString ';%d'];
    end

    dataCtr = 1;

    while(dataLine ~= -1)
      tmpData = textscan(dataLine,formatString);
      data.time(dataCtr,1)     = tmpData{1};
      data.id(dataCtr,1)       = tmpData{2};
      data.idParent(dataCtr,1) = tmpData{3};

      data.x1(dataCtr,1)        = tmpData{4};
      data.y1(dataCtr,1)        = tmpData{5};
      data.z1(dataCtr,1)        = tmpData{6};

      data.x2(dataCtr,1)        = tmpData{7};
      data.y2(dataCtr,1)        = tmpData{8};
      data.z2(dataCtr,1)        = tmpData{9};


      data.radie(dataCtr,1)    = tmpData{10};

      for i=1:nConc
        data.conc(dataCtr,i) = tmpData{10+i};
      end

      dataLine = fgetl(fid);      
      dataCtr = dataCtr + 1;
    end

    fclose(fid);

  end

  function writePovRay(renderPath, outFile, t)

    disp(sprintf('Writing %s%s (time=%d)', renderPath, outFile, t))

    fid = fopen([renderPath outFile],'w');

    fprintf(fid,'#include "colors.inc"\n\n');

    % Place camera
    cameraStr = 'camera {\n  location <%f, %f, %f>\n  look_at <%f,%f,%f>\n\n}';
    cameraCoord = [0 0 -1200/cameraZoomIn]; 
    cameraTarg = [0 0 0];                                                      

    fprintf(fid, cameraStr, [cameraCoord cameraTarg]);


    lightStr = 'light_source { <%f, %f, %f> color White}\n\n';

    fprintf(fid, lightStr, [300 500 -500]);
    fprintf(fid, lightStr, [-100 200 150]);

    tmp = abs(data.time - t);
    idx = find(tmp == min(tmp));

    fsSomaStr = [ 'object{ \n' ...
		  '  sphere { <%f, %f, %f>, %f\n' ...
                  '    pigment { color rgb<%f, %f, %f> }\n\n' ...
                  '    finish {\n      phong 1\n    }\n  }\n}\n\n'];


    fsNeurite = ['object {\n' ...
              '  cylinder { <%f, %f, %f>, <%f, %f, %f>, %f\n' ...
              '    pigment { color rgb<%f, %f, %f> }\n\n' ...    
              '    finish {\n      phong 1\n    }\n  }\n}\n'];   

    fsSomaColour = [0.8 0.4 0.4];

    for i=1:length(idx)
      if(data.idParent(idx(i)) == -1)
        disp('Drawing a sphere')
        % Draw the sphere, this is soma without a parent

        % Place soma

        fsSomaCoord = [data.x1(idx(i)) data.y1(idx(i)) data.z1(idx(i))]*1e6;

        somaRadie = data.radie(idx(i))*1e6;

        fprintf(fid, fsSomaStr, [fsSomaCoord somaRadie fsSomaColour]);

      else
        disp('Drawing a cylinder')   
        % Draw a cylinder

        xp1 = data.x1(idx(i))*1e6;
        yp1 = data.y1(idx(i))*1e6;
        zp1 = data.z1(idx(i))*1e6;

        xp2 = data.x2(idx(i))*1e6;
        yp2 = data.y2(idx(i))*1e6;
        zp2 = data.z2(idx(i))*1e6;

        neuriteRadie = data.radie(idx(i))*1e6;

        fprintf(fid, fsNeurite, ...
	        [xp1 yp1 zp1 xp2 yp2 zp2 radieScale*neuriteRadie 0.5*fsSomaColour]);

      end
    end

    fclose(fid);

  end

  function renderMovie(filePath, fileName, speedup)

    renderPath = strcat(filePath, 'pics/');

    if(~exist(renderPath))
      mkdir(renderPath)
    else
      system(sprintf('rm %s*pov', renderPath))
      system(sprintf('rm %s*png', renderPath))
    end

    if(isempty(timePoint))
      timePoint = 1:speedup:length(data.time)
		end

    for i= timePoint
      t = data.time(i);
      outFile = sprintf('%s-%0.4d.pov', strrep(fileName,'.txt',''), t)
      pngFile = sprintf('%s-%0.4d.png', strrep(fileName,'.txt',''), t)

      writePovRay(renderPath, outFile, t);

      curPwd = pwd;
      cd(renderPath);

      % system(sprintf('povray -geometry 1024x768 Output_Alpha=on %s', outFile))
      %system(sprintf('povray -geometry 1024x768 %s', outFile))
      system(sprintf('povray -D -geometry 640x480 %s', outFile))

      convertCmd = sprintf('convert %s -depth 8 %s', pngFile, pngFile)
%      convertCmd = sprintf('convert %s -depth 8 %s', ...
%			   pngFile, strrep(pngFile,'.png','.jpg'))
      system(convertCmd)

      cd(curPwd)

    end

    renderMask = strrep(fileName,'.txt','');

    curPwd = pwd;
    cd(renderPath);
    % -r gives frame rate
    % -y overwrites output file
    % -an no audio
    ffmpegCommand = sprintf('ffmpeg -r 1 -an -i %s-%%04d.png %s.avi', ...
			    renderMask,renderMask)
%    ffmpegCommand = sprintf('ffmpeg -r 1 -an -i %s-%%04d.jpg %s.avi', ...
%			    renderMask,renderMask)


  	if(numel(timePoint) > 1)
      system(ffmpegCommand)
		end

    cd(curPwd);    

  end


  importFile(filePath, fileName);
  renderMovie(filePath, fileName, speedup);

end
