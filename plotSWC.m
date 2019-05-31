% plotSWC, use morphData read with readSWC.m
%
% Send bug reports to:
% johannes.hjorth@cncr.vu.nl

% print -dpng -r300 cell4-testfig.png  


function plotSWC(morphData, labelGCflag)

  if(~exist('labelGCflag'))
    labelGCflag = true;
  end

  % close all

  function [cx,cy,cz] = makeCyl(xp,yp,zp,r, col)

    [ucx,ucy,ucz] = cylinder();
    dx = xp(2)-xp(1);
    dy = yp(2)-yp(1);
    dz = zp(2)-zp(1);

    d = sqrt(dx^2+dy^2+dz^2);

    theta = acos(dz/sqrt(dx^2+dy^2+dz^2));
    fi = atan(dy/dx);

    if(dx < 0)
      fi = fi + pi;
    end

    rotMat = [cos(fi) -sin(fi) 0; sin(fi) cos(fi) 0; 0 0 1] ...
      *[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];

    cx = ucx * r;
    cy = ucy * r;
    cz = ucz * d;

    for i = 1:numel(cx);

      tmp = rotMat*[cx(i);cy(i);cz(i)];
      cx(i) = tmp(1) + xp(1);
      cy(i) = tmp(2) + yp(1);
      cz(i) = tmp(3) + zp(1);

    end

    surf(cx,cy,cz,col*ones(size(cx)));

  end

  % First handle the soma separately. The neurolucida file only has
  % 2D information for the soma, so we need to make a guess for 3D.

  somaIdx = find(morphData.type == 1);

  if(length(somaIdx) > 1)
    sx = morphData.x(somaIdx);
    sy = morphData.y(somaIdx);
    sz = morphData.z(somaIdx);

    somaX = mean(sx);
    somaY = mean(sy);
    somaZ = mean(sz);

    somaR = max(sqrt((sx-somaX).^2 + (sy-somaY).^2 + (sz-somaZ).^2));

    [x,y,z] = sphere(20);
    x = x*somaR;
    y = y*somaR;
    z = z*somaR;

    surf(x+somaX,y+somaY,z+somaZ,0*ones(size(x)));
    hold on

  elseif(length(somaIdx) == 1)
    [x,y,z] = sphere(20);
     x = x*morphData.r(somaIdx);
          y = y*morphData.r(somaIdx);
          z = z*morphData.r(somaIdx);

     surf(x+morphData.x(somaIdx), ....
	  y+morphData.y(somaIdx), ...
	  z+morphData.z(somaIdx), ...
	  0*ones(size(x)));
          hold on
  end

  for i = 1:length(morphData.id)
    % Due to planar somas in the SWC files from neurolucida we
    % have to handle them in a special way, see above
    %
    % if(morphData.parent(i) == -1)
    %   switch(morphData.type(i))
    %     case 1
    %	    [x,y,z] = sphere(20);
    %       x = x*morphData.r(i);
    %       y = y*morphData.r(i);
    %       z = z*morphData.r(i);
    %
    %       surf(x+morphData.x(i),y+morphData.y(i),z+morphData.z(i), ...
    % 	         0*ones(size(x)));
    %       hold on
    %
    %     otherwise
    %	    fprintf('Unknown type: %s', morphData.type(i))
    %   end
    % end

    % Locate the parent
    idx = find(morphData.id == morphData.parent(i));
    p = [idx i];

    switch(morphData.type(i))
      case 1
        % Already handled above, ignore.
      case 2
        %plot3(morphData.x(p), morphData.y(p), morphData.z(p),'r');
        makeCyl(morphData.x(p), morphData.y(p), morphData.z(p), ...
		morphData.r(i),1);
      case 3
        %plot3(morphData.x(p), morphData.y(p), morphData.z(p),'b');
        makeCyl(morphData.x(p), morphData.y(p), morphData.z(p), ...
		morphData.r(i),0.5);
      case 4
        %plot3(morphData.x(p), morphData.y(p), morphData.z(p),'k');
        makeCyl(morphData.x(p), morphData.y(p), morphData.z(p), ...
		morphData.r(i),0);


      otherwise
      fprintf('Unknown type: %d\n', morphData.type(i))
    end

    hold on

  end

  if(labelGCflag)
    % Label the end points

    endPoints = setdiff(morphData.id,morphData.parent);
    for i = 1:length(endPoints)
      text(morphData.x(endPoints(i)), ...
  	   morphData.y(endPoints(i)), ...
  	   morphData.z(endPoints(i)), ...
	   num2str(i-1),'fontsize',15)
    end

    xlabel('x')
    ylabel('y')
    zlabel('z')

  else
    box off
  end

  axis equal

  shading flat
  % lighting phong

end
