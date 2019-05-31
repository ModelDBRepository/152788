% Plots data from iGrow read with readData.m
% plotTime is a time, like 0, or 'end' for last time point.
% lineStyle and lineWidth specifies how the plot should look
% markGC marks a growth cone with a star, numbered 1,2,3...

function makeDendrogram(data,plotTime,lineStyle,lineWidth,markGC,lineColor)

  %close all
  format short g

  simplifyTree = true;

  if(~exist('lineStyle'))
    lineStyle = 'k-';
  end

  if(~exist('lineWidth'))
    lineWidth = 1;
  end

  if(~exist('plotTime'))
    plotTime = 0;
  elseif(strcmp(plotTime,'end'))
    plotTime = max(data.time);
  else
    % Find closest time point
    t = unique(data.time);
    at = abs(t - plotTime);
    plotTime = t(find(at == min(at),1));
  end

  fprintf('Showing plot time %d\n', plotTime)

  nodeSpacing = 1;

  idx = find(data.time == plotTime);

  % Create a new tree representation

  % 1. Find soma
  somaIdx = idx(find(data.parentID(idx) == -1));

  if(somaIdx ~= idx(1))
    disp('We assumed soma compartment was first... was not the case.')
    keyboard
  end

  reducedTree.parent(1)    = -1;
  reducedTree.children{1}  = []; 
  reducedTree.startDist(1) = 0;

  somaChildrenIdx = idx(find(data.parentID(idx) == data.ID(somaIdx)));
  nextIdx = 2; 

  % This assumes soma is first one
  reducedTree.oldID = data.ID(idx);
  reducedTree.endDist = data.dist(idx);
  reducedTree.plotY = NaN*zeros(size(idx));

  if(strcmp(lineStyle,'gradient'))
    reducedTree.conc = data.tubulinConc(idx);

    reducedTree.colour = zeros(length(idx),3);
    colMap = colormap('winter');

    minConc = min(reducedTree.conc);
    maxConc = max(reducedTree.conc);

    for j = 1:3 
      reducedTree.colour(:,j) = interp1(linspace(minConc, maxConc,64), ...
					colMap(:,j), ...
					reducedTree.conc);
      % We got rounding errors leading to one colourvalue being 
      % 1e-15 larger than 1
      reducedTree.colour(:,j) = min(1,reducedTree.colour(:,j));
    end

    if(nnz(isnan(reducedTree.colour)))
      disp('Found NaN colour')
      keyboard
    end

    simplifyTree = false;

  end


  if(lineWidth > 0)
    reducedTree.lineWidth = lineWidth*ones(size(idx));
  else
    % Negative values are used to scale (exclude soma radius)
    reducedTree.lineWidth = data.r(idx)*abs(lineWidth)/max(data.r(idx(2:end)));
  end

  for i = 1:length(idx)
    reducedTree.children{i} = [];
  end

  % Loop through all branch compartments

  for i = 2:length(idx)
    % Find parent compartment in the reduced tree
    pIdx = find(reducedTree.oldID == data.parentID(idx(i)));

    reducedTree.parent(i) = pIdx;

    reducedTree.children{pIdx}(end+1) = i;
    reducedTree.startDist(i) = reducedTree.endDist(pIdx);

  end

  GCidx = [];

  % This is a lookup table to mark the right growth cone
  for i = 1:length(reducedTree.children)
    if(isempty(reducedTree.children{i}))
      GCidx(end+1) = i;
    end
  end

  if(simplifyTree)

    % Reduce the tree, merging compartments that are part of the same
    % branch segment, and  marking redundant compartments for deletion

    for i = 1:length(reducedTree.children)
      if(numel(reducedTree.children{i}) == 1)
        % Only one child, merge this compartment with the child
        pIdx = reducedTree.parent(i);

        if(pIdx == -1)
          % The soma was the parent, do not merge
          continue
        end

        cIdx = reducedTree.children{i};

        % keyboard

        % Update the line width
        cLen = reducedTree.endDist(cIdx)-reducedTree.startDist(cIdx);
        pLen = reducedTree.endDist(pIdx)-reducedTree.startDist(pIdx);

        reducedTree.lineWidth(cIdx) = ...
            (reducedTree.lineWidth(cIdx)*cLen+reducedTree.lineWidth(pIdx)*pLen)...
            / (cLen + pLen);

        % Bypass the current compartment
        j = find(reducedTree.children{pIdx} == i);
        reducedTree.children{pIdx}(j) = cIdx;

        reducedTree.parent(cIdx) = pIdx;
        reducedTree.startDist(cIdx) = reducedTree.startDist(i);


        % Mark the old compartment as removed
        reducedTree.parent(i) = NaN;
        reducedTree.children{i} = NaN;
        reducedTree.startDist(i) = NaN;
        reducedTree.endDist(i) = NaN;
        reducedTree.oldID(i) = NaN;
        reducedTree.plotY(i) = NaN;
      end
    end
  end

  % Do a depth first search to number the leafs

  neuronStack = 1;
  leafCtr = 1;

  while(~isempty(neuronStack))

    curSegment = neuronStack(end);

    if(~isnan(reducedTree.plotY(curSegment)))
      % This one is done, go up to parent
      parSegment = reducedTree.parent(curSegment);

      if(parSegment == -1)
        neuronStack(end) = [];
        continue
      end

      childIdx = find(reducedTree.children{parSegment} == curSegment);

      if(childIdx < length(reducedTree.children{parSegment}))
        % Go to the next child
        neuronStack(end) = reducedTree.children{parSegment}(childIdx+1);
      else
        % Last child
        if(~isnan(reducedTree.plotY(parSegment)))
          disp('aargh')
          keyboard
        end

        reducedTree.plotY(parSegment) = ...
            mean(reducedTree.plotY(reducedTree.children{parSegment}));

        neuronStack(end) = [];        
      end
    elseif(isempty(reducedTree.children{curSegment}))
      % No children, assign plotY and increment leaf counter

      if(~isnan(reducedTree.plotY(curSegment)))
        disp('aargh2')
        keyboard
      end

      reducedTree.plotY(curSegment) = leafCtr;
      leafCtr = leafCtr + 1;
    else
      neuronStack(end+1) = reducedTree.children{curSegment}(1);
      % Go through the children
    end

  end


  disp('Plot time')

  for i = 1:length(reducedTree.children)
    if(~isnan(reducedTree.oldID(i)))

      if(strcmp(lineStyle,'gradient'))

        % Plot the concentration gradient in the histogram
        plot([reducedTree.startDist(i) reducedTree.endDist(i)]*1e6, ...
             reducedTree.plotY(i)*[1 1], ...
             'color',reducedTree.colour(i,:), ...
             'linewidth',reducedTree.lineWidth(i));

        hold on

        pIdx = reducedTree.parent(i);

        if(pIdx > 0)
          plot([reducedTree.endDist(pIdx) reducedTree.startDist(i)]*1e6, ...
	       [reducedTree.plotY(pIdx) reducedTree.plotY(i)], ...
	       'color', reducedTree.colour(i,:), ...
	       'linewidth',reducedTree.lineWidth(i));
	    
        end

      elseif(strcmp(lineStyle,'growthcones'))

        hold on
        
        if(~exist('lineColor'))
          lineColor = [0 0 0];
        end

        for j = 1:numel(GCidx)
          x = reducedTree.endDist(GCidx(j))*[1 1]*1e6;
          y = reducedTree.plotY(GCidx(j)) + [-0.3 0.3];
          plot(x,y,'color',lineColor);
        end
        
      else
        if(exist('lineColor'))
          % User defined colour
          plot([reducedTree.startDist(i) reducedTree.endDist(i)]*1e6, ...
	       reducedTree.plotY(i)*[1 1], ...
	       lineStyle, ...
	       'linewidth',reducedTree.lineWidth(i), ...
	       'color',lineColor);
          hold on

          pIdx = reducedTree.parent(i);

          if(pIdx > 0)
            plot([reducedTree.endDist(pIdx) reducedTree.startDist(i)]*1e6, ...
 	         [reducedTree.plotY(pIdx) reducedTree.plotY(i)], ...
                 lineStyle, ...
                 'linewidth',reducedTree.lineWidth(i), ...
                 'color',lineColor)
            
          end
        else
          % Use simple colours
          plot([reducedTree.startDist(i) reducedTree.endDist(i)]*1e6, ...
               reducedTree.plotY(i)*[1 1], ...
               lineStyle,'linewidth',reducedTree.lineWidth(i));
          hold on

          pIdx = reducedTree.parent(i);

          if(pIdx > 0)
            plot([reducedTree.endDist(pIdx) reducedTree.startDist(i)]*1e6, ...
 	         [reducedTree.plotY(pIdx) reducedTree.plotY(i)], ...
	         lineStyle,'linewidth',reducedTree.lineWidth(i))

          end
        end
      end

    end
  end

% keyboard

  if(exist('markGC') & ~isnan(markGC))
    plot(reducedTree.endDist(GCidx(markGC))*1e6 + 15, ...
         reducedTree.plotY(GCidx(markGC)),'k*')
  end

  set(gca,'ytick',[])
  axis tight
  a = axis;
  a(1)=-10;
  a(2) = a(2)+10;
  a(3) = -1;
  a(4) = a(4)+1;
  axis(a);
  box off

  if(strcmp(lineStyle,'gradient'))
    c = colorbar;
    yTick = linspace(1,64,5);
    yTickLabel = {};
    yTickLabelVal = linspace(minConc,maxConc,numel(yTick));
    for i = 1:length(yTick)
      % We get it in microMolar
      yTickLabel{i} = sprintf('%.2f',1e3*yTickLabelVal(i));
    end

    set(c,'ytick',yTick,'yticklabel',yTickLabel);

  end

  xlabel('Distance from soma (\mum)','fontsize',20)
  set(gca,'fontsize',20)

  %keyboard


end
