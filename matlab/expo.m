% For Henrik: start with e..g. expo('AFALK_duo.dat')

function expo(filename)
    % EXPO(FILENAME) fit and exponential to decimated PD and pressure data
    
    
    if ~exist('filename', 'var')
        [filename, path] = uigetfile('*.dat', 'Select a PD/pressure file')
%         path = 'C:\Users\evan\Documents\neurosim\ibs-tools\test-data\';
%         filename = 'acast1_jej.dat';
%         path = '/Users/thomas.e/neurosim/ibs-tools/test-data/';
        [~, filename, ext] = fileparts(filename);
    else
        [path, filename, ext] = fileparts(filename);
    end
    
    [pdAxes, pressAxes] = setup;
    
    ud.path = path;
    ud.filename = filename;
    ud.ext = ext;
    set(gcf, 'UserData', ud);
    
    junk = load(strcat(path, filename, ext));
    
    PD = -junk(:,3);
    startAnalysis(PD, pdAxes, 'PD');
    
    press = junk(:,2);
    startAnalysis(press, pressAxes, 'pressure');

    % Set up mouse movement tracker for the mouse hover
    setDefaultHandlers()
    
    function startAnalysis(signal, theseAxes, sigName)
        signal = signal - min(signal);
        tm = 0:length(signal)-1;
        [maxSignal, indx] = max(signal);
        
        % Time to signal peak
        t2pSignal = (indx-24)*5;        
        
        % Plot original data
        line(theseAxes, tm*5, signal, 'Color', 'blue');
        meanSignal = mean(signal(24:end-24));
        
        % Add the analysis lines.
        % Note: This needs to be done after the data are drawn
        %       and before text us added, so that the axes are 
        %       fully initialised.
        setRiseLine(theseAxes, signal)
        setExpLine(theseAxes, signal)
        
        ylim = get(theseAxes, 'YLim');
        xt = 20;
        dy = ylim(2)/12;
        yt = ylim(2)*0.6 + 1*dy;
        
        % Display time to peak
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('time to %s peak = %g', sigName, t2pSignal), ...
            'fontsize',12,'FontWeight','bold');
        saveDatum(theseAxes, t2pSignal, 'timeToPeak');
        
        % Display max PD
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('Max %s = %g', sigName, maxSignal), ...
            'fontsize',12,'FontWeight','bold');        
        saveDatum(theseAxes, maxSignal, 'max');
        
        % Display mean pressure
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('mean %s = %g', sigName, meanSignal), ...
            'fontsize',12,'FontWeight','bold');
        saveDatum(theseAxes, meanSignal, 'mean');
    end
    
    function [t, expFit, thalf] = makeExpFit(istart, iend, signal)
        maxSignal = max(signal(istart:iend));
        % Initial estimate for lambda
        y = signal(istart:iend)/maxSignal + 1e-3;
        t = (1:length(y))-1;
        logstart = exp(polyfit(t',log(y),1));
        lam = logstart(2);
        % Refine fit
        Options = optimset('TolX', 1e-3, 'Display', 'off');
        lambda = fminsearch(@expfit, lam, Options, t', y);
        thalf = log(2)/lambda*5;
        
        % exponential fit
        expFit = exp(-lambda*t)*maxSignal;
    end
    
    function setExpLine(theseAxes, signal)
        [~, indx] = max(signal);
        [t, expFit, thalf] = makeExpFit(indx, length(signal), signal);
        
        line(theseAxes, ...
            (t+indx-1)*5, expFit, ...
            'Color', 'red', ...
            'LineWidth', 2, ...
            'tag', 'exp line');
        
        ylim = get(theseAxes, 'YLim');
        xt = 20;
        yt = ylim(2)*0.6;
        
        text(theseAxes, xt, yt, ...
            sprintf('t_1_/_2 = %g', thalf), ...
            'fontsize',12, ...
            'FontWeight','bold', ...
            'tag', 'halflife text');                
        saveDatum(theseAxes, thalf, 'decayHalflife');

        % Add the begin and end drag lines
        tstart = (t(1)+indx-1)*5;
        tend   = (t(end)+indx-1)*5;
        l = line(theseAxes, ...
            [tstart tstart], ylim, ...
            'color', 'black', ...
            'LineWidth', 2, ...
            'tag', 'drag line');
        ud1.hoverFunction = {@expLineHoverHandler, 2, l, signal};
        set(l, 'UserData', ud1);
        
        l = line(theseAxes, ...
            [tend tend], ylim, ...
            'color', 'black', ...
            'LineWidth', 2, ...
            'tag', 'drag line');
        ud2.hoverFunction = {@expLineHoverHandler, 2, l, signal};
        set(l, 'UserData', ud2);
    end
    
    function setRiseLine(theseAxes, signal)
        [~, indx] = max(signal);
        
        % Linear regression on initial slope
        tl = (24:indx)';
        signall = signal(24:indx);
        X = [ones(size(tl)) tl];
        a = X\signall;
        
        % Linear fit
        tm = [24 indx]*5 - 1*5; % Subtract time zero offset
        signalLine = a(1) + a(2)*[24 indx];
        drawRiseLine(theseAxes, signalLine, tm)
    end
    
    function drawRiseLine(theseAxes, signalLine, tm)
        RoR = (signalLine(2)-signalLine(1))/(tm(2)-tm(1));
        
        line(theseAxes, ...
            tm, signalLine, ...
            'Color', 'red', ...
            'LineWidth', 2, ...
            'tag', 'rise line data');
        
        x = tm(1);
        y = signalLine(1);
        l = line(theseAxes, ...
            x, y, ...
            'Color', 'red', ...
            'MarkerFaceColor', 'red', ...
            'Marker', 'square', ...
            'MarkerSize', 8, ...
            'tag', 'rise line start');
        ud1.hoverFunction = {@riseLineHoverHandler, 0, l};
        set(l, 'UserData', ud1)
        
        x = tm(end);
        y = signalLine(end);
        l = line(theseAxes, ...
            x, y, ...
            'Color', 'red', ...
            'MarkerFaceColor', 'red', ...
            'Marker', 'square', ...
            'MarkerSize', 8, ...
            'tag', 'rise line end');
        ud2.hoverFunction = {@riseLineHoverHandler, 0, l};
        set(l, 'UserData', ud2)
        
        
        % Display rate of rise
        ylim = get(theseAxes, 'YLim');
        xt = 20;
        dy = ylim(2)/12;
        yt = ylim(2)*0.6 + dy;
        text(theseAxes, ...
            xt, yt, ...
            sprintf('Rate of rise = %g', RoR), ...
            'fontsize', 12, ...
            'FontWeight','bold', ...
            'tag', 'rate of rise');
        
        saveDatum(theseAxes, RoR, 'rateOfRise');
        
    end
    
    function saveDatum(axes, val, field)
        ud = get(gcf, 'UserData');
        type = get(axes, 'Tag');
        ud.(strcat(field, type)) = val;
        set(gcf, 'UserData', ud);
    end
    
    function riseLineClickHandler(~, ~, l)
        set(gcf, ...
            'WindowButtonMotionFcn', {@riseLineDragHandler, l}, ...
            'WindowButtonUpFcn', @setDefaultHandlers ...
            );
    end
    
    function riseLineDragHandler(~, ~, lhandle)
        cp = get(gca, 'CurrentPoint');
        x = cp(1, 1);
        y = cp(1, 2);
        
        ldata = findobj(gca, 'tag', 'rise line data');
        xdata = get(ldata, 'XData');
        ydata = get(ldata, 'YData');
        
        detail = get(lhandle, 'tag');
        switch detail
            case 'rise line start'
                xstart = x;
                ystart = y;
                xend = xdata(end);
                yend = ydata(end);
            case 'rise line end'
                xstart = xdata(1);
                ystart = ydata(1);
                xend = x;
                yend = y;
        end
                
       set(lhandle, ...
           'XData', x, ...
           'YData', y);
       set(ldata, ...
           'XData', [xstart, xend], ...
           'YData', [ystart, yend]);
       
       RoR = (yend-ystart)/(xend-xstart);
       th = findobj(gca, 'tag', 'rate of rise');
       set(th, 'String', sprintf('Rate of rise = %g', RoR));
       saveDatum(gca, RoR, 'rateOfRise');
    end
    
    function setDefaultHandlers(varargin)
        set(gcf, ...
            'WindowButtonMotionFcn', @mousemover, ...
            'CloseRequestFcn', {@expoQuit, gcf}, ...
            'WindowButtonDownFcn', '', ...
            'WindowButtonUpFcn', '' ...
            );
        setptr(gcf, 'arrow');
    end
    
    function accepted = riseLineHoverHandler(args)            
        delta = 1.5;
        l = args{1};
        ptr = args{2};
        pnt = [get(l, 'XData'), get(l, 'YData')];
        
        x = ptr(1);
        y = ptr(2);
        xlo = pnt(1) - delta;
        xhi = pnt(1) + delta;
        ylo = pnt(2) - delta;
        yhi = pnt(2) + delta;
        
        if xlo<=x && x<=xhi && ylo<=y && y<=yhi
            set(gcf, 'Pointer', 'hand')
            set(gcf, 'WindowButtonDownFcn', {@riseLineClickHandler, l});
            accepted = true;
        else
            setDefaultHandlers();
            accepted = false;
        end            
        
    end
    
    function accepted = expLineHoverHandler(args)            
        delta  = 1.5;
        l      = args{1};
        signal = args{2};
        ptr    = args{3};
        pnt    = [get(l, 'XData'), get(l, 'YData')];
        
        x   = ptr(1);
        xlo = pnt(1) - delta;
        xhi = pnt(1) + delta;
        
        if xlo<=x && x<=xhi
            setptr(gcf, 'lrdrag')
            set(gcf, 'WindowButtonDownFcn', {@expLineClickHandler, l, signal});
            accepted = true;
        else
            setDefaultHandlers();
            accepted = false;
        end            
        
    end
    
    function expLineClickHandler(~, ~, l, signal)
        set(gcf, ...
            'WindowButtonMotionFcn', {@expLineDragHandler, l, signal}, ...
            'WindowButtonUpFcn', @setDefaultHandlers ...
            );        
    end
    
    function expLineDragHandler(~, ~, dragLine, signal)
        cp = get(gca, 'CurrentPoint');
        x = cp(1, 1);
        tmax = (length(signal)-1)*5;
        if x>tmax
            x = tmax;
        end
        if x<0
            x = 0;
        end
        
        set(dragLine, 'XData', [x, x]);
        
        dragLines = findobj(gca, 'tag', 'drag line');
        x = get(dragLines(1), 'XData');
        xlo = x(1);
        x = get(dragLines(2), 'XData');
        xhi = x(1);
        
        if xlo>xhi
            x = xlo;
            xlo = xhi;
            xhi = x;
        end
                
        istart = round(xlo/5)+1;
        iend   = round(xhi/5)+1;
        [t, expFit, thalf] = makeExpFit(istart, iend, signal);
        expLine = findobj(gca, 'tag', 'exp line');
        set(expLine, ...
            'XData', (t+istart-1)*5, ...
            'YData', expFit);
        
        th = findobj(gca, 'tag', 'halflife text');
        set(th, 'String', sprintf('t_1_/_2 = %g', thalf));
        saveDatum(gca, thalf, 'decayHalflife');
    end
    
    function mousemover(o, ~)
        
        % Determine which (if any) axes the mouse is over
        axs = findobj(o, 'type', 'axes')';
        inaxes = 0;
        for ax=axs
            %subplot(ax)
            p = get(ax, 'CurrentPoint');
            p = p(1, 1:2);
            xlim = get(ax, 'xlim');
            ylim = get(ax, 'ylim');
            if pinrect(p, [xlim ylim])  % inside an axes
                inaxes = 1;
                break
            end
        end

        % Not in any axes
        if ~inaxes
            setptr(gcf, 'arrow');
            return
        end

        
        gs = findall(ax);
        n = length(gs);
        priorities = zeros(1, n);
        handlers = cell(1, n);
        args = cell(1, n);
        indx = 1;
        for i=1:n
            g = gs(i);
            ud = get(g, 'UserData');
            if ~isempty(ud)
                if isfield(ud, 'hoverFunction')
                    handlers{indx} = ud.hoverFunction{1};
                    priorities(indx) = ud.hoverFunction{2};
                    args{indx} = {ud.hoverFunction{3:end}, p};
                    indx = indx + 1;
                end
            end
        end

        indx = indx - 1;
        [~, order] = sort(priorities(1:indx));
        for i=1:indx
            ii = order(i);
            f = handlers{ii};
            arg = args{ii};
            if f(arg)
                return
            end            
        end
        
        setptr(gcf, 'arrow');
        
    end
    
    function expoQuit(~, ~, ~)
        closereq
    end
    
    function saveToExcel(~, ~)
        ud = get(gcf, 'UserData');
        fh = fopen(strcat(ud.path, ud.filename, '.csv'), 'w');
        
        try
            fprintf(fh, '"mean pressure", "max pressure", "time to peak pressure", "rate of pressure rise", "pressure decay half life", "mean PD", "max PD", "time to peak PD", "rate of PD rise", "PD decay half life"\n');
            fprintf(fh, '%g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n', ...
            ud.meanPressure, ...
            ud.maxPressure, ...
            ud.timeToPeakPressure, ...
            ud.rateOfRisePressure, ...
            ud.decayHalflifePressure, ...
            ud.meanPD, ...
            ud.maxPD, ...
            ud.timeToPeakPD, ...
            ud.rateOfRisePD, ...
            ud.decayHalflifePD ...
            );
        catch ME
            disp('Error writing CSV data:')
            disp(ME)
        end

        fclose(fh);
    end
    
    function bool = pinrect(pts,rect)
        %PINRECT Determine if points lie in or on rectangle.
        %   Inputs:
        %     pts - n-by-2 array of [x,y] data
        %     rect - 1-by-4 vector of [xlim ylim] for the rectangle
        %   Outputs:
        %     bool - length n binary vector
        
        %   Copyright 1988-2002 The MathWorks, Inc.
        % $Revision: 1.8 $
        
        [i,~] = find(isnan(pts));
        bool = (pts(:,1)<rect(1))|(pts(:,1)>rect(2))|...
            (pts(:,2)<rect(3))|(pts(:,2)>rect(4));
        bool = ~bool;
        bool(i) = 0;
    end
    
    function [pdAxes, pressAxes] = setup
        UI_NAME = 'Pulse Analyzer';
        try
            delete(findobj('Name', UI_NAME))
        catch
            % Nothing
        end
        
        h1 = figure( ...
            'Units', 'characters', ...
            'Color',[0.8 1.0 0.8], ...
            'MenuBar', 'none', ...
            'Name', UI_NAME, ...
            'NumberTitle', 'off', ...
            'Position', [36.6667   15.5000  140   50], ...
            'toolbar',  'figure');
        
        pdAxes = axes( ...
            'Position', [0.03 0.084016393442623 0.96 0.379098360655738], ...
            'XminorTick', 'on', ...
            'XminorGrid', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Tag', 'PD');
        
        pressAxes = axes( ...
            'Position', [0.03 0.532786885245902 0.96 0.381147540983607], ...
            'XminorTick', 'on', ...
            'XminorGrid', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Tag', 'Pressure');
        
        uicontrol( ...
            'Style', 'text', ...
            'Parent', h1, ...
            'Units', 'normalized', ...
            'BackgroundColor', [0.8 1 0.8], ...
            'FontWeight', 'bold', ...
            'FontSize', 14, ...
            'Position', [0.32761087267525 0.47 0.353361945636624 0.03], ...
            'String', 'Potential Difference');
        
        uicontrol( ...
            'Style', 'text', ...
            'Parent', h1, ...
            'Units', 'normalized', ...
            'BackgroundColor', [0.8 1 0.8], ...
            'FontSize', 14, ...
            'FontWeight', 'bold', ...
            'Position', [0.313304721030043 0.92 0.386266094420601 0.04], ...
            'String', 'Pressure');
        
        
        uicontrol( ...
            'Parent', h1, ...
            'Units', 'normalized', ...
            'BackgroundColor', [0.996108949416342 0.996108949416342 0], ...
            'FontWeight', 'bold', ...
            'FontSize', 12, ...
            'ForegroundColor', 'black', ...
            'Position', [0.89 0.003 0.1 0.04], ...
            'String', 'save to Excel', ...
            'callback', @saveToExcel);
        
    end
    
    function err = expfit(lambda, t, y)
        
        A = zeros(length(t),length(lambda));
        for j = 1:length(lambda)
            A(:,j) = exp(-lambda(j)*t);
        end
        c = A\y;
        z = A*c;
        %set(Plothandle,'ydata',z)
        %drawnow
        err = norm(z-y);
    end
    
end