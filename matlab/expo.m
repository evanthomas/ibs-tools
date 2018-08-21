% For Henrik: start with e..g. expo('AFALK_duo.dat')

function expo(filename)
    % EXPO(FILENAME) fit and exponential to decimated PD data
    
    
    if ~exist('filename', 'var')
        %         [filename, path] = uigetfile('*.dat', 'Select a PD/pressure file')
        filename = 'acast1_jej.dat';
        path = '/Users/thomas.e/neurosim/ibs-tools/test-data/';
%         path = 'C:\Users\evan\Documents\neurosim\ibs-tools\test-data\';
    end
    
    [pdAxes, pressAxes] = setup;
    
    junk = load(strcat(path, filename));
    
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
        
        setRiseLine(theseAxes, signal, sigName)
        setExpLine(theseAxes, signal, sigName)
        
        
        % Plot original data
        line(theseAxes, tm*5, signal, 'Color', 'blue');
        meanSignal = mean(signal(24:end-24));
        
        xt = 20;
        dy = maxSignal/12;
        yt = maxSignal*0.7 + dy;
        
        % Display time to peak
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('time to %s peak = %g', sigName, t2pSignal), ...
            'fontsize',12,'FontWeight','bold');
        
        % Display max PD
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('Max %s = %g', sigName, maxSignal), ...
            'fontsize',12,'FontWeight','bold');
        
        % Display mean pressure
        yt = yt + dy;
        text(theseAxes, xt, yt, ...
            sprintf('mean %s = %g', sigName, meanSignal), ...
            'fontsize',12,'FontWeight','bold');
        
    end
    
    function setExpLine(theseAxes, signal, sigName)
        [maxSignal, indx] = max(signal);
        % Initial estimate for lambda
        y = signal(indx:end)/maxSignal + 1e-3;
        t = 1:length(y);
        logstart = exp(polyfit(t',log(y),1));
        lam = logstart(2);
        % Refine fit
        Options = optimset('TolX', 1e-3, 'Display', 'off');
        lambda = fminsearch(@expfit, lam, Options, t', y);
        thalf = log(2)/lambda*5;
        
        % exponential fit
        t = [0 t];
        expFit = exp(-lambda*t)*maxSignal;
        
        line(theseAxes, (t+indx)*5, expFit, 'Color', 'red');
        
        xt = 20;
        dy = maxSignal/12;
        yt = maxSignal*0.7;% + 2*dy;
        
        text(theseAxes, xt, yt, ...
            sprintf('t_1_/_2 = %g', thalf), ...
            'fontsize',12,'FontWeight','bold');
        
        % Add the begin and end drag lines
        ylim = get(theseAxes, 'YLim');
        tstart = (t(1)+indx)*5;
        tend   = (t(end)+indx)*5;
        line(theseAxes, ...
            [tstart tstart], ylim, ...
            'color', 'black', ...
            'LineWidth', 2 ...
            );
        line(theseAxes, ...
            [tend tend], ylim, ...
            'color', 'black', ...
            'LineWidth', 2 ...
            );
    end
    
    function setRiseLine(theseAxes, signal, sigName)
        [maxSignal, indx] = max(signal);
        
        % Linear regression on initial slope
        tl = (24:indx)';
        signall = signal(24:indx);
        X = [ones(size(tl)) tl];
        a = X\signall;
        
        % Time for signal to return with 10% of basline
%         x = signal(indx:end);
%         signalRet = min(find(x<0.1*maxSignal));
%         if isempty(signalRet)
%             signalRet = (length(signal)-indx)*5;
%         else
%             signalRet = signalRet*5;
%         end
        
        % Linear fit
        tm = [24 indx]*5;
        signalLine = a(1) + a(2)*[24 indx];
        drawRiseLine(theseAxes, signalLine, tm, sigName)
    end
    
    function drawRiseLine(theseAxes, signalLine, tm, sigName)
        RoR = (signalLine(2)-signalLine(1))/(tm(2)-tm(1));
        
        line(theseAxes, ...
            tm, signalLine, ...
            'Color', 'red');
        
        x = tm(1);
        y = signalLine(1);
        ud1.hoverFunction = {@riseLineHoverHandler, [x, y], signalLine, 'end'};
        line(theseAxes, ...
            x, y, ...
            'Color', 'red', ...
            'MarkerFaceColor', 'red', ...
            'Marker', 'square', ...
            'MarkerSize', 8, ...
            'UserData', ud1);
        
        x = tm(end);
        y = signalLine(end);
        ud2.hoverFunction = {@riseLineHoverHandler, [x, y], signalLine, 'start'};
        line(theseAxes, ...
            x, y, ...
            'Color', 'red', ...
            'MarkerFaceColor', 'red', ...
            'Marker', 'square', ...
            'MarkerSize', 8, ...
            'UserData', ud2);
        
        
        % Display rate of rise
        ylim = get(theseAxes, 'YLim');
        xt = 20;
        dy = ylim(2)/12;
        yt = ylim(2)*0.7 + dy;
        text(theseAxes, xt, yt, ...
            sprintf('Rate of rise = %g', RoR), ...
            'fontsize',12,'FontWeight','bold');
        
    end
    
    function riseLineClickHandler(o, evt, signal, detail)
        set(gcf, 'WindowButtonMotionFcn', {@riseLineDragHandler, signal, detail});        
    end
    
    function riseLineDragHandler(o, evt, signal, detail)
        get(gca, 'CurrentPoint')
    end
    
    function setDefaultHandlers
        set(gcf, ...
            'WindowButtonMotionFcn', @mousemover, ...
            'CloseRequestFcn', {@expoQuit, gcf}, ...
            'WindowButtonDownFcn', '' ...
            );
    end
    
    function accepted = riseLineHoverHandler(args)            
        delta = 1.5;        
        pnt = args{1};
        signal= args{2};
        detail = args{3};
        ptr = args{4};
        
        x = ptr(1);
        y = ptr(2);
        xlo = pnt(1) - delta;
        xhi = pnt(1) + delta;
        ylo = pnt(2) - delta;
        yhi = pnt(2) + delta;
        
        if xlo<=x & x<=xhi & ylo<=y & y<=yhi
            set(gcf, 'Pointer', 'hand')
            set(gcf, 'WindowButtonDownFcn', {@riseLineClickHandler, signal, detail});
            accepted = true;
        else
            setDefaultHandlers();
            accepted = false;
        end            
        
    end
    
    function mousemover(o, evt, hoverTargets)
        
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
        for i=1:length(gs)
            g = gs(i);
            ud = get(g, 'UserData');
            if ~isempty(ud)
                if isfield(ud, 'hoverFunction')
                    f = ud.hoverFunction{1};
                    args = {ud.hoverFunction{2:end}};
                    args{end+1} = p;
                    if f(args)
                        return
                    end
                end
            end
        end
        
        setptr(gcf, 'arrow');
        
    end
    
    function expoQuit(o, crap, h)
        closereq
    end
    
    function saveToExcel(o, crap, h)
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
        
        [i,j] = find(isnan(pts));
        bool = (pts(:,1)<rect(1))|(pts(:,1)>rect(2))|...
            (pts(:,2)<rect(3))|(pts(:,2)>rect(4));
        bool = ~bool;
        bool(i) = 0;
    end
    
    function [pdAxes, pressAxes] = setup(hoverTargets)
        UI_NAME = 'Pulse Analyzer';
        try
            delete(findobj('Name', UI_NAME))
        catch
            % Nothing
        end
        
        ud.mode = 'move';
        h1 = figure( ...
            'Units', 'characters', ...
            'Color',[0.8 1.0 0.8], ...
            'MenuBar', 'none', ...
            'Name', UI_NAME, ...
            'NumberTitle', 'off', ...
            'Position', [36.6667   15.5000  140   50], ...
            'UserData', ud, ...
            'toolbar',  'figure');
        
        pdAxes = axes( ...
            'Position', [0.03 0.084016393442623 0.96 0.379098360655738], ...
            'XminorTick', 'on', ...
            'XminorGrid', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Tag', 'pd axes');
        
        pressAxes = axes( ...
            'Position', [0.03 0.532786885245902 0.96 0.381147540983607], ...
            'XminorTick', 'on', ...
            'XminorGrid', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Tag', 'press axes');
        
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
            'callback', {@saveToExcel, h1});
        
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