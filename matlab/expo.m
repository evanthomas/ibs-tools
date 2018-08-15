% For Henrik: start with e..g. expo('AFALK_duo.dat')

function expo(filename)
    % EXPO(FILENAME) fit and exponential to decimated PD data
    
    
    if ~exist('filename', 'var')
        %     [filename, path] = uigetfile('*.dat', 'Select a PD/pressure file')
        filename = 'acast1_jej.dat';
        path = '/Users/thomas.e/neurosim/ibs-tools/test-data/';
    end
    
    [pdAxes, pressAxes] = setup;
    
    junk = load(strcat(path, filename));
    
    PD = -junk(:,3);
    startAnalysis(PD, pdAxes, 'PD');
    
    
%     % Pressure
%     press = junk(:,2);
%     press = press - min(press);
%     press = press/max(press)*maxPD;
%     meanP = mean(press(24:end-24));
    
    
    
%     line(tm*5, press, 'Color', 'blue');
    
    %  legend('PD', 'pressure', 'exp fit')
    %title(sprintf('%s', texescape(n)))
        
end

function startAnalysis(signal, theseAxes, sigName)
    signal = signal - min(signal);
    meanSignal = mean(signal(24:end-24));
    [maxSignal, indx] = max(signal);
    y = signal(indx:end)/maxSignal + 1e-3;
    
    % Time to signal peak
    t2pSignal = (indx-24)*5;
    % Initial estimate for lambda
    t = 1:length(y);
    logstart = exp(polyfit(t',log(y),1));
    lam = logstart(2);
    % Refine fit
    Options = optimset('TolX', 1e-3, 'Display', 'off');
    lambda = fminsearch(@expfit, lam, Options, t', y);
    thalf = log(2)/lambda*5;
    
    % exponential fit
    t = [0 t];
    fit = exp(-lambda*t)*maxSignal;
    tm = 0:length(signal)-1;
    
    % Linear regression on initial slope
    tl = (24:indx)';
    signall = signal(24:indx);
    X = [ones(size(tl)) tl];
    a = X\signall;
    
    % Phase III duration
    PIIIdur = (length(signal)-48)*5;
    
    % Time for signal to return with 10% of basline
    x = signal(indx:end);
    signalRet = min(find(x<0.1*maxSignal));
    if isempty(signalRet)
        signalRet = (length(signal)-indx)*5;
    else
        signalRet = signalRet*5;
    end
    
    
    % Linear fit
    signalLength = a(1) + a(2)*tl;
    RoR = a(2)/5;
    
    % Plot results
    line(theseAxes, tm*5, signal, 'Color', 'blue');
    line(theseAxes, (t+indx)*5, fit, 'Color', 'red');
    line(theseAxes, tl*5, signalLength, 'Color', 'red');
    
    % Display t_1/2
    xt = 20;
    yt = maxSignal*0.7;
    dy = maxSignal/12;
    text(theseAxes, xt, yt, ...
        sprintf('t_1_/_2 = %g', thalf), ...
        'fontsize',12,'FontWeight','bold');
    
    % Display time to peak
    yt = yt + dy;
    text(theseAxes, xt, yt, ...
        sprintf('time to %s peak = %g', sigName, t2pSignal), ...
        'fontsize',12,'FontWeight','bold');
    
    % Display rate of rise
    yt = yt + dy;
    text(theseAxes, xt, yt, ...
        sprintf('Rate of rise = %g', RoR), ...
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


function mousemover(o, crap)
end

function expoQuit(o, crap, h)
    closereq
end

function saveToExcel(o, crap, h)
end


function [pdAxes, pressAxes] = setup
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
    set(h1, ...
        'WindowButtonMotionFcn', @mousemover, ...
        'CloseRequestFcn', {@expoQuit, h1} ...
        );
    
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
