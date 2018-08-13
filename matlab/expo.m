% For Henrik: start with e..g. expo('AFALK_duo.dat')

function expo(filename)
% EXPO(FILENAME) fit and exponential to decimated PD and pressure data

%
% Constants
% 
  DT = 0.2; % sample period seconds
  ROR_START = 24; % Samples from peak to measure rate of rise
  EPS = 1e-3; % Small number to avoid divide/3 issues

  [pd_axes, press_axes] = setup;
  
  if ~exist('filename', 'var')
      %[filename, path] = uigetfile('*.dat', 'Select PD/pressure file')
      filename = 'acast1_jej.dat';
      path = 'C:\Users\evan\Documents\neurosim\ibs-tools\test-data\';
  end
  
  
  junk = load(strcat(path, filename));

  PD = -junk(:,3);
  PD = PD - min(PD);
  [maxPD, indx] = max(PD);
  y = PD(indx:end)/maxPD + EPS;
 
  % Time to PD peak
  t2pPD = (indx-ROR_START)/DT;
%   
  % Initial estimate for lambda
  t = 1:length(y);
  logstart = exp(polyfit(t',log(y),1));
  lam = logstart(2);
  % Refine fit
  Options = optimset('TolX', 1e-3, 'Display', 'off');
  lambda = fminsearch(@expfit, lam, Options, t', y);
  thalf = log(2)/lambda/DT;

  % exponential fit
  t = [0 t];
  fit = exp(-lambda*t)*maxPD;

  % Linear regression on initial slope
  tl = (ROR_START:indx)';
  pdl = PD(ROR_START:indx);
  X = [ones(size(tl)) tl];
  a = X\pdl;
  
  % Phase III duration
  PIIIdur = (length(PD)-ROR_START*2)/DT;
  
  % Time for PD to return with 10% of basline
  x = PD(indx:end);
  PDret = min(find(x<0.1*maxPD));
  if isempty(PDret)
    PDret = (length(PD)-indx)/DT;
  else
    PDret = PDret/DT;
  end
  
   
  % Pressure
  press = junk(:,2);
  press = press - min(press);
  press = press/max(press)*maxPD;
  meanP = mean(press(24:end-24));
  
  tm = 0:length(press)-1;
  line(press_axes, tm*5, press, 'Color', 'blue');

  % Linear fit
  pdl = a(1) + a(2)*tl;
  RoR = a(2)/5;


  % Plot results
  line(pd_axes, tm*5, PD, 'Color', 'blue');
  line(pd_axes, (t+indx)*5, fit, 'Color', 'red');
  line(pd_axes, tl*5, pdl, 'Color', 'red');
  
  % Display t_1/2
  xt = 100;
  yt = 1;
  text(pd_axes, xt, yt, ...
       sprintf('t_1_/_2 = %g', thalf), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display time to peak
  text(pd_axes, xt, yt+1*(maxPD/5), ...
       sprintf('time to PD peak = %g', t2pPD), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display rate of rise
  text(pd_axes, xt, yt+2*(maxPD/5), ...
       sprintf('Rate of rise = %g', RoR), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display max PD
  text(pd_axes, xt, yt+3*(maxPD/5), ...
       sprintf('Max PD = %g', maxPD), ...
       'fontsize',8,'FontWeight','bold');

  % Display mean pressure
  text(pd_axes, xt, yt+4*(maxPD/5), ...
       sprintf('mean pressure = %g', meanP), ...
       'fontsize',8,'FontWeight','bold');

 
end

function [pd_axes, press_axes] = setup
  UI_NAME = 'Pulse analyzer';
  try
      delete(findobj('Name', UI_NAME))
  end

  ud.mode = 'move';
  h1 = figure( ...
      'Units', 'characters', ...
      'Color',[0.8 1.0 0.8], ...
      'MenuBar', 'none', ...
      'Name', UI_NAME, ...
      'NumberTitle', 'off', ...
      'Position', [36.6667   15.5000  150   45], ...
      'UserData', ud, ...
      'toolbar',  'figure');

  set(h1, ...
      'WindowButtonMotionFcn', @mousemover, ...
      'CloseRequestFcn', {@expQuit, h1} ...
      );

  pd_axes = axes( ...
      'Position', [0.03 0.084016393442623 0.96 0.379098360655738], ...
      'XminorTick', 'on', ...
      'XminorGrid', 'on', ...
      'XGrid', 'on', ...
      'YGrid', 'on', ...
      'Tag', 'pd axes');

  press_axes = axes( ...
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
      'FontSize', 10, ...
      'Position', [0.32761087267525 0.47 0.353361945636624 0.03], ...
      'String', 'Potential Difference');
  
  uicontrol( ...
      'Style', 'text', ...
      'Parent', h1, ...
      'Units', 'normalized', ...
      'BackgroundColor', [0.8 1 0.8], ...
      'FontSize', 10, ...
      'FontWeight', 'bold', ...
      'Position', [0.313304721030043 0.92 0.386266094420601 0.04], ...
      'String', 'Pressure');
  
  uicontrol( ...
      'Style', 'push', ...
      'Parent', h1, ...
      'Units', 'normalized', ...
      'BackgroundColor', [0.996108949416342 0.644541084916457 0], ...
      'FontSize', 8, ...
      'Position', [0.9 0.003 0.09 0.04], ...
      'String', 'save to excel', ...
      'Enable', 'on', ...
      'callback', @save_to_excel);
 
end
  
function mousemover(o, crap)
end

function save_to_excel(o,crap)
end

function expQuit(o, crap, h)
closereq
end

function err = expfit(lambda, t, y)

  A = zeros(length(t),length(lambda));
  for j = 1:length(lambda)
    A(:,j) = exp(-lambda(j)*t);
  end
  c = A\y;
  z = A*c;
  err = norm(z-y);
end
