% For Henrik: start with e..g. expo('AFALK_duo.dat')

function [thalf, meanP, t2pPD, maxPD, RoR, PDret, PIIIdur] = ...
    expo(filename, display)
% EXPO(FILENAME) fit and exponential to decimated PD data


  if ischar(filename)
    junk = load(filename);
    [p, n, e] = fileparts(filename);
  else
    junk = filename;
    n = 'raw input';
  end

  PD = -junk(:,3);
  PD = PD - min(PD);
  [maxPD, indx] = max(PD);
  y = PD(indx:end)/maxPD + 1e-3;

  % Time to PD peak
  t2pPD = (indx-24)*5;

  % Pressure
  press = junk(:,2);
  press = press - min(press);
  press = press/max(press)*maxPD;
  meanP = mean(press(24:end-24));
  
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
  fit = exp(-lambda*t)*maxPD;
  tm = 0:length(PD)-1;

  % Linear regression on initial slope
  tl = (24:indx)';
  pdl = PD(24:indx);
  X = [ones(size(tl)) tl];
  a = X\pdl;
  
  % Phase III duration
  PIIIdur = (length(PD)-48)*5;
  
  % Time for PD to return with 10% of basline
  x = PD(indx:end);
  PDret = min(find(x<0.1*maxPD));
  if isempty(PDret)
    PDret = (length(PD)-indx)*5;
  else
    PDret = PDret*5;
  end
  
  
  % Linear fit
  pdl = a(1) + a(2)*tl;
  RoR = a(2)/5;

  if nargin==2 & strcmp(upper(display), 'NODISPLAY')==1
    return
  end

  newplot

  % Plot results
  line(tm*5, PD, 'Color', 'blue');
  line(tm*5, press, 'Color', 8.5*[.1 .1 .1]);
  line((t+indx)*5, fit, 'Color', 'red');
  line(tl*5, pdl, 'Color', 'red');
  
  % Display t_1/2
  xt = 100;
  yt = 1;
  text(xt, yt, ...
       sprintf('t_1_/_2 = %g', thalf), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display time to peak
  text(xt, yt+1*(maxPD/5), ...
       sprintf('time to PD peak = %g', t2pPD), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display rate of rise
  text(xt, yt+2*(maxPD/5), ...
       sprintf('Rate of rise = %g', RoR), ...
       'fontsize',8,'FontWeight','bold');
  
  % Display max PD
  text(xt, yt+3*(maxPD/5), ...
       sprintf('Max PD = %g', maxPD), ...
       'fontsize',8,'FontWeight','bold');

  % Display mean pressure
  text(xt, yt+4*(maxPD/5), ...
       sprintf('mean pressure = %g', meanP), ...
       'fontsize',8,'FontWeight','bold');
  
  
%  legend('PD', 'pressure', 'exp fit')
  %title(sprintf('%s', texescape(n)))
  
  
 [meanP,maxPD,RoR,t2pPD,PIIIdur,thalf]'
  
  
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
  
  
% Changelog
% EAT: 2018-08-09
%      1. fileparts only return 3 values


