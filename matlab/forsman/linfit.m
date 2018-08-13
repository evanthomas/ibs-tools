

function [slope, intercept, r] = linfit(filename, display)
% EXPO(FILENAME) fit and exponential to decimated PD data

  if ischar(filename)
    junk = load(filename);
    [p, n, e, v] = fileparts(filename);
  else
    junk = filename;
    n = 'raw input';
  end
  PD = -junk(:,3);
  PD = PD - min(PD);
  
  % Pressure
  press = junk(:,2);
  
  % Linear regression
  X = [ones(size(press)) press];
  a = X\PD;
  y = a(1) + a(2)*press;

  slope = a(2);
  intercept = a(1);
  r = corrcoef(press, PD);
  r = r(1,2);

  if nargin==2 & strcmp(upper(display), 'NODISPLAY')==1
    return
  end
  
  % Plot
  plot(press, PD, '.', press, y, 'r')
  
  xlim = get(gca, 'XLim');
  xt = xlim(1);
  yt = (max(PD)+min(PD))/2;
  
  % Regression parameters
  text(xt, yt, ...
       sprintf('PD = %g*press + %g', a(2), a(1)), ...
       'fontsize',8,'FontWeight','bold');
  
  % Correlation coefficient
  text(xt, yt+3, ...
       sprintf('R = %g', r), ...
       'fontsize',8,'FontWeight','bold');
  
  
  xlabel('press'), ylabel('PD')
  %title(texescape(n));

  [a(2);r]