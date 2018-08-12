lam = [1];,
trace = 0;,
tol = .01;,


t=ExpData(:,1);
y=ExpData(:,2);
% Sometimes the program does not correctly identify the percentile
% points
if length(t)>10
  logstart=exp(polyfit(t,log(max(.1,abs(y))),1));
  lam=logstart(2);
  
  % Old style function call - produces warning under matlab 6.5
  %lambda = fmins('expfit',lam,[trace tol])
  Options = optimset('TolX', tol, 'Display', 'off');
  lambda = fminsearch('expfit', lam, Options);
  
  A = zeros(length(t),length(lambda));
  for j = 1:length(lambda)
    A(:,j) = exp(-lambda(j)*t);
  end
  c = A\y;
  expcurve = A*c;
  H=line((tt+ii1+t)/60,yy2+expcurve);
  set(H, 'color', 'g', ...
	 'linewidth', 2, ...
	 'ButtonDownFcn', ...
	 ['trigg = ', num2str(trigg) '; typ=EXP; FixaPar']);

  Hexp(trigg)=H;
else
  lambda = NaN;
end
