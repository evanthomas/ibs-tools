function err = expfit(lambda)
%   FITFUN(lambda) returns the error between the data and the
%   values computed by the current function of lambda.
%   FITFUN assumes a function of the form
%
%     y =  c(1)*exp(-lambda(1)*t) + ... + c(n)*exp(-lambda(n)*t)
%
%   with n linear parameters and n nonlinear parameters.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.3 $  $Date: 1997/11/21 23:25:46 $

global ExpData

t=ExpData(:,1);
y=ExpData(:,2);
A = zeros(length(t),length(lambda));
for j = 1:length(lambda)
   A(:,j) = exp(-lambda(j)*t);
end
c = A\y;
z = A*c;
%set(Plothandle,'ydata',z)
%drawnow
err = norm(z-y);

