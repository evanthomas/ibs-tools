%figure(3); clf;
%plot(PD1(1:N1),'.');
%hold on, plot(10*starts(1:N1),'r'), hold off

%trigg=0;
%for i=ups,
%   trigg=trigg+1;
%   PlotRes       
%end
save(['trigg_' namn '.mat'], 'ups', 'downs','ResF');
