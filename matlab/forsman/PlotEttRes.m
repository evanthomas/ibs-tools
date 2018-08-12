tt=round(ups(trigg)/4);
tlow  = max(1, tt-300);
thigh = min(tt+700, length(PD1));
ty=PD1(tlow:thigh);
tP=P1(tlow:thigh);
slask=mean(tP(320:500));   %för att få lagom skala för
			   %tryckpuls (get suitable scale for
			   %the presssure pulse)
tP=.7*mean(ty(320:500))*tP/slask;


clf, H=plot((tlow:thigh)/60, tP,   (tlow:thigh)/60, ty); 
set(H(1),'color',8.5*[.1 .1 .1],'linestyle','-','linewidth',2)
set(H(2),'color','b','linestyle','-','marker','*')

set(gcf,'Name',['Puls ', num2str(trigg), ' - ', namn], ...
	'NumberTitle','off', ...
	'toolbar', 'figure')


axis([(tlow)/60 (thigh)/60 floor(min(ty)) 2+ceil(max(ty))]);

line([tt tt]/60, [floor(min(ty)) 2+ceil(max(ty))],'color','r')

% Editing mode for FixaPar
mode = MODE_DO;
typ=HIGH;
FixaPar

mode = MODE_DO;
typ=LOW;
FixaPar

mode = MODE_DO;
typ=PD_REBOUND;
FixaPar

mode = MODE_DO;
typ=SLOPE;
FixaPar

mode = MODE_DO;
typ=HIGH2;
FixaPar
