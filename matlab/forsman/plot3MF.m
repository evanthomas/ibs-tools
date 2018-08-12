% plot3MF
global ExpData

  % Indexes into results matix
  STARTS       =  1;
  ENDS         =  2;
  LOW          =  3;
  II_LOW       =  4;
  HIGH         =  5;
  HIGH2        =  6;
  II_HIGH      =  7;
  II_HIGH2     =  8;
  II_P10       =  9;
  II_P90       = 10;
  P10          = 11;
  P90          = 12;
  SLOPE        = 13;
  EXPAMP       = 14;
  T_12         = 15;
  PMEAN        = 16;
  AUC          = 17;
  RISE_X_START = 18;
  RISE_X_END   = 19;
  RISE_Y_START = 20;
  RISE_Y_END   = 21;
  PDMEAN       = 22;
  PD_REBOUND   = 23;
  PD_REBOUND_X = 24;
  PRESSURE     = 25;
  EXP          = 26;
  
  % Editing modes
  MODE_DO    = 1;
  MODE_UNDO  = 2;
  MODE_DEL   = 3;
  MODE_MOVE  = 4;
  MODE_BOGUS = 5;
  
  %[B,A] = butter(3,120/120);   %LP vid frekvensen 3/min 
  close all

  try
    [namn, fil] = filegui;
  catch
    return
  end
  
  data=load(fil);
  if min(size(data))==2,
    P = [0*(1:800)'; data(:,1); 0*(1:1200)'];
    PD = [0*(1:800)'; data(:,2); 0*(1:1200)'];
  else
    P = data(:,2);
    PD = data(:,3);
  end
  P = P(1:length(P)-rem(length(P),16));
  PD = PD(1:1*length(P));
  
  figure(1)
  set(gcf,'Name',['Trigger edit window - ' namn], ...
	'NumberTitle','off', ...
	'toolbar', 'figure')
  for i=1:2,
    %subplot(4,1,1); 
    %plot(a(:,[1])); 
    %axis([1,length(a), 0 150])
    %title('Duodenal distension, rat, 010703')
    %ylabel('Blood pressure, mmHg')
    tt=(1:length(P)/2)+((i-1)*length(P)/2);
    subplot(4,1,1+2*(i-1)); 
    plot(tt/240,P(tt)); 
    axis([min(tt/240) max(tt/240) -10 100])
    ylabel('Intestinal pressure, mmHg')
    if i==1,
      title(namn)
    end
    
    subplot(4,1,2+2*(i-1)); 
    plot(tt/240,PD(tt)); 
    axis([min(tt/240) max(tt/240) 0 12])
    ylabel('Intestinal PD, mV')
    if i==2,
      xlabel('Time, minutes')
    end
  end
  
  
  %Hitta triggers
  P=P';
  N=length(P);
  dP=[0 diff(P)];
  d2P=[0 diff(dP)];
  d2P2=[0 0 diff(dP(3:N)-dP(1:N-2)) 0];
  d2P3=[0 0 diff(dP(4:N)-dP(1:N-3)) 0 0];  
  d2P4=[0 0 0 diff(dP(5:N)-dP(1:N-4)) 0 0];
  dmax=(max(abs([d2P; d2P2; d2P3; d2P4])));
  
  
  bp=(1/240)*[ones(1,240) zeros(1,241) -ones(1,240)]; 
  ap=1;
  ups=0*P;
  ups(1:N-360)=filter(bp,ap,P(361:N));
  u=find((dmax>5) & (ups>2.5));
  starts=0*P;
  starts(u)=1+starts(u);
  u=1+find(diff(starts)>.5);
  
  starts=0*P; 
  starts(u)=1+starts(u);
  
  bp2=[1 -ones(1,240)];
  ap2=1;
  ups=filter(bp2,ap2,starts);
  ups=find(ups==1);
  ups=round(ups);   %Trigg i sample;

  % ups now contains candidate triggers based on a rate of change
  % criteria.  We also require pressure to be in a +/- 15% band of
  % one of the target pressures, 20s after the start of the pulse.
  % The target pressures are:
  Ptarget = [5 10 20 40 60 80];
  kindx = [];
  for i=1:length(ups)
    keep = 0;
    u = ups(i);
    pcand = P(u+4*80);
    for p=Ptarget
      if p*0.85<=pcand & pcand<=p*1.15
	keep = 1;
	break
      end
    end
    if ~keep, kindx(end+1)=i; end
  end
  ups(kindx) = [];      
      
  downs=ups+4*300;
  ResF = [];
  previous_results = 0;
  %Ta istället sparade triggers? (take instead saved triggers)
  D=dir(['trigg_' namn '.mat']);
  if length(D)>0,
    load(['trigg_' namn '.mat']);
    if ~isempty(ResF)
      previous_results = 1;
    end
  end
  
  starts=0*P;
  starts(ups)=1+starts(ups);
  ends=0*P;
  ends(downs)=1+ends(downs);
  
  bx=.1;
  bwL=.10;
  bw=.13;
  % Fixa triggers
  figure(1), hold on
  mode=1; %Flytta (move)
  CommandSpara=['save trigg_' namn '.mat ups downs'];
  
  %   for i=ups,
  % take away
  Ebuttons(1) = uicontrol('style','push', ...
			  'string','Ta Bort', ...
			  'FontSize',8,  ...
			  'FontWeight', 'normal', ...
			  'units','normal', ...
			  'pos',[bx .01 bwL .05], ...
			  'backgroundcolor', [255 177 113]/255, ...
			  'call','mode=0;uiresume');
  % move
  Ebuttons(1) = uicontrol('style','push', ...
			  'string','Flytta', ...
			  'FontSize',8, ...
			  'FontWeight', 'normal', ...
			  'units','normal', ...
			  'pos',[bx+bwL+.008 0.01 bwL .05], ...
			  'backgroundcolor', [255 177 113]/255, ...
			  'call','mode=1;uiresume');
  % add
  Ebuttons(2) = uicontrol('style','push', ...
			  'string',['Lägg till'], ...
			  'FontSize',8, ...
			  'FontWeight', 'normal', ...
			  'units','normal', ...
			  'pos',[bx+2*(bwL+.008) .01 bwL .05], ...
			  'backgroundcolor', [255 177 113]/255, ...
			  'call','AddTrigger;uiresume');
  
  % save starting points (red triggers)
  Ebuttons(3) = uicontrol('style','push', ...
			  'string','Spara starter', ...
			  'FontSize', 8, ...
			  'FontWeight', 'normal', ...
			  'units','normal', ...
			  'backgroundcolor', [231 229 178]/255, ...
			  'pos',[bx+3*(bw+.008) .01 bw .05], ...
			  'call','eval(CommandSpara);uiresume');
  
  Ebuttons(4) = uicontrol('style','push', ...
			  'string','Starter OK', ...
			  'FontSize',8, ...
			  'FontWeight', 'normal', ...
			  'units','normal', ...
			  'pos',[bx+4*(bw+.008) .01 bw .05], ...
			  'backgroundcolor', [231 229 178]/255, ...
			  'call','mode=2;uiresume');  %Gå vidare (happy)
  if previous_results
    Ebuttons(5) = uicontrol('Style', 'Checkbox', ...
			    'string', 'use saved results', ...
			    'FontSize',8, ...
			    'FontWeight', 'normal', ...
			    'units','normal', ...
			    'backgroundcolor', [231 229 178]/255, ...
			    'pos',[bx+5*(bw+.008) .01 bw+0.05 .05]);
  end
			    
  
  
  %if 0,
  %   if round(mod(j-.1,4))==1
  %      HES(ceil(j/4),1)=subplot('position',[0.01    0.57    0.48    0.37]);
  %   elseif round(mod(j-.1,4))==2
  %      HES(ceil(j/4),2)=subplot('position',[0.51    0.57    0.48    0.37]);
  %   elseif round(mod(j-.1,4))==3
  %      HES(ceil(j/4),3)=subplot('position',[0.01    0.13    0.48    0.37]);
  %   elseif round(mod(j-.1,4))==4
  %      HES(ceil(j/4),4)=subplot('position',[0.51    0.13    0.48    0.37]);
  %   end   
  %end   
  y=[-10 100];
  x=[1 1]/(60*4);
  UPS = [];
  DOWNS = [];
  for l=1:length(ups),
    if ups(l)<=length(P)/2,
      subplot(4,1,1),
    else 
      subplot(4,1,3),
    end
    UPS(l)=line(ups(l)*x,y,'color','r');
    if downs(l)<=length(P)/2,
      subplot(4,1,1),
    else 
      subplot(4,1,3),
    end
    DOWNS(l)=line(downs(l)*x,y,'color','y');
  end
  
  for l=1:length(ups);
    % remove/move trigger
    s = ['spec = ', num2str(l), ...
	 '; if mode==0, TabortTrigg, elseif mode==1, ', ...
	 'FlyttaTrigg, end;uiresume'];
    set(UPS(l), 'ButtonDownFcn', s);
    % move end
    s = ['spec = ', num2str(l), '; FlyttaSlut;uiresume'];
    set(DOWNS(l), 'ButtonDownFcn', s);
  end
  
  % Matlab bug?
  drawnow
  while (mode<2)
    uiwait
  end

  if previous_results
    previous_results = get(Ebuttons(5), 'Value');
  end
  
  close(1) % Close the trigger edit window
  
  ResF(:,STARTS)=ups';
  ResF(:,ENDS)=downs';
  

  % Resample to 1Hz
  % Note: type of filter used will be different depending on the
  % toolbox.
  try
    % First try signal processing toolbox
    PD1=decimate(PD(1:N-mod(N,4)),4);
    P1=decimate(P(1:N-mod(N,4)),4);
  catch
    % Try system identification toolbox
    PD1 = idresamp(PD(1:N-mod(N,4)), 4);
    P1 =  idresamp(P(1:N-mod(N,4))', 4)';
  end
  N1=length(PD1);

  
  ups=round(ups);
  downs=round(downs);
  %   starts=0*PD1;
  %   starts(ups)=1+starts(ups);
  
  figure(3), clf, plot(PD1(1:N1),'.');
  set(gcf,'Name',['Pulse selector - ' namn], ...
	  'NumberTitle','off', ...
	  'toolbar', 'figure')
    
  uicontrol('Style', 'Push','string','Excel..',...
	    'FontSize',8, 'FontWeight', 'normal', ...
	    'units','normal','pos',[bx+5*(bw+.01) .01 bw .05], ...
	    'call','excelgui');
  uicontrol('style','push', ...
	    'string','Spara', ...
	    'FontSize',8, ...
	    'FontWeight', 'normal', ...
	    'units','normal', ...
	    'pos',[bx+4*(bw+.01) .01 bw .05], ...
	    'call','SavePars');  %Gå vidare
  uicontrol('Style', 'Push','string','New file',...
	    'FontSize',8, 'FontWeight', 'normal', ...
	    'units','normal','pos',[0.05 .01 bw .05], ...
	    'call','plot3MF');
  movegui(gcf, 'south')
  
  
  NR_LOW=40;
  
  slask=axis;
  ymax=slask(4);
  
  trigg=0;
  for i=round(ups/4),
    trigg=trigg+1;
    %och trigglinjer
    PlotRes
    %Fixa lutningar och toppar
    figure(3+trigg);
    PlotEttRes

    %   for i=ups,
    %      Ebuttons(1) = uicontrol('style','push','string','Spara starter','FontSize',8, 'FontWeight', 'normal', ...
    %     'units','normal','pos',[bx+3*(bw+.01) .01 bw .05], ...
    %     'call','CommandSpara');
    
    uicontrol('style','push','string','Delete',...
	      'FontSize',8, 'FontWeight', 'normal', ...
	      'units','normal','pos',[bx+2*(bw+.01) .01 bw .05], ...
	      'call','mode=MODE_DEL;');    
    
    movegui(gcf, 'north')

  end
  
% Bring the selector window to the front
figure(3)

mode = MODE_MOVE;
