% FixaPar
tt=round(ups(trigg)/4);   %ii i sekunder
textyy = get(gca, 'Ylim');
txdlt = (textyy(2)-textyy(1))/15;
textyy = (textyy(2)+textyy(1))/2;

if (typ==HIGH),
  try
    H=HyyEtt(trigg);
    set(H,'visible','off')
    delete(H)
  
    Ht=HtmaxEtt(trigg);
    set(Ht,'visible','off')
    delete(Ht)
  end
  
  switch mode
   case MODE_UNDO
    yy=yy_old;
    ii=ii_old;
    yy_old=ResF(trigg,HIGH);
    ii_old=ResF(trigg,II_HIGH);
    ResF(trigg,HIGH)=yy;
    ResF(trigg,II_HIGH)=ii;

   case MODE_MOVE
    yy_old=ResF(trigg,HIGH);
    ii_old=ResF(trigg,II_HIGH);  
    [ii,yy] = ginput(1);
    ii=60*(ii-tt/60);
    ii=round(ii);
    ResF(trigg,HIGH)=yy;
    ResF(trigg,II_HIGH)=ii;
    
   case MODE_DEL
    yy_old=ResF(trigg,HIGH);
    ii_old=ResF(trigg,II_HIGH);
    ResF(trigg,HIGH)=NaN;
    ResF(trigg,II_HIGH)=NaN;
    
   case MODE_DO
    yy=ResF(trigg,HIGH);
    ii=ResF(trigg,II_HIGH);
    
   otherwise
    error('logic error in FixaPar')
  end
    

  if ~isnan(ResF(trigg,HIGH))
    H=line(tt/60+[ii-NR_LOW/2 ii+NR_LOW/2]/60,ResF(trigg,HIGH)*[1 1]);
    set(H,'color','k','linewidth',2);
    set(H,'ButtonDownFcn', ...
	['trigg = ', num2str(trigg) '; typ=HIGH; FixaPar']);
    HyyEtt(trigg)=H;
    
    try delete(HtmaxEtt(trigg)), end
    Ht=text((tt-200)/60,textyy+5*txdlt,...
	    ['max intrapulse PD = ',num2str(round(10*ResF(trigg,HIGH))/10),' mV'],...
	    'fontsize',8,'FontWeight','bold'); 
    HtmaxEtt(trigg)=Ht;
  end
end


if (typ==HIGH2),
  try
    H=Hyy2Ett(trigg);
    set(H,'visible','off')
    delete(H)
    
    Ht=Htmax2Ett(trigg);
    set(Ht,'visible','off')
    delete(Ht)
  end
  
  switch mode
   case MODE_UNDO
    yy=yy_old;
    ii=ii_old;
    yy_old=ResF(trigg,HIGH2);
    ii_old=ResF(trigg,II_HIGH2);
    ResF(trigg,HIGH2)=yy;
    ResF(trigg,II_HIGH2)=ii;    
    
   case MODE_MOVE    
    [ii,yy] = ginput(1);
    ii=60*(ii-tt/60);
    ii=round(ii);
    yy_old=ResF(trigg,HIGH2);
    ii_old=ResF(trigg,II_HIGH2);
    ResF(trigg,HIGH2)=yy;
    ResF(trigg,II_HIGH2)=ii;
    
   case MODE_DEL
    yy_old=ResF(trigg,HIGH2);
    ii_old=ResF(trigg,II_HIGH2);
    ResF(trigg,HIGH2)=NaN;
    ResF(trigg,II_HIGH2)=NaN;
    
   case MODE_DO
    yy=ResF(trigg,HIGH2);
    ii=ResF(trigg,II_HIGH2);
    
   otherwise
    error('logic error in FixaPar')
  end
  
  
  if ~isnan(ResF(trigg,HIGH2))
    H=line(tt/60+[ii-NR_LOW/2 ii+NR_LOW/2]/60,yy*[1 1]);
    set(H,'color','k','linewidth',2);
    set(H,'ButtonDownFcn',...
	  ['trigg = ', num2str(trigg) '; typ=HIGH2; FixaPar']);
    Hyy2Ett(trigg)=H;
    
    Ht=text((tt-200)/60,textyy-1*txdlt,...
	    ['PD at end = ',sprintf('%0.2f',ResF(trigg,HIGH2)),' mV'],...
	    'fontsize',8,'FontWeight','bold'); 
    Htmax2Ett(trigg)=Ht;
  end
  
end

if (typ==PD_REBOUND)
  try
    H=HpdrEtt(trigg);
    set(H,'visible','off')
    delete(H)
    
    H=HtpdrEtt(trigg);
    set(H,'visible','off')
    delete(H)
  end

  switch mode
   case MODE_DO
    yy=ResF(trigg,PD_REBOUND);
    ii=ResF(trigg,PD_REBOUND_X);
    
   case MODE_DEL
    yy_old=ResF(trigg,PD_REBOUND);
    ii_old=ResF(trigg,PD_REBOUND_X);
    ResF(trigg,PD_REBOUND)=NaN;
    ResF(trigg,PD_REBOUND_X)=NaN;
    
   case MODE_UNDO
    yy=yy_old;
    ii=ii_old;
    yy_old=ResF(trigg,PD_REBOUND);
    ii_old=ResF(trigg,PD_REBOUND_X);
    ResF(trigg,PD_REBOUND)=yy;
    ResF(trigg,PD_REBOUND_X)=ii;
    
   case MODE_MOVE
    yy_old=ResF(trigg,PD_REBOUND);
    ii_old=ResF(trigg,PD_REBOUND_X);
    [ii,yy] = ginput(1);
    ii=60*(ii-tt/60);
    ii=round(ii);
    ResF(trigg,PD_REBOUND)=yy;
    ResF(trigg,PD_REBOUND_X)=ii;
    
   otherwise
    error('Logic error in FixaPar')
  end
    
  if ~isnan(ResF(trigg, PD_REBOUND_X))
    xboobie = (tt + ResF(trigg, PD_REBOUND_X) + [-1 1]*NR_LOW/2)/60;
    yboobie = ResF(trigg,PD_REBOUND)*[1 1];
    H=line(xboobie, yboobie);
    set(H,'color','k','linewidth',2);
    set(H,'ButtonDownFcn',...
	  ['trigg=', num2str(trigg) '; typ=PD_REBOUND;FixaPar']);
    HpdrEtt(trigg)=H;
    
    str = sprintf('PD rebound = %0.2f mV',ResF(trigg,PD_REBOUND));
    HtpdrEtt(trigg)=text((tt-200)/60,textyy, str, ...	
			 'fontsize',8,'FontWeight','bold');
  end
end

if (typ==LOW),
  try
    H=HlowEtt(trigg);
    set(H,'visible','off')
    delete(H)
    
    H=HtlowEtt(trigg);
    set(H,'visible','off')
    delete(H)    
  end

  switch mode
   case MODE_DO
    low=ResF(trigg,LOW);
    ii=ResF(trigg,II_LOW);
    
   case MODE_DEL
    low_old=ResF(trigg,LOW);
    ii_old=ResF(trigg,II_LOW);
    ResF(trigg,LOW)=NaN;
    ResF(trigg,II_LOW)=NaN;
    
   case MODE_UNDO
    low=low_old;
    ii=ii_old;
    low_old=ResF(trigg,LOW);   
    ii_old=ResF(trigg,II_LOW);   
    ResF(trigg,LOW)=low;
    ResF(trigg,II_LOW)=ii;
    
   case MODE_MOVE
    [ii,low] = ginput(1);
    ii=60*(ii-tt/60);
    ii=round(ii);
    low_old=ResF(trigg,LOW);
    ii_old=ResF(trigg,II_LOW);
    ResF(trigg,LOW)=low;
    ResF(trigg,II_LOW)=ii;
    
   otherwise
    error('Logic error in FixaPar')
  end
    
  if ~isnan(ResF(trigg,LOW))
    H=line(tt/60+[ii-NR_LOW/2 ii+NR_LOW/2]/60,ResF(trigg,LOW)*[1 1]);
    set(H,'color','k','linewidth',2);
    set(H,'ButtonDownFcn',...
	  ['trigg = ', num2str(trigg) '; typ=LOW; FixaPar']);
    HlowEtt(trigg)=H;
    
    str = num2str(round(10*ResF(trigg,LOW))/10);
    Ht=text((tt-200)/60,textyy+1*txdlt, ['PD at start = ',str,' mV'],...
	    'fontsize',8,'FontWeight','bold');
    HtlowEtt(trigg)=Ht;
  end
end


if (typ==SLOPE),
  try
    H=HslopeEtt(trigg);
    set(H,'visible','off')
    set(H,'color','r')
    
    Ht=HtkEtt(trigg);
    set(Ht,'visible','off')
    delete(Ht)
  end

  switch mode
   case MODE_DO 
    yy=ResF(trigg,[RISE_Y_START RISE_Y_END]);
    ii=ResF(trigg,[RISE_X_START RISE_X_END]);
    
   case MODE_DEL
    yy_old=ResF(trigg,[RISE_Y_START RISE_Y_END]);
    ii_old=ResF(trigg,[RISE_X_START RISE_X_END]);
    ResF(trigg,[RISE_Y_START RISE_Y_END]) = [NaN NaN];
    ResF(trigg,[RISE_X_START RISE_X_END]) = [NaN NaN];
    delete(H)
    
   case MODE_UNDO
    yy_p = yy_old;
    ii_p = ii_old;
    yy_old=ResF(trigg,[RISE_Y_START RISE_Y_END]);
    ii_old=ResF(trigg,[RISE_X_START RISE_X_END]);
    ResF(trigg,[RISE_Y_START RISE_Y_END]) = yy_p;
    ResF(trigg,[RISE_X_START RISE_X_END]) = ii_p;
    
   case MODE_MOVE
    yy_old=ResF(trigg,[RISE_Y_START RISE_Y_END]);
    ii_old=ResF(trigg,[RISE_X_START RISE_X_END]);
    slask=get(gca,'CurrentPoint');
    x=slask(1,1);
    y=slask(1,2);
    if abs(y-ResF(trigg,RISE_Y_END))<abs(y-ResF(trigg,RISE_Y_START))
      [ii,yy] = ginput(1);
      ii=60*(ii-tt/60);
      
      ResF(trigg,RISE_Y_END)=yy;
      ResF(trigg,RISE_X_END)=ii;
    else   
      [ii,yy] = ginput(1);
      ii=60*(ii-tt/60);

      ResF(trigg,RISE_Y_START)=yy;
      ResF(trigg,RISE_X_START)=ii;
    end
   otherwise
    error('Logic error in FixaPar')
  end

  
  if ~isnan(ResF(trigg,RISE_X_START))
    ResF(trigg, SLOPE) = (ResF(trigg,RISE_Y_END)-ResF(trigg,RISE_Y_START)) ...
	/ (ResF(trigg,RISE_X_END)-ResF(trigg,RISE_X_START));
    slask=ResF(trigg, SLOPE);
    Ht=text((tt-200)/60,textyy+2*txdlt,['k = ',num2str(slask),' mV/min'], ...
	    'fontsize',8,'FontWeight','bold'); 
    HtkEtt(trigg)=Ht;
    
    
    % Regression line for PD rise
    xboobie = (tt+ResF(trigg,[RISE_X_START RISE_X_END]))/60;
    yboobie = ResF(trigg,[RISE_Y_START RISE_Y_END]);
    H=line(xboobie, yboobie,'color','g','linewidth',2);
    set(H,'ButtonDownFcn',...
	  ['trigg = ', num2str(trigg) '; typ=SLOPE; FixaPar']);
    set(H,'color','k','linewidth',2);
    HslopeEtt(trigg)=H;
  end
end

if typ==EXP 
  switch mode
   case MODE_MOVE
    return
    
   case MODE_DEL
    t_12_old = ResF(trigg, T_12);
    ResF(trigg, T_12) = NaN;
    
   case MODE_UNDO
    t_p = ResF(trigg, T_12);
    ResF(trigg, T_12) = t_12_old;
    t_12_old = t_p;
  end
end

if mode==MODE_DO
  text((tt-200)/60,textyy+4*txdlt,...
	    sprintf('mean press = %0.2f mmHg',ResF(trigg,PRESSURE)),...
	    'fontsize',8,'FontWeight','bold');
end

% Fix the exponential time course
try
  H=Hexp(trigg);
  set(H,'visible','off')
  delete(H)
  H = HtThalf(trigg);
  set(H,'visible','off')
  delete(H)
end
tt=round(ups(trigg)/4);
ii1=ResF(trigg,II_HIGH);
ii2=round(ResF(trigg,II_HIGH2));
yy2=ResF(trigg,HIGH2);
if ~isnan(ii1) & ~isnan(ii2) & ~isnan(ResF(trigg, T_12))
  yexp=PD1(tt+ii1:tt+ii2)-yy2;
  texp=(0:length(yexp)-1)';
  ExpData=[texp yexp];
  PlotExp
  ResF(trigg,T_12) = 0.693/lambda;
  if ~isnan(lambda)
    Ht=text((tt-200)/60,textyy+2.8*txdlt,...
	    ['t_1_/_2 = ',sprintf('%0.2f',.693/lambda),' s'],...
	    'fontsize',8,'FontWeight','bold'); 
    HtThalf(trigg)=Ht;
  end
end

% Put up the Undo button (or Ångra button as it appears to be
% Swedish)
if mode~=MODE_DO
  str = ['undo' num2str(trigg)];
  H = findobj('tag', str);
  if isempty(H)
    Butt = uicontrol('style','push','string','Ångra',...
		     'FontSize',8, 'FontWeight', 'normal', ...
		     'units','normal','pos',[bx+3*(bw+.01) .01 bw .05], ...
		     'call','mode=MODE_UNDO;FixaPar', ...
		     'Tag', str);
  end
end

% Reset the move for moving
mode = MODE_MOVE;

return  
