% FlyttaTrigg


  set(UPS(spec),'color','c')

  [x,y]=ginput(1);
  
  ups(spec)=round(x*60*4);
  ups=sort(ups);
  
  for l=1:length(UPS)
    set(UPS(l),'visible','off')
  end
  delete(UPS);
  clear UPS
  
  y=[0 100];
  x=[1 1]/(60*4);
  
  for l=1:length(ups),
    if ups(l)<=length(P)/2,
      subplot(4,1,1),
    else 
      subplot(4,1,3),
    end
    UPS(l)=line(ups(l)*x,y,'color','r');
    s = ['spec = ', num2str(l), ...
	 '; if mode==0, TabortTrigg, elseif mode==1, FlyttaTrigg, end'];
    set(UPS(l), 'ButtonDownFcn', s);
  end

  
if previous_results
  delete(Ebuttons(5))
  ResF = [];
end
previous_results = 0;
