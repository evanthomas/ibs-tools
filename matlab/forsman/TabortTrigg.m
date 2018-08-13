% TabortTrigg
  
  for l=1:length(UPS)
    set(UPS(l),  'visible','off')
    set(DOWNS(l),'visible','off')
  end
  delete(UPS);
  delete(DOWNS);
  
  clear UPS
  clear DOWNS
  %delete(UPS),
  ups=ups(find((1:length(ups))~=spec));
  downs=downs(find((1:length(downs))~=spec));
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
    set(UPS(l),'ButtonDownFcn',['spec = ', num2str(l) '; if mode==0, TabortTrigg, elseif mode==1, FlyttaTrigg, end']);
    set(DOWNS(l),'ButtonDownFcn',['spec = ', num2str(l) '; FlyttaSlut']);
  end
  
if previous_results
  delete(Ebuttons(5))
  ResF = [];
end
previous_results = 0;
