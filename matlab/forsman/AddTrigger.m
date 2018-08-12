
[x,y]=ginput(1);
for l=1:length(UPS)
   set(UPS(l),'visible','off')
   set(DOWNS(l),'visible','off')
end
delete(UPS);
delete(DOWNS);

clear UPS
clear DOWNS

ups(length(ups)+1)=round(x*60*4);
downs(length(ups))=ups(length(ups))+280*4;
ups=sort(ups);
downs=sort(downs);
y=[0 100];
x=[1 1]/(60*4);
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
