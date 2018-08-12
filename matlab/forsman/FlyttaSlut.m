%TabortTrigg

set(DOWNS(spec),'color','c')

[x,y]=ginput(1);

downs(spec)=round(x*60*4);
downs=sort(downs);

for l=1:length(UPS)
   set(DOWNS(l),'visible','off')
end
delete(DOWNS);
clear DOWNS

y=[0 100];
x=[1 1]/(60*4);

for l=1:length(ups),
   if downs(l)<=length(P)/2,
      subplot(4,1,1),
   else 
      subplot(4,1,3),
   end
   DOWNS(l)=line(downs(l)*x,y,'color','y');
   set(DOWNS(l),'ButtonDownFcn',['spec = ', num2str(l) '; FlyttaSlut']);
end

if previous_results
  delete(Ebuttons(5))
  ResF = [];
end
previous_results = 0;
