
% data mellan 4 och 5 minuter

file=input('PD file?');
X =load(file);
x=[X(480+960:480+1200,1)];
for i=1:5
    y=[X(i*2160+480+960:i*2160+480+1200,i)];
    x=[x,y];
end



%mean(x(1,:))