
% data vid 4 och 5 minuter


x=[ContrinphasePD(480+1200-240,:);ContrinphasePD(480+1200,:)];
for i=1:5
    y=[ContrinphasePD(i*2160+480+1200-240,:);ContrinphasePD(i*2160+480+1200,:)];
    x=[x;y];
end

%mean(x(1,:))