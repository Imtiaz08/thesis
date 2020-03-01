function [localrange] = Getlocalrange(car1,car2)
epoch1 = car1(:,1:2);
epoch2 = car2(:,1:2);
if length(epoch1)<=length(epoch2)
    epoch = epoch1;
else
    epoch = epoch2;
end
for t = 1:length(epoch)
    index1 = find(epoch1(:,2) == epoch(t,2));
    index2 = find(epoch2(:,2) == epoch(t,2));
    if isempty(index1)||isempty(index2)
        continue
    else
        localrange(t,1:2) = epoch(t,1:2);
        localrange(t,3) = norm(car1(index1,3:5)-car2(index2,3:5));
    end
end