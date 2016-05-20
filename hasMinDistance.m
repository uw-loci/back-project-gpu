function id=hasMinDistance(point, list, distance)
%returns 0 if pointid from the list has a distance to point < distance

id=0;
distance=distance^2;

for i=1:size(list, 2)
    if (list(1,i)-point(1))^2+(list(2,i)-point(2))^2<distance
        id=i;
    end
end