function wc = wcfromccv2(wpts, cpts, cc)
%find the three points in cam space closest to cc

%distlist(1,:)=1:size(wpts,2);
distlist=sqrt((cpts(1,:)-cc(1)).^2+(cpts(2,:)-cc(2)).^2);

[distlist, order]=sort(distlist);

dlist=zeros(2,size(wpts,2));

dlist(1,:)=order;
dlist(2,:)=distlist;

%distances=dlist(2, 1:3);
%points=dlist(1,1:3);

%  cc
%  cpts(:,points(1))
%  cpts(:,points(2))
%  cpts(:,points(3))
%  wpts(:,points(1))
%  wpts(:,points(2))
%  wpts(:,points(3))

%interpolate to get wc
%coord 1
pind=2;
while (abs(wpts(1,order(1))-wpts(1,order(pind))))<0.05
    pind=pind+1;
end

%pind
a=(cpts(1,order(1))-cpts(1,order(pind)))/(wpts(1,order(1))-wpts(1,order(pind)));
wc(1)=wpts(1,order(1))-(cpts(1,order(1))-cc(1))/a;

pind=2;

while (abs(wpts(2,order(1))-wpts(2,order(pind))))<0.05
    pind=pind+1;
end

%pind
a=(cpts(2,order(1))-cpts(2,order(pind)))/(wpts(2,order(1))-wpts(2,order(pind)));
wc(2)=wpts(2,order(1))-(cpts(2,order(1))-cc(2))/a;

% wpts(1,points(1))
% wpts(1,points(pind))
% cpts(1,points(1))
% cpts(1,points(pind))
% cc(1)
% (cpts(1,points(1))-cc(1))
% a
