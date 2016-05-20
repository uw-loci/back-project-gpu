function calib=generateCameraCalibration(img)
img = imread(img);
figure; imshow(img);
img=img(:, end:-1:1, :);
%permute img here ...
r=double(squeeze(img(:,:,1)));
g=double(squeeze(img(:,:,2)));
b=double(squeeze(img(:,:,3)));

%img(1:35, 1:35, :)=0;

v=double(max(img,[],3)-min(img,[],3));

dr=zeros(size(r));
dg=zeros(size(g));
db=zeros(size(b));

for i=1:size(img, 1)
    for j=1:size(img,2)
        m=max(img(i,j,:));
        if(r(i,j)==m) 
            dr(i,j)=m;
        end
        if(g(i,j)==m) 
            dg(i,j)=m;
        end
        if(b(i,j)==m) 
            db(i,j)=m;
        end
    end
end

% dr=abs(r-g./2-b./2);
% dg=abs(g-r./2-b./2);
% db=abs(b-g./2-r./2);

s=squeeze(sum(img,3));

%r=r./g./b;
%g=g./r./b;
%b=b./g./r;
rp=dr.*v;
gp=dg.*v;
bp=db.*v;

rp=rp./max(rp(:));
gp=gp./max(gp(:));
bp=bp./max(bp(:));

rp(rp<0.3)=0;
gp(gp<0.3)=0;
bp(bp<0.3)=0;

v=v./max(v(:));

blackp=(max(s(:))-s);
blackp=blackp./max(blackp(:));
blackp(v>0.1)=0;
blackp(blackp<0.9)=0;
blackp(rp~=0)=0;
blackp(gp~=0)=0;
blackp(bp~=0)=0;

figure; imagesc(v);
figure; imagesc(s);
figure; imagesc(rp);
figure; imagesc(gp);
figure; imagesc(bp);
%figure; imagesc(blackp);
%return;

ir=0;
ig=0;
ib=0;
iblack=0;

ptlistr=zeros(2,1);
ptlistg=zeros(2,1);
ptlistb=zeros(2,1);
ptlistblack=zeros(2,1);

for windowx=1:20:(size(img,1)-40)
    for windowy=1:20:(size(img,2)-40)
        section=rp(windowx:(windowx+40), windowy:(windowy+40));
        if sum(section(:))~=0
            if sum(squeeze(section(1,:)))+sum(squeeze(section(end,:)))+sum(squeeze(section(:,1)))+sum(squeeze(section(:,end)))==0
                ir=ir+1;
                %ptlistr(:,ir)=extractPointfromRegion(section);%+[windowx windowy];
                [y,x]=extractPointfromRegion(section);
                ptlistr(1,ir)=x+windowx;
                ptlistr(2,ir)=y+windowy;
            end
        end
        
        section=gp(windowx:(windowx+40), windowy:(windowy+40));
        if sum(section(:))~=0
            if sum(squeeze(section(1,:)))+sum(squeeze(section(end,:)))+sum(squeeze(section(:,1)))+sum(squeeze(section(:,end)))==0
                ig=ig+1;
                %ptlistr(:,ir)=extractPointfromRegion(section);%+[windowx windowy];
                [y,x]=extractPointfromRegion(section);
                ptlistg(1,ig)=x+windowx;
                ptlistg(2,ig)=y+windowy;
            end
        end
        
        section=bp(windowx:(windowx+40), windowy:(windowy+40));
        if sum(section(:))~=0
            if sum(squeeze(section(1,:)))+sum(squeeze(section(end,:)))+sum(squeeze(section(:,1)))+sum(squeeze(section(:,end)))==0
                ib=ib+1;
                %ptlistr(:,ir)=extractPointfromRegion(section);%+[windowx windowy];
                [y,x]=extractPointfromRegion(section);
                ptlistb(1,ib)=x+windowx;
                ptlistb(2,ib)=y+windowy;
            end
        end
        
        section=blackp(windowx:(windowx+40), windowy:(windowy+40));
        if sum(section(:))~=0
            if sum(squeeze(section(1,:)))+sum(squeeze(section(end,:)))+sum(squeeze(section(:,1)))+sum(squeeze(section(:,end)))==0
                iblack=iblack+1;
                %ptlistr(:,ir)=extractPointfromRegion(section);%+[windowx windowy];
                [y,x]=extractPointfromRegion(section);
                ptlistblack(1,iblack)=x+windowx;
                ptlistblack(2,iblack)=y+windowy;
            end
        end
    end
end

figure; 
hold on;
plot(ptlistr(2,:), ptlistr(1,:), 'r.');
plot(ptlistg(2,:), ptlistg(1,:), 'g.');
plot(ptlistb(2,:), ptlistb(1,:), 'b.');
%plot(ptlistblack(2,:), ptlistblack(1,:), 'black.');
hold off;
%return
%in world coordinates all the coortinates of a point set have to be
%multiples of twenty (if one of the points of the set is chosen as the origin)

worldptr=double(ptlistr);
worldptg=double(ptlistg);
worldptb=double(ptlistb);

worldptr(1,:)=worldptr(1,:)./size(img,1).*1.5./0.8129;
worldptr(2,:)=worldptr(2,:)./size(img,2).*1.5./0.6082;
worldptg(1,:)=worldptg(1,:)./size(img,1).*1.5./0.8129;
worldptg(2,:)=worldptg(2,:)./size(img,2).*1.5./0.6082;
worldptb(1,:)=worldptb(1,:)./size(img,1).*1.5./0.8129;
worldptb(2,:)=worldptb(2,:)./size(img,2).*1.5./0.6082;

refpt=10;

worldptr(1,:)=worldptr(1,:)-worldptr(1,1);
worldptr(2,:)=worldptr(2,:)-worldptr(2,1);
worldptg(1,:)=worldptg(1,:)-worldptg(1,1);
worldptg(2,:)=worldptg(2,:)-worldptg(2,1);
worldptb(1,:)=worldptb(1,:)-worldptb(1,1);
worldptb(2,:)=worldptb(2,:)-worldptb(2,1);

roimin=[-0.1 -0.01];
roimax=[2 2];

cwptr=zeros(2, 1);
ccptr=zeros(2, 1);
j=0;

for i=1:size(worldptr,2)
    if (worldptr(1,i)>roimin(1)) && (worldptr(1,i)<roimax(1)) && (worldptr(2,i)>roimin(2)) && (worldptr(2,i)<roimax(2))
        if j==0 || hasMinDistance(worldptr(:,i), cwptr, 0.1)==0
            j=j+1;
            cwptr(1,j)=worldptr(1,i);
            cwptr(2,j)=worldptr(2,i);
            ccptr(1,j)=ptlistr(1,i);
            ccptr(2,j)=ptlistr(2,i);
        end
    end
end

roimin=[-0.1 -0.01];
roimax=[2 2];

cwptg=zeros(2, 1);
ccptg=zeros(2, 1);
j=0;

for i=1:size(worldptg,2)
    if (worldptg(1,i)>roimin(1)) && (worldptg(1,i)<roimax(1)) && (worldptg(2,i)>roimin(2)) && (worldptg(2,i)<roimax(2))
        if j==0 || hasMinDistance(worldptg(:,i), cwptg, 0.1)==0
            j=j+1;
            cwptg(1,j)=worldptg(1,i);
            cwptg(2,j)=worldptg(2,i);
            ccptg(1,j)=ptlistg(1,i);
            ccptg(2,j)=ptlistg(2,i);
        end
    end
end

roimin=[-0.1 -0.01];
roimax=[2 2];

cwptb=zeros(2, 1);
ccptb=zeros(2, 1);
j=0;

for i=1:size(worldptb,2)
    if (worldptb(1,i)>roimin(1)) && (worldptb(1,i)<roimax(1)) && (worldptb(2,i)>roimin(2)) && (worldptb(2,i)<roimax(2))
        if j==0 || hasMinDistance(worldptb(:,i), cwptb, 0.1)==0
            j=j+1;
            cwptb(1,j)=worldptb(1,i);
            cwptb(2,j)=worldptb(2,i);
            ccptb(1,j)=ptlistb(1,i);
            ccptb(2,j)=ptlistb(2,i);
        end
    end
end

% cwptr=cwptr([2 1], :);
% cwptg=cwptg([2 1], :);
% cwptb=cwptb([2 1], :);
% ccptr=ccptr([2 1], :);
% ccptg=ccptg([2 1], :);
% ccptb=ccptb([2 1], :);

figure; plot(cwptr(1,:), cwptr(2,:), 'r.');
figure; plot(cwptg(1,:), cwptg(2,:), 'g.');
figure; plot(cwptb(1,:), cwptb(2,:), 'b.');

offset=0.05;%this has to be tuned for each dataset ... it determines where the cutoff betweeen different pixel rows is made

cwptr=round((offset+cwptr)./0.2).*0.2-offset;
cwptg=round((offset+cwptg)./0.2).*0.2-offset;
cwptb=round((offset+cwptb)./0.2).*0.2-offset;

cwptr(1,:)=cwptr(1,:)+0.1;
cwptg(2,:)=cwptg(2,:)-0.1;
cwptb(2,:)=cwptb(2,:)+0.2;

cwpt=[cwptr cwptg cwptb];
ccpt=[ccptr ccptg ccptb];

% cwpt=-cwpt;
% ccpt=-ccpt;
% cwptr=-cwptr;
% ccptr=-ccptr;
% cwptg=-cwptg;
% ccptg=-ccptg;
% cwptb=-cwptb;
% ccptb=-ccptb;

testpoints=[
    100 50;
    545 336;
    224 99;
    400 400;
    ]';

wc=zeros(2,size(testpoints,2));
for i=1:size(testpoints,2)
    wc(:,i)=wcfromccv2(cwpt, ccpt, testpoints(:,i));
end

figure; 
hold on;
plot(cwptr(1,:), cwptr(2,:), 'r.');
plot(cwptg(1,:), cwptg(2,:), 'g.');
plot(cwptb(1,:), cwptb(2,:), 'b.');
plot(wc(1,:), wc(2,:), 'mo');
hold off;

figure;
hold on;
plot(ccptr(1,:), ccptr(2,:), 'r.');
plot(ccptg(1,:), ccptg(2,:), 'g.');
plot(ccptb(1,:), ccptb(2,:), 'b.');
plot(testpoints(1,:), testpoints(2,:), 'mo');
hold off;

calib.world=cwpt./2;
calib.cam=ccpt;
