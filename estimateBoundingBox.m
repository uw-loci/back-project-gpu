function bb = estimateBoundingBox(ds)
    cpcom=sum(ds.cameraPos,1)./size(ds.cameraPos,1);
    lpcom=sum(ds.cameraPos,1)./size(ds.cameraPos,1);
    com=(cpcom+lpcom)./2;
    range=sqrt(sum((ds.laserOrigin-lpcom).^2));
    range=range+sqrt(sum((ds.cameraOrigin-cpcom).^2));
    range=ds.t0+ds.deltat*ds.t-range;
    bb.rightbottomback=[com(1)-range com(2)-range com(3)];
    bb.lefttopfront=[com(1)+range com(2)+range com(3)+range];
end
