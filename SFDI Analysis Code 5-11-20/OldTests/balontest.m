%t=d.z/v.z
%tx=v.x * d.z/v.z
%ty=v.y * d.z/v.z
%supply x and y coordinates for center of rectangular area bxc and byc
%supply width in bx and by
%bx0=bxc - bx/2;
%bx1=bxc + bx/2;
%by0=byc - by/2;
%by1=byc + by/2;
%if(tx>bx0) & (tx<bx1) & if(ty>by0) & if(ty<by1)

bz=600/scale;
bxc=xscaled/2;
byc=yscaled/2;
bx=20/scale;
by=20/scale;
bx0=bxc - bx/2;
bx1=bxc + bx/2;
by0=byc - by/2;
by1=byc + by/2;
fn = fieldnames(detphoton);
propf=find(strcmp('prop',fn));
if(~isempty(propf))
    fn(propf,:)=[];
end
unitf=find(strcmp('unitinmm',fn));
if(~isempty(unitf))
    fn(unitf,:)=[];
end
dataf=find(strcmp('data',fn));
if(~isempty(dataf))
    fn(dataf,:)=[];
end

%%
for i = 1:length(detphoton)
    walon(i).t=bz./detphoton(i).v(:,3);
    walon(i).tx=detphoton(i).p(:,1)+detphoton(i).v(:,1).*walon(i).t;
    walon(i).ty=detphoton(i).p(:,2)+detphoton(i).v(:,2).*walon(i).t;
    walon(i).b = (walon(i).tx>bx0) & (walon(i).tx<bx1) & (walon(i).ty>by0) & (walon(i).ty<by1);
    for k=1:numel(fn)
        detset(i).(fn{k})=detphoton(i).(fn{k})(walon(i).b,:);
    end
    detset(i).fp=[walon(i).tx(walon(i).b) walon(i).ty(walon(i).b) repmat(bz,sum(walon(i).b),1)];
    
end