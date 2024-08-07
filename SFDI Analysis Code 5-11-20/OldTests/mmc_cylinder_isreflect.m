clear;

load 'cylinder.mat'

[x,y,z] = size(volume);
x_mm = unitinmm * x;
y_mm = unitinmm * y;
z_mm = unitinmm * z;
srcdir = [0 0 -1];
srcdir = srcdir/norm(srcdir);
srcpos = [x_mm/2 y_mm/2 z_mm-4];
srcparam1 = [8.5 0 0];
srcdef=struct('srctype','disk',...
              'srcpos',srcpos,...
              'srcdir',srcdir,...
              'srcparam1',srcparam1);  
detsize = 12;
detpos = [x_mm/2-detsize/2 y_mm/2-detsize/2 4];
detdef =struct('srctype','planar',...
              'srcpos',detpos,...
              'srcdir',[1 0 0],...
              'srcparam1',[detsize 0 0],...
              'srcparam2',[0 detsize 0]);
     
%% create mesh
volume(volume == 0) = 1;
volume = volume - 1;
opt.distbound=1;    % set max distance that deviates from the level-set
opt.radbound=4;      % set surface triangle maximum size
opt.autoregion=0;     % don't save interior points
opt.A = diag([unitinmm,unitinmm,unitinmm]); % include voxel size in mm as scaling matrix
opt.B = zeros(3,1); % no translation
[node, elem, face] = v2m(volume, 1:max(volume(:)), opt, 1000, 'cgalmesh');
tooth = figure('name','cylinder');
plotmesh(node, elem);
colorbar;
node = node(:,1:3);
[node,elem]=mmcaddsrc(node,elem,srcdef);
tooth_src = figure;
plotmesh(node, elem);
colorbar;
[node,elem]=mmcadddet(node,elem,detdef);
tooth_det = figure;
plotmesh(node, elem);
colorbar;
%% create cfg for simulation  
cfg.nphoton = 1e8;
cfg.prop = [0 0 1 1;
           0.04 0.5 0.9 1.633;
           0.15 6.6 0.96 1.54]; 
cfg.maxdetphoton = cfg.nphoton;
cfg.node = node;
cfg.elem = elem;
cfg.srctype=srcdef.srctype;
cfg.srcpos=srcdef.srcpos;
cfg.srcdir=srcdef.srcdir;
cfg.srcparam1=srcdef.srcparam1;
cfg.detpos=[detdef.srcpos 1];
cfg.elemprop=cfg.elem(:,5);
cfg.elem=cfg.elem(:,1:4);
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.detpos = [detdef.srcpos 0];
cfg.detparam1 = [12 0 0 128];
cfg.detparam2 = [0 12 0 128];
cfg.issaveexit = 2; 
        
for isreflect = [0:3] 
    cfg.isreflect = isreflect;
    [fluence, detphoton, cfg] = mmclab(cfg);
    im = log(sum(detphoton.data,3)');
    tooth_figure = figure('name',strcat('MMC cylinder reflect:', int2str(cfg.isreflect)));
    imagesc(im);
    colorbar;
end