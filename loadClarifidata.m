%procdir=Left_Hand_... folder
%name=Left_hand_2022.03.03...... file name
procdir=uigetdir();
folders=regexp(procdir,'\','split');
name=char(folders(6)); %Might need to change this number depending on how deep directory is.
r=1;
out.LEDwaves=[470,525,625,730,850,850];
out.ctitles={'Mel','HbT1','HbT2','A','StO2','HbT','HbO','HbR'};
out.factor=[1,1,1,1,1];
out.limits=[10,75,150,4,100,250,150,150];
out.name=name;
factor=[100,2300,2300,1,100];
pixelWidth=1024; % Allraw dim 3
pixelHeight=109; % Allraw dim 4
outname='DaemonTestPhantom';

% Assign data structures
exes=[ceil(pixelHeight./2)-5:ceil(pixelHeight./2)+5];
whys=[ceil(pixelWidth./2)-5:ceil(pixelWidth./2)+5];
RdMapsPlanar=zeros(pixelHeight,pixelWidth,5);
RdMapsSFD=zeros(pixelHeight,pixelWidth,1);
ChromMaps=zeros(pixelHeight,pixelWidth,5);
MaskMaps=zeros(pixelHeight,pixelWidth,1);
HeightMaps=zeros(pixelHeight,pixelWidth,1);

% Load data
disp(['Loading Clarifi data: ' name]);
filestring = [procdir '\' name];
disp(['Loading Planar Rd Rep ' int2str(r)])
filename_raw = [filestring '-' int2str(r-1) '_' outname '_planar_calibrated_rd'];
if exist(filename_raw) % Only keep going if data exists   
    
    % Load planar data
    [fid_raw, message] = fopen(filename_raw, 'rb');
    for i=1:5
        RdMapsPlanar(:,:,i) = fread(fid_raw,[pixelWidth pixelHeight],'float')';
    end
    fclose(fid_raw);
    
    % Load SFD data    
    disp(['Loading SFD Rd Rep ' int2str(r)])
    filename_raw = [filestring '-' int2str(r-1) '_' outname  '_calibrated_rd'];
    [fid_raw, message] = fopen(filename_raw, 'rb');
    for i=1:1
        RdMapsSFD(:,:,i) = fread(fid_raw,[pixelWidth pixelHeight],'float')';
    end
    fclose(fid_raw);
    
    % Load Chromophore data
    disp(['Loading Chromophore ' int2str(r)]);
    filename_raw = [filestring '-' int2str(r-1) '_' outname  '_chrom'];
    [fid_raw, message] = fopen(filename_raw, 'rb');
    for i=1:5
        ChromMaps(:,:,i) = fread(fid_raw,[pixelWidth pixelHeight],'float')';
    end
    fclose(fid_raw);
    
    % Load Mask data    
    disp(['Loading Mask ' int2str(r)]);
    filename_raw = [filestring '-' int2str(r-1) '_' outname  '_mask'];
    [fid_raw, message] = fopen(filename_raw, 'rb');
    MaskMaps(:,:) = fread(fid_raw,[pixelWidth pixelHeight],'float')';
    fclose(fid_raw);
    
    
    % Load Height data    
    disp(['Loading Height ' int2str(r)]);
    filename_raw = [filestring '-' int2str(r-1) '_' outname  '_height'];
    [fid_raw, message] = fopen(filename_raw, 'rb');
    for i=1:1
    HeightMaps(:,:,i) = fread(fid_raw,[pixelWidth pixelHeight],'float')';
    end
    fclose(fid_raw);
    
    % Load color data
    cname=[filestring '-' int2str(r-1) '_' outname  '_reg.bmp'];
    out.color=imread(cname);
        
    % Assign maps to "out" structure and    
    out.RdMaps(:,:,[1:5])=RdMapsPlanar(:,:,[1:5]);
    out.RdMaps(:,:,[6])=RdMapsSFD(:,:,1);
    
    for c=1:5
        out.ChromMaps(:,:,c)=out.factor(c).*ChromMaps(:,:,c);
    end
    out.ChromMaps(:,:,6)=(2/3).*out.ChromMaps(:,:,2)+(1/3).*out.ChromMaps(:,:,3);
    out.ChromMaps(:,:,7)=ChromMaps(:,:,5).*out.ChromMaps(:,:,6);
    out.ChromMaps(:,:,8)=(1-ChromMaps(:,:,5)).*out.ChromMaps(:,:,6);
    
    out.HeightMaps=HeightMaps;
    out.MaskMaps=MaskMaps;
    
    % Get ROI data
    out.Rd_roi=squeeze(mean(mean(out.RdMaps(exes,whys,:))));
    out.Chrom_roi=squeeze(mean(mean(out.ChromMaps(exes,whys,:))));
    out.Height_roi=squeeze(mean(mean(HeightMaps(exes,whys,:))));
    
    %f=plotSkinDataDefault(out)
else
    disp(['Filename ' name ' missing...moving on'])
end

