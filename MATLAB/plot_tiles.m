clear; close; clc;
%%PLOT_TILES plots HYCOM tile bounds
% 
%  MCB, USM, 2019-12-18

%% Options
IEEE = 'ieee-be'; % binary format (big-endian)
ieee='b';

% use tflg=5, for M2 tides only
% also change the scalefactor sclf
RUNNM = '221'; 
bandp = 'M2_JGAU/';      
vernum = '_v12'; 
tflg=5; 
sclf=0.3; % M2 tides; FINISHED

% Database
dirbase = ['/data/mbui/projects/HYCOM_chunk/GLBc0.04/expt_' RUNNM(1:2) '.' RUNNM(3) '/']   

% Dimensions
nx=9000; ny=7055;
lenrec = nx*ny+2;
idm=nx; jdm=ny; nz=41;

% Open files
dirin3 = [dirbase 'regional/'];

fnameg = 'regional.grid.a';
fidg=fopen([dirin3 fnameg],'r',ieee);

fnamed = 'regional.depth.a';
fidd=fopen([dirin3 fnamed],'r',ieee);

% Read coordinates 
plon  = readH_vars(fidg,idm,jdm,'plon',1);
plat  = readH_vars(fidg,idm,jdm,'plat',1);
depth = readH_vars(fidd,idm,jdm,'depth',1);

depth(isnan(depth)) = 0;

fclose(fidg);
fclose(fidd);

%% plot bins boxes tiles
stp = 4;
figA4L
contour(plon(1:stp:end,1:stp:end),plat(1:stp:end,1:stp:end),depth(1:stp:end,1:stp:end),[0.1 0.1],'k-')
holder
contour(plon(1:stp:end,1:stp:end),plat(1:stp:end,1:stp:end),depth(1:stp:end,1:stp:end),[2000 2000],'k-')

blki=150, blkj=200; %explodes on shepard
holder
ntx = floor(idm/blki);
nty = floor(jdm/blkj);

is = 1;
%for i=1:26
for i=1:ntx
    js = 1;   
    ib = [is is+blki-1 is+blki-1 is is];
    for j=1:nty
        jb = [js js js+blkj-1 js+blkj-1 js];
        for ii=1:4
            xx(ii) = plon(jb(ii),ib(ii));
            yy(ii) = plat(jb(ii),ib(ii));           
        end
        plot(xx,yy)
        text(mean(xx(1:4)),mean(yy(1:4)),[num2str(j) ',' num2str(i)]);
        %drawnow
        js = jb(3)+1;            
    end
    is = ib(2)+1;
end
