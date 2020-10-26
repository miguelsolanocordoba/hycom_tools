function hycom = read_hycoma_grd(runnum,blki,blkj)
%%READ_HYCOMA_GRD reads in HYCOM's grid (.BinF)
% HYCOM = READ_HYCOM reads HYCOM binaries (.BinF) and saves output 
% into a Matlab structure:
%
% hycom.lon   % longitude 
% hycom.lat   % latitude 
% hycom.h     % depth 
% 
% Created: October 25, 2020 by M. Solano 

% Format 
IEEE = 'ieee-be';
addpath /data/msolano/Matlab

dirin = '/data2/msolano/hycom/ATLc0.02/expt_04.3/'; % EXPT_19.0 

%% Experiment and tile number 
%runnum=043;  blki=1; blkj=1;
runnumstr = sprintf('%3.3d',runnum);
blkistr = sprintf('%2.2d',blki);
blkjstr = sprintf('%2.2d',blkj);

% Directories
fprintf('\nReading HYCOM files (read_hycoma)\n')
fprintf('Input directory: %s\n',dirin)
fprintf('iTile = %d\n',blki)
fprintf('jTile = %d\n',blkj)

% Grid file data
depfile = [dirin 'griddata/depth_' runnumstr '_blk_' ...
           blkistr '_' blkjstr '.BinF'];
lonfile = [dirin 'griddata/plon_' runnumstr '_blk_' ...
           blkistr '_' blkjstr '.BinF'];
latfile = [dirin 'griddata/plat_' runnumstr '_blk_' ...
           blkistr '_' blkjstr '.BinF'];

% Dimensions
nbf = 3;  % halo/padding 

nx=52; ny=38; nz=32; nt=624;
nxb=nx+nbf*2;
nyb=ny+nbf*2;

lenrec2 = nxb*nyb+2;
t = datenum(2019,1,1):datenum(0,0,0,1,0,0):datenum(2019,1,15);
nt = numel(t); 

% Load grid
fiddep = fopen(depfile,'r',IEEE);
fidlon = fopen(lonfile,'r',IEEE);
fidlat = fopen(latfile,'r',IEEE);

depdata = fread(fiddep,lenrec2,'single');
londata = fread(fidlon,lenrec2,'single');
latdata = fread(fidlat,lenrec2,'single');

lon = []; lat = []; depth = [];

depth(:,:) = permute(reshape(depdata(2:end-1),[nxb nyb]),[2 1]);
lon(:,:) = permute(reshape(londata(2:end-1),[nxb nyb]),[2 1]);
lat(:,:) = permute(reshape(latdata(2:end-1),[nxb nyb]),[2 1]);

fclose(fiddep);
fclose(fidlon);
fclose(fidlat);

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf]; 
a = [nbf+1:ny+nbf]; 

%% Save output to hycom (structure) 
hycom.lon = lon(a,b);         % longitude 
hycom.lat = lat(a,b);         % latitude 
hycom.h   = depth(a,b);       % depth 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
