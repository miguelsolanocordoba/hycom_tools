function hycom = read_hycomcoord(model,runnum,iblk,jblk)
%%READ_HYCOM reads in HYCOM's tiled output (lon,lat,depth only).
%  HYCOM = READ_HYCOMCOORD(MODEL,RUNNUM,BLKI,BLKJ) reads hycom tiled output
%  in binary format (.BinF). HYCOM simulations are stored according to:
% 
%  MODEL = Simulation case name. ('GLBc0.04' or 'ATLc0.02')
%  RUNNUM = Experiment number. (190 or 221)
%
%  The output is partitioned into 'tiles', indicated in the x-direction 
%  by BLKI and in the y-direction by BLKJ. The output is saved in Matlab 
%  structure format HYCOM. 
% 
%  hycom.lon   % longitude 
%  hycom.lat   % latitude 
%  hycom.h     % depth 
% 
% Example: 
%
% hycom = read_hycom('GLBc0.04',190,40,15)
% 
% Created: December 18, 2020 by M. Solano 

% Format 
IEEE = 'ieee-be';

% Dimensions (depends on model) 
switch model 
  case 'GLBc0.04'
     nx=150; ny=200; nz=41; nt=624;
  case 'ATLc0.02'
     nx=129; ny=194; nz=41; nt=730;
  otherwise 
    disp('Error: MODEL must be either GLBc0.04 or ATLc0.02')
end 

% Dimensions (continued)
nbf = 3;  % halo/padding 
nxb=nx+nbf*2;
nyb=ny+nbf*2;

lenrec2 = nxb*nyb+2;

%% Experiment and tile number 
% North Atlantic > runnum=221; jblk=27; iblk=45;
% South Pacific > runnum=190;  jblk=15; iblk=25;
% Amazon (1) > runnum=190;     jblk=19; iblk=40;
% Amazon (2) > runnum=190;     jblk=19; iblk=41;
% Amazon (3) > runnum=190;     jblk=18; iblk=40;
% Amazon (4) > runnum=190;     jblk=18; iblk=41;

% Directories
expt = num2str(runnum); 
dirin = ['/data2/msolano/hycom/' model '/expt_' expt(1:2) '.' expt(3) '/']; % 
runnumstr = num2str(runnum);
iblkstr = sprintf('%.2d',iblk); 
jblkstr = sprintf('%.2d',jblk); 

fprintf('\nReading HYCOM coordinates (lon,lat) and bathymetry (depth)\n')
fprintf('Input directory: %s/%s\n',dirin,'griddata')
fprintf('iTile = %s\n',iblkstr)
fprintf('jTile = %s\n',jblkstr)

% Grid file data
depfile = [dirin 'griddata/depth_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
lonfile = [dirin 'griddata/plon_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
latfile = [dirin 'griddata/plat_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];

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

fprintf('Done reading coordinates!\n')
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
