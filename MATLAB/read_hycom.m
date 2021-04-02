function hycom = read_hycom(model,runnum,iblk,jblk)
%%READ_HYCOM reads in HYCOM's tiled output.
%  HYCOM = READ_HYCOM(MODEL,RUNNUM,BLKI,BLKJ) reads hycom tiled output
%  in binary format (.BinF). HYCOM simulations are stored according to:
% 
%  MODEL = Simulation case name. ('GLBc0.04' or 'ATLc0.02')
%  RUNNUM = Experiment number. (190 or 221)
%
%  The output is partitioned into 'tiles', indicated in the x-direction 
%  by BLKI and in the y-direction by BLKJ. The output is saved in Matlab 
%  structure format HYCOM. 
% 
%  hycom.time  % time (in datenum format)
%  hycom.lon   % longitude 
%  hycom.lat   % latitude 
%  hycom.h     % depth 
%  hycom.dz    % layer thickness
%  hycom.uiso  % baroclinic velocity (u) 
%  hycom.viso  % baroclinic velocity (v) 
%  hycom.rho   % density  
% 
% Created: July 10, 2020 by M. Solano 

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

switch runnum
  case 190
     t = datenum(2019,5,21):datenum(0,0,0,1,0,0):datenum(2019,6,19,23,0,0);
  case 221
     t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,30,23,0,0);
end

% Dimensions (continued)
nt = numel(t); 
nbf = 3;  % halo/padding 
nxb=nx+nbf*2;
nyb=ny+nbf*2;

lenrec2 = nxb*nyb+2; % record length

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

fprintf('\nReading HYCOM variables: lon, lat, depth, dz, uiso, viso, rho\n')
fprintf('Input directory: %s\n',dirin)
fprintf('iTile = %s\n',iblkstr)
fprintf('jTile = %s\n',jblkstr)

% Grid file data
depfile = [dirin 'griddata/depth_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
lonfile = [dirin 'griddata/plon_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
latfile = [dirin 'griddata/plat_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];

% Variables 
fname1 = [dirin 'u_iso/u_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
fname2 = [dirin 'v_iso/v_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
fname3 = [dirin 'thknss/thknss_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
fname4 = [dirin 'sig/sig_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
%fname5 = [dirin 'temp/T_' num2str(runnum) '_blk_' ...
%           jblkstr '_' iblkstr '.BinF'];
%fname6 = [dirin 'sal/S_' num2str(runnum) '_blk_' ...
%           jblkstr '_' iblkstr '.BinF'];

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

%% Read variables
% Open files
fid1 = fopen(fname1,'r',IEEE);
fid2 = fopen(fname2,'r',IEEE);
fid3 = fopen(fname3,'r',IEEE);
fid4 = fopen(fname4,'r',IEEE);
%fid5 = fopen(fname5,'r',IEEE);
%fid6 = fopen(fname6,'r',IEEE);

%% Read variables: u_iso, v_iso, sig, thknss
uiso = [];
viso = [];
thknss = [];
sig = [];
%sal = [];
%temp = [];

% extract layer thickness, u,v in space and time
fprintf('\nReading HYCOM output: \n') 
for i=1:nt
    
    fprintf('%d/%d\n',i,nt)	
    
    for k=1:nz
        alldata1 = fread(fid1,lenrec2,'single');
        alldata2 = fread(fid2,lenrec2,'single');
        alldata3 = fread(fid3,lenrec2,'single');
        alldata4 = fread(fid4,lenrec2,'single');
%        alldata5 = fread(fid5,lenrec2,'single');
%        alldata6 = fread(fid6,lenrec2,'single');
        
        uiso(:,:,k,i)   = permute(reshape(alldata1(2:end-1),[nxb nyb]),[2 1]);
        viso(:,:,k,i)   = permute(reshape(alldata2(2:end-1),[nxb nyb]),[2 1]);
        thknss(:,:,k,i) = permute(reshape(alldata3(2:end-1),[nxb nyb]),[2 1]);
        sig(:,:,k,i)    = permute(reshape(alldata4(2:end-1),[nxb nyb]),[2 1]);
%        temp(:,:,k,i)    = permute(reshape(alldata5(2:end-1),[nxb nyb]),[2 1]);
%        sal(:,:,k,i)    = permute(reshape(alldata6(2:end-1),[nxb nyb]),[2 1]);
        
    end
end
    
fprintf('\nDone reading variables!\n')
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
%fclose(fid5);
%fclose(fid6);

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf]; 
a = [nbf+1:ny+nbf]; 

%% Save output to hycom (structure) 
hycom.time = t;               % time (in datenum format)
hycom.lon = lon(a,b);         % longitude 
hycom.lat = lat(a,b);         % latitude 
hycom.h   = depth(a,b);       % depth 
hycom.dz  = thknss(a,b,:,:);  % layer thickness
hycom.uiso = uiso(a,b,:,:);   % baroclinic velocity (u) 
hycom.viso = viso(a,b,:,:);   % baroclinic velocity (v) 
hycom.rho = sig(a,b,:,:);     % density  
%hycom.salt = sal(a,b,:,:);    % salinity 
%hycom.temp = temp(a,b,:,:);   % temperature 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
