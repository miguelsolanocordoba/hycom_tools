function hycom = read_hycomvar(model,runnum,iblk,jblk,var)
%%READ_HYCOMVAR reads in a variable from HYCOM's output (.BinF)
% HYCOM = READ_HYCOMVAR(MODEL,RUNNUM,BLKI,BLKJ,VAR) reads a single 
% variable VAR from HYCOM tile (BLKI,BLKJ) under input directory 
% /data2/msolano/hycom/MODEL/expt_RUNNUM/. 
%
% VAR can be:  uiso, viso, sig, temp, salt, ssh, steric or srfhgt. 
% The variable and time vector are loaded into a structure format:
%
% hycom.uiso  % baroclinic velocity (u) 
% hycom.viso  % baroclinic velocity (v) 
% hycom.sig   % density  
% hycom.temp  % temperature 
% hycom.salt  % salinity
% hycom.ubar  % eastward barotropic velocity
% hycom.vbar  % northward barotropic velocity
% hycom.ssh   % sea surface height 
% 
% hycom.time  % time (in datenum format)
%
% Created: October 3, 2020 by M. Solano 

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

%% Directories
expt = num2str(runnum);
dirin = ['/data2/msolano/hycom/' model '/expt_' expt(1:2) '.' expt(3) '/']; %
runnumstr = num2str(runnum);
iblkstr = sprintf('%.2d',iblk);
jblkstr = sprintf('%.2d',jblk);

fprintf('\nReading HYCOM variable: %s\n',var)
fprintf('Input directory: %s%s\n',dirin,var)
fprintf('iTile = %s\n',iblkstr)
fprintf('jTile = %s\n',jblkstr)

% Variables 
switch var
   case 'uiso'
fname = [dirin 'u_iso/u_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 3; 
   case 'viso'
fname = [dirin 'v_iso/v_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 3; 
   case 'sig'
fname = [dirin 'sig/sig_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 3; 
   case 'temp'
fname = [dirin 'temp/T_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 3; 
   case 'salt'
fname = [dirin 'sal/S_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 3; 
   case 'ssh'
fname = [dirin 'srfhgt/srfhgt_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 2; 
   case 'steric'
fname = [dirin 'steric/steric_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 2; 
   case 'srfhgt'
fname = [dirin 'srfhgt/srfhgt_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 2; 
   case 'ubar'
fname = [dirin 'ubaro/ubaro_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 2; 
   case 'vbar'
fname = [dirin 'vbaro/vbaro_' num2str(runnum) '_blk_' ...
           jblkstr '_' iblkstr '.BinF'];
      vardim = 2; 
   otherwise 
      disp('VAR must be: uiso, viso, ubar, vbar, sig, temp, salt, ssh, steric, srfhgt')
      return
end


%% Dimensions
nbf = 3;  % halo/padding 

nx=150; ny=200; nz=41; 
nxb=nx+nbf*2;
nyb=ny+nbf*2;

% time vector
lenrec2 = nxb*nyb+2;
switch runnum
  case 190
     t = datenum(2019,5,21):datenum(0,0,0,1,0,0):datenum(2019,6,19,23,0,0);
  case 221
     t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,30,23,0,0);
end
nt = numel(t); 

%% Read variables
% Open files
fid = fopen(fname,'r',IEEE);

%% Read selected variable
vars = [];

% extract layer thickness, u,v in space and time
if vardim==3
   for i=1:nt
       for k=1:nz
           alldata = fread(fid,lenrec2,'single');
           vars(:,:,k,i)   = permute(reshape(alldata(2:end-1),[nxb nyb]),[2 1]);
       end
   end
elseif vardim==2
   for i=1:nt
       alldata = fread(fid,lenrec2,'single');
       vars(:,:,i)   = permute(reshape(alldata(2:end-1),[nxb nyb]),[2 1]);
   end
end


fprintf('Done reading variable!\n')
fclose(fid);

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf]; 
a = [nbf+1:ny+nbf]; 

%% Save output to hycom (structure) 
hycom.time = t;               % time (in datenum format)
switch vardim
   case 2  % 2D variable
      hycom.vars = vars(a,b,:);    
   case 3  % 3D variable
      hycom.vars = vars(a,b,:,:);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
