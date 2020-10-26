function [lon,lat,depth] = read_hgrida(xtiles,ytiles)
%%READ_HGRIDA reads HYCOM's grid for the ATLc0.02 simulations 
% [LON,LAT,DEPTH] = READ_HGRIDA(VARNAME,XTILES,YTILES) reads HYCOM grid
% LON/LAT/DEPTH on select tiles XTILES [xstart xend], YILES [ystart yend]
% from binaries (.BinF).
%
% INPUT (for Atlc_0.02) 
% XTILES-> range (1-52)
% YTILES-> range (1-38)
%
% OUTPUT:
% 'lon'   % longitude
% 'lat'   % latitude
% 'depth' % depth (positive down)
% 
% Created: October 25, 2020 by M. Solano 

clc; close all; 

% Format 
IEEE = 'ieee-be';
addpath /data/msolano/Matlab

dirin = '/data2/msolano/hycom/ATLc0.02/expt_04.3/'; % EXPT_19.0 

%% Experiment and tile number 
% Expt and tile numbers
runnum=043; 
%blkis=1; blkie=52; blki = blkie-blkis+1; 
%blkjs=1; blkje=38; blkj = blkje-blkjs+1;
blkis = xtiles(1); blkie = xtiles(2); blki = blkie-blkis+1;
blkjs = ytiles(1); blkje = ytiles(2); blkj = blkje-blkjs+1;
numfiles = blki*blkj;

runnumstr = sprintf('%3.3d',runnum);

% Directories
fprintf('\nReading HYCOM grid (Atlc_0.02; expt_04.3\n')
fprintf('Input directory: %s\n',dirin)
fprintf('Tiles in x-direction = %d\n',blki)
fprintf('Tiles in y-direction = %d\n',blkj)

% Dimensions
nbf = 3;        % halo/padding 
nx=129; ny=194; % tile dimensions 
nxb=nx+nbf*2;   % tile with padding
nyb=ny+nbf*2;   % tiel with padding

lenrec2 = nxb*nyb+2; % length of record

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf];
a = [nbf+1:ny+nbf];

% Initilize empty variables
depth = [];
lon = [];
lat = [];

ii = 0; 
for i = blkis:blkie
    ii = ii + 1; 
    jj = 0;
    for j = blkjs:blkje
       jj = jj + 1;

       blkistr = sprintf('%2.2d',i);
       blkjstr = sprintf('%2.2d',j);
       
       % Grid file data
       depfile = [dirin 'griddata/depth_' runnumstr '_blk_' ...
                  blkjstr '_' blkistr '.BinF'];
       lonfile = [dirin 'griddata/plon_' runnumstr '_blk_' ...
                  blkjstr '_' blkistr '.BinF'];
       latfile = [dirin 'griddata/plat_' runnumstr '_blk_' ...
                  blkjstr '_' blkistr '.BinF'];
       
       % Open files
       fiddep = fopen(depfile,'r',IEEE);
       fidlon = fopen(lonfile,'r',IEEE);
       fidlat = fopen(latfile,'r',IEEE);
       
       % Read in grid
       depdata = fread(fiddep,lenrec2,'single');
       londata = fread(fidlon,lenrec2,'single');
       latdata = fread(fidlat,lenrec2,'single');
       
       % Load and reshape
       lon1 = []; lat1 = []; depth1 = [];
       depth1 = permute(reshape(depdata(2:end-1),[nxb nyb]),[2 1]);
       lon1 = permute(reshape(londata(2:end-1),[nxb nyb]),[2 1]);
       lat1 = permute(reshape(latdata(2:end-1),[nxb nyb]),[2 1]);

       % Discard padding
       depth2(:,:,ii,jj) = depth1(a,b); 
       lon2(:,:,ii,jj) = lon1(a,b); 
       lat2(:,:,ii,jj) = lat1(a,b); 

       % Clear and close
       clear lon1 lat1 depth1
       fclose(fiddep);
       fclose(fidlon);
       fclose(fidlat);

   end
end


% Check size and stitch tiles 
[ny,nx,~,~] = size(lon2); 
lon = zeros(ny*blkj,nx*blki); 
lat = zeros(ny*blkj,nx*blki); 
depth = zeros(ny*blkj,nx*blki); 

% Stitch tiles
for i = 1:blki
   indsx = 1 + (i-1)*nx;
   index = indsx + nx - 1;
   for j =1:blkj
      indsy = 1 + (j-1)*ny;
      indey = indsy + ny - 1;

      lon(indsy:indey,indsx:index) = lon2(:,:,i,j); 
      lat(indsy:indey,indsx:index) = lat2(:,:,i,j); 
      depth(indsy:indey,indsx:index) = depth2(:,:,i,j); 

   end
end

save('/data/msolano/atlantic_coord.mat','lon','lat','depth');


%% Save output to hycom (structure) 
%hycom.lon = lon(a,b);         % longitude 
%hycom.lat = lat(a,b);         % latitude 
%hycom.h   = depth(a,b);       % depth 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
