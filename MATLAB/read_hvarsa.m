function [varout,t] = read_hvarsa(varname,xtiles,ytiles)
%%READ_HVARSA reads HYCOM's vars for the ATLc0.02 simulations 
% VAROUT = READ_HVARSA(VARNAME,XTILES,YTILES) reads HYCOM variable 
% VARNAME on select tiles XTILES [xstart xend], YILES [ystart yend]
% from binaries (.BinF). 
%
% XTILES: range (1-52)
% YTILES: range (1-38)
%
% VARNAME:
% 'srfhgt'   % sea surface height 
% 
% Created: October 25, 2020 by M. Solano 

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
fprintf('\nReading HYCOM vars (Atlc_0.02; expt_04.3)\n')
fprintf('Input directory: %s\n',dirin)
fprintf('Tiles in x-direction = %d\n',blki)
fprintf('Tiles in y-direction = %d\n',blkj)

% Dimensions
t = datenum(2019,7,1):datenum(0,0,0,1,0,0):datenum(2019,7,31,23,0,0);
nt = numel(t);  % time (length)
nbf = 3;        % halo/padding 
nx=129; ny=194; % tile dimensions 
kdm = 32;       % vertical levels
nxb=nx+nbf*2;   % tile with padding
nyb=ny+nbf*2;   % tiel with padding

lenrec2 = nxb*nyb+2; % length of record

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf];
a = [nbf+1:ny+nbf];

% Initilize empty variables
ii = 0; 
count = 0; 
for i = blkis:blkie
    ii = ii + 1; 
    jj = 0;
    for j = blkjs:blkje
       jj = jj + 1;
       count = count + 1;

       fprintf('Reading (%d/%d)\n',count,numfiles);

       blkistr = sprintf('%2.2d',i);
       blkjstr = sprintf('%2.2d',j);
       
       % Grid file data
       varfile = [dirin varname '/' varname '_' runnumstr '_blk_' ...
                  blkjstr '_' blkistr '.BinF'];
       
       % Open files
       var1 = [];
       fidvar = fopen(varfile,'r',IEEE);
       for k = 1:nt
       
          % Read in grid
          vardata = fread(fidvar,lenrec2,'single');
          
          % Load and reshape
          var1(:,:,k) = permute(reshape(vardata(2:end-1),[nxb nyb]),[2 1]);

          % Discard padding
          var2(:,:,ii,jj,k) = var1(a,b,k); 

       end

       % Clear and close
       clear var1
       fclose(fidvar);

   end
end


% Check size and stitch tiles 
[ny,nx,~,~,~] = size(var2); 
varout = zeros(ny*blkj,nx*blki,k); 

% Stitch tiles
for i = 1:blki
   indsx = 1 + (i-1)*nx;
   index = indsx + nx - 1;
   for j =1:blkj
      indsy = 1 + (j-1)*ny;
      indey = indsy + ny - 1;

      for k =1:kdm
         varout(indsy:indey,indsx:index,k) = var2(:,:,i,j,k); 
      end

   end
end


%% Save output to hycom (structure) 
%hycom.lon = lon(a,b);         % longitude 
%hycom.lat = lat(a,b);         % latitude 
%hycom.h   = depth(a,b);       % depth 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
