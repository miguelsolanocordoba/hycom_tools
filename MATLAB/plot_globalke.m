clear; close all; clc
%%PLOT_GLOBALKE plots global modal kinetic energy
%
% 
% Created: January 11, 2021 by M. Solano

%% Input
stp = 2;     % step size for plotting 
fbin = 2;    % frequency band (1=sub, 2=D1, 3=D2, 4=D2+)
fbinstr = num2str(fbin);

cmin = -3; cmax = 2;     % colorbar limits

switch fbin
  case 0
    fnm = 'D0'; 
  case 1
    fnm = 'D1'; 
  case 2
    fnm = 'D2'; 
  case 3
    fnm = 'D3'; 
  case 4
    fnm = 'D4'; 
  case 5
    fnm = 'D5'; 
  case 6
    fnm = 'D6'; 
  case 7
    fnm = 'HH'; 
  otherwise
    disp('fbin has a range [1-4], try again')
    return
end

% Global HYCOM (GLBc0.04) 
iblk = 60;          % number of tiles in x-dir
jblk = 35;          % number of tiles in y-dir
tilestr = num2str(iblk*jblk); % number of tiles in string

% Tile dimensions
nxt=150; nyt=200;   % Horizontal TILE dimensions (without padding) 
idm = iblk*nxt;     % number of grid cells in x-dir
jdm = jblk*nyt;     % number of grid cells in y-dir
nbf = 3;            % halo/padding
nxb = nxt + nbf*2;  % Horizontal TILE dimensions (with padding) 
nyb = nyt + nbf*2;
nmodes = 5;         % number of modes 

lenrec = nxb*nyb+2; % *Length of record* 

% Input directory 
dirin = '/data2/msolano/hycom/GLBc0.04/expt_19.0/ke_rotspec/';

%% Main loop: open and read spectral energy binaries (BinF) 
kecw = zeros(idm,jdm,nmodes); 
keccw = zeros(idm,jdm,nmodes); 

fprintf('\nReading files:')
count = 0; 
for i = 1:iblk
   for j = 1:jblk
      count = count + 1; 
      fprintf('\n%s/%s',num2str(count),tilestr);

      iblkstr = sprintf('%.2d',i);   
      jblkstr = sprintf('%.2d',j); 
     
      % open file
      fname1 = [dirin fnm '_cw_190_' jblkstr '_' iblkstr '.BinF'];
      fname2 = [dirin fnm '_ccw_190_' jblkstr '_' iblkstr '.BinF'];
      fid1 = fopen(fname1,'r','ieee-be'); 
      fid2 = fopen(fname2,'r','ieee-be'); 

      % pre-allocate for speed-up
      datam1 = zeros(lenrec,1);
      datam2 = zeros(lenrec,1);
      var1 = zeros(nxt,nyt); 
      var2 = zeros(nxt,nyt); 
      ke1 = zeros(nxt,nyt); 
      ke2 = zeros(nxt,nyt); 
       
      % loop over nmodes
      for k = 1:nmodes

         modstr = num2str(k); 
      
         % Read the entire record
         datam1 = fread(fid1,lenrec,'single');
         var1 = squeeze(reshape(datam1(2:end-1),[nxb,nyb]));
         ke1 = var1(nbf+1:end-nbf,nbf+1:end-nbf); 

         datam2 = fread(fid2,lenrec,'single');
         var2 = squeeze(reshape(datam2(2:end-1),[nxb,nyb]));
         ke2 = var2(nbf+1:end-nbf,nbf+1:end-nbf); 

	 kecw(1+nxt*(i-1):nxt*i,1+nyt*(j-1):nyt*j,k) = ke1;
	 keccw(1+nxt*(i-1):nxt*i,1+nyt*(j-1):nyt*j,k) = ke2;

      end 

      fclose(fid1);
      fclose(fid2);

   end
end

fprintf('\nReading done!\n')
fprintf('Closing files\n')

%% Pre-processing
% Load grid 
load global_coord.mat 
[~,ind] = min(abs(lon(100,:)-90)); 

% Plot sub-set only (see stp on Input) 
lon = lon(1:stp:end,1:stp:end); 
lat = lat(1:stp:end,1:stp:end); 
depth = depth(1:stp:end,1:stp:end); 
kecw = kecw(1:stp:end,1:stp:end,:); 
keccw = keccw(1:stp:end,1:stp:end,:); 

% masking
kecw(kecw==0)=nan;
keccw(keccw==0)=nan;
kecw = permute(kecw,[2 1 3]); 
keccw = permute(keccw,[2 1 3]); 

depth2 = zeros(size(kecw));
for i=1:nmodes
   depth2(:,:,i) = depth;
end
kecw(depth2<250)=nan;
keccw(depth2<250)=nan;

% Figure options
%contvec = [-4000 -1000 -50 -0.1];
contvec = [-2000 -0.1];  % isobaths for contours
stpb = 2;                % step size for boundaries
fntsz = 5;               % fontsize for gca
xmax = max(max(lon)); xmin = xmax-360; 
ymin = -60; ymax = 60; 

numph = 1;     numpv = 5;
hors  = 0.125; hore =0.15;
vers  = 0.125; vere =0.125;
Dsh   = 0.0; Dsv  =0.0;
pos = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);

%% Figure
fprintf('Plotting...\n')

%% 1) Decomposed by mode (CW+CCW)
numph = 1;     numpv = 3;
hors  = 0.125; hore =0.15;
vers  = 0.125; vere =0.08;
Dsh   = 0.0; Dsv  =0.0;
pos = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);

figure; colormap('jet');
pcolor(lon,lat,squeeze(log10(sum(kecw,3)+sum(keccw,3)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]);
title(cb,'log10 [W/m]')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); %text(85,50,'1','fontsize',8)

print('-dpng','-r600',[fnm '_KE_global.png'])

figure; colormap('jet');
subplot('position',pos(1,:)); % Mode 1
pcolor(lon,lat,log10(squeeze(kecw(:,:,1)+keccw(:,:,1)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]); 
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,30,'1','fontsize',8)

subplot('position',pos(2,:)); % Mode 2
pcolor(lon,lat,log10(squeeze(kecw(:,:,2)+keccw(:,:,2)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,30,'2','fontsize',8)

subplot('position',pos(3,:)); % Mode 3
pcolor(lon,lat,log10(squeeze(kecw(:,:,3)+keccw(:,:,3)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]); 
shading flat; hold on; cb=colorbar; caxis([cmin cmax]);
set(cb,'position',[0.81 0.31 0.02 0.4],'tickdir','out','ticks',[cmin:cmax])
title(cb,'log_{10}[m^2/s^2]')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90E','135E','180W','135W','90W','45W','0','45E'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,30,'3','fontsize',8)

print('-dpng','-r600',[ fnm '_KE_global_modes.png'])

return

%% 2) Decomposed by direction CW/CCW (sum of all modes)
%numph = 1;     numpv = 3;
%hors  = 0.125; hore =0.15;
%vers  = 0.125; vere =0.08;
%Dsh   = 0.0; Dsv  =0.0;
%pos = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);
%
%figure; colormap('jet');
%subplot('position',pos(1,:)); % CW+CCW
%pcolor(lon,lat,log10(squeeze(sum(kecw+keccw,3)))); hold on
%pbaspect([xmax-xmin ymax-ymin 1]);
%shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
%contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
%       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
%axis([xmin xmax ymin ymax]); yticks([-50:25:50])
%set(gca,'fontsize',fntsz,'tickdir','out'); text(82,50,'CW+CCW','fontsize',8)
%
%subplot('position',pos(2,:)); % CW
%pcolor(lon,lat,log10(squeeze(sum(kecw,3)))); hold on
%pbaspect([xmax-xmin ymax-ymin 1]);
%shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
%contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
%       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
%axis([xmin xmax ymin ymax]); yticks([-50:25:50])
%set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'CW','fontsize',8)
%
%subplot('position',pos(3,:)); % CCW
%pcolor(lon,lat,log10(squeeze(sum(keccw,3)))); hold on
%pbaspect([xmax-xmin ymax-ymin 1]);
%shading flat; hold on; cb=colorbar; caxis([cmin cmax]);
%set(cb,'position',[0.78 0.31 0.02 0.4],'tickdir','out','ticks',[cmin:cmax])
%title(cb,'log_{10}[m^2/s^2]')
%contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
%       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
%axis([xmin xmax ymin ymax]); yticks([-50:25:50])
%xticks([90:45:xmax]); xticklabels({'90E','135E','180W','135W','90W','45W','0','45E'})
%set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'CCW','fontsize',8)
%
%print('-dpng','-r600',[ fnm '_global_all.png'])

fprintf('Done!\n')

