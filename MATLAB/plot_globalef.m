clear; close all; clc
%%PLOT_GLOBALEF plots global modal energy flux
%
% 
% Created: January 11, 2021 by M. Solano

%% Input
stp = 2;     % step size for plotting 
fbin = 2;    % frequency band (1=sub, 2=D1, 3=D2, 4=D2+)
fbinstr = num2str(fbin);

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
dirin = '/data2/msolano/hycom/GLBc0.04/expt_19.0/energy_flux_test/';

%% Main loop: open and read spectral energy binaries (BinF) 
efx = zeros(idm,jdm,nmodes); 
efy = zeros(idm,jdm,nmodes); 

fprintf('\nReading files:')
count = 0; 
for i = 1:iblk
   for j = 1:jblk
      count = count + 1; 
      fprintf('\n%s/%s',num2str(count),tilestr);

      iblkstr = sprintf('%.2d',i);   
      jblkstr = sprintf('%.2d',j); 
     
      % open file
      fname1 = [dirin 'FEx_190_' fnm '_' jblkstr '_' iblkstr '.BinF'];
      fname2 = [dirin 'FEy_190_' fnm '_' jblkstr '_' iblkstr '.BinF'];
      fid1 = fopen(fname1,'r','ieee-be'); 
      fid2 = fopen(fname2,'r','ieee-be'); 

      % pre-allocate for speed-up
      datam1 = zeros(lenrec,1);
      datam2 = zeros(lenrec,1);
      var1 = zeros(nxt,nyt); 
      var2 = zeros(nxt,nyt); 
      efx1 = zeros(nxt,nyt); 
      efy1 = zeros(nxt,nyt); 
       
      % loop over nmodes
      for k = 1:nmodes

         modstr = num2str(k); 
      
         % Read the entire record
         datam1 = fread(fid1,lenrec,'single');
         var1 = squeeze(reshape(datam1(2:end-1),[nxb,nyb]));
         efx1 = var1(nbf+1:end-nbf,nbf+1:end-nbf); 

         datam2 = fread(fid2,lenrec,'single');
         var2 = squeeze(reshape(datam2(2:end-1),[nxb,nyb]));
         efy1 = var2(nbf+1:end-nbf,nbf+1:end-nbf); 

	 efx(1+nxt*(i-1):nxt*i,1+nyt*(j-1):nyt*j,k) = efx1;
	 efy(1+nxt*(i-1):nxt*i,1+nyt*(j-1):nyt*j,k) = efy1;

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
efx = efx(1:stp:end,1:stp:end,:); 
efy = efy(1:stp:end,1:stp:end,:); 

% masking
efx(efx==0)=nan;
efy(efy==0)=nan;
efx = permute(efx,[2 1 3]); 
efy = permute(efy,[2 1 3]); 

depth2 = zeros(size(efx));
for i=1:nmodes
   depth2(:,:,i) = depth;
end
efx(depth2<250)=nan;
efy(depth2<250)=nan;

% Figure options
%contvec = [-4000 -1000 -50 -0.1];
contvec = [-3000 -250 -0.1];  % isobaths for contours
stpb = 2;                % step size for boundaries
fntsz = 5;               % fontsize for gca
xmax = max(max(lon)); xmin = xmax-360; 
ymin = -60; ymax = 60; 

%% Figure                                                                   
fprintf('Plotting...\n')

%% EF magnitude
cmin = -1; cmax = 5;     % colorbar limits

numph = 1;     numpv = 3;
hors  = 0.125; hore =0.15;
vers  = 0.125; vere =0.08;
Dsh   = 0.0; Dsv  =0.0;
pos = subplot_hor_vertpos(numph,numpv,hors,hore,vers,vere,Dsh,Dsv);

figure; colormap('jet');
pcolor(lon,lat,squeeze(log10(sqrt(sum(efx,3).^2+sum(efy,3).^2)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); 
title(cb,'log10 [W/m]')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); %text(85,50,'1','fontsize',8)

print('-dpng','-r600',[fnm '_EF_global.png'])

figure; colormap('jet');
subplot('position',pos(1,:)); % Mode 1
pcolor(lon,lat,squeeze(log10(sqrt(efx(:,:,1).^2+efy(:,:,1).^2)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'1','fontsize',8)

subplot('position',pos(2,:)); % Mode 2
pcolor(lon,lat,squeeze(log10(sqrt(efx(:,:,2).^2+efy(:,:,2).^2)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'2','fontsize',8)

subplot('position',pos(3,:)); % Mode 3
pcolor(lon,lat,squeeze(log10(sqrt(efx(:,:,3).^2+efy(:,:,3).^2)))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); %colorbar('off')
set(cb,'position',[0.81 0.31 0.02 0.4],'tickdir','out','ticks',[cmin:cmax])
title(cb,'log10 [W/m]')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'3','fontsize',8)

print('-dpng','-r600',[fnm '_EF_global_modes.png'])

return

%% EFx
cmin = -5; cmax = 5;     % colorbar limits

figure; colormap('jet');
subplot('position',pos(1,:)); % Mode 1
pcolor(lon,lat,squeeze(efx(:,:,1))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'1','fontsize',8)

subplot('position',pos(2,:)); % Mode 2
pcolor(lon,lat,squeeze(efx(:,:,2))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'2','fontsize',8)

subplot('position',pos(3,:)); % Mode 3
pcolor(lon,lat,squeeze(efx(:,:,3))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); %colorbar('off')
set(cb,'position',[0.78 0.31 0.02 0.4],'tickdir','out','ticks',[cmin:cmax])
title(cb,'[W/m]')
subplot('position',pos(3,:)); % Mode 3
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'3','fontsize',8)

print('-dpng','-r300',[fnm '_EFx.png'])


%% EFy
figure; colormap('jet');
subplot('position',pos(1,:)); % Mode 1
pcolor(lon,lat,squeeze(efy(:,:,1))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'1','fontsize',8)

subplot('position',pos(2,:)); % Mode 2
pcolor(lon,lat,squeeze(efy(:,:,2))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); colorbar('off')
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'2','fontsize',8)

subplot('position',pos(3,:)); % Mode 3
pcolor(lon,lat,squeeze(efy(:,:,3))); hold on
pbaspect([xmax-xmin ymax-ymin 1]);
shading flat; cb=colorbar; caxis([cmin cmax]); %colorbar('off')
set(cb,'position',[0.78 0.31 0.02 0.4],'tickdir','out','ticks',[cmin:cmax])
title(cb,'[W/m]')
subplot('position',pos(3,:)); % Mode 3
contour(lon(1:stpb:end,1:stpb:end),lat(1:stpb:end,1:stpb:end),...
       -depth(1:stpb:end,1:stpb:end),[contvec],'k','LineWidth',0.25);
axis([xmin xmax ymin ymax]); yticks([-50:25:50])
xticks([90:45:xmax]); xticklabels({'90','135','-180','-135','-90','-45','0','45'})
set(gca,'fontsize',fntsz,'tickdir','out'); text(85,50,'3','fontsize',8)

print('-dpng','-r300',[fnm '_EFy.png'])

fprintf('Done!\n')

