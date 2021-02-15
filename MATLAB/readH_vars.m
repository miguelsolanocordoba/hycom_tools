function aout = readH_vars(fid,idm,jdm,varname,LH);
% MCB, NRL, 2012/10/19
% readH_vars reads variables varname at a certain layer LH in hycom
% assumes the file is already opened with identifier FID
% the variables are listed in the .b file at line LS
% idm and jdm are lengths of x and y axis
%
% list of variables ====================================
%
% layered variables --------------------------------------    
% montg1
% srfhgt
% steric
% surflx
% salflx
% bl_dpth
% mix_dpth
% covice
% thkice
% temice
% u_btrop
% v_btrop
% u-vel
% v-vel
% thknss
% temp
% salin
% 
% grid --------------------------------------    
% plon
% plat
% qlon
% qlat
% ulon
% ulat
% vlon
% vlat
% pang
% pscx
% pscy
% qscx
% qscy
% uscx
% uscy
% vscx
% vscy
% cori
% pasp
%     
% topography --------------------------------------    
% depth
%
% sealevel --------------------------------------    
% ssh
%
% drag --------------------------------------    
% wdrg
% qdrg
%
% equilibrium tide ----------------------------------
% M2 S2 k1 o1 n2 p1 k2 q1

% % test
% varname = 'u-vel';
% LH = 2;
% 
% dirin = '/u/duck/scratch/richman/hycom_tide/';
% ftimeseries = 'Sept_2004_hourly';
% fnamev = [ftimeseries '/185_archv.2004_268_00.a'];
% fid=fopen([dirin fnamev],'r',ieee); 

% test test test test
%fid=3; idm=4500; jdm=3298; varname='plon'; LH=1;



%% fields --------------------------------------
if     strcmp(varname,'montg1');   LS=1;  k=1;
elseif strcmp(varname,'srfhgt');   LS=2;  k=1;
elseif strcmp(varname,'steric');   LS=3;  k=1;
elseif strcmp(varname,'surflx');   LS=4;  k=1;
elseif strcmp(varname,'salflx');   LS=5;  k=1;
elseif strcmp(varname,'bl_dpth');  LS=6;  k=1;
elseif strcmp(varname,'mix_dpth'); LS=7;  k=1;
elseif strcmp(varname,'covice');   LS=8;  k=1;
elseif strcmp(varname,'thkice');   LS=9;  k=1;
elseif strcmp(varname,'temice');   LS=10; k=1;
elseif strcmp(varname,'u_btrop');  LS=11; k=1;
elseif strcmp(varname,'v_btrop');  LS=12; k=1;
elseif strcmp(varname,'u-vel');    LS=13; k=LH;
elseif strcmp(varname,'v-vel');    LS=14; k=LH;
elseif strcmp(varname,'thknss');   LS=15; k=LH;
elseif strcmp(varname,'temp');     LS=16; k=LH;
elseif strcmp(varname,'salin');    LS=17; k=LH;

% tides 011.2013l 
% TIDEM2Re, TIDEM2Im, TIDES2Re, TIDES2Im, TIDEK1Re, TIDEK1Im
% TIDEO1Re, TIDEO1Im, TIDEN2Re, TIDEN2Im, TIDEM4Re, TIDEM4Im
elseif strcmp(varname,'TIDEM2Re'); LS=1;  k=1;
elseif strcmp(varname,'TIDEM2Im'); LS=2;  k=1;
elseif strcmp(varname,'TIDES2Re'); LS=3;  k=1;    
elseif strcmp(varname,'TIDES2Im'); LS=4;  k=1;        
elseif strcmp(varname,'TIDEK1Re'); LS=5;  k=1;    
elseif strcmp(varname,'TIDEK1Im'); LS=6;  k=1;        
elseif strcmp(varname,'TIDEO1Re'); LS=7;  k=1;    
elseif strcmp(varname,'TIDEO1Im'); LS=8;  k=1;        
elseif strcmp(varname,'TIDEN2Re'); LS=9;  k=1;    
elseif strcmp(varname,'TIDEN2Im'); LS=10; k=1;        
elseif strcmp(varname,'TIDEM4Re'); LS=11; k=1;    
elseif strcmp(varname,'TIDEM4Im'); LS=12; k=1;        
    
% grid --------------------------------------    
elseif strcmp(varname,'plon');   LS=1;   k=1;
elseif strcmp(varname,'plat');   LS=2;   k=1;
elseif strcmp(varname,'qlon');   LS=3;   k=1;
elseif strcmp(varname,'qlat');   LS=4;   k=1;
elseif strcmp(varname,'ulon');   LS=5;   k=1;
elseif strcmp(varname,'ulat');   LS=6;   k=1;
elseif strcmp(varname,'vlon');   LS=7;   k=1;
elseif strcmp(varname,'vlat');   LS=8;   k=1;
elseif strcmp(varname,'pang');   LS=9;   k=1;
elseif strcmp(varname,'pscx');   LS=10;  k=1;
elseif strcmp(varname,'pscy');   LS=11;  k=1;
elseif strcmp(varname,'qscx');   LS=12;  k=1;
elseif strcmp(varname,'qscy');   LS=13;  k=1;
elseif strcmp(varname,'uscx');   LS=14;  k=1;
elseif strcmp(varname,'uscy');   LS=15;  k=1;
elseif strcmp(varname,'vscx');   LS=16;  k=1;
elseif strcmp(varname,'vscy');   LS=17;  k=1;
elseif strcmp(varname,'cori');   LS=18;  k=1;
elseif strcmp(varname,'pasp');   LS=19;  k=1;
    
% topography --------------------------------------    
elseif strcmp(varname,'depth');   LS=1;  k=1;

% sealevel --------------------------------------    
elseif strcmp(varname,'ssh');   LS=1;  k=1;

% drag --------------------------------------    
elseif strcmp(varname,'wdrg');   LS=1;  k=1; % linear wave drag
elseif strcmp(varname,'qdrg');   LS=1;  k=1; % quadratic bottom drag    
    
% bodytide: M2 S2 k1 o1 n2 p1 k2 q1    
elseif strcmp(varname,'TIDEM2Re'); LS=1;  k=1;
elseif strcmp(varname,'TIDEM2Im'); LS=2;  k=1;
elseif strcmp(varname,'TIDES2Re'); LS=3;  k=1;
elseif strcmp(varname,'TIDES2Im'); LS=4;  k=1;
elseif strcmp(varname,'TIDEK1Re'); LS=5;  k=1;
elseif strcmp(varname,'TIDEK1Im'); LS=6;  k=1;
elseif strcmp(varname,'TIDEO1Re'); LS=7;  k=1;
elseif strcmp(varname,'TIDEO1Im'); LS=8;  k=1;
elseif strcmp(varname,'TIDEN2Re'); LS=9;  k=1;
elseif strcmp(varname,'TIDEN2Im'); LS=10;  k=1;
elseif strcmp(varname,'TIDEP1Re'); LS=11;  k=1;
elseif strcmp(varname,'TIDEP1Im'); LS=12;  k=1;
elseif strcmp(varname,'TIDEK2Re'); LS=13;  k=1;
elseif strcmp(varname,'TIDEK2Im'); LS=14;  k=1;
elseif strcmp(varname,'TIDEQ1Re'); LS=15;  k=1;
elseif strcmp(varname,'TIDEQ1Im'); LS=16;  k=1;
end 

%disp([num2str([LS k])])

%% extract data --------------------------------------
nt=idm*jdm;   
nta=floor((nt+4095)/4096)*4096;  %total length of line

% -- search from the beginning
% LS-1 is starting position
start_after_numbytes = ( (LS-1)+(k-1)*5 )*nta*4;
status = fseek(fid,start_after_numbytes,-1); % multiply times 4 bytes !
aout   = fread(fid,[1,    nt],'real*4');
aout   = reshape(aout,idm,jdm);

% --- so for now reshape works I need to transpose the matrix
aout=aout';
 
% --- remove land mask huge=2^100=1.2677*10^30 
aout(aout>1e30)=NaN; 

return
