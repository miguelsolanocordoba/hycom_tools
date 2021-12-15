      PROGRAM PRES_CALC

! pres_calc_c004_v13.f, ! Maarten Buijsman, USM, 2020/01/08
!  fix fix dx dy error for conversion
!  nx*ny = 150*200 = 60*35 = 9000*7000
! pres_calc_c004_v12.f, ! Maarten Buijsman, USM, 2016/08/04
!  improved filtering with zero padding on right and 
!  constant density and layers are excluded from filtering
!  to speed up process
! pres_calc_c004_v11.f, ! Maarten Buijsman, USM, 2016/08/01
!  uses improved band-pass filter with NS=2 and zero padding!
! pres_calc_c004_v01.f, ! Maarten Buijsman, USM, 2016/07/22
!  density perturbation computed in each layer 
!  using drho/dz * eta instead of montgomery potential
!  alternative method: compute perturbation pressure -eta*drho_dz at 
!  density interfaces and then average to layer centers ....
!  !! MAKE SURE thin layers do not cause super big fluxes 
! montg_calc_c008_v13.f, ! Maarten Buijsman, USM, 2016/05/25
!  for 1/25 grid
!  include density time series
!  do we get fluxes on the shelf?
! montg_calc_c008_v11.f, ! Maarten Buijsman, NRL, 2015/04/23
!  do analysis for M2 and expt_05.1, 
! montg_calc_c008_v2.f, ! Maarten Buijsman, NRL, 2015/02/11
! v11. import: DO_FILT, maxobs, FC1, FC2, ftnm, tidcon; make sure numcon
! is adjusted in program
! v2. adjust for 41L

      IMPLICIT NONE
      
      integer idm,jdm,kdm,blki,blkj ! block sizes
      integer maxobs                ! number of observations
      integer mbdy                  ! size of halo
      real grav                     ! gravitational acceleration
      real RHO_0                    ! reference density
      real flag                     ! missing value flag
      real thktol                   ! tolerance for layer thickness. To preven blowup
      real thkdrg                   ! bottom boundary layer of linear drag
      real drgscl                   ! scale factor
      real cbmin,vonk,z0,cbp        ! for quadratic drag computation
      real thkbot                   ! thickness of bottom quadratic boundary layer (m)
      real onemm                    ! one mm 
      real FC1,FC2,dt               ! filter params  
      integer NS                    ! filter params
      integer TIVAR                 ! (=1) then time varying fields are stored
      integer RILE                  ! number of time steps to omit to reduce effect of ringing
      integer DO_FILT               ! switch for filtering (if 0 mean is still removed)
      integer CB_VAR                ! if 1, use spatially varying drag 
      integer jblks,jblke !1,17
      integer iblks,iblke !1,10
      character runid*3
     
      parameter(TIVAR=0)              ! (=1) then time varying fields are stored
      parameter(CB_VAR=1)              ! in case of 0, use constant drag  
      parameter(RILE=24)              ! exclude #time steps affected by ringing 
c      parameter(RILE=0)               ! exclude #time steps affected by ringing 
c      parameter(DO_FILT=1)             ! 1=filtering,0=no filtering,2=harmonic analysis  
                                       ! 3=do both 
c      parameter(maxobs=25)            ! length of time series (hours)
c      parameter(maxobs=169)           ! length of time series (hours)
c      parameter(maxobs=120)           ! length of time series (hours)                        
c      parameter(maxobs=744)           ! length of time series (hours)                        

      parameter(idm=9000,jdm=7000,kdm=41)
      parameter(blki=150,blkj=200)    ! GLBc0.04 => 60*35
      parameter(mbdy=3)               ! halo
c      parameter(blki=450,blkj=194)    ! GLBa0.08
c      parameter(blki=900,blkj=194)    ! GLBa0.08
c      parameter(jblks=1,jblke=17,iblks=1,iblke=5) !block idices
      
c      parameter(drgscl=0.000001)
      parameter(grav=9.806)
      parameter(thktol=0.0)       
      parameter(thkdrg=500.0)       
c      parameter(thkdrg=10000.0)       
      parameter(RHO_0=1034.0)
      parameter(onemm=0.001)     
      parameter(cbmin=0.0025,vonk=0.4,z0=0.01,thkbot=10.0) !cbmin=0.0025, vonk=0.4, z0=10mm (0.01) 
      parameter(flag=2.0**100)

! filter parameters
      parameter(NS=4)             !new; default number of sections (2*NS=order of filter)
      parameter(dt=1.0)           !Time step in hr
      REAL*8 A(NS),B(NS),C(NS),D(NS),E(NS)
     
      character ftnm*10 !DO NOT FORGET TO CHANGE THE FILT STRING


      character infil10*240,infil11*240,infillat*240
      character infil12*240,infil13*240,infil14*240
      character infil15*240,infil16*240,infil17*240,infil18*240
      character outfile1*60,outfile2*60,outfile3*60
      character outfile4*60,outfile5*60,outfile6*60,outfile7*60
      character outfile8*60,outfile9*60,outfile10*60,outfile11*60
      character outfile12*60,outfile13*60,outfile14*60
      character outfile15*60,outfile16*60,outfile17*60,outfile18*60
      character outfile19*60,outfile20*60,outfile21*60,outfile22*60
      character outfile23*60,outfile24*60,outfile25*60
      character outfile26*60,outfile27*60,outfile28*60,outfile29*60
      character outfile30*60,outfile31*60
      character*2 numstr(60)
      character*2 constr(8)

      integer i,j,k,l,iblk,jblk,joe,joe2,kbot
      real grav_r,Mbar,blk_mean_ub,blk_mean_vb,hatv,hatu,
     & grad_mprime_u,grad_mprime_v,depthu,depthv,
     & blk_mean_u,blk_mean_v,ptopl,pbotl,qdpu,dragu,dragv,pbop,
     & thkbop1,ubot,vbot,phi,plo,uu,vv,vmag,pbot 

c 3D variables
      real blk_sigavg,blk_thknss,blk_mean_t,layer_depth,init_depth,
     & eta_anom,M2,Mprime,blk_u,blk_v,blk_sigt,blk_mean_s,
     & init_dep_rho,drhodz,rhop1,rhop2,pres,pave,psum
c      dimension blk_sigavg(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_thknss(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_sigt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_mean_t(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm) !time-mean thickness
      dimension blk_mean_s(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm) !time-mean density
      dimension layer_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm+1) !kdm+1
      dimension init_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm+1)  !kdm+1 
      dimension init_dep_rho(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)  !kdm 
      dimension drhodz(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)        !kdm 
      dimension rhop1(kdm)        !kdm 
      dimension rhop2(kdm)        !kdm 
      dimension pres(kdm+1)       !kdm + 1
      dimension pave(kdm)        !kdm 
      dimension eta_anom(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
c      dimension M2(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
c      dimension Mprime(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_u(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_v(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)

c 3D variables (x,y,t)
      real, allocatable :: drag(:,:,:) ! 1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs

c 2D variables
      integer blk_coast
      dimension blk_coast(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      real blk_depth,blk_ubaro,blk_vbaro,blk_uscx,blk_vscy,
     & conv_u,conv_v,conv,u_at_p,v_at_p,duflux,dvflux,
     & dufluxa,dvfluxa,conva,dufluxs,dvfluxs,convs,
     & dufluxaver,dvfluxaver,convaver,dev,blk_drgten,blk_cb,
     & util5,util6,Dlin_ubc,Dlin_vbc,Dlin_ubt,Dlin_vbt,
     & Dlin_bc,Dlin_bt,Dlin_bta,Dlin_bca,thkbop,  
     & Dqdr_ubc,Dqdr_vbc,Dqdr_ubt,Dqdr_vbt,
     & Dqdr_bc,Dqdr_bt,Dqdr_bta,Dqdr_bca,
     & Dlin_uu,Dlin_vv,Dlin_uq,Dlin_vq,Dlin_ux,Dlin_vx,
     & Dqdr_uu,Dqdr_vv,Dqdr_uq,Dqdr_vq,Dqdr_ux,Dqdr_vx,
     & Dlin_bcq,Dlin_btq,Dlin_x,Dqdr_bcq,Dqdr_btq,Dqdr_x,
     & Dlin_bcqa,Dlin_btqa,Dlin_xa,Dqdr_bcqa,Dqdr_btqa,Dqdr_xa,
     & udhdx,vdhdy

      dimension blk_cb(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_drgten(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_ubaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_vbaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_uscx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_vscy(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
c      dimension conv_u(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
c      dimension conv_v(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension conv(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension duflux(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension dvflux(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension dufluxa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension dvfluxa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension conva(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension util5(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension util6(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_ubc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_vbc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_ubt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_vbt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      dimension Dlin_bt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_bc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_bta(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_bca(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_ubc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_vbc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_ubt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_vbt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bc(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bta(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bca(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension thkbop(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      dimension Dlin_uu(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_vv(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_uq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_vq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_ux(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_vx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_uu(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_vv(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_uq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_vq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_ux(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_vx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      dimension Dlin_bcq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_btq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_x(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bcq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_btq(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_x(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_bcqa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_btqa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dlin_xa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_bcqa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_btqa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension Dqdr_xa(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension udhdx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension vdhdy(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)


c--------filter parameters-------------------
c 3D
      real, allocatable :: blk_t_filt(:,:,:,:),blk_u_filt(:,:,:,:),
     & blk_v_filt(:,:,:,:),
     & blk_s_filt(:,:,:,:) ! s=sig,t=thickness 1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs

c 2D
      real, allocatable :: blk_ubarof(:,:,:),blk_vbarof(:,:,:) !1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60'/

! tides information -------------------------------------------------    
      integer tidcon1,numcon,numcon2,ncon,ipt
      integer tidcon
      parameter(ncon=8)
      CHARACTER*8   tidcons
      CHARACTER*2   TideMode(ncon)
      CHARACTER*24  Tides
      DATA TideMode/'M2','S2','K1','O1','N2','P1','K2','Q1'/
      LOGICAL tide_on(ncon)

c      integer,parameter :: tidcon=00010011 !Q1K2P1N2O1K1S2M2, 0=off,1=on
                                           !      N2    S2M2 selected
c      integer,parameter :: tidcon=00000001 !Q1K2P1N2O1K1S2M2, 0=off,1=on
                                           !              M2 selected
c      parameter(numcon=1) ! ADJUST!!!!!!!!!!!!!

! least square fit parameters-----------------------------
! maxobs = NDATA is number of data time steps
! MA is number of functions (2x number of constituents) 
      INTEGER MA !MA=numcon*2
c      PARAMETER(MA=numcon*2)

      REAL, allocatable :: Yh(:), Eh(:) !maxobs
      REAL, allocatable :: VECTh(:,:) !maxobs,MA
      REAL, allocatable :: Ah(:) !MA
      DOUBLE PRECISION, allocatable :: omsel(:) !numcon
      INTEGER, allocatable :: LISTA(:) !MA
c      REAL Uh(maxobs,MA),Vh(MA,MA)
c      REAL Wh(MA),Bh(maxobs)
c      REAL THRESH,CHISQh
c      PARAMETER(Eh=1.0)
c      Eh(:)=1.0
c      do i=1,MA
c       LISTA(i)=i
c      enddo

c      write(6,*) maxobs,MA,Eh

c ------------------------------------------------------------
! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid
      read(*,*) drgscl
!DO_FILT, maxobs, FC1, FC2, ftnm, tidcon; make sure numcon
! fix maxobs and numcon
      read(*,*) DO_FILT  !integer
      read(*,*) maxobs   !integer
      read(*,*) FC1      !real
      read(*,*) FC2      !real
      read(*,*) ftnm     !string
      read(*,*) tidcons  !char
  
      read( tidcons, '(i8.8)' ) tidcon
      write (6,*) tidcons,tidcon
      call flush(6)


      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,f10.7)') ' drgscl: ',drgscl
      call flush(6)

      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      write (6,'(a,2i4,2f5.1,i10)') ' DO_FILT maxobs FC1 FC2 tidcon: ',
     & DO_FILT, maxobs, FC1,FC2,tidcon
      call flush(6)

      write (6,'(a,a)') ' ftnm:',ftnm
      call flush(6)

      if (DO_FILT.eq.1) then
        write(6,'(a,a)') 'bandpass filtering applied: ',ftnm
      elseif (DO_FILT.eq.2) then
        write(6,'(a,a)') 'harmonic analysis applied: ',ftnm
      elseif (DO_FILT.eq.3) then
      write(6,'(a,a)') 'filtering and harmonic analysis applied: ',ftnm
      else 
        write(6,'(a,a)') 'NO filtering: ',ftnm
      endif  
      call flush(6)

      write(6,'(a,i2)') 'number of time steps excluded: ',RILE
      call flush(6)

c      WRITE(6,'(a,a)')'Tidal Modes included: ',trim(TIDES)
c      call flush(6)

c ------------------------------------------------------------
c invert band-pass hours
      FC1 = 1/FC1
      FC2 = 1/FC2

c ------------------------------------------------------------
c determine what tidal constituents are inlcuded
c tidcon 1 digit per constituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
      tidcon1 = tidcon
      do i=1,ncon
        tide_on(i) = mod(tidcon1,10) .eq. 1
c        write(6,*) i,tidcon1,tide_on(i),mod(tidcon1,10)
        tidcon1    = tidcon1/10      ! shift by one decimal digit
      enddo

      TIDES='                        '
      ipt=1
      numcon2=0
      do i=1,ncon
        if(tide_on(i))then
           TIDES(ipt:ipt+1)=TideMode(i)
c           write(6,*) ipt,TideMode(i)
           ipt=ipt+3
           numcon2=numcon2+1
        endif
      end do

      numcon = numcon2
      MA     = numcon*2

      WRITE(6,'(a,a)')'Tidal Modes included: ',trim(TIDES)
      call flush(6)

      WRITE(6,'(a,i2,i2)')'numcon MA is: ',numcon,MA
      call flush(6)


c -----------------------------------------------------------
c allocate variables
      ALLOCATE(drag(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs)) 
      ALLOCATE(blk_t_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_u_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_v_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_s_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs)) 
      ALLOCATE(blk_ubarof(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_vbarof(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs))
      ALLOCATE(Yh(maxobs),Eh(maxobs))
      ALLOCATE(VECTh(maxobs,MA)) 
      ALLOCATE(Ah(MA))
      ALLOCATE(omsel(numcon))
      ALLOCATE(LISTA(MA))

c some other freq stuff
      Eh(:)=1.0
      do i=1,MA
        LISTA(i)=i
      enddo
  
!! -----------------------------------------------------------
! prepare coeffiecent matrix a cos omt + b sin omt
! to be used in SVDFIT
!! -----------------------------------------------------------
      call AFUNC(VECTh,maxobs,MA,tide_on,dt,omsel)
c      write(6,*) tide_on,omsel 

c test --------------------------------------
c        do i=1,maxobs
c          Yh(i)=2.5*cos(omsel(1)*i*dt/24.0D0)
c     &         +1.5*cos(omsel(2)*i*dt/24.0D0-0.5d0) 
c     &         +0.5*cos(omsel(3)*i*dt/24.0D0-1.0d0) 
cc     &         +1.0*cos(omsel(4)*i*dt/24.0D0-2.0d0) 
cc     &         +0.15*cos(omsel(5)*i*dt/24.0D0-3.0d0) 
cc          write(6,*) Yh(i),Eh(i)
c        end do 
cc           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
c            CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
cc
c        do i=1,MA,2
c          write(6,*) sqrt(Ah(i)**2+Ah(i+1)**2),atan2(Ah(i+1),Ah(i))
c        enddo 
c        stop

!! loop over blocks ----------------------------------------------------
      do jblk=jblks,jblke
        do iblk=iblks,iblke  
!! loop over blocks ----------------------------------------------------
          write(6,*) 'blk_',jblk,iblk
          call flush(6)

          write(6,*) 'opening files for reading and writing ...'
          call flush(6)
          
!!--------------------------------------------------------------------------
!! OPEN 2D GRIDDATA FILES            
!!--------------------------------------------------------------------------

! open the time invariant data files and read data these files 
! are written as 3-dimesional arrays
! (lon,lat,depth) - (i,j,k)

! coast
         infil10='griddata/coast_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_coast

          close(10)

c          write(6,*) blk_coast(1:10,1:10)

! depth
         infil11='griddata/'
     &  //'depth_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'	  

          write(6,*) 'read',infil11	
          open(unit=11,file=infil11,status='old',
     &                           form='unformatted')

          read(11) blk_depth

          write(6,*) 'max depth is ',maxval(maxval(blk_depth(:,:),1))
c          write(6,*) 'depth is ',blk_depth(1,:)

          close(11)

! uscx
         infil10='griddata/uscx_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_uscx
c          write(6,*) blk_uscx(1:10,1)

          close(10)

! vscy
         infil10='griddata/vscy_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_vscy
c          write(6,*) blk_vscy(1,1:10)

          close(10)

!! linear wave drag
!         infil10='griddata/drgten_'//runid//
!     &'_blk_'
!     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!          write(6,*) 'read',infil10
!          open(unit=10,file=infil10,status='old',
!     &                           form='unformatted')
!
!          read(10) blk_drgten
!c          write(6,*) blk_vscy(1,1:10)
!
!          close(10)
!
!! quadratic bottom drag
!         infil10='griddata/cb_'//runid//
!     &'_blk_'
!     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!          write(6,*) 'read',infil10
!          open(unit=10,file=infil10,status='old',
!     &                           form='unformatted')
!
!          read(10) blk_cb
!c          write(6,*) blk_cb(1,1:10)
!
!          close(10)


!!--------------------------------------------------------------------------
!! OPEN 3D average files         
!!--------------------------------------------------------------------------

c         infil10='sigavg/sigavg_'//runid//
c     &'_blk_'
c     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
c
c          write(6,*) 'read',infil10
c          open(unit=10,file=infil10,status='old',
c     &                           form='unformatted')
c          do k=1,kdm
c            read(10) blk_sigavg(:,:,k)
c            write(6,*) k,blk_sigavg(100,100,k)
c          enddo   
c
c          close(10)
!
c          stop	  

!!--------------------------------------------------------------------------
!! OPEN 3D time files         
!!--------------------------------------------------------------------------
          infil12='sig/'
     &  //'sig_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=12,file=infil12,status='old',form='unformatted')

          infil13='thknss/'
     &  //'thknss_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=13,file=infil13,status='old',form='unformatted')

          infil14='u_iso/'
     &  //'u_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=14,file=infil14,status='old',form='unformatted')

          infil15='v_iso/'
     &  //'v_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=15,file=infil15,status='old',form='unformatted')


!!--------------------------------------------------------------------------
!! OPEN 2D time files         
!!--------------------------------------------------------------------------

        infil16='ubaro/'     
     &  //'ubaro_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=16,file=infil16,status='old',form='unformatted')

        infil17='vbaro/'     
     &  //'vbaro_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
     
        open(unit=17,file=infil17,status='old',form='unformatted')

!!--------------------------------------------------------------------------
!! OPEN 2D output files         
!!--------------------------------------------------------------------------

        if (TIVAR.EQ.1) then

!barotropic to baroclinic conversion rate
        outfile7='out/' 
     &  //'conv_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 

        open(unit=107,file=outfile7,status='replace',
     &     form='unformatted',access='sequential')

!depth-integrated x flux at p
        outfile10='out/' 
     &  //'duflux_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 

        open(unit=110,file=outfile10,status='replace',
     &      form='unformatted',access='sequential')

!depth-integrated y flux at p
        outfile11='out/' 
     &  //'dvflux_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'  

        open(unit=111,file=outfile11,status='replace',
     &      form='unformatted',access='sequential')

! lin dissipation as a function of time
!        outfile20='out/' 
!     &  //'Dlinbc_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=120,file=outfile20,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile21='out/' 
!     &  //'Dlinbt_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=121,file=outfile21,status='replace',
!     &     form='unformatted',access='sequential')

! quad dissipation as a function of time
!        outfile22='out/' 
!     &  //'Dqdrbc_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=122,file=outfile22,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile23='out/' 
!     &  //'Dqdrbt_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=123,file=outfile23,status='replace',
!     &     form='unformatted',access='sequential')

        end if
  
! always save time-mean fields------------------------------------------------
!barotropic to baroclinic conversion rate
        outfile12='out/' 
     &  //'conva_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 

        open(unit=112,file=outfile12,status='replace',
     &     form='unformatted',access='sequential')

!depth-integrated x flux at p
        outfile13='out/' 
     &  //'dufluxa_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 

        open(unit=113,file=outfile13,status='replace',
     &      form='unformatted',access='sequential')

!depth-integrated y flux at p
        outfile14='out/' 
     &  //'dvfluxa_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'  

        open(unit=114,file=outfile14,status='replace',
     &      form='unformatted',access='sequential')

!time-mean lin dissipation
!        outfile18='out/' 
!     &  //'Dlinbca_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=118,file=outfile18,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile19='out/' 
!     &  //'Dlinbta_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=119,file=outfile19,status='replace',
!     &     form='unformatted',access='sequential')

!time-mean qdr dissipation
!        outfile24='out/' 
!     &  //'Dqdrbca_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=124,file=outfile24,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile25='out/' 
!     &  //'Dqdrbta_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
!
!        open(unit=125,file=outfile25,status='replace',
!     &     form='unformatted',access='sequential')

!time-mean qdr dissipation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        outfile26='out/'
!     &  //'Dlinbcqa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=126,file=outfile26,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile27='out/'
!     &  //'Dlinbtqa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=127,file=outfile27,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile28='out/'
!     &  //'Dlinxa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=128,file=outfile28,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile29='out/'
!     &  //'Dqdrbcqa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=129,file=outfile29,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile30='out/'
!     &  //'Dqdrbtqa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=130,file=outfile30,status='replace',
!     &     form='unformatted',access='sequential')
!
!        outfile31='out/'
!     &  //'Dqdrxa_'//runid//
!     &'_'//trim(ftnm)//'_'
!     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=131,file=outfile31,status='replace',
!     &     form='unformatted',access='sequential')


!!--------------------------------------------------------------------------
!! Departure from resting depth zprime
!!--------------------------------------------------------------------------
!zprime start at k=2

! filter layer thickness, get mean thickness and M2 thickness
! then compute layer depth

        joe=1
	if (joe.EQ.1) then 
			  
	  blk_ubarof(:,:,:)=0.0
	  blk_vbarof(:,:,:)=0.0
	  
	  blk_s_filt(:,:,:,:)=0.0 
	  blk_t_filt(:,:,:,:)=0.0 
	  blk_u_filt(:,:,:,:)=0.0  
	  blk_v_filt(:,:,:,:)=0.0  

          do l=1,maxobs

            read(16) blk_ubaro
            read(17) blk_vbaro

            blk_ubarof(:,:,l)=blk_ubaro(:,:)
            blk_vbarof(:,:,l)=blk_vbaro(:,:)

 	    do k=1,kdm

	      read(12) blk_sigt(:,:,k)	     
	      read(13) blk_thknss(:,:,k)	     
              read(14) blk_u(:,:,k)
              read(15) blk_v(:,:,k)
 	   
              blk_s_filt(:,:,k,l)=blk_sigt(:,:,k)
              blk_t_filt(:,:,k,l)=blk_thknss(:,:,k)
              blk_u_filt(:,:,k,l)=blk_u(:,:,k)
              blk_v_filt(:,:,k,l)=blk_v(:,:,k)

c      write(6,*)l,k,blk_thknss(100,100,k),blk_t_filt(100,100,k,l)

c      write(6,*)'max=',k,maxval(blk_u(:,:,k)),maxval(blk_v(:,:,k)),
c     & maxval(blk_thknss(:,:,k))
	   
            enddo
          enddo !----first maxobs

	 close(12)
	 close(13)
	 close(14)
	 close(15)
	 close(16)
	 close(17)
  
!-----------set flag values to zero else filtering code doesn't work----	    
	 do l=1,maxobs
	   do j=1-mbdy,blkj+mbdy
             do i=1-mbdy,blki+mbdy  
	      
               if (l.eq.1.AND.blk_depth(i,j).EQ.flag) then
                 blk_depth(i,j)=0.0
!                 write(6,*) i,j
               endif
 
               if (blk_ubarof(i,j,l).eq.flag) blk_ubarof(i,j,l)=0.0
               if (blk_vbarof(i,j,l).eq.flag) blk_vbarof(i,j,l)=0.0

               do k=1,kdm
		   
	     if (blk_s_filt(i,j,k,l).eq.flag) blk_s_filt(i,j,k,l)=0.0
	     if (blk_t_filt(i,j,k,l).eq.flag) blk_t_filt(i,j,k,l)=0.0
             if (blk_u_filt(i,j,k,l).eq.flag) blk_u_filt(i,j,k,l)=0.0
             if (blk_v_filt(i,j,k,l).eq.flag) blk_v_filt(i,j,k,l)=0.0	

               enddo !k
	     enddo   !i
	   enddo     !j 
	            
	 end do!-----second maxobs! 	
	     	    	     	     		  
         write(6,*) 'finished reading time series'
         call flush(6)

       end if ! joe-----to run reading time series

!! =========================================================================
! compute Cd*|u|, where |u| is the abs. tot. velocity averaged over the bbl
! the u,v and thickness are unfiltered !!
!! =========================================================================

! I guess this adds some time ....
! per my silly convention layer_depth is negative .....
! blk_depth is positive
c       joe2=2 
       joe2=1 
       if (joe2.EQ.1) then 
       do l=1,maxobs

         layer_depth = 0.0 	   
         drag(:,:,l) = 0.0

         do j=1-mbdy,blkj+mbdy-1

           do i=1-mbdy,blki+mbdy-1

             if (blk_coast(i,j).eq.1) then

! compute layer depth
               layer_depth(i,j,kdm+1)=-blk_depth(i,j)  

               do k=kdm,2,-1   ! surface is zero

                 layer_depth(i,j,k) = layer_depth(i,j,k+1)
     &                             + blk_t_filt(i,j,k,l)
               enddo

! before filtering, within loop 
! quadratic bottom stress is applied over bbl or total depth is less
               thkbop1 = min(thkbot,blk_depth(i,j)) 
               ubot=0.0 
               vbot=0.0
               pbop = blk_depth(i,j)-thkbop1  !top of bot. b.b.l.
               phi  = max(-layer_depth(i,j,1),pbop)     !lower layer

               do k=1,kdm
		   
! -----------------------------------------------------------------------	    
! compute mean velocity in quadratic drag bbl
                 plo = phi  ! max(p(i,j,k),pbop)       !if k is above bbl, take pbop
                 phi = max(-layer_depth(i,j,k+1),pbop) !if k+1 is above bbl, take pbop

! mean velocities at p
! 0.5 is added below the loop when drag is computed
                 uu = blk_u_filt(i,j,k,l) + blk_u_filt(i+1,j,k,l) 
                 vv = blk_v_filt(i,j,k,l) + blk_v_filt(i,j+1,k,l) 

                 ubot = ubot + uu*(phi-plo) !transports
                 vbot = vbot + vv*(phi-plo)

!               if (thkbop1.LE.thkbot.AND.thkbop1.GE.0.0) then
!               if (j.EQ.24-3.AND.i.EQ.151-3) then
c               if (j.EQ.60-3.AND.i.EQ.117-3) then
c                   write(6,*) k,blk_t_filt(i,j,k,l),(phi-plo),ubot,vbot
c                   call flush(6)
c               endif 

               enddo !k
	      
! correct with total bbl thickness, or less if shallower
! ad barotropic velocity

               ubot=ubot/thkbop1
     &            + (blk_ubarof(i,j,l) + blk_ubarof(i+1,j,l)) 
               vbot=vbot/thkbop1
     &            + (blk_vbarof(i,j,l) + blk_vbarof(i,j+1,l)) 

! drag = Cb * |v| [m/s]
! include 1/thkbop for the fraction of layer calculation

               if (CB_VAR.eq.1) then   !spatially varying
c old                 cbp = max( cbmin,(vonk/log(0.5*blk_depth(i,j)/z0))**2 ) 
                 cbp = blk_cb(i,j)
               else
                 cbp = cbmin
               end if 
 
               vmag = 0.5*sqrt(ubot**2+vbot**2)   ! m/s
               drag(i,j,l)= cbp*vmag/thkbop1      ! m/s*1/m => 1/s

c               if (j.EQ.24-3.AND.i.EQ.151-3) then
c               if (j.EQ.60-3.AND.i.EQ.117-3) then
c      write(6,*)i,j,ubot,vbot,blk_depth(i,j),drag(i,j,l),vmag
c                 call flush(6)
c               endif 


             end if   !coast
    
           enddo !i

         enddo   !j

! test ===================================
!         i=151-3
!         j = 24-3
!         write(6,*) drag(i-1,j,l),drag(i,j,l),drag(i+1,j,l),
!     &     drag(i+2,j,l)
!         write(6,*) drag(i,j-1,l),drag(i,j,l),drag(i,j+1,l),
!     &     drag(i,j+2,l)
!       call flush(6)   
!       stop
! test ===================================

       end do!-----second maxobs! 
       endif !joe2

       write(6,*) 'finished computing unfiltered bottom drag'
       call flush(6)

!!--------------------------------------------------------------------------                        
!START FILTERING AT EACH GRID POINT
!!--------------------------------------------------------------------------	       	

       joe=1 !!------------To start filtering
       if (joe.EQ.1) then         
         write(6,*) 'Starting filtering ...'

!----------------------------------------------------------------
! find filter constants
! used in SUBROUTINE BPFILT(X,LX,NS,A,B,C,D,E,RILE)
         call BPDES(FC1,FC2,dt,NS,A,B,C,D,E)
c         write(6,*) A,B,C,D,E   
	      	
         blk_mean_t(:,:,:) = 0.0
         blk_mean_s(:,:,:) = 0.0

         do j=1-mbdy,blkj+mbdy     !j->latiude
          
c           write(6,*)'j = ',j,' of ',blkj 
        
	   do i=1-mbdy,blki+mbdy   !i->longitude

             if (blk_coast(i,j).eq.1) then

!----------------------------------------------------------------	  
         !ubar vbar filtering
         blk_mean_ub = sum(blk_ubarof(i,j,:))/  !mean of time series 
     &                 (max(1,size(blk_ubarof(i,j,:))))
         blk_ubarof(i,j,:)= blk_ubarof(i,j,:)-blk_mean_ub  !subtract mean       
      
         if (DO_FILT.eq.1) then

           CALL BPFILT(blk_ubarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)

         elseif (DO_FILT.eq.2) then
! harmonic analysis
! Yh is the data to be fitted, Eh is uncertainty = 1,
           Yh = blk_ubarof(i,j,:)

c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
            CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
c           do k=1,MA,2
c             write(6,*) sqrt(Ah(k)**2+Ah(k+1)**2),atan2(Ah(k+1),Ah(k))
c           enddo 

           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
c           do k=1,maxobs
c             write(6,*) k,Yh(k),blk_ubarof(i,j,k) 
c           enddo

           blk_ubarof(i,j,:) = Yh

         elseif (DO_FILT.eq.3) then
           CALL BPFILT(blk_ubarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)
           Yh = blk_ubarof(i,j,:) 
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_ubarof(i,j,:) = Yh
         endif


         blk_mean_vb = sum(blk_vbarof(i,j,:))/             !mean of time series
     &                 (max(1,size(blk_vbarof(i,j,:))))
         blk_vbarof(i,j,:)= blk_vbarof(i,j,:)-blk_mean_vb  !subtract mean       

         if (DO_FILT.eq.1) then

           CALL BPFILT(blk_vbarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)

         elseif (DO_FILT.eq.2) then
           Yh = blk_vbarof(i,j,:)

c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)

           blk_vbarof(i,j,:) = Yh
         elseif (DO_FILT.eq.3) then
           CALL BPFILT(blk_vbarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)
           Yh=blk_vbarof(i,j,:)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_vbarof(i,j,:) = Yh
         endif

               do k=1,kdm

!----------------------------------------------------------------        
! only filter when the layer is above the bottom


!----------------------------------------------------------------	  
         !thickness filtering   
 
         blk_mean_t(i,j,k) = sum(blk_t_filt(i,j,k,:))/    !mean of time series
     &                 (max(1,size(blk_t_filt(i,j,k,:))))     

         blk_mean_s(i,j,k) = sum(blk_s_filt(i,j,k,:))/ !mean of timeseries
     &                 (max(1,size(blk_s_filt(i,j,k,:))))

!----------------------------------------------------------------        
! only filter when the layer has a tickness
         if(blk_mean_t(i,j,k).gt.0.0) then 

! speed it up by omitting constant layer thicknesses
         if(abs(blk_mean_t(i,j,k)-blk_t_filt(i,j,k,1)).GT.0.001.OR.
     &      abs(blk_mean_t(i,j,k)-blk_t_filt(i,j,k,maxobs)).GT.0.001)
     &   then

         if (DO_FILT.eq.1) then
           blk_t_filt(i,j,k,:)= blk_t_filt(i,j,k,:)-blk_mean_t(i,j,k) !subtract mean 

           CALL BPFILT(blk_t_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

           blk_t_filt(i,j,k,:)=blk_t_filt(i,j,k,:)+blk_mean_t(i,j,k)

         elseif (DO_FILT.eq.2) then
           Yh = blk_t_filt(i,j,k,:)-blk_mean_t(i,j,k)!subtract mean 

c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)

           blk_t_filt(i,j,k,:) = Yh+blk_mean_t(i,j,k)
         elseif (DO_FILT.eq.3) then
           blk_t_filt(i,j,k,:)= blk_t_filt(i,j,k,:)-blk_mean_t(i,j,k)
           CALL BPFILT(blk_t_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

           Yh = blk_t_filt(i,j,k,:) 
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_t_filt(i,j,k,:) = Yh+blk_mean_t(i,j,k)
         endif
         endif !0.001

c         write(6,*)i,j,k,blk_mean_t(i,j,k)

!----------------------------------------------------------------        
         !density filtering   

c      write(6,*)i,j,k, blk_mean_s(i,j,k)

! speed it up by omitting layers of constant density
         if(abs(blk_mean_s(i,j,k)-blk_s_filt(i,j,k,1)).GT.0.001.OR.
     &      abs(blk_mean_s(i,j,k)-blk_s_filt(i,j,k,maxobs)).GT.0.001)
     &   then

         if (DO_FILT.eq.1) then
           blk_s_filt(i,j,k,:)= blk_s_filt(i,j,k,:)-blk_mean_s(i,j,k) !subtract mean 

           CALL BPFILT(blk_s_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

           blk_s_filt(i,j,k,:)=blk_s_filt(i,j,k,:)+blk_mean_s(i,j,k) 

         elseif (DO_FILT.eq.2) then
           Yh = blk_s_filt(i,j,k,:)-blk_mean_s(i,j,k)!subtract mean 

c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)

           blk_s_filt(i,j,k,:) = Yh+blk_mean_s(i,j,k)
         elseif (DO_FILT.eq.3) then
           blk_s_filt(i,j,k,:)= blk_s_filt(i,j,k,:)-blk_mean_s(i,j,k)
           CALL BPFILT(blk_s_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)
           Yh = blk_s_filt(i,j,k,:)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_s_filt(i,j,k,:) = Yh+blk_mean_s(i,j,k)
         endif
         endif !0.001

!----------------------------------------------------------------	  
         !baroclinic velocity filtering   

         blk_mean_u = sum(blk_u_filt(i,j,k,:))/              !mean of time series
     &                 (max(1,size(blk_u_filt(i,j,k,:))))     
         blk_u_filt(i,j,k,:)= blk_u_filt(i,j,k,:)-blk_mean_u  !subtract mean	

         if (DO_FILT.eq.1) then
           CALL BPFILT(blk_u_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

         elseif (DO_FILT.eq.2) then
           Yh = blk_u_filt(i,j,k,:)      !subtract mean 

c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)

           blk_u_filt(i,j,k,:) = Yh
         elseif (DO_FILT.eq.3) then
           CALL BPFILT(blk_u_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)
           Yh = blk_u_filt(i,j,k,:)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_u_filt(i,j,k,:) = Yh
         endif

         blk_mean_v = sum(blk_v_filt(i,j,k,:))/              !mean of time series
     &                 (max(1,size(blk_v_filt(i,j,k,:))))     
         blk_v_filt(i,j,k,:)= blk_v_filt(i,j,k,:)-blk_mean_v  !subtract mean	
 
         if (DO_FILT.eq.1) then
         
           CALL BPFILT(blk_v_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

         elseif (DO_FILT.eq.2) then
           Yh = blk_v_filt(i,j,k,:)      !subtract mean 
           
c           CALL SVDFIT(Yh,Eh,maxobs,Ah,MA,VECTh)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           
           blk_v_filt(i,j,k,:) = Yh
         elseif (DO_FILT.eq.3) then
           CALL BPFILT(blk_v_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)
           Yh = blk_v_filt(i,j,k,:)
           CALL LFIT(Yh,Eh,maxobs,Ah,MA,LISTA,MA,MA,VECTh)
           CALL TSERCON(Yh,maxobs,MA,Ah,dt,omsel)
           blk_v_filt(i,j,k,:) = Yh
         endif

         else
c           if(i.eq.150.and.j.eq.150) write(6,*)k,blk_mean_t(i,j,k)
           blk_s_filt(i,j,k,:)=0.0
           blk_u_filt(i,j,k,:)=0.0
           blk_v_filt(i,j,k,:)=0.0
           blk_t_filt(i,j,k,:)=0.0
         endif ! thickness criterium

               end do !k
  
             else !----coast
  
               blk_s_filt(i,j,:,:)=0.0
               blk_u_filt(i,j,:,:)=0.0
               blk_v_filt(i,j,:,:)=0.0
               blk_t_filt(i,j,:,:)=0.0
               blk_ubarof(i,j,:)=0.0
               blk_vbarof(i,j,:)=0.0  

             end if !coast	 	  
           end do ! i       
         end do  ! j
        
         write(6,*) 'Finished filtering...'  
         call flush(6)
      
       end if   !end joe-------------------------------
     

      joe=1         !used to start/stop running main program joe NOT 1 ->don't run
      if (joe.EQ.1) then	  	   
        write(6,*) 'Starting main program...'
        call flush(6)

!!--------------------------------------------------------------------------
!! MAIN LOOP FOR COMPUTING CONVERSION AND FLUXES      
!!--------------------------------------------------------------------------

! calculate the mean layer depths by summing from the bottom up
! init_depth should be positive

! use land switch here?
! what about layers with a small thickness?

! initialize variables
        dufluxa(:,:) = 0.0  !"a" is mean
        dvfluxa(:,:) = 0.0
        conva(:,:)   = 0.0
!        Dlin_bca(:,:) = 0.0
!        Dlin_bta(:,:) = 0.0
!        Dqdr_bca(:,:) = 0.0
!        Dqdr_bta(:,:) = 0.0
c cross-terms
!        Dlin_bcqa(:,:) = 0.0
!        Dlin_btqa(:,:) = 0.0        
!        Dlin_xa(:,:)   = 0.0
!        Dqdr_bcqa(:,:) = 0.0
!        Dqdr_btqa(:,:) = 0.0
!        Dqdr_xa(:,:)   = 0.0

        init_depth(:,:,:) = 0.0 
        init_dep_rho(:,:,:) = 0.0 
        drhodz(:,:,:) = 0.0
        
        do j=1-mbdy,blkj+mbdy

          do i=1-mbdy,blki+mbdy
     
            init_depth(i,j,kdm+1)=blk_depth(i,j)            

! -----------------------------------------------------------------------
! do as in HYCOM - momtum.f
! normalize drag with the minimum of blk_depth and thkdrg
! do not forget to include the scale factor !!
! blk_drgten in linear drag in [m/s]
! util5 = C/h = rh/h in [1/s]
!            util6(i,j) = max(0.1,thkdrg)  ![m] max of thkdrg or depth
!            util5(i,j) = blk_drgten(i,j)*drgscl/
!     &  min(util6(i,j),blk_depth(i,j))

! quadratic bottom stress is applied over bbl or total depth is less
            thkbop(i,j) = min(thkbot,blk_depth(i,j)) 

! init_depth should be positive
            do k=kdm,2,-1
              init_depth(i,j,k) = init_depth(i,j,k+1)
     &                                - blk_mean_t(i,j,k)

            enddo
c              write(6,*)i,j,1,blk_mean_t(i,j,1),init_depth(i,j,1)

! compute the mean (resting) depth of the layer centers
! init_dep_rho is also positive
! length is kdm
            do k=1,kdm
              init_dep_rho(i,j,k) = 0.5*
     &           (init_depth(i,j,k) + init_depth(i,j,k+1))
            enddo

! compute mean density gradient drhodz
! drhodz is only computed at interfaces, surface = 0
! length is kdm
! init_dep_rho is positive => k's are flipped!

            do k=1,kdm-1
            
! prevent inf values             
      if(init_dep_rho(i,j,k+1) - init_dep_rho(i,j,k)>0.0)then
              drhodz(i,j,k+1) = 
     &            (blk_mean_s(i,j,k)     - blk_mean_s(i,j,k+1))
     &           /(init_dep_rho(i,j,k+1) - init_dep_rho(i,j,k))
      endif

            enddo


          enddo
	    
        enddo

! ===================================================================
! loop over time
! ===================================================================
        do l=1,maxobs

! compute u*dhdx and v*dhdy
! old and bad: 0.5*(blk_uscx(i-1,j)  + blk_uscx(i,j))
! new:              blk_uscx(i,j)
          udhdx(:,:)=0.0
          vdhdy(:,:)=0.0

          do j=2-mbdy,blkj+mbdy

            do i=2-mbdy,blki+mbdy
              !flip because of minus sign
              udhdx(i,j) = (blk_depth(i-1,j) - blk_depth(i,j))
     &                   /  blk_uscx(i,j) !dx
     &                   *  blk_ubarof(i,j,l)          

              vdhdy(i,j) = (blk_depth(i,j-1) - blk_depth(i,j))
     &                   /  blk_vscy(i,j) !dy
     &                   *  blk_vbarof(i,j,l)

c         write(6,*)l,udhdx(1-mbdy:2,50),udhdx(blki-1:blki+mbdy,50)
c         write(6,*)l,vdhdy(1-mbdy:2,50),
c     &               vdhdy(blki-1:blki+mbdy,50)

            enddo
          
           enddo  

c        write(6,*)'udhdx',l,maxval(abs(udhdx),1)
c        write(6,*)'vdhdy',l,maxval(abs(vdhdy),1)


! read the time dependent data        
! initialize layer_depth at the bottom
          layer_depth(:,:,:) = 0.0 	   
          eta_anom(:,:,:) = 0.0
          duflux(:,:)=0.0
          dvflux(:,:)=0.0
          conv(:,:)=0.0

          do j=1-mbdy,blkj+mbdy

            do i=1-mbdy,blki+mbdy

              do k=1,kdm

! why do we do this mapping?  
! time-consuming? 
                blk_thknss(i,j,k)= blk_t_filt(i,j,k,l)
                blk_u(i,j,k)     = blk_u_filt(i,j,k,l)
                blk_v(i,j,k)     = blk_v_filt(i,j,k,l)
                blk_sigt(i,j,k)  = blk_s_filt(i,j,k,l)

              enddo

            enddo

          enddo  

          do j=1-mbdy,blkj+mbdy

            do i=1-mbdy,blki+mbdy

              if (blk_coast(i,j).EQ.1) then ! apply everywhere to speed up?

! determine position of layer interface
! they are negative!
! for each time step
  
             layer_depth(i,j,kdm+1)=-blk_depth(i,j)
 
              do k=kdm,2,-1   ! surface is zero

                layer_depth(i,j,k) = layer_depth(i,j,k+1)
     &                                + blk_thknss(i,j,k)

              enddo

!! ============================================================
!! calculate pressure fluxes and conversion
!! ============================================================

! eta_anom is positive upward
! eta_anom(1) is surface and is zero
! eat_anom(kdm) is last value ABOVE bottom
! layer_depth is negative, init_depth is positive
! rhop1 and rhop2 have length kdm

             rhop1(:)=0.0
             rhop2(:)=0.0
             kbot = 1 
             do k=2,kdm   !eta_anom at bottom is always zero and excluded

                if (blk_depth(i,j).gt.init_depth(i,j,k)) then  !when above bottom

                    eta_anom(i,j,k)=layer_depth(i,j,k)+init_depth(i,j,k)
                    kbot = k

! density perturbation at layer interface
! rhop1 has length kdm
                    rhop1(k) = -1.0*eta_anom(i,j,k)*drhodz(i,j,k)

                endif
              enddo

! mean perturbation denisty at layer center
! remove the time-mean density
              do k=1,kdm-1 
                 rhop2(k) = 0.5*(rhop1(k)+rhop1(k+1))
     &           + blk_sigt(i,j,k) - blk_mean_s(i,j,k)
              enddo
              rhop2(kdm) = 0.5*rhop1(kdm) 
     &           + blk_sigt(i,j,kdm) - blk_mean_s(i,j,kdm)

! integrate the perturbation pressure
! pres is kdm+1 
! pave is at layer centers
! pave is kdm 

              pres(:) = 0.0
              pave(:) = 0.0
              psum = 0.0 
              pbot = 0.0  
              do k=1,kdm

c                  if((i.eq.150).AND.(j.eq.150)) then
c       write(6,*)k,layer_depth(i,j,k),blk_depth(i,j),init_depth(i,j,k),
c     & rhop2(k),eta_anom(i,j,k)
c                  endif

! integrate the perturbation pressure
! using time-varying thickness is better
                pres(k+1) = pres(k) + grav*rhop2(k)*blk_thknss(i,j,k) 
                pave(k)   = (pres(k+1) + pres(k))*0.5                
                psum = psum + pave(k)*blk_thknss(i,j,k)

c                  if((i.eq.150).AND.(j.eq.150)) then
c       write(6,*)k,pave(k),psum,blk_thknss(i,j,k)
c                  endif

              enddo 

! remove depth-mean pressure
              pave(:) = pave(:) - psum/blk_depth(i,j)

! get bottom pressure
              pbot = pres(kbot+1) - psum/blk_depth(i,j)


c                  if((i.eq.150).AND.(j.eq.150)) then
c       write(6,*)l,pave,pbot
c                  endif

! compute pressure fluxes
              do k=1,kdm

! multiply with mean velocity
! correct for thickness at u and v locations??
                if(j.GT.1-mbdy.AND.j.LT.blkj+mbdy.AND.
     &           i.GT.1-mbdy.AND.i.LT.blki+mbdy) then 

                  duflux(i,j) = duflux(i,j) 
     &                  + 0.5*(blk_u(i,j,k) + blk_u(i+1,j,k))
     &                  *pave(k)*blk_thknss(i,j,k)
                  dvflux(i,j) = dvflux(i,j) 
     &                  + 0.5*(blk_v(i,j,k) + blk_v(i,j+1,k))
     &                  *pave(k)*blk_thknss(i,j,k)

c                  if (duflux(i,j).GT.1E+20) then
c      write(6,*)i,j,k,blk_depth(i,j),pave(k),pbot,
c     &      blk_u(i,j,k),blk_u(i+1,j,k),
c     &      blk_coast(i,j),blk_coast(i+1,j)         
c                  endif



                  if(k.eq.1) then
! conversion
                    conv(i,j) = pbot
     &                      *(0.5*(udhdx(i,j)+udhdx(i+1,j))            
     &                      + 0.5*(vdhdy(i,j)+vdhdy(i,j+1)))            
                  endif 

                endif

              enddo

              endif !blk_coast

            enddo
    
          enddo

c        write(6,*)'conv',l,maxval(abs(conv),1)
c        write(6,*)'fx',l,maxval(abs(duflux),1)
c        write(6,*)'fy',l,maxval(abs(dvflux),1)
c          write(6,*)'fx',l,duflux(100,:)
c          write(6,*)'de',l,blk_depth(100,:)
c          stop

!! ============================================================
! uscx vscy

!
!            Dlin_ubc(:,:) = 0.0
!	    Dlin_vbc(:,:) = 0.0
!	    Dlin_bc(:,:)  = 0.0
!
!            Dlin_ubt(:,:) = 0.0
!	    Dlin_vbt(:,:) = 0.0
!	    Dlin_bt(:,:)  = 0.0
!
!            Dqdr_ubc(:,:) = 0.0
!	    Dqdr_vbc(:,:) = 0.0
!	    Dqdr_bc(:,:)  = 0.0
!
!            Dqdr_ubt(:,:) = 0.0
!	    Dqdr_vbt(:,:) = 0.0
!	    Dqdr_bt(:,:)  = 0.0

c cross-terms
!            Dqdr_uu(:,:) = 0.0
!            Dqdr_vv(:,:) = 0.0
!            Dqdr_uq(:,:) = 0.0
!            Dqdr_vq(:,:) = 0.0
!            Dqdr_ux(:,:) = 0.0
!            Dqdr_vx(:,:) = 0.0
!            Dlin_uu(:,:) = 0.0
!            Dlin_vv(:,:) = 0.0
!            Dlin_uq(:,:) = 0.0
!            Dlin_vq(:,:) = 0.0
!            Dlin_ux(:,:) = 0.0
!            Dlin_vx(:,:) = 0.0
!
!            Dqdr_btq(:,:) = 0.0
!            Dqdr_bcq(:,:) = 0.0
!            Dqdr_x(:,:)   = 0.0
!            Dlin_btq(:,:) = 0.0
!            Dlin_bcq(:,:) = 0.0
!            Dlin_x(:,:)   = 0.0


!            do j=1-mbdy,blkj+mbdy-1
!
!              do i=1-mbdy,blki+mbdy-1
!
!! use bottom at u and v to compute layer thickness
!                depthu=min(blk_depth(i+1,j),blk_depth(i,j)) ! at second u
!                depthv=min(blk_depth(i,j+1),blk_depth(i,j)) 
!
!                do k=1,kdm ! k loop k - k - k - k - k - k - k - k - k - k - k - k - k 1
!
!! ============================================================================
!! quadratic drag 
!! ============================================================================
!
!! D = Cd*|u|*u*u
!! HYCOM: D = Cd*|u|/hbbl*1/hl*hl*u*u
!!          = Cd*|u|/hbbl*u*u  
!
!                  if (thkbot.GT.0.0) then
!! along u ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! top and bottom of layer in [m]
!                    ptopl=min(depthu,-0.5*(layer_depth(i+1,j,k) !layer_depth.<0
!     &                                    +layer_depth(i,j,k))) ! at second u
!                    pbotl=min(depthu,-0.5*(layer_depth(i+1,j,k+1  )
!     &                                    +layer_depth(i,j,k+1  )))
!
!! fraction of drag over the layer
!                    pbop=depthu-0.5*(thkbop(i,j)+thkbop(i+1,j)) ! top of bot. b.b.l.
!                    dragu = 0.5*( drag(i,j,l)+drag(i+1,j,l) )*  ! omit *qdpu!
!     &                      ( max(pbop,pbotl)-
!     &                        max(pbop,min(ptopl,pbotl-onemm)) )
!
!! baroclinic dissipation due to u*u'
!                    Dqdr_ubc(i+1,j) = Dqdr_ubc(i+1,j) + dragu*
!     &        (blk_u(i+1,j,k)+blk_ubarof(i+1,j,l))*blk_u(i+1,j,k)
!
!! barotropic dissipation due to u*U
!                    Dqdr_ubt(i+1,j) = Dqdr_ubt(i+1,j) + dragu*
!     &        (blk_u(i+1,j,k)+blk_ubarof(i+1,j,l))*blk_ubarof(i+1,j,l)
!
!
!! pure BT, BC, and cross terms =======================================
!                    Dqdr_uu(i+1,j) = Dqdr_uu(i+1,j) + dragu*
!     &        blk_u(i+1,j,k)*blk_u(i+1,j,k)
!
!                    Dqdr_uq(i+1,j) = Dqdr_uq(i+1,j) + dragu*
!     &        blk_ubarof(i+1,j,l)*blk_ubarof(i+1,j,l)
!
!                    Dqdr_ux(i+1,j) = Dqdr_ux(i+1,j) + dragu*
!     &        2*blk_u(i+1,j,k)*blk_ubarof(i+1,j,l)
!
!
!! test output
!!               if (j.EQ.63-3.AND.i.EQ.116-3) then
!!      write(6,*)k,hatu,max(pbop,pbotl),max(pbop,min(ptopl,pbotl-onemm))
!!     &  ,dragu,dragv
!!               call flush(6)
!!               end if
!
!! along v ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! top and bottom of layer in [m]
!                    ptopl=min(depthv,-0.5*(layer_depth(i,j+1,k) !layer_depth.<0
!     &                                    +layer_depth(i,j,k))) ! at second u
!                    pbotl=min(depthv,-0.5*(layer_depth(i,j+1,k+1  )
!     &                                    +layer_depth(i,j,k+1  )))
!
!! fraction of drag over the layer
!                    pbop=depthv-0.5*(thkbop(i,j)+thkbop(i,j+1)) ! top of bot. b.b.l.
!                    dragv = 0.5*( drag(i,j,l)+drag(i,j+1,l) )*  ! omit *qdpu!
!     &                      ( max(pbop,pbotl)-
!     &                        max(pbop,min(ptopl,pbotl-onemm)) )
!
!! baroclinic dissipation due to u*u'
!                    Dqdr_vbc(i,j+1) = Dqdr_vbc(i,j+1) + dragv*
!     &        (blk_v(i,j+1,k)+blk_vbarof(i,j+1,l))*blk_v(i,j+1,k)
!
!! barotropic dissipation due to u*U
!                    Dqdr_vbt(i,j+1) = Dqdr_vbt(i,j+1) + dragv*
!     &        (blk_v(i,j+1,k)+blk_vbarof(i,j+1,l))*blk_vbarof(i,j+1,l)
!
!
!! pure BT, BC, and cross terms =======================================
!                    Dqdr_vv(i,j+1) = Dqdr_vv(i,j+1) + dragv*
!     &        blk_v(i,j+1,k)*blk_v(i,j+1,k)
!
!                    Dqdr_vq(i,j+1) = Dqdr_vq(i,j+1) + dragv*
!     &        blk_vbarof(i,j+1,l)*blk_vbarof(i,j+1,l)
!
!                    Dqdr_vx(i,j+1) = Dqdr_vx(i,j+1) + dragv*
!     &        2*blk_v(i,j+1,k)*blk_vbarof(i,j+1,l)
!
!                  end if 
!     
!! ============================================================================           
!! linear drag 
!! ============================================================================
!
!! integrate over layer depths to get something like
!! C = rh [m/s]
!! D = rho0 * rh * ub^2  (ub is velocity in bbl)
!! D = rho0 * sum ul^2 * rh/500 * dh * hf/hatu * hatu   (hf=fraction of bbl above layer) 
!!   = rho0 * sum ul^2 * rh/500 * dh * hf               (hatu cancels!)
!! not in pressure coordinates, but in depth coordinates (confusing 1)
!! computation is for u at i+1 (noy i-1)                 (confusing 2)
!! layer_depth<0                                         (confusing 3)
!
!                  if (thkdrg.GT.0.0) then
!! along u ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!! top and bottom of layer in [m]
!                    ptopl=min(depthu,-0.5*(layer_depth(i+1,j,k) !layer_depth.<0
!     &                                    +layer_depth(i,j,k))) ! at second u
!                    pbotl=min(depthu,-0.5*(layer_depth(i+1,j,k+1  )
!     &                                    +layer_depth(i,j,k+1  )))
!
!! fraction of drag over the layer
!! do not need this  qdpu  = 1.0/max(hatu,onemm)
!! util5 = C/h
!! depthu is depth at u points
!! if util6 > depthu, pbop is negative: drag applied to entire layer
!! dragu = C/h*h in [m/s]
!                    pbop  = depthu - 0.5*(util6(i,j)+util6(i+1,j)) 
!                    dragu = 0.5*( util5(i,j)+util5(i+1,j) )*   ! omit *qdpu!
!     &                      ( max(pbop,pbotl)-
!     &                        max(pbop,min(ptopl,pbotl-onemm)) )
!
!! baroclinic dissipation due to u*u'
!                    Dlin_ubc(i+1,j) = Dlin_ubc(i+1,j) + dragu*
!     &        (blk_u(i+1,j,k)+blk_ubarof(i+1,j,l))*blk_u(i+1,j,k)
!
!! barotropic dissipation due to u*U
!                    Dlin_ubt(i+1,j) = Dlin_ubt(i+1,j) + dragu*
!     &        (blk_u(i+1,j,k)+blk_ubarof(i+1,j,l))*blk_ubarof(i+1,j,l)
!
!! pure BT, BC, and cross terms =======================================
!                    Dlin_uu(i+1,j) = Dlin_uu(i+1,j) + dragu*
!     &        blk_u(i+1,j,k)*blk_u(i+1,j,k)
!
!                    Dlin_uq(i+1,j) = Dlin_uq(i+1,j) + dragu*
!     &        blk_ubarof(i+1,j,l)*blk_ubarof(i+1,j,l)
!
!                    Dlin_ux(i+1,j) = Dlin_ux(i+1,j) + dragu*
!     &        2*blk_u(i+1,j,k)*blk_ubarof(i+1,j,l)
!
!
!! test output
!!                    if (j.eq.148.AND.i.eq.60) then
!!      write(6,*) k,hatu,max(pbop,pbotl),max(pbop,min(ptopl,pbotl-onemm))
!!     &  ,dragu
!!                    end if
!
!! along v ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    ptopl=min(depthv,-0.5*(layer_depth(i,j+1,k) !layer_depth.<0
!     &                                    +layer_depth(i,j,k))) ! at second u
!                    pbotl=min(depthv,-0.5*(layer_depth(i,j+1,k+1  )
!     &                                    +layer_depth(i,j,k+1  )))
!
!                    pbop  = depthv - 0.5*(util6(i,j)+util6(i,j+1)) 
!
!                    dragv = 0.5*( util5(i,j)+util5(i,j+1) )*   ! omit *qdpu!
!     &                      ( max(pbop,pbotl)-
!     &                        max(pbop,min(ptopl,pbotl-onemm)) )
!
!! baroclinic dissipation due to u*u' (Kang and Fringer 2011, JPO)
!                    Dlin_vbc(i,j+1) = Dlin_vbc(i,j+1) + dragv*
!     &        (blk_v(i,j+1,k)+blk_vbarof(i,j+1,l))*blk_v(i,j+1,k)
!
!! barotropic dissipation due to u*U (Kang and Fringer 2011, JPO) 
!                    Dlin_vbt(i,j+1) = Dlin_vbt(i,j+1) + dragv*
!     &        (blk_v(i,j+1,k)+blk_vbarof(i,j+1,l))*blk_vbarof(i,j+1,l)
!
!! pure BT, BC, and cross terms =======================================
!                    Dlin_vv(i,j+1) = Dlin_vv(i,j+1) + dragv*
!     &        blk_v(i,j+1,k)*blk_v(i,j+1,k)
!
!                    Dlin_vq(i,j+1) = Dlin_vq(i,j+1) + dragv*
!     &        blk_vbarof(i,j+1,l)*blk_vbarof(i,j+1,l)
!
!                    Dlin_vx(i,j+1) = Dlin_vx(i,j+1) + dragv*
!     &        2*blk_v(i,j+1,k)*blk_vbarof(i,j+1,l)
!
!                  end if 
!
!                enddo ! k loop k - k - k - k - k - k - k - k - k - k - k - k - k 1
!
!!               if (j.EQ.63-3.AND.i.EQ.116-3) then
!!       write(6,*) i,j,Dqdr_vbc(i,j+1),Dqdr_vbt(i,j+1),
!!     & Dqdr_vbc(i,j),Dqdr_vbt(i,j)
!!              call flush(6)
!!              stop
!!               end if
!
!
!
!              enddo ! i
!
!            enddo  ! j
!
!! now that we have the gradients at u and v points we use the 4-point average
!! to calculate the conversion rates at the pressure points for each time step
!! the same is true for dissipation D = C*u*u + C*v*v !!!
!
!            do j=1-mbdy,blkj+mbdy-1
!
!              do i=1-mbdy,blki+mbdy-1	      
!
!                Dlin_bc(i,j) = 0.5 * RHO_0
!     &                          * ( Dlin_ubc(i,j) + Dlin_ubc(i+1,j)
!     &                            + Dlin_vbc(i,j) + Dlin_vbc(i,j+1) )       
!
!                Dlin_bt(i,j) = 0.5 * RHO_0
!     &                          * ( Dlin_ubt(i,j) + Dlin_ubt(i+1,j)
!     &                            + Dlin_vbt(i,j) + Dlin_vbt(i,j+1) )       
!
!                Dlin_bcq(i,j) = 0.5 * RHO_0
!     &                          * ( Dlin_uu(i,j) + Dlin_uu(i+1,j)
!     &                            + Dlin_vv(i,j) + Dlin_vv(i,j+1) )       
!
!                Dlin_btq(i,j) = 0.5 * RHO_0
!     &                          * ( Dlin_uq(i,j) + Dlin_uq(i+1,j)
!     &                            + Dlin_vq(i,j) + Dlin_vq(i,j+1) )       
!
!                Dlin_x(i,j) = 0.5 * RHO_0
!     &                          * ( Dlin_ux(i,j) + Dlin_ux(i+1,j)
!     &                            + Dlin_vx(i,j) + Dlin_vx(i,j+1) )       
!
!
!                Dqdr_bc(i,j) = 0.5 * RHO_0
!     &                          * ( Dqdr_ubc(i,j) + Dqdr_ubc(i+1,j)
!     &                            + Dqdr_vbc(i,j) + Dqdr_vbc(i,j+1) )       
!
!                Dqdr_bt(i,j) = 0.5 * RHO_0
!     &                          * ( Dqdr_ubt(i,j) + Dqdr_ubt(i+1,j)
!     &                            + Dqdr_vbt(i,j) + Dqdr_vbt(i,j+1) )       
!    
!                Dqdr_bcq(i,j) = 0.5 * RHO_0
!     &                          * ( Dqdr_uu(i,j) + Dqdr_uu(i+1,j)
!     &                            + Dqdr_vv(i,j) + Dqdr_vv(i,j+1) )       
!
!                Dqdr_btq(i,j) = 0.5 * RHO_0
!     &                          * ( Dqdr_uq(i,j) + Dqdr_uq(i+1,j)
!     &                            + Dqdr_vq(i,j) + Dqdr_vq(i,j+1) )       
!
!                Dqdr_x(i,j) = 0.5 * RHO_0
!     &                          * ( Dqdr_ux(i,j) + Dqdr_ux(i+1,j)
!     &                            + Dqdr_vx(i,j) + Dqdr_vx(i,j+1) )       
!
!
!              enddo
!
!            enddo

!! ============================================================
!! calculate the time-mean values
!! ============================================================

            do j=1-mbdy,blkj+mbdy-1

    	      do i=1-mbdy,blki+mbdy-1


! keep track of the time-mean fields, avoid polluted parts due to ringing
                 if ((l.GT.RILE).AND.(l.LE.maxobs-RILE)) then
	           dufluxa(i,j) = dufluxa(i,j)+duflux(i,j)
                   dvfluxa(i,j) = dvfluxa(i,j)+dvflux(i,j)
                   conva(i,j)   = conva(i,j)  +conv(i,j)

!                   Dlin_bca(i,j)  = Dlin_bca(i,j)  + Dlin_bc(i,j)
!                   Dlin_bta(i,j)  = Dlin_bta(i,j)  + Dlin_bt(i,j)
!                   Dlin_bcqa(i,j) = Dlin_bcqa(i,j) + Dlin_bcq(i,j)
!                   Dlin_btqa(i,j) = Dlin_btqa(i,j) + Dlin_btq(i,j)
!                   Dlin_xa(i,j)   = Dlin_xa(i,j)   + Dlin_x(i,j)
!
!                   Dqdr_bca(i,j)  = Dqdr_bca(i,j)  + Dqdr_bc(i,j)
!                   Dqdr_bta(i,j)  = Dqdr_bta(i,j)  + Dqdr_bt(i,j)
!                   Dqdr_bcqa(i,j) = Dqdr_bcqa(i,j) + Dqdr_bcq(i,j)
!                   Dqdr_btqa(i,j) = Dqdr_btqa(i,j) + Dqdr_btq(i,j)
!                   Dqdr_xa(i,j)   = Dqdr_xa(i,j)   + Dqdr_x(i,j)

c          write(6,*)l,conva(1-mbdy:2,50),conva(blki-1:blki+mbdy,50) 
c          write(6,*)l,blk_depth(1-mbdy:2,50),
c     &               blk_depth(blki-1:blki+mbdy,50) 

                 endif 
   
 	       enddo
 
	     enddo   

            if (TIVAR.EQ.1) then
              write(107) conv	
	      write(110) duflux
	      write(111) dvflux
!              write(120) Dlin_bc
!              write(121) Dlin_bt
!              write(122) Dqdr_bc
!              write(123) Dqdr_bt
            endif  

        enddo ! maxobs

! finish computing mean
        conva     = conva/float(maxobs-2*RILE)
        dufluxa   = dufluxa/float(maxobs-2*RILE)
	dvfluxa   = dvfluxa/float(maxobs-2*RILE)
!        Dlin_bca  = Dlin_bca/float(maxobs-2*RILE)
!        Dlin_bta  = Dlin_bta/float(maxobs-2*RILE)
!        Dlin_bcqa = Dlin_bcqa/float(maxobs-2*RILE)
!        Dlin_btqa = Dlin_btqa/float(maxobs-2*RILE)
!        Dlin_xa   = Dlin_xa/float(maxobs-2*RILE)
!
!        Dqdr_bca  = Dqdr_bca/float(maxobs-2*RILE)
!        Dqdr_bta  = Dqdr_bta/float(maxobs-2*RILE)
!        Dqdr_bcqa = Dqdr_bcqa/float(maxobs-2*RILE)
!        Dqdr_btqa = Dqdr_btqa/float(maxobs-2*RILE)
!        Dqdr_xa   = Dqdr_xa/float(maxobs-2*RILE)

        if (TIVAR.EQ.1) then
          close(107)
          close(110)
          close(111)
!          close(120)
!          close(121)
!          close(122)
!          close(123)
        endif

        write(112) conva
        write(113) dufluxa
	write(114) dvfluxa

!        write(118) Dlin_bca
!	write(119) Dlin_bta
!
!        write(124) Dqdr_bca
!	write(125) Dqdr_bta
!
!        write(126) Dlin_bcqa
!        write(127) Dlin_btqa
!        write(128) Dlin_xa
!
!        write(129) Dqdr_bcqa
!        write(130) Dqdr_btqa
!        write(131) Dqdr_xa
!
!        write(6,*)'Dlin_bta max=',maxval(Dlin_bta)
!        write(6,*)'Dqdr_bta max=',maxval(Dqdr_bta)
        write(6,*)'maxval fxa',maxval(abs(dufluxa))
        write(6,*)'maxval fya',maxval(abs(dvfluxa))
        write(6,*)'maxval conva',maxval(abs(conva)),
     &                maxobs-2*RILE

        close(112)
        close(113)
        close(114)

!        close(118)
!        close(119)
!
!        close(124)
!        close(125)
!
!        close(126)
!        close(127)
!        close(128)
!        close(129)
!        close(130)
!        close(131)

      endif !for joe check  

!! loop over blocks ----------------------------------------------------
      enddo ! iblk
      enddo ! jblk
!! loop over blocks ----------------------------------------------------

      write(6,*) 'Finished...'

      stop

      end program PRES_CALC

C ---------------------------------------------------
      SUBROUTINE BPFILT(X,LX,NS,A,B,C,D,E,RI)
C APPLIES A BANDPASS FILTER TO ARRAY X
C CALLS SBR REVERS 
C BASED ON INFO IN STEARNS (1975)
C DON ALBERT 12/21/82

      REAL*8 A(NS),B(NS),C(NS),D(NS),E(NS)
      REAL*4 X(LX)
      REAL*4 Y(11,5)
      INTEGER LX2,RI,aa
      REAL*4 X2(LX+RI*4)
      PI=3.1415926536
C HANDLES FILTERS UP TO NS=10
C APPLIES THE FILTER TWICE, ONCE IN EACH DIRECTION, 
C SO THE PHASE IS ZERO.

c      write(6,*)'XB ',X
c      write(6,*)'As ',A,B,C,D,E

C padd with zeros on RIGHT side
C to avoid bad right side ringing
      aa=RI*4
      LX2 = LX+aa
      X2(:)    = 0.0
      X2(1:LX) = X(1:LX)

      IFLAG=0
5     CONTINUE
      IFLAG=IFLAG+1
c      write(6,*)'IFLAG',IFLAG
      IF(IFLAG.EQ.3) then
        X(1:LX)=X2(1:LX)
        RETURN
      ENDIF
      NS1=NS+1
      DO 10 I = 1,NS1
        DO 10 J=1,4
          Y(I,J)=0.0
10    CONTINUE
      DO 40 M=1,LX2
        Y(1,5)=X2(M)
        DO 20 N=1,NS
          Y(N+1,5)=A(N)*(Y(N,5)-2.*Y(N,3)+Y(N,1))
     &             -B(N)*Y(N+1,4)-C(N)*Y(N+1,3)
     &             -D(N)*Y(N+1,2)-E(N)*Y(N+1,1)
20       CONTINUE
         DO 30 N=1,NS1
           DO 30 MM=1,4
           Y(N,MM)=Y(N,MM+1)
30       CONTINUE
         X2(M)=Y(NS1,5)
40     CONTINUE
c       write(6,*)'XR',X 
       CALL   REVERS(X2,LX2)
       GOTO 5
       END

C ---------------------------------------------------
      SUBROUTINE REVERS(X,LX)
C REVERSES A TIME SERIES
C REFERENCE: ROBINSON (1967), p28
      DIMENSION X(LX)
c      write(6,*)'X2',LX,X

      L=LX/2
      DO 10 I=1,L
        J=LX-I
        TEMP=X(I)
        X(I)=X(J+1)
        X(J+1)=TEMP
10    CONTINUE
      RETURN
      END

C ---------------------------------------------------
      SUBROUTINE BPDES(F1,F2,DT,NS,A,B,C,D,E)
C BANDPASS BUTTERWORTH DIGITAL FILTER DESIGN SBR

      REAL*8 A(NS),B(NS),C(NS),D(NS),E(NS)
      PI=3.1415926536
      W1=SIN(F1*PI*DT)/COS(F1*PI*DT)
      W2=SIN(F2*PI*DT)/COS(F2*PI*DT)
c      write(6,*)'W1,W2',W1,W2
      WC=W2-W1
      Q=WC*WC+2.*W1*W2
      S=W1*W1*W2*W2
      DO 10 K=1,NS
        CS=COS(FLOAT(2*(K+NS)-1)*PI/FLOAT(4*NS))
        P=-2.*WC*CS
        R=P*W1*W2
        X=1.+ P + Q + R + S
        A(K)=WC*WC/X
        B(K)=( - 4. - 2.*P + 2.*R + 4.*S)/X
        C(K)=(6. - 2.*Q + 6.*S)/X
        D(K)=( - 4. + 2.*P - 2.*R + 4.*S)/X
        E(K)=(1. - P + Q - R + S)/X

c        write(6,*)'filt',K,CS,P,R,X,A(K)  

10    CONTINUE
      RETURN
      END



!-------------------------------------------------------
      SUBROUTINE TSERCON(TOUT,NDATA,MA,Ah,dt,omsel)
!-------------------------------------------------------
! construct time series for selected frequencies
      IMPLICIT NONE
      INTEGER i,j,NDATA,MA
      REAL TOUT(NDATA),dt,Ah(MA)
      DOUBLE PRECISION tt,omsel(MA/2)

      do i=1,NDATA
        TOUT(i)=0.0
        tt = i*dt/24.0d0  !in days
        do j=1,MA/2
          TOUT(i) = TOUT(i) + Ah(2*j-1)*cos(omsel(j)*tt)
     &                      + Ah(2*j)*sin(omsel(j)*tt)
        end do
      end do 
      RETURN
      END !-------------END OF HPDES--------
      
!-------------------------------------------------------
      SUBROUTINE AFUNC(VECT,NDATA,MA,tide_on,dt,omsel)
!-------------------------------------------------------
!     used in SVDFIT
!     makes the function matrix a*cos om*t + b*sin om*t 
!     horiz: MA funcs, vert: NDATA time
!     physical = logical matrix dimensions
!     conmat='11101100'  ! M2 S2 N2 K2 K1 O1 P1 Q1

      IMPLICIT NONE
      INTEGER i,j,k,NDATA,MA,numom
      REAL dt
      REAL VECT(NDATA,MA)
      CHARACTER conmat*8
      PARAMETER(numom=8)
      LOGICAL tide_on(numom)
      DOUBLE PRECISION tt,omsel(MA/2)
! Q1K2P1N2O1K1S2M2      
      DOUBLE PRECISION,parameter,dimension(numom) :: om = [
     & 12.1408331767d0,   !1  M2 rad/day
     & 12.5663706144d0,   !2  S2
     & 6.3003880821d0,    !3  K1 
     & 5.8404450946d0,    !4  O1
     & 11.9128060356d0,   !5  N2
     & 6.2659825322d0,    !6  P1 
     & 12.6007762061d0,   !7  K2
     & 5.6124179535d0]    !8  Q1

!      write(6,*)'here?',om(1),om(8),NDATA,MA,conmat,dt 

      do i=1,NDATA
         tt = i*dt/24.0d0  !in days
c         write(6,*) tt
         k=0
         do j=1,numom
            if (tide_on(j)) then
               k=k+1
!           a*cos omt + b*sin omt
               VECT(i,2*k-1) = cos(om(j)*tt) 
               VECT(i,2*k)   = sin(om(j)*tt)
            endif
         enddo
c         write(6,*)i,k,VECT(i,:)
      enddo

! return selected frequencies in rad/day
      k=0
      do j=1,numom
        if (tide_on(j)) then
          k=k+1
          omsel(k)=om(j)
        endif
      enddo
      END

!-------------------------------------------------------
      SUBROUTINE LFIT(Y,SIG,NDATA,A,MA,LISTA,MFIT,NCVM,VECT)
c      SUBROUTINE LFIT(X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,NCVM,CHISQ,FUN
c     *CS)
!-------------------------------------------------------
! from numerical recipes
      PARAMETER (MMAX=50)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),LISTA(MFIT),
     *    COVAR(NCVM,NCVM),BETA(MMAX),AFUNC(MMAX)
      REAL VECT(NDATA, MA)

      KK=MFIT+1
      DO 12 J=1,MA
        IHIT=0
        DO 11 K=1,MFIT
          IF (LISTA(K).EQ.J) IHIT=IHIT+1
11      CONTINUE
        IF (IHIT.EQ.0) THEN
          LISTA(KK)=J
          KK=KK+1
        ELSE IF (IHIT.GT.1) THEN
c          PAUSE 'Improper set in LISTA'
          write(6,*) 'Improper set in LISTA'
        ENDIF
12    CONTINUE
c      IF (KK.NE.(MA+1)) PAUSE 'Improper set in LISTA'
      IF (KK.NE.(MA+1)) write(6,*) 'Improper set in LISTA'
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE

! extract data from function
      DO 18 I=1,NDATA
!        CALL FUNCS(X(I),AFUNC,MA)
        AFUNC(1:MA) = VECT(I,:)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO 15 J=MFIT+1,MA
            YM=YM-A(LISTA(J))*AFUNC(LISTA(J))
15        CONTINUE
        ENDIF
        SIG2I=1./SIG(I)**2
        DO 17 J=1,MFIT
          WT=AFUNC(LISTA(J))*SIG2I
          DO 16 K=1,J
            COVAR(J,K)=COVAR(J,K)+WT*AFUNC(LISTA(K))
16        CONTINUE
          BETA(J)=BETA(J)+YM*WT
17      CONTINUE
18    CONTINUE
      IF (MFIT.GT.1) THEN
        DO 21 J=2,MFIT
          DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
19        CONTINUE
21      CONTINUE
      ENDIF
      CALL GAUSSJ(COVAR,MFIT,NCVM,BETA,1,1)
      DO 22 J=1,MFIT
        A(LISTA(J))=BETA(J)
22    CONTINUE
      CHISQ=0.
c      DO 24 I=1,NDATA
c        CALL FUNCS(X(I),AFUNC,MA)
c        SUM=0.
c        DO 23 J=1,MA
c          SUM=SUM+A(J)*AFUNC(J)
c23      CONTINUE
c        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
c24    CONTINUE
c      CALL COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
      RETURN
      END

!-------------------------------------------------------
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
!-------------------------------------------------------
! from numerical recipes
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
c                PAUSE 'Singular matrix'
                write(6,*) 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
c        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        IF (A(ICOL,ICOL).EQ.0.) write(6,*) 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

!-------------------------------------------------------
      SUBROUTINE COVSRT(COVAR,NCVM,MA,LISTA,MFIT)
!-------------------------------------------------------
! from numerical recipes
      DIMENSION COVAR(NCVM,NCVM),LISTA(MFIT)
      DO 12 J=1,MA-1
        DO 11 I=J+1,MA
          COVAR(I,J)=0.
11      CONTINUE
12    CONTINUE
      DO 14 I=1,MFIT-1
        DO 13 J=I+1,MFIT
          IF(LISTA(J).GT.LISTA(I)) THEN
            COVAR(LISTA(J),LISTA(I))=COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J))=COVAR(I,J)
          ENDIF
13      CONTINUE
14    CONTINUE
      SWAP=COVAR(1,1)
      DO 15 J=1,MA
        COVAR(1,J)=COVAR(J,J)
        COVAR(J,J)=0.
15    CONTINUE
      COVAR(LISTA(1),LISTA(1))=SWAP
      DO 16 J=2,MFIT
        COVAR(LISTA(J),LISTA(J))=COVAR(1,J)
16    CONTINUE
      DO 18 J=2,MA
        DO 17 I=1,J-1
          COVAR(I,J)=COVAR(J,I)
17      CONTINUE
18    CONTINUE
      RETURN
      END
