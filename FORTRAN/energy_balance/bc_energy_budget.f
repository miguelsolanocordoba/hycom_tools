      PROGRAM BC_ENERGY_BUDGET

! bc_energy_budget.f - Computes the baroclinic energy budget in HYCOM 
! Created by: Miguel Solano, USM, 2021/08/03
! Modified from pres_calc.f and KE_APE.f by Maarten Buijsman
! 
! This program reads in tiled output from HYCOM and computes the
! time-mean and depth-averaged baroclinic (BC) energy budget (W/m):  
! 
!                 div(F') = C - epsilon' - D'
! 
! Where F is the energy flux in (W/m^2), C is the barotropic (BT) to 
! BC conversion, epsilon is the turbulent dissipation, D is the 
! dissipation due to bottom drag and the primes(') denote BC 
! quantities due to wave motion. 
! For details see Kang and Fringer (2012). 


!=====================================================================!
! --------------------- Variable Declaration ------------------------ ! 
!=====================================================================!
      IMPLICIT NONE

      integer idm,jdm,kdm           ! Global grid tiles (idm,jdm,kdm)
      integer jblks,jblke           ! y-tile start(jblks) and end (jblke)
      integer iblks,iblke           ! x-tile start(iblks) and end (iblke)
      integer blki,blkj             ! Tile size (blki,blkj)
      integer maxobs                ! Number of observations
      integer mbdy                  ! Halo size
      integer NS                    ! Filter order parameter
      integer RILE                  ! Number of time steps to omit to reduce effect of ringing
      integer DO_FILT               ! Filtering (1= band-pass filter, 0= no filter)

      integer i,j,k,l,iblk,jblk,joe,joe2,kbot
      integer blk_coast

      real grav                     ! gravitational acceleration
      real RHO_0                    ! reference density
      real flag                     ! missing value flag
      real thktol                   ! tolerance for layer thickness.
      real thkdrg                   ! bottom boundary layer of lin drag
      real drgscl                   ! scale factor
      real cbmin,vonk,z0,cbp        ! for quadratic drag computation
      real thkbot                   ! thickness of bottom quadraticboundary layer (m)
      real onemm                    ! one mm 
      real FC1,FC2,dt               ! filter params  

      character runid*3
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

! tides information -------------------------------------------------    
      integer tidcon1,numcon,numcon2,ncon,ipt
      integer tidcon
      parameter(ncon=8)
      CHARACTER*8   tidcons
      CHARACTER*2   TideMode(ncon)
      CHARACTER*24  Tides
      DATA TideMode/'M2','S2','K1','O1','N2','P1','K2','Q1'/
      LOGICAL tide_on(ncon)


c--------filter parameters-------------------
c 3D
      real, allocatable :: blk_t_filt(:,:,:,:),
     & blk_u_filt(:,:,:,:),
     & blk_v_filt(:,:,:,:),
     & blk_s_filt(:,:,:,:) 

c 2D
      real, allocatable :: blk_ubarof(:,:,:),blk_vbarof(:,:,:)
!1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60'/

      INTEGER MA !MA=numcon*2
c      PARAMETER(MA=numcon*2)

      REAL, allocatable :: Yh(:), Eh(:) !maxobs
      REAL, allocatable :: VECTh(:,:) !maxobs,MA
      REAL, allocatable :: Ah(:) !MA
      DOUBLE PRECISION, allocatable :: omsel(:) !numcon
      INTEGER, allocatable :: LISTA(:) !MA
!=====================================================================!
! ------------------ Memory Allocation (static) --------------------- ! 
!=====================================================================!



!=====================================================================!
! ----------------------- Read Input File --------------------------- ! 
!=====================================================================!
! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid

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
      else
        write(6,'(a,a)') 'NO filtering: ',ftnm
      endif
      call flush(6)

      write(6,'(a,i2)') 'number of time steps excluded: ',RILE
      call flush(6)

! Invert band-pass hours
      FC1 = 1/FC1
      FC2 = 1/FC2


!=====================================================================!
! ----------------------- Tide Constituents ------------------------- ! 
!=====================================================================!
! Determine what tidal constituents are included. 
! tidecon 1 digit per consituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
      tidcon1 = tidcon
      do i=1,ncon
        tide_on(i) = mod(tidcon1,10) .eq. 1
        tidcon1    = tidcon1/10      ! shift by one decimal digit
      enddo

      TIDES='                        '
      ipt=1
      numcon2=0
      do i=1,ncon
        if(tide_on(i))then
           TIDES(ipt:ipt+1)=TideMode(i)
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

!=====================================================================!
! ------------------ Memory Allocation (dynamic) -------------------- ! 
!=====================================================================!
      allocate(blk_t_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_u_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_v_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs),
     & blk_s_filt(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,maxobs))
      allocate(blk_ubarof(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_vbarof(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs))
      allocate(Yh(maxobs),Eh(maxobs))
      allocate(VECTh(maxobs,MA))
      allocate(Ah(MA))
      allocate(omsel(numcon))
      allocate(LISTA(MA))

      Eh(:)=1.0
      do i=1,MA
        LISTA(i)=i
      enddo

!=====================================================================!
! ------------------ Read Grid and Output files --------------------- ! 
!=====================================================================!
! Read HYCOM tiled grid and output files.  

! MAIN LOOP over tiles specified in input file (iblks-iblke,jblks-jblke)
      do jblk=jblks,jblke
        do iblk=iblks,iblke

          write(6,*) 'blk_',jblk,iblk
          call flush(6)

          write(6,*) 'Reading grid...'
          call flush(6)


!! Grid files (2D):(x,y) 
! coast
         infil10='griddata/coast_'//runid//
     &           '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          write(6,*) 'Reading: ',infil10
          open(unit=10,file=infil10,status='old',form='unformatted')
          read(10) blk_coast
          close(10)

! depth
         infil11='griddata/'//'depth_'//runid//
     &           '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF' 
          write(6,*) 'Reading: ',infil11     
          open(unit=11,file=infil11,status='old',form='unformatted')
          read(11) blk_depth
          write(6,*) 'Max depth = ',maxval(maxval(blk_depth(:,:),1))
          close(11)

! uscx
         infil12='griddata/uscx_'//runid//
     &           '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          write(6,*) 'Reading: ',infil12
          open(unit=12,file=infil12,status='old',form='unformatted')
          read(12) blk_uscx
          close(12)

! vscy
         infil13='griddata/vscy_'//runid//
     &           '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          write(6,*) 'read',infil13
          open(unit=13,file=infil13,status='old',form='unformatted')
          read(13) blk_vscy


!! Output variables (3D):(x,y,t)
! ubaro (barotropic velocity in x-dir)
          infil18='ubaro/ubaro_'//runid//
     &          '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=18,file=infil18,status='old',form='unformatted')

! vbaro (barotropic velocity in y-dir)
          infil19='vbaro/vbaro_'//runid//
     &          '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=19,file=infil19,status='old',form='unformatted')

!! Output variables (4D):(x,y,z,t)
! sig (sigma2 density) 
          infil14='sig/sig_'//runid//
     &            '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=14,file=infil14,status='old',form='unformatted')

! thknss (layer thickness)
          infil15='thknss/thknss_'//runid//
     &            '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=15,file=infil15,status='old',form='unformatted')

! uiso (baroclinic velocity in x-dir) 
          infil16='u_iso/u_'//runid//
     &            '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=16,file=infil16,status='old',form='unformatted')

! viso (baroclinic velocity in y-dir) 
          infil17='v_iso/v_'//runid//
     &            '_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
          open(unit=17,file=infil17,status='old',form='unformatted')


!=====================================================================!
! ----------------------- Bandpass Filtering ------------------------ ! 
!=====================================================================!
! Initialize all variables to 0 

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

            enddo ! kdm
          enddo ! maxobs

         close(12)
         close(13)
         close(14)
         close(15)
         close(16)
         close(17)

! Set land values (blk_depth=flag) to 0. Needed for filtering. 
         do l=1,maxobs
           do j=1-mbdy,blkj+mbdy
             do i=1-mbdy,blki+mbdy
        
               if (l.eq.1.AND.blk_depth(i,j).EQ.flag) then
                 blk_depth(i,j)=0.0
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

         write(6,*) 'Starting filtering ...'

! find filter constants
! used in SUBROUTINE BPFILT(X,LX,NS,A,B,C,D,E,RILE)
         call BPDES(FC1,FC2,dt,NS,A,B,C,D,E)

         blk_mean_t(:,:,:) = 0.0
         blk_mean_s(:,:,:) = 0.0

         do j=1-mbdy,blkj+mbdy     !j->latiude
           do i=1-mbdy,blki+mbdy   !i->longitude
             if (blk_coast(i,j).eq.1) then

!! Ubar 
               blk_mean_ub = sum(blk_ubarof(i,j,:))/  !mean of time series 
     &                        (max(1,size(blk_ubarof(i,j,:))))
               blk_ubarof(i,j,:)= blk_ubarof(i,j,:)-blk_mean_ub  !subtractmean       
               CALL BPFILT(blk_ubarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)


!! Vbar 
               blk_mean_vb = sum(blk_vbarof(i,j,:))/             !mean of time series
     &                        (max(1,size(blk_vbarof(i,j,:))))
               blk_vbarof(i,j,:)= blk_vbarof(i,j,:)-blk_mean_vb  !subtract mean       
               CALL BPFILT(blk_vbarof(i,j,:),maxobs,NS,A,B,C,D,E,RILE)


! Loop over depth to filter 4D variables (uiso,viso,sigma2,thknss)
               do k=1,kdm
! only filter when the layer is above the bottom

                 blk_mean_t(i,j,k) = sum(blk_t_filt(i,j,k,:))/    
     &                      (max(1,size(blk_t_filt(i,j,k,:))))

                 blk_mean_s(i,j,k) = sum(blk_s_filt(i,j,k,:))/ 
     &                      (max(1,size(blk_s_filt(i,j,k,:))))

!----------------------------------------------------------------        
!! Thnkss
! only filter when the layer has a nonzero thickness
                 if(blk_mean_t(i,j,k).gt.0.0) then

! speed it up by omitting constant layer thicknesses
         if(abs(blk_mean_t(i,j,k)-blk_t_filt(i,j,k,1)).GT.0.001.OR.
     &      abs(blk_mean_t(i,j,k)-blk_t_filt(i,j,k,maxobs)).GT.0.001)
     &   then

           blk_t_filt(i,j,k,:)= blk_t_filt(i,j,k,:)-blk_mean_t(i,j,k)
           CALL BPFILT(blk_t_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)
           blk_t_filt(i,j,k,:)=blk_t_filt(i,j,k,:)+blk_mean_t(i,j,k)
         endif ! thknss>0.0

!! Sigma2
! speed it up by omitting layers of constant density
         if(abs(blk_mean_s(i,j,k)-blk_s_filt(i,j,k,1)).GT.0.001.OR.
     &      abs(blk_mean_s(i,j,k)-blk_s_filt(i,j,k,maxobs)).GT.0.001)
     &   then

           blk_s_filt(i,j,k,:)= blk_s_filt(i,j,k,:)-blk_mean_s(i,j,k)
           CALL BPFILT(blk_s_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)
           blk_s_filt(i,j,k,:)=blk_s_filt(i,j,k,:)+blk_mean_s(i,j,k)

         endif ! thknss>0.0

         !baroclinic velocity filtering   

!! Uiso
         blk_mean_u = sum(blk_u_filt(i,j,k,:))/           
     &                 (max(1,size(blk_u_filt(i,j,k,:))))
         blk_u_filt(i,j,k,:)= blk_u_filt(i,j,k,:)-blk_mean_u 

           CALL BPFILT(blk_u_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

!! Viso
         blk_mean_v = sum(blk_v_filt(i,j,k,:))/             
     &                 (max(1,size(blk_v_filt(i,j,k,:))))
         blk_v_filt(i,j,k,:)= blk_v_filt(i,j,k,:)-blk_mean_v 
           CALL BPFILT(blk_v_filt(i,j,k,:),maxobs,NS,A,B,C,D,E,RILE)

         else
           blk_s_filt(i,j,k,:)=0.0
           blk_u_filt(i,j,k,:)=0.0
           blk_v_filt(i,j,k,:)=0.0
           blk_t_filt(i,j,k,:)=0.0
         endif ! thickness criterium

               enddo !k

               else !----coast

                 blk_s_filt(i,j,:,:)=0.0
                 blk_u_filt(i,j,:,:)=0.0
                 blk_v_filt(i,j,:,:)=0.0
                 blk_t_filt(i,j,:,:)=0.0
                 blk_ubarof(i,j,:)=0.0
                 blk_vbarof(i,j,:)=0.0

             endif !coast                
           enddo ! i       
         enddo  ! j

         write(6,*) 'Finished filtering...'
         call flush(6)


!=====================================================================!
! ---------------------- Compute Pressure Flux ---------------------- ! 
!=====================================================================!
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

        init_depth(:,:,:) = 0.0
        init_dep_rho(:,:,:) = 0.0
        drhodz(:,:,:) = 0.0

        do j=1-mbdy,blkj+mbdy
          do i=1-mbdy,blki+mbdy

            init_depth(i,j,kdm+1)=blk_depth(i,j)

! init_depth should be positive
            do k=kdm,2,-1
              init_depth(i,j,k) = init_depth(i,j,k+1)-blk_mean_t(i,j,k)
            enddo

! compute the mean (resting) depth of the layer centers
! init_dep_rho is also positive. length is kdm
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

            enddo !kdm
          enddo !blki
        enddo !blkj

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

              if (blk_coast(i,j).EQ.1) then ! apply everywhere to speed
up?

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
             do k=2,kdm   !eta_anom at bottom is always zero and
excluded

                if (blk_depth(i,j).gt.init_depth(i,j,k)) then  !when
above bottom

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
     ! using time-varying thickness is better
     ! pres is kdm+1 
     ! pave is at layer centers
     ! pave is kdm 

              pres(:) = 0.0
              pave(:) = 0.0
              psum = 0.0
              pbot = 0.0

              do k=1,kdm
                pres(k+1) = pres(k) + grav*rhop2(k)*blk_thknss(i,j,k)
                pave(k)   = (pres(k+1) + pres(k))*0.5
                psum = psum + pave(k)*blk_thknss(i,j,k)

              enddo

     ! remove depth-mean pressure
              pave(:) = pave(:) - psum/blk_depth(i,j)

     ! get bottom pressure
              pbot = pres(kbot+1) - psum/blk_depth(i,j)

     ! compue conversion
                    conv(i,j) = pbot
     &                      *(0.5*(udhdx(i,j)+udhdx(i+1,j))
     &                      + 0.5*(vdhdy(i,j)+vdhdy(i,j+1)))

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

                  if(k.eq.1) then
                  endif
                endif

              enddo

              endif !blk_coast

            enddo
          enddo





!=====================================================================!
! ----------------------- Compute KE and APE ------------------------ ! 
!=====================================================================!

!=====================================================================!
! ---------------------- Compute Advective Flux --------------------- ! 
!=====================================================================!

!=====================================================================!
! ----------------------- Compute Conversion ------------------------ ! 
!=====================================================================!

      END PROGRAM BC_ENERGY_BUDGET




!=====================================================================!
! ----------------------- SUBROUTINES (BPF) ------------------------- ! 
!=====================================================================!

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

C padd with zeros on RIGHT side
C to avoid bad right side ringing
      aa=RI*4
      LX2 = LX+aa
      X2(:)    = 0.0
      X2(1:LX) = X(1:LX)

      IFLAG=0
5     CONTINUE
      IFLAG=IFLAG+1
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
       CALL   REVERS(X2,LX2)
       GOTO 5
       END

C ---------------------------------------------------
      SUBROUTINE REVERS(X,LX)
C REVERSES A TIME SERIES
C REFERENCE: ROBINSON (1967), p28
      DIMENSION X(LX)

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

10    CONTINUE
      RETURN
      END
