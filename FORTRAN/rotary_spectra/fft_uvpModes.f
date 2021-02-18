      PROGRAM FFT_UVPMODES
! This program computes the fft of modal amplitudes u, v and p. 
!
! Run script: fft_190.com
! This script should be modified for each run. It contains the model 
! used (GLBc_0.04), experiment/runid (19.0/190) and tiles over which 
! the fft is computed (iblks-iblke,jblks-jblke). It also sets the 
! path to input/output, jobid and other options.  
! 
! Input files:
! These are the files needed to run the program. They include the grid's 
! sea/land mask and modal amplitudes. 
!    griddata/coast_*        - Grid sea/land mask
!    modeslayeramp/modAU_*   - Modal amplitude U 
!    modeslayeramp/modAV_*   - Modal amplitude V 
!    modeslayeramp/modAP_*   - Modal amplitude P 
!
! Output files: 
!    out_spectral/fft_*_umod - fft of modal amplitude U 
!    out_spectral/fft_*_vmod - fft of modal amplitude V 
!    out_spectral/fft_*_pmod - fft of modal amplitude P 
!
! fft_uvpmodes.f; Miguel S Solano, USM, 2021/02/16

!---- Variable declaration
      implicit none

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
      real omf                      ! angular velocity
      real fcor,pi                  ! coriolis frequency, pi
      integer NS                    ! filter params
      integer TIVAR                 ! (=1) then time varying fields are stored
      integer CB_VAR                ! if 1, use spatially varying drag 
      integer jblks,jblke !1,17
      integer iblks,iblke !1,10
      integer MEIG
      character runid*3

      parameter(MEIG=5)
      parameter(TIVAR=0)              ! (=1) then time varying fields are stored
      parameter(CB_VAR=1)              ! in case of 0, use constant drag  

      ! Global HYCOM parameters (GLBc0.04)
      parameter(idm=9000,jdm=7000,kdm=41)
      parameter(blki=150,blkj=200)    ! GLBc0.04
      parameter(mbdy=3)               ! halo
      parameter(grav=9.806)
      parameter(thktol=0.0)
      parameter(thkdrg=500.0)
      parameter(RHO_0=1034.0)
      parameter(onemm=0.001)
      parameter(cbmin=0.0025,vonk=0.4,z0=0.01,thkbot=10.0) !cbmin=0.0025, vonk=0.4, z0=10mm (0.01) 
      parameter(flag=2.0**100)
      parameter(omf=7.292115900231276/100000.0)

      character ftnm*10 
      character infil10*240,infil11*240,infil12*240
      character infil13*240,infil14*240,infil15*240,infil16*240
      character outfile11*60,outfile12*60,outfile13*60
      character*2 numstr(90)
      character*2 constr(8)

      integer i,j,k,m,l,iblk,jblk,joe,joe2,kbot
      integer*8 plan

! 2D variables (grid)
      integer blk_coast
      dimension blk_coast(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

! 4D variables (uvp modal amplitudes and their fft)
      real zz,blk_mean_u,blk_mean_v, blk_mean_p
      real, allocatable :: fft_modau(:,:,:,:),fft_modav(:,:,:,:)
      real, allocatable :: fft_modap(:,:,:,:)

      real, allocatable :: modAU(:,:,:,:)  ! (x,y,mod,time)
      real, allocatable :: modAV(:,:,:,:)
      real, allocatable :: modAP(:,:,:,:)
      real, allocatable :: var(:)

! FFTW subroutine
      double precision, allocatable :: t1(:),f1(:),pp1(:) 
      double precision, allocatable :: ufft(:),vfft(:),pfft(:) 
      double complex, allocatable :: u1(:),v1(:),p1(:)
      parameter(pi=3.141592654d0)

! Numstr
      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/


!---- READ INPUT FILE
! Tile numbers are read from input file created from submit script(.com)
! Make sure input file matches the read(*,*) order below. 

! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid

! fix maxobs and numcon
      read(*,*) maxobs   !integer
      read(*,*) ftnm     !string

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      write (6,'(a,i5)') 'maxobs: ',  maxobs
      call flush(6)

      write (6,'(a,a)') ' ftnm:',ftnm
      call flush(6)


!----------------------------------------------------------------------
! Allocate variables

      allocate (t1(maxobs),f1(maxobs),pp1(maxobs))
      allocate (u1(maxobs),v1(maxobs),p1(maxobs))
      allocate (ufft(maxobs),vfft(maxobs),pfft(maxobs))

      allocate(var(1:maxobs), 
     & fft_modau(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG),
     & fft_modav(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG),
     & fft_modap(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG),
     & modAU(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & modAV(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & modAP(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs))


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


!!--------------------------------------------------------------------------
!! OPEN and load time varying 3D files          
!!--------------------------------------------------------------------------
! Modal amplitude (U) 
          infil12='modeslayramp/'
     &  //'modAU_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'reading file: ',infil12
        open(unit=12,file=infil12,status='old',form='unformatted')

! Modal amplitude (V) 
          infil13='modeslayramp/'
     &  //'modAV_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'reading file: ',infil13
        open(unit=13,file=infil13,status='old',form='unformatted')

! Modal amplitude (P) 
          infil14='modeslayramp/'
     &  //'modAP_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'reading file: ',infil14
        open(unit=14,file=infil14,status='old',form='unformatted')


!!--------------------------------------------------------------------------
! load time-varying modal amplitudes  
!!--------------------------------------------------------------------------
      modAU(:,:,:,:) = 0.0
      modAV(:,:,:,:) = 0.0
      modAP(:,:,:,:) = 0.0

      do l=1,maxobs

          read(12) modAU(:,:,:,l)
          read(13) modAV(:,:,:,l)
          read(14) modAP(:,:,:,l)

      enddo ! maxobs

      close(12)
      close(13)
      close(14)

      write(6,*) 'finished reading mode amps'
      call flush(6)

      fft_modau(:,:,:,:)=0.0
      fft_modav(:,:,:,:)=0.0
      fft_modap(:,:,:,:)=0.0

!---------Calling rot_spec_fftw-------------------------------------
      ! loop over modes
      do m=1,MEIG

         do l=1,maxobs
           t1(l)=dble(l) 
         enddo!-----second maxobs!     
            
         do j=1-mbdy,blkj+mbdy !j->latitude
             zz=0d0

             do i=1-mbdy,blki+mbdy !i->longitude
         blk_mean_u = sum(modAU(i,j,m,:))/  !mean of time series 
     &                 (max(1,size(modAU(i,j,m,:))))
         modAU(i,j,m,:)= modAU(i,j,m,:)-blk_mean_u  !subtract mean       

         blk_mean_v = sum(modAV(i,j,m,:))/  !mean of time series 
     &                 (max(1,size(modAV(i,j,m,:))))
         modAV(i,j,m,:)= modAV(i,j,m,:)-blk_mean_v  !subtract mean       

         blk_mean_p = sum(modAP(i,j,m,:))/  !mean of time series 
     &                 (max(1,size(modAP(i,j,m,:))))
         modAP(i,j,m,:)= modAP(i,j,m,:)-blk_mean_p  !subtract mean       

                if (blk_coast(i,j).eq.1) then
                u1 = dble(modAU(i,j,m,:))
                v1 = dble(modAV(i,j,m,:))
                p1 = dble(modAP(i,j,m,:))

                call COMPUTE_FFTW(u1,t1,maxobs,ufft,f1,pp1)
                call COMPUTE_FFTW(v1,t1,maxobs,vfft,f1,pp1)
                call COMPUTE_FFTW(p1,t1,maxobs,pfft,f1,pp1)
                 
                fft_modau(i,j,:,m) = ufft
                fft_modav(i,j,:,m) = vfft
                fft_modap(i,j,:,m) = pfft

                zz = zz+1d0

                endif ! blk_coast
             enddo   !i
           enddo     !j 

         enddo ! modes


!!--------------------------------------------------------------------------
!! OPEN 2D output files         
!!--------------------------------------------------------------------------
        outfile11='out_spectral/'//'fft_'//runid//
     &'_umod_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        outfile12='out_spectral/'//'fft_'//runid//
     &'_vmod_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        outfile13='out_spectral/'//'fft_'//runid//
     &'_pmod_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=111,file=outfile11,status='replace',
     &      form='unformatted',access='sequential')
        open(unit=112,file=outfile12,status='replace',
     &      form='unformatted',access='sequential')
        open(unit=113,file=outfile13,status='replace',
     &      form='unformatted',access='sequential')


!-----------write output--------------------------------------------------
      do m=1,MEIG 
         write(111) fft_modau(:,:,:,m)
         write(112) fft_modav(:,:,:,m)
         write(113) fft_modap(:,:,:,m)
      enddo ! MEIG

      close(111)
      close(112)
      close(113)


!! loop over blocks ----------------------------------------------------
      enddo ! iblk
      enddo ! jblk

      write(6,*) 'Finished...'

      stop

      END ! PROGRAM FFT_UVPMODES



!-------------------------------------------------------------------------------
      SUBROUTINE COMPUTE_FFTW(u1,t1,n,ufft,f1,pp1)
!-------------------------------------------------------------------------------
! COMPUTE_FFTW computes the fast fourier transform (ufft) of variable
! u1. 
! Miguel Solano, 2021/02/16

        implicit none
        include 'fftw3.f'
        integer n,k
        integer*8 plan
        double precision del,t1(n),f1(n),pp1(n),ufft(n)
        double complex u1f(n), u1(n)

        call dfftw_plan_dft_1d(plan,n,u1,u1f,1,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan,u1,u1f)
        call dfftw_destroy_plan(plan)

        ufft = dreal(u1f)
        
        del = t1(2) - t1(1)  
        do k = 1,n/2+1
           f1(k) = dble(k-1)/dble(n)/del
        enddo
        do k = n/2+2,n
           f1(k) = dble(k-1-n)/dble(n)/del
        enddo

        pp1 = 1/f1


      END ! SUBROUTINE COMPUTE_FFTW
