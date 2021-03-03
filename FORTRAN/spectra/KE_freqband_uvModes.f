      PROGRAM ROTARY_SPECTRA_UVMODES
! rotary_spectra_v01.f; Keshav J Raja, USM, 2020/08/25
! rotary_spectra_uvmodes.f; Miguel S Solano, USM, 2021/02/17
!
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
      real FC1,FC2,FC11,FC22,dt     ! filter params  
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
      character*2 numstr(90)
      character*2 constr(8)

      integer i,j,k,m,l,iblk,jblk,joe,joe2,kbot
      integer*8 plan

! 2D variables
      integer blk_coast
      dimension blk_coast(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      real blk_depth,blk_ubaro,blk_vbaro
      dimension blk_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_vbaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_ubaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

! ---2d+t variables-----
      real, allocatable :: blk_u_surf(:,:,:),blk_v_surf(:,:,:)

      real zz,blk_mean_u,blk_mean_v
      real,allocatable :: blk_cw(:,:,:,:),blk_ccw(:,:,:,:)
      real,allocatable :: cw1(:,:,:,:),ccw1(:,:,:,:)

! 3D variables
      real blk_u,blk_v
      dimension blk_u(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_v(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)

c 4D variables
      real blk_Uvec
      dimension blk_Uvec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,MEIG)

      real, allocatable :: modAU(:,:,:,:)  ! modal amplitudes x,y,mod,time
      real, allocatable :: modAV(:,:,:,:)
      real, allocatable :: var(:)

! Frequency bins and indices
      integer find1s,find1e,find2s,find2e
      integer find3s,find3e,find4s,find4e
!      parameter(find1s=1, find1e=21,find2s=24,find2e=36)
!      parameter(find3s=48,find3e=72,find4s=75,find4e=360)

      integer finds,finde
!      parameter(finds=1, finde=25) !D0 [0.033 - 0.833]
!      parameter(finds=26, finde=33) !D1 [0.8667 - 1.1]
!      parameter(finds=54, finde=63) !D2 [1.8 - 2.1]
!      parameter(finds=81, finde=360) !D3 [2.7 - 12]
      parameter(finds=105, finde=135) !D [3.5 - 4.5]

! Spectral energy 
      real sef_cw1,sef_cw2,sef_cw3,sef_cw4
      real sef_ccw1,sef_ccw2,sef_ccw3,sef_ccw4
      dimension sef_cw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_cw2(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_cw3(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_cw4(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_ccw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_ccw2(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_ccw3(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)
      dimension sef_ccw4(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG)


!--------filter parameters-------------------

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/


! rot_spec information-----------------------------------------------
      double precision,allocatable :: t1(:),f1(:),p1(:) ! maxobs
      double complex,allocatable :: u1(:),v1(:)
      double precision,allocatable :: cw(:),ccw(:),puv(:),quv(:)
      parameter(pi=3.141592654d0)

! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid

!DO_FILT, maxobs, FC1, FC2, ftnm, tidcon; make sure numcon
! fix maxobs and numcon
      read(*,*) maxobs   !integer
      read(*,*) ftnm     !string

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)


      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      write (6,'(a,i5)') 'maxobs: ', maxobs
      call flush(6)

      write (6,'(a,a)') ' ftnm:',ftnm
      call flush(6)


!----------------------------------------------------------------------
! allocate variables
      allocate (blk_u_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs),
     & blk_v_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs),
     & cw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs/2,1:MEIG),
     & ccw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs/2,1:MEIG),
     & blk_cw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG),
     & blk_ccw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG))

      allocate (t1(maxobs),f1(maxobs),p1(maxobs),u1(maxobs),v1(maxobs),
     & cw(maxobs),ccw(maxobs),puv(maxobs),quv(maxobs)) 

      ! allocate variables
      allocate(
     & modAU(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & modAV(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & var(1:maxobs))


!! loop over blocks ----------------------------------------------------
      do jblk=jblks,jblke
        do iblk=iblks,iblke
!! loop over blocks ----------------------------------------------------
          write(6,*) 'blk_',jblk,iblk
          call flush(6)

          write(6,*) 'opening files for reading and writing ...'
          call flush(6)


!!--------------------------------------------------------------------------
!! OPEN and load time varying 3D files          
!!--------------------------------------------------------------------------
          infil12='rotary_spectra/'
     &  //'rotspec_cw_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'reading: ',infil12
        open(unit=12,file=infil12,status='old',form='unformatted')


          infil13='rotary_spectra/'
     &  //'rotspec_ccw_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'reading: ',infil13
        open(unit=13,file=infil13,status='old',form='unformatted')


!!--------------------------------------------------------------------------
! load rotary spectra coefficients 
!!--------------------------------------------------------------------------
      cw1(:,:,:,:) = 0.0
      ccw1(:,:,:,:) = 0.0
      do m=1,MEIG
          read(12) cw1(:,:,:,m)
          read(13) ccw1(:,:,:,m)
      enddo ! maxobs

      close(12)
      close(13)

      write(6,*) 'finished reading rotary spectra coefficients'
      call flush(6)


!!--------------------------------------------------------------------------
! Compute spectral energy  
!!--------------------------------------------------------------------------
        sef_cw1(:,:,:) = 0.0
        sef_ccw1(:,:,:) = 0.0

! Loop over grid points in 1 tile (150x200)
        do i=1-mbdy,blki+mbdy
          do j=1-mbdy,blkj+mbdy
       
            do m=1,MEIG
!               sef_cw1(i,j,m) = sum(cw1(i,j,find1s:find1e,m)) 

               do l=finds,finde
                 sef_cw1(i,j,m) = sef_cw1(i,j,m) + cw1(i,j,l,m)
                 sef_ccw1(i,j,m) = sef_ccw1(i,j,m) + ccw1(i,j,l,m)
               enddo

             enddo

          enddo !blkj
        enddo !blki 


!!--------------------------------------------------------------------------
!! OPEN 2D output files         
!!--------------------------------------------------------------------------
        outfile11='ke_rotspec/'//'D4_cw_'//runid//
     &'_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        outfile12='ke_rotspec/'//'D4_ccw_'//runid//
     &'_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=111,file=outfile11,status='replace',
     &      form='unformatted',access='sequential')
        open(unit=112,file=outfile12,status='replace',
     &      form='unformatted',access='sequential')

!-----------write output--------------------------------------------------
        do m=1,MEIG
            write(111) sef_cw1(:,:,m)
            write(112) sef_ccw1(:,:,m)
        enddo
        
!! loop over blocks ----------------------------------------------------
        enddo ! iblk
      enddo ! jblk

      write(6,*) 'Finished...'

      stop

      END ! PROGRAM ROTARY_SPECTRA
