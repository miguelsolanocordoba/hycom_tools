      PROGRAM ROTARY_SPECTRA_UVMODES
! This program computes the rotary spectra of modal amplitudes u/v. 
!
! Rotary spectra is computed as (see subroutine ROT_SPEC_FFTW): 
!
!   u1f = fft(u)
!   v1f = fft(v) 
!
!   pu = dreal(u1f)**2 + dimag(u1f)**2
!   pv = dreal(v1f)**2 + dimag(v1f)**2
!
!   puv = dreal(u1f)*dreal(v1f) + dimag(u1f)*dimag(v1f)
!   quv =-dreal(u1f)*dimag(v1f) + dreal(v1f)*dimag(u1f)
!
!   cw = (pu+pv-2d0*quv)/8.d0
!   ccw= (pu+pv+2d0*quv)/8.d0
! 
! rotary_spectra_v01.f; Keshav J Raja, USM, 2020/08/25
! rotary_spectra_uvmodes.f; Miguel S Solano, USM, 2021/02/17

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
      character infil10*240,infil11*240,infil12*240,infil13*240
      character outfile1*60,outfile2*60,outfile3*60
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
      real,allocatable :: blk_rotspec(:,:,:,:)
      real,allocatable :: blk_cw(:,:,:,:),blk_ccw(:,:,:,:)
      real,allocatable :: cw1(:,:,:,:),ccw1(:,:,:,:)

! 3D variables
      real blk_u,blk_v
      dimension blk_u(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_v(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)

! 4D variables
      real blk_Uvec
      dimension blk_Uvec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm,MEIG)

      real, allocatable :: modAU(:,:,:,:)  ! modal amplitudes x,y,mod,time
      real, allocatable :: modAV(:,:,:,:)
      real, allocatable :: var(:)

! rot_spec information-----------------------------------------------
      double precision,allocatable :: t1(:),f1(:),p1(:) ! maxobs
      double complex,allocatable :: u1(:),v1(:)
      double precision,allocatable :: cw(:),ccw(:),puv(:),quv(:)
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

      read(*,*) jblks    ! starting tile in y-dir
      read(*,*) jblke    ! ending tile in y-dir
      read(*,*) iblks    ! starting tile in x-dir
      read(*,*) iblke    ! ending tile in x-dir
      read(*,*) runid    ! Experiment number (e.g. 190)
      read(*,*) maxobs   ! length of time array (integer)
      read(*,*) ftnm     ! file identifier (string)

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      write (6,'(a,i5)') 'maxobs: ', maxobs
      call flush(6)

      write (6,'(a,a)') ' ftnm:',ftnm
      call flush(6)


!---- Dynamic memory allocation
! Note, blk_cw/ccw contain the double-sided spectra and cw1/ccw1 only
! the single sided spectra

      allocate (blk_u_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs),
     & blk_v_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs),
     & cw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs/2,1:MEIG),
     & ccw1(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs/2,1:MEIG),
     & blk_cw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG),
     & blk_ccw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:maxobs,1:MEIG))

      allocate (t1(maxobs),f1(maxobs),p1(maxobs),u1(maxobs),v1(maxobs),
     & cw(maxobs),ccw(maxobs),puv(maxobs),quv(maxobs)) 

      allocate(
     & modAU(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & modAV(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG,1:maxobs),
     & var(1:maxobs))


!-----------------------------------------------!
!---- MAIN LOOP: Loops over tiles iblk and jblk !
!-----------------------------------------------!
      do jblk=jblks,jblke
        do iblk=iblks,iblke

          write(6,*) 'blk_',jblk,iblk
          call flush(6)

          write(6,*) 'opening files for reading and writing ...'
          call flush(6)

! Load coast
          infil10='griddata/coast_'//runid//'_blk_'
     &            //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',form='unformatted')

          read(10) blk_coast
          close(10)

! Open modal amplitude files
          infil12='modeslayramp/modAU_'//runid//'_'
     &            //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'reading file: ',infil12
          open(unit=12,file=infil12,status='old',form='unformatted')

          infil13='modeslayramp/modAV_'//runid//'_'
     &            //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'reading file: ',infil13
          open(unit=13,file=infil13,status='old',form='unformatted')


! Load modal amplitudes modAU = f(x,y,mod,time) 
          modAU(:,:,:,:) = 0.0
          modAV(:,:,:,:) = 0.0

          do l=1,maxobs
            read(12) modAU(:,:,:,l)
            read(13) modAV(:,:,:,l)
          enddo ! maxobs

          close(12)
          close(13)

          write(6,*) 'finished reading mode amps'
          call flush(6)

          blk_cw(:,:,:,:)=0.0
!          blk_ccw(:,:,:,:)=0.0

! Start loop over modAU(i,j,m,:) [3 do loops]
! Loop over modes (MEIG) first
          do m=1,MEIG

            do l=1,maxobs ! define time (t1)
             t1(l)=dble(l) 
            enddo!

            
            do j=1-mbdy,blkj+mbdy !j->latitude
              zz=0d0
              do i=1-mbdy,blki+mbdy !i->longitude
                blk_mean_u = sum(modAU(i,j,m,:))/  !mean of time series 
     &                       (max(1,size(modAU(i,j,m,:))))
                modAU(i,j,m,:)= modAU(i,j,m,:)-blk_mean_u  !subtract mean       
                blk_mean_v = sum(modAV(i,j,m,:))/  !mean of time series 
     &                       (max(1,size(modAV(i,j,m,:))))
                modAV(i,j,m,:)= modAV(i,j,m,:)-blk_mean_v  !subtract mean       

                if (blk_coast(i,j).eq.1) then
                  u1 = dble(modAU(i,j,m,:))
                  v1 = dble(modAV(i,j,m,:))
  
                  call ROT_SPEC_FFTW(u1,v1,t1,maxobs,cw,ccw,f1,p1)
                   
                  blk_cw(i,j,:,m) = cw
!                  blk_ccw(i,j,:,m) = ccw

                  zz = zz+1d0
                endif ! blk_coast

              enddo   ! iblk (within tile)
            enddo     ! jblk (within tile) 

          enddo ! modes (MEIG)


!--------- 
! The double sided spectra (CW) contains both CW and CCW
! (CW) has size maxobs. CW is from 2 to 361 and CCW from 361 to 720. 
      do l=2,maxobs/2+1
          ccw1(:,:,l-1,:) = blk_cw(:,:,l,:)
      enddo

! The second half of the double-sided spectra is reversed. Therefore,
! 
      do l=maxobs,maxobs/2+1,-1
          cw1(:,:,l-1,:) = blk_cw(:,:,l,:)
      enddo


!--- OUTPUT
! output files
        outfile1='rotary_spectra/rotspec_cw_'//runid//'_'
     &            //numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        outfile2='rotary_spectra/rotspec_ccw_'//runid//'_'
     &            //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=111,file=outfile1,status='replace',
     &      form='unformatted',access='sequential')
        open(unit=112,file=outfile2,status='replace',
     &      form='unformatted',access='sequential')

! Write output
      do m=1,MEIG 
       write(111) cw1(:,:,:,m)
       write(112) ccw1(:,:,:,m)
      enddo ! MEIG

      close(111)
      close(112)


      enddo ! iblk
      enddo ! jblk

      write(6,*) 'Finished...'

      stop

      END ! PROGRAM ROTARY_SPECTRA


!-------------------------------------------------------------------------------
      SUBROUTINE ROT_SPEC_FFTW(u1,v1,t1,n,cw,ccw,f1,p1)
!-------------------------------------------------------------------------------
        implicit none
        include 'fftw3.f'
        integer n,k
        integer*8 plan
        double precision del,t1(n),f1(n),p1(n)
        double complex u1f(n),v1f(n),u1(n),v1(n)
        double precision pu(n),pv(n),cw(n),ccw(n),puv(n),quv(n)

        call dfftw_plan_dft_1d(plan,n,u1,u1f,1,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, u1, u1f)
        call dfftw_destroy_plan(plan)

        call dfftw_plan_dft_1d(plan,n,v1,v1f,1,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, v1, v1f)
        call dfftw_destroy_plan(plan)

        pu = dreal(u1f)**2 + dimag(u1f)**2
        pv = dreal(v1f)**2 + dimag(v1f)**2

        puv = dreal(u1f)*dreal(v1f) + dimag(u1f)*dimag(v1f)
        quv =-dreal(u1f)*dimag(v1f) + dreal(v1f)*dimag(u1f)

        cw = (pu+pv-2d0*quv)/8.d0
        ccw= (pu+pv+2d0*quv)/8.d0
        
        del = t1(2) - t1(1)  

        do k = 1,n/2+1
           f1(k) = dble(k-1)/dble(n)/del
        enddo
        do k = n/2+2,n
           f1(k) = dble(k-1-n)/dble(n)/del
        enddo

        p1 = 1/f1


      END ! SUBROUTINE ROT_SPEC_FFTW
