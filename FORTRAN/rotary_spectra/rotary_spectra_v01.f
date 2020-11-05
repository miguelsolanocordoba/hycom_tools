      PROGRAM ROTARY_SPECTRA
! rotary_spectra_v01.f; Keshav J Raja, USM, 2020/08/25
!   plot rotary spectra from u and v velocities
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
      integer RILE                  ! number of time steps to omit to reduce effect of ringing
      integer DO_FILT               ! switch for filtering (if 0 mean is still removed)
      integer CB_VAR                ! if 1, use spatially varying drag 
      integer jblks,jblke !1,17
      integer iblks,iblke !1,10
      character runid*3

      parameter(TIVAR=0)              ! (=1) then time varying fields are stored
      parameter(CB_VAR=1)              ! in case of 0, use constant drag  
      parameter(RILE=24)              ! exclude #time steps affected by ringing 

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
! filter parameters
!      parameter(NS=4)             !old; 
      parameter(NS=2)             !new; default number of sections (2*NS=order of filter)
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
      character*2 numstr(90)
      character*2 constr(8)

      integer i,j,k,l,iblk,jblk,joe,joe2,kbot
      integer*8 plan

! 2D variables
      integer blk_coast
      dimension blk_coast(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

      real blk_depth,blk_ubaro,blk_vbaro,blk_uscx,blk_vscy,blk_plat
      dimension blk_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_vbaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_ubaro(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_uscx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_vscy(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_plat(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
! ---2d+t variables-----
      real, allocatable :: blk_u_surf(:,:,:),blk_v_surf(:,:,:)

      real zz,blk_mean_u,blk_mean_v
      real,allocatable :: blk_cw(:,:,:),blk_ccw(:,:,:),
     &            blk_rot_spec(:,:,:)
! 3D variables
      real blk_u,blk_v
      dimension blk_u(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)
      dimension blk_v(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm)


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
      read(*,*) DO_FILT  !integer
      read(*,*) maxobs   !integer
!      read(*,*) FC11     !real
!      read(*,*) FC22     !real
      read(*,*) ftnm     !string

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)


      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      write (6,'(a,1i4,i10)') ' DO_FILT maxobs: ',
     & DO_FILT, maxobs
      call flush(6)

!     write (6,'(a,2F3.1)') ' FC11 FC22: ',
!     & FC11, FC22
!      call flush(6)

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

!----------------------------------------------------------------------
! allocate variables
      allocate (blk_u_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_v_surf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_cw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_ccw(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,maxobs),
     & blk_rot_spec(1-mbdy:blkj+mbdy,maxobs,2))

      allocate (t1(maxobs),f1(maxobs),p1(maxobs),u1(maxobs),v1(maxobs),
     & cw(maxobs),ccw(maxobs),puv(maxobs),quv(maxobs)) 

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

! lat and lon
         infil10='griddata/plat_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_plat

          write(6,*) 'max lat is ',maxval(maxval(blk_plat(:,:),1))
          close(10)

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

!!--------------------------------------------------------------------------
!! OPEN 3D time files         
!!--------------------------------------------------------------------------

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

!        infil16='ubaro/'
!     &  //'ubaro_'//runid//
!     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=16,file=infil16,status='old',form='unformatted')
!
!        infil17='vbaro/'
!     &  //'vbaro_'//runid//
!     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
!
!        open(unit=17,file=infil17,status='old',form='unformatted')

!!--------------------------------------------------------------------------
!! OPEN 2D output files         
!!--------------------------------------------------------------------------

        outfile10='out_rotspec/'
     &  //'rot_spec_'//runid//
     &'_'//trim(ftnm)//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=110,file=outfile10,status='replace',
     &      form='unformatted',access='sequential')


!

!----------------------------------------------------------------------------
          blk_u_surf(:,:,:)=0.0
          blk_v_surf(:,:,:)=0.0

          do l=1,maxobs

!            read(16) blk_ubaro
!            read(17) blk_vbaro
!
!            blk_ubarof(:,:,l)=blk_ubaro(:,:)
!            blk_vbarof(:,:,l)=blk_vbaro(:,:)

            do k=1,kdm
              read(14) blk_u(:,:,k)
              read(15) blk_v(:,:,k)
            enddo

            blk_u_surf(:,:,l)=blk_u(:,:,1)
            blk_v_surf(:,:,l)=blk_v(:,:,1)

          enddo !----first maxobs

         close(14)
         close(15)
!         close(16)
!         close(17)

!         blk_u_surfu = blk_u_filt + blk_ubarof
!         blk_v_surfu = blk_v_filt + blk_vbarof

!-----------set flag values to zero else filtering code doesn't work----            
         do l=1,maxobs
           do j=1-mbdy,blkj+mbdy
             do i=1-mbdy,blki+mbdy
        
               if (l.eq.1.AND.blk_depth(i,j).EQ.flag) then
                 blk_depth(i,j)=0.0
!                 write(6,*) i,j
               endif

               if (blk_u_surf(i,j,l).eq.flag) blk_u_surf(i,j,l)=0.0
               if (blk_v_surf(i,j,l).eq.flag) blk_v_surf(i,j,l)=0.0


             enddo   !i
           enddo     !j 

           t1(l)=dble(l) 
         end do!-----second maxobs!     


!---------Calling rot_spec_fftw------------------------------------------

          blk_cw(:,:,:)=0.0
          blk_ccw(:,:,:)=0.0
            
           do j=1-mbdy,blkj+mbdy !j->latitude
             zz=0.0
             do i=1-mbdy,blki+mbdy !i->longitude
         blk_mean_u = sum(blk_u_surf(i,j,:))/  !mean of time series 
     &                 (max(1,size(blk_u_surf(i,j,:))))
         blk_u_surf(i,j,:)= blk_u_surf(i,j,:)-blk_mean_u  !subtract mean       
         blk_mean_v = sum(blk_v_surf(i,j,:))/  !mean of time series 
     &                 (max(1,size(blk_v_surf(i,j,:))))
         blk_v_surf(i,j,:)= blk_v_surf(i,j,:)-blk_mean_v  !subtract mean       


!                if (blk_coast(i,j).eq.1) then
                u1 = dble(blk_u_surf(i,j,:))
                v1 = dble(blk_v_surf(i,j,:))

                call ROT_SPEC_FFTW(u1,v1,t1,maxobs,cw,ccw,f1,p1)
         if (i.eq.50.AND.j.eq.50) then
           write(*,*)cw
         endif
                 
!                blk_cw(j,:) = blk_cw(j,:) + cw
!                blk_ccw(j,:) = blk_ccw(j,:) + ccw 
                blk_cw(i,j,:) = cw
                blk_ccw(i,j,:) = ccw
                zz = zz+1
!                endif ! blk_coast
             enddo   !i
           enddo     !j 
          blk_rot_spec(:,:,1) = sum(blk_cw,1)/zz
          blk_rot_spec(:,:,2) = sum(blk_ccw,1)/zz


!-----------write output--------------------------------------------------
      
       write(110) blk_rot_spec

      close(110)

!! loop over blocks ----------------------------------------------------
      enddo ! iblk
      enddo ! jblk
!! loop over blocks ----------------------------------------------------

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

!       call dfftw_plan_dft_1d(plan,n,psi,psif,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan,n,u1,u1f,1,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, u1, u1f)
        call dfftw_destroy_plan(plan)

!       call dfftw_plan_dft_1d(plan,n,psi,psif,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_plan_dft_1d(plan,n,v1,v1f,1,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, v1, v1f)
        call dfftw_destroy_plan(plan)

!       pu = fu.*conj(fu); pv = fv.*conj(fv);
        pu = dreal(u1f)**2 + dimag(u1f)**2
        pv = dreal(v1f)**2 + dimag(v1f)**2

!        puv =  real(fu).*real(fv) + imag(fu).*imag(fv);
!        quv = -real(fu).*imag(fv) + real(fv).*imag(fu);
        puv = dreal(u1f)*dreal(v1f) + dimag(u1f)*dimag(v1f)
        quv =-dreal(u1f)*dimag(v1f) + dreal(v1f)*dimag(u1f)

!        cw   = (pu+pv-2*quv)./8;
!        ccw  = (pu+pv+2*quv)./8;
        cw = (pu+pv-2d0*quv)/8.d0
        ccw= (pu+pv+2d0*quv)/8.d0
        
        del = t1(2) - t1(1)  
!       freq  = 1/dt*(0:(n/2))/n;
        do k = 1,n/2+1
           f1(k) = dble(k-1)/dble(n)/del
        enddo
        do k = n/2+2,n
           f1(k) = dble(k-1-n)/dble(n)/del
        enddo

        p1 = 1/f1


      END ! SUBROUTINE ROT_SPEC_FFTW
