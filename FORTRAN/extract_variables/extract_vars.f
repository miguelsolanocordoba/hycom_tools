      PROGRAM READ_ARCHIVE
       
      IMPLICIT NONE
      
c---USED TO TEST THE EXTRACTION OF A SINGLE FIELD: THKNSS---------      
c extract_vars_TS_c004_v12.f, MCB, 2019-1-8
c   extracts T and S and saves it to a file for 1/25
c   COMMENT OUT also extract mixed layer depth and normalizes it by 9806.0
c extract_vars_c004_v12b.f, MCB, USM, 2018-1-27
c   for 1/25 and made able to run on shepard
c   allowing for fewer variables to be stored in memory
c   only these 3D global variables: u,v,thknss 
c extract_vars_c004_v12.f, MCB, USM, 2017-11-14
c   for 1/25 grid
c extract_vars_c008_v12.f, USM, 2017-05-12
c   added DO3D to exclude 3D variables extraction
c extract_vars_c008_v11.f, USM, 2016-08-02
c   prepared for gordon
c extract_vars_c008_v4.f, USM, 2016-06-23
c   also extract density f(t)
c extract_vars_c008_v3.f, USM, 2015-05-08
c extract_vars_c008_v2.f, nrl, 2015-02-05
c extract_vars_c008_v1.f, nrl, 2014-06-04
c extract_vars_41L_v1.f, nrl, 2014-04-29
c extract_vars_v6.f mcb, nrl, 2013-11-15
c extract_vars_v5.f mcb, nrl, 2013-10-30
c v2 now 30 tiles along x
c v1 new grid
c v1: redo for 41 layers sims, added bottom drag, cb, extraction 
c v6: extract linear drag fields r*h at p
c global extraction
c loads a global *.a file, adds halo, and saves chunks
       
      integer idm,jdm,kdm,blki,blkj,halocd,halo
      integer jblks,jblke !1,35
      integer iblks,iblke !1,60
      integer DO3D
      real flag,smxl
      real salin,c1,c2,c3,c4,c5,c6,c7,c8,c9,rc6
      real maxvar 

c      parameter(DO3D=1) ! if 0 than 3D fields are NOT stored
      parameter(DO3D=0) ! if 0 than 3D fields are NOT stored

c glbc004
      parameter(idm=9000,jdm=7055,kdm=41) !150 x 60 = 9000 and 200 x 35 = 7000
      parameter (blki=150,blkj=200,halo=3) ! modes

      parameter(flag=2.0**100)
      parameter(smxl=9806.0)

      character flnm*240,outfil*240,flnm2*240,flnm3*240,flnm4*240
      character outfil1*240,outfil2*240
      character flnm5*240
      character blank*5,runid*3
      character*2 numstr(90)
      integer i,j,k,n2drec,nrecl,numfls,iblk,jblk,n,ic,flip_mcb
      
c      integer, allocatable, dimension(:,:) :: coast
      real, allocatable, dimension(:) :: w
      real, allocatable, dimension(:) :: w1
      real, allocatable, dimension(:) :: w2
      real, allocatable, dimension(:,:) :: surflx
      real, allocatable, dimension(:,:) :: salflx
c      real, allocatable, dimension(:,:) :: ubaro
c      real, allocatable, dimension(:,:) :: depth,drgten,cbdrag,
c     &                                     plon,plat,qlon,qlat,
c     &                                     ulon,ulat,vlon,vlat,
c     &                                     pang,pscx,pscy,
c     &                                     qscx,qscy,uscx,uscy,
c     &                                     vscx,vscy,cori,pasp,
c     &                                     srfhgt,ubaro,vbaro,steric
      real, allocatable, dimension(:,:,:) :: u,v,thknss
c      real, allocatable, dimension(:,:,:) :: temp,sig,siga
c      real srfhgt(idm,jdm),ubaro(idm,jdm),vbaro(idm,jdm)
c      real u(idm,jdm,kdm),v(idm,jdm,kdm),thknss(idm,jdm,kdm)

      data blank/'     '/
      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/


c == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE ==

c old eq of state
c      parameter(c1= 9.903308E+00,c2=-1.618075E-02,c3= 7.819166E-01,
c     & c4=-6.593939E-03,c5=-2.896464E-03,c6= 3.038697E-05,rc6=1.0/c6,
c     & c7= 3.266933E-05,c8= 1.180109E-04,c9= 3.399511E-06)

      integer, parameter ::
     &  sigver=6  !17-term sigma-2
csig0&  sigver=5  !17-term sigma-0
c
c      real    sig,dsigdt,dsigds,tofsig,sofsig,kappaf,
c     &        sigloc,dsiglocdt,dsiglocds
c
      real    sig_n,sig_d,sig_q, dsigdt_n,dsigdt_d, dsigds_n,dsigds_d
      real    kappaf1
      real    sigloc_n,sigloc_d,sigloc_q,
     &        dsiglocdt_n,dsiglocdt_d, dsiglocds_n,dsiglocds_d
c
c      real    r,s,t,pdb,prs
      real    r,pdb,prs
      integer kkf
c
      real, parameter ::
     &   aone =1.0,
     &   ahalf=1.0/2.0,
     &   a3rd =1.0/3.0, athird =a3rd,
     &   a4th =1.0/4.0, afourth=a4th
c
cc      real, parameter ::
cc     &   c1= 1.0e-01,      !not used, required to compile mxkrtm.f
cc     &   c2= 1.0e-02,      !not used, required to compile mxkrtm.f
cc     &   c3= 1.0e-03,      !not used, required to compile mxkrtm.f
cc     &   c4= 1.0e-04,      !not used, required to compile mxkrtm.f
cc     &   c5= 1.0e-05,      !not used, required to compile mxkrtm.f
cc     &   c6= 1.0e-06,      !not used, required to compile mxkrtm.f
cc     &   c7= 1.0e-07       !not used, required to compile mxkrtm.f
c
c --- Jackett, McDougall, Feistel, Wright and Griffies (2006), 
c --- Algorithms for Density, Potential Temperature, Conservative
c --- Temperature, and the Freezing Temperature of Seawater, JAOT
c
c --- coefficients for 25-term rational function sigloc().
      real, parameter ::
     &   c001= 9.9984085444849347d+02,     !num. constant    coefficent
     &   c002= 7.3471625860981584d+00,     !num.    T        coefficent
     &   c003=-5.3211231792841769d-02,     !num.    T^2      coefficent
     &   c004= 3.6492439109814549d-04,     !num.    T^3      coefficent
     &   c005= 2.5880571023991390d+00,     !num.       S     coefficent
     &   c006= 6.7168282786692355d-03,     !num.    T  S     coefficent
     &   c007= 1.9203202055760151d-03,     !num.       S^2   coefficent
     &   c008= 1.0000000000000000d+00,     !den. constant    coefficent
     &   c009= 7.2815210113327091d-03,     !den.    T        coefficent
     &   c010=-4.4787265461983921d-05,     !den.    T^2      coefficent
     &   c011= 3.3851002965802430d-07,     !den.    T^3      coefficent
     &   c012= 1.3651202389758572d-10,     !den.    T^4      coefficent
     &   c013= 1.7632126669040377d-03,     !den.       S     coefficent
     &   c014= 8.8066583251206474d-06,     !den.    T  S     coefficent
     &   c015= 1.8832689434804897d-10,     !den.    T^3S     coefficent
     &   c016= 5.7463776745432097d-06,     !den.    T  S^1.5 coefficent
     &   c017= 1.4716275472242334d-09      !den.    T^3S^1.5 coefficent
      real, parameter ::
     &   c018= 1.1798263740430364d-02,     !num. P           coefficent
     &   c019= 9.8920219266399117d-08,     !num. P  T^2      coefficent
     &   c020= 4.6996642771754730d-06,     !num. P     S     coefficent
     &   c021= 2.5862187075154352d-08,     !num. P^2         coefficent
     &   c022= 3.2921414007960662d-12,     !num. P^2T^2      coefficent
     &   c023= 6.7103246285651894d-06,     !den. P           coefficent
     &   c024= 2.4461698007024582d-17,     !den. P^2T^3      coefficent
     &   c025= 9.1534417604289062d-18      !den. P^3T        coefficent
c --- additional coefficients for dsiglocdt().
      real, parameter ::
     &   c031= 7.3471625860981580d+00,     !num. constant    coefficent
     &   c032=-1.0642246358568354d-01,     !num.    T        coefficent
     &   c033= 1.0947731732944364d-03,     !num.    T^2      coefficent
     &   c034= 6.7168282786692355d-03,     !num.       S     coefficent
     &   c035= 7.2815210113327090d-03,     !den. constant    coefficent
     &   c036=-8.9574530923967840d-05,     !den.    T        coefficent
     &   c037= 1.0155300889740728d-06,     !den.    T^2      coefficent
     &   c038= 5.4604809559034290d-10,     !den.    T^3      coefficent
     &   c039=-8.8066583251206470d-06,     !den.       S     coefficent
     &   c040= 5.6498068304414700d-10,     !den.    T^2S     coefficent
     &   c041= 2.9432550944484670d-09,     !den.    T  S^1.5 coefficent
     &   c042= 1.9784043853279823d-07,     !num. P  T        coefficent
     &   c043= 6.5842828015921320d-12,     !num. P^2T        coefficent
     &   c044= 7.3385094021073750d-17,     !den. P^2T^2      coefficent
     &   c045= 9.1534417604289060d-18      !den. P^3         coefficent
c --- additional coefficients for dsiglocds().
      real, parameter ::
     &   c051= 2.5880571023991390d+00,     !num. constant    coefficent
     &   c052= 6.7168282786692355d-03,     !num.    T        coefficent
     &   c053= 3.8406404111520300d-03,     !num.       S     coefficent
     &   c054= 1.7632126669040377d-03,     !den. constant    coefficent
     &   c055=-8.8066583251206470d-06,     !den.    T        coefficent
     &   c056= 1.8832689434804897d-10,     !den.    T^3      coefficent
     &   c057= 8.6195665118148150d-06,     !den.       S^0.5 coefficent
     &   c058= 2.2074413208363504d-09,     !den.    T^2S^0.5 coefficent
     &   c059= 4.6996642771754730d-06      !num. P           coefficent
c
      real, parameter :: sqrmin=0.d0       !sqrt arg can't be negative
c --- reference pressure.
      real, parameter :: prs2pdb=1.d-4     !Pascals to dbar
csig0 real, parameter :: pref=   0.d0      !ref. pressure in Pascals, sigma0
      real, parameter :: pref=2000.d4      !ref. pressure in Pascals, sigma2
      real, parameter :: rpdb=pref*prs2pdb !ref. pressure in dbar
c --- coefficients for 17-term rational function sig() at rpdb.
      real, parameter ::
     &   c101=c001+(c018-c021*rpdb)*rpdb, !num. constant    coefficent
     &   c103=c003+(c019-c022*rpdb)*rpdb, !num.    T^2      coefficent
     &   c105=c005+c020*rpdb,             !num.       S     coefficent
     &   c108=c008+c023*rpdb,             !den. constant    coefficent
     &   c109=c009-c025*rpdb**3,          !den.    T        coefficent
     &   c111=c011-c024*rpdb**2           !den.    T^3      coefficent
c --- additional coefficients for dsigdt().
      real, parameter ::
     &   c132=c032+(c042-c043*rpdb)*rpdb, !num.    T        coefficent
     &   c135=c035-c045*rpdb**3,          !den. constant    coefficent
     &   c137=c037-c044*rpdb**2           !den.    T^2      coefficent
c --- additional coefficients for dsigds().
      real, parameter ::
     &   c151=c051+c059*rpdb              !num. constant    coefficent
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter ::
     &   rhoref=1.d3  !rhoref=qthref kg/m^3
      real, parameter ::
     &   sclkap=1.e-11
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)

c == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE ==
      
c       READ(5,*) FILENAME

c      n2drec = size of output 2-d array, multiple of 4096
c
       n2drec = ((idm*jdm+4095)/4096)*4096
c
c      initialize I/O buffer
c
      allocate( w1(n2drec) )
      allocate( w2(n2drec) )
      allocate( surflx(1-halo:idm+halo,1-halo:jdm+halo) ) !surface heat flux 
      allocate( salflx(1-halo:idm+halo,1-halo:jdm+halo) ) !surface salt flux 
c      allocate( ubaro(1-halo:idm+halo,1-halo:jdm+halo) ) !mixed layer depth 
c      allocate( u(1-halo:idm+halo,1-halo:jdm+halo,kdm) ) !T 
c      allocate( v(1-halo:idm+halo,1-halo:jdm+halo,kdm) ) !S
c      allocate( thknss(1-halo:idm+halo,1-halo:jdm+halo,kdm) )
c      allocate( temp(1-halo:idm+halo,1-halo:jdm+halo,kdm) )
c      allocate( sig(1-halo:idm+halo,1-halo:jdm+halo,kdm) )
c      allocate( siga(1-halo:idm+halo,1-halo:jdm+halo,kdm) )

c      siga(:,:,:) = 0.0
c
      inquire(iolength=nrecl) w1
      inquire(iolength=nrecl) w2
c
! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write(6,'(a)') blank 
      call flush(6)
      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      read(*,*) numfls
      numfls = numfls - 4 !exclude grid, depth, and 2 drag files

      write(6,'(a)') blank 
      call flush(6)
      write (6,'(a,i5)') ' number of time steps: ',numfls
      call flush(6)

c read regional.depth.a
	read (*,'(a)') flnm2
	write(6,'(a)') blank 
	call flush(6)
	write (6,'(2a)') ' depth file: ',trim(flnm2)

c read regional.grid.a
	read (*,'(a)') flnm3
	write(6,'(a)') blank 
	call flush(6)
	write (6,'(2a)') ' grid file: ',trim(flnm3)

c read linear drag
	read (*,'(a)') flnm4
	write(6,'(a)') blank 
	call flush(6)
	write (6,'(2a)') ' drag file: ',trim(flnm4)

c read quadratic drag
        read (*,'(a)') flnm5
        write(6,'(a)') blank
        call flush(6)
        write (6,'(2a)') ' drag file: ',trim(flnm5)


      do n=1,numfls !time
c      do n=1,1

        read (*,'(a)') flnm
	write(6,'(a)') blank 
        call flush(6)
        write (6,'(i5,a,a)') n,') input file: ',trim(flnm)
        call flush(6)
	write(6,'(a)') blank
	call flush(6)

	open(unit=9, file=flnm,form='unformatted', status='OLD',
     &     access='direct', recl=nrecl, action='READ')



c layered variables --------------------------------------    
c montg1 1
c srfhgt 2
c steric 3
c surflx 4
c salflx 5
c bl_dpth 6
c mix_dpth 7
c covice 8
c thkice 9
c temice 10
c u_btrop 11
c v_btrop 12
c u-vel 13
c v-vel 14
c thknss 15
c temp 16
c salin 17

! -----------------------------------------------------------------
! surface heat/salt flux - 2D var

        read(9,rec=4) w1
        read(9,rec=5) w2

        do j=1,jdm

          do i=1,idm

            surflx(i,j)=w1(i+(j-1)*idm)
            salflx(i,j)=w2(i+(j-1)*idm)
c            ubaro(i,j)=w(i+(j-1)*idm)/smxl
c            write(*,*) i,j,vbaro(i,j) ! mcb
          enddo

        enddo



! -----------------------------------------------------------------
! only do 3D variables if necesarry 

        if(DO3D.EQ.1) then   !DO3D

! read T, every 5th line

        do k=1,kdm ! loop over depth

	  read(9,rec=16+5*(k-1)) w

	  do j=1,jdm

	    do i=1,idm

	      u(i,j,k)=w(i+(j-1)*idm)

	    enddo

	  enddo

! read baroclinic S

	  read(9,rec=17+5*(k-1)) w

	  do j=1,jdm

	    do i=1,idm

	      v(i,j,k)=w(i+(j-1)*idm)

c              if (k.eq.41) write(*,*) i,j,k,v(i,j,k) ! mcb

	    enddo

	  enddo


        enddo ! loop over depth

        endif !DO3D
        close(9)

c        stop

	do j=1,halo
          do i=1,idm

	    ic=idm+1-i
c take a value that is "reflected at the northern boundary"
c why use ic, ic is reverse of idm index??
c	      srfhgt(i,jdm+j)   =srfhgt(ic,jdm-j-1)
c	      steric(i,jdm+j)   =steric(ic,jdm-j-1)
c	      ubaro(i,jdm+j)    =ubaro(i,jdm-j-1)
	      surflx(i,jdm+j)    =surflx(i,jdm-j-1)
	      salflx(i,jdm+j)    =salflx(i,jdm-j-1)
c	      vbaro(i,jdm+j)    =vbaro(i,jdm-j-1)
c	      u(i,jdm+j,:)      =u(i,jdm-j-1,:)
c	      v(i,jdm+j,:)      =v(i,jdm-j-1,:)
c             thknss(i,jdm+j,:) =thknss(i,jdm-j-1,:)
c	      sig(i,jdm+j,:)    =sig(i,jdm-j-1,:)
c	      siga(i,jdm+j,:)   =siga(i,jdm-j-1,:)

	  enddo !j

        enddo !i

c do the x-s on both sides
c	srfhgt(-2:0,:)   =srfhgt(idm-halo+1:idm,:)
c	steric(-2:0,:)   =steric(idm-halo+1:idm,:)
c	ubaro(-2:0,:)    =ubaro(idm-halo+1:idm,:)
	surflx(-2:0,:)    =surflx(idm-halo+1:idm,:)
	salflx(-2:0,:)    =salflx(idm-halo+1:idm,:)
c	vbaro(-2:0,:)    =vbaro(idm-halo+1:idm,:)
c	u(-2:0,:,:)      =u(idm-halo+1:idm,:,:)
c	v(-2:0,:,:)      =v(idm-halo+1:idm,:,:)
c	thknss(-2:0,:,:) =thknss(idm-halo+1:idm,:,:)
c	sig(-2:0,:,:)    =sig(idm-halo+1:idm,:,:)
c	siga(-2:0,:,:)   =siga(idm-halo+1:idm,:,:)

c	srfhgt(idm+1:idm+halo,:)   =srfhgt(1:3,:)
c	steric(idm+1:idm+halo,:)   =steric(1:3,:)
c	ubaro(idm+1:idm+halo,:)    =ubaro(1:3,:)
	surflx(idm+1:idm+halo,:)    =surflx(1:3,:)
	salflx(idm+1:idm+halo,:)    =salflx(1:3,:)
c	vbaro(idm+1:idm+halo,:)    =vbaro(1:3,:)
c	u(idm+1:idm+halo,:,:)      =u(1:3,:,:)
c	v(idm+1:idm+halo,:,:)      =v(1:3,:,:)
c	thknss(idm+1:idm+halo,:,:) =thknss(1:3,:,:)
c	sig(idm+1:idm+halo,:,:)    =sig(1:3,:,:)
c	siga(idm+1:idm+halo,:,:)   =siga(1:3,:,:)

c padd the south boundary with zeros
c	srfhgt(:,1:3)=0.0
c	steric(:,1:3)=0.0
c	ubaro(:,1:3)=0.0
	surflx(:,1:3)=0.0
	salflx(:,1:3)=0.0
c	vbaro(:,1:3)=0.0
c	u(:,1:3,:)=0.0
c	v(:,1:3,:)=0.0
c	thknss(:,1:3,:)=0.0
c	sig(:,1:3,:)=0.0
c	siga(:,1:3,:)=0.0

        write(*,*) '======== finished reading global data =========='


        do jblk=jblks,jblke !1,17
	
         do iblk=iblks,iblke !1,5

c 2D mixed-layer depth -------------------------------------------------------------
c         outfil='mix_dpth/mix_dpth_'//runid//
c     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
c     &'.BinF'
c
c         open(unit=99,file=outfil,status='unknown',form='unformatted',
c     &                access='append')
c
c         write(99) ubaro(1+blki*(iblk-1)-halo:blki*iblk+halo,
c     &                    1+blkj*(jblk-1)-halo:blkj*jblk+halo)
c
c         close(99)
c
c         write(6,'(2a)') 'saved to: ',  outfil
c         write(6,'(a)') blank
c         call flush(6)  

c 2D surface fluxes -------------------------------------------------------------
         outfil1='surface/surflx_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
         open(unit=98,file=outfil1,status='unknown',form='unformatted',
     &                access='append')

         outfil2='surface/salflx_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
         open(unit=99,file=outfil2,status='unknown',form='unformatted',
     &                access='append')


         write(98) surflx(1+blki*(iblk-1)-halo:blki*iblk+halo,
     &                    1+blkj*(jblk-1)-halo:blkj*jblk+halo)

         write(99) salflx(1+blki*(iblk-1)-halo:blki*iblk+halo,
     &                    1+blkj*(jblk-1)-halo:blkj*jblk+halo)

         close(98)
         close(99)

         write(6,'(2a)') 'saved to: ',  outfil1
         write(6,'(2a)') 'saved to: ',  outfil2
         write(6,'(a)') blank
         call flush(6)  

! -----------------------------------------------------------------
! only do 3D variables if necesarry 

         if (DO3D.EQ.1) then !DO3D

c u velocity -------------------------------------------------------------
	   outfil='temp/T_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'

	   open(unit=99,file=outfil,status='unknown',form='unformatted',
     & 		      access='append')

	 do k=1,kdm
	       
	   write(99) u(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)

	 enddo
        
         close(99)

	 write(6,'(2a)') 'saved to: ',  outfil
         write(6,'(a)') blank 
         call flush(6)	


c v velocity -------------------------------------------------------------
	  outfil='sal/S_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'

	   open(unit=99,file=outfil,status='unknown',form='unformatted',
     & 		      access='append')

	   do k=1,kdm
	  
	   write(99) v(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			 1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)

	   enddo
	  
	   close(99)

	 write(6,'(2a)') 'saved to: ',  outfil
         write(6,'(a)') blank 
         call flush(6)	


        endif !DO3D


        enddo

        enddo
	
      enddo !time
      
      stop
      end
