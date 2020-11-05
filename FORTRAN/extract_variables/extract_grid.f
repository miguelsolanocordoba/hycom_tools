      PROGRAM READ_ARCHIVE
       
      IMPLICIT NONE
      
c---USED TO TEST THE EXTRACTION OF A SINGLE FIELD: THKNSS---------      
c extract_grid_c004_v1.f, USM, MCB, 2018-05-08
c   extract all grid parameters only
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
      integer jblks,jblke 
      integer iblks,iblke
      real flag
      real salin,c1,c2,c3,c4,c5,c6,c7,c8,c9,rc6
      real maxvar 
c glbc008
c      parameter(idm=4500,jdm=3528,kdm=41) 
c      parameter (blki=150,blkj=196,halo=3)
c glbc004
c      parameter(idm=9000,jdm=7055,kdm=41) !180 x 50 = 9000 and 200 x 35 = 7000
c      parameter (blki=150,blkj=200,halo=3)
c atlc002
c      parameter(idm=6708,jdm=7372,kdm=32) !
      parameter(idm=6709,jdm=7373,kdm=32) !
      parameter (blki=129,blkj=194,halo=3) ! modes

      parameter(flag=2.0**100)

      character flnm*240,outfil*240,flnm2*240,flnm3*240,flnm4*240
      character flnm5*240
      character blank*5,runid*3
      character*2 numstr(60)
      integer i,j,k,n2drec,nrecl,numfls,iblk,jblk,n,ic,flip_mcb
      
      integer, allocatable, dimension(:,:) :: coast
      real, allocatable, dimension(:) :: w
      real, allocatable, dimension(:,:) :: depth,drgten,cbdrag,
     &                                     plon,plat,qlon,qlat,
     &                                     ulon,ulat,vlon,vlat,
     &                                     pang,pscx,pscy,
     &                                     qscx,qscy,uscx,uscy,
     &                                     vscx,vscy,cori,pasp,
     &                                     srfhgt,ubaro,vbaro,steric
      real, allocatable, dimension(:,:,:) :: u,v,thknss
      real, allocatable, dimension(:,:,:) :: temp,sig,siga
c      real srfhgt(idm,jdm),ubaro(idm,jdm),vbaro(idm,jdm)
c      real u(idm,jdm,kdm),v(idm,jdm,kdm),thknss(idm,jdm,kdm)

      data blank/'     '/
      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60'/


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
c       n2drec = ((idm*jdm+4095)/4096)*4096
       n2drec = 49467392
       write(*,*) n2drec
c
c      initialize I/O buffer
c
      allocate( w(n2drec) )

      allocate( depth(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( coast(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( plon(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( plat(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( pscx(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( pscy(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( uscx(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( uscy(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( vscx(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( vscy(1-halo:idm+halo,1-halo:jdm+halo) )
      allocate( pang(1-halo:idm+halo,1-halo:jdm+halo) )

c
      inquire(iolength=nrecl) w
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

c read regional.depth.a
	read (*,'(a)') flnm2
        write(*,*) flnm2
	write(6,'(a)') blank 
	call flush(6)
	write (6,'(2a)') ' depth file: ',trim(flnm2)

	open(unit=10, file=flnm2,form='unformatted', status='OLD',
     &     access='direct', recl=nrecl, action='READ')


c read regional.grid.a
	read (*,'(a)') flnm3
        write(*,*) flnm3
	write(6,'(a)') blank 
	call flush(6)
	write (6,'(2a)') ' grid file: ',trim(flnm3)

	open(unit=11, file=flnm3,form='unformatted', status='OLD',
     &     access='direct', recl=nrecl, action='READ')


c map griddata variables -----------------------------------------------

c depth
	read(10,rec=1) w   
	do j=1,jdm
	  do i=1,idm
	    depth(i,j)=w(i+(j-1)*idm)

c            write(*,*) i,j,depth(i,j)

C coast is land mask
            if     (depth(i,j).ne.flag) then
              coast(i,j) = 1  ! sea
            else
              coast(i,j) = 0  ! land
c              write(*,*) i,j,coast(i,j),depth(i,j)
            endif

	  enddo
	enddo
	close(10)

c pang
	read(11,rec=9) w   
	do j=1,jdm
	  do i=1,idm
	    pang(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo
c plon
	read(11,rec=1) w   
	do j=1,jdm
	  do i=1,idm
	    plon(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c plat
	read(11,rec=2) w   
	do j=1,jdm
	  do i=1,idm
	    plat(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c pscx
	read(11,rec=10) w   
	do j=1,jdm
	  do i=1,idm
	    pscx(i,j)=w(i+(j-1)*idm)
c              write(*,*) i,j,pang(i,j),plon(i,j),plat(i,j),pscx(i,j)
	  enddo
	enddo

c pscy
	read(11,rec=11) w   
	do j=1,jdm
	  do i=1,idm
	    pscy(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c uscx
	read(11,rec=14) w   
	do j=1,jdm
	  do i=1,idm
	    uscx(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c uscy
	read(11,rec=15) w   
	do j=1,jdm
	  do i=1,idm
	    uscy(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c vscx
	read(11,rec=16) w   
	do j=1,jdm
	  do i=1,idm
	    vscx(i,j)=w(i+(j-1)*idm)
	  enddo
	enddo

c vscy
	read(11,rec=17) w   
	do j=1,jdm
	  do i=1,idm
	    vscy(i,j)=w(i+(j-1)*idm)
c              write(*,*) i,j,uscx(i,j),uscy(i,j),vscx(i,j),vscy(i,j)
c     &   ,pscy(i,j)
	  enddo
	enddo

        close(11)

c        stop


	do j=1,halo
          do i=1,idm

	    ic=idm+1-i
c take a value that is "reflected at the northern boundary"
c why use ic, ic is reverse of idm index??

         	      coast(i,jdm+j)=coast(ic,jdm-j-1)
         	      depth(i,jdm+j)=depth(ic,jdm-j-1)
         	      uscx(i,jdm+j)=uscx(ic,jdm-j-1)
         	      vscx(i,jdm+j)=vscx(ic,jdm-j-1)
         	      pscx(i,jdm+j)=pscx(ic,jdm-j-1)
         	      uscy(i,jdm+j)=uscy(ic,jdm-j-1)
         	      vscy(i,jdm+j)=vscy(ic,jdm-j-1)
         	      pscy(i,jdm+j)=pscy(ic,jdm-j-1)
         	      plon(i,jdm+j)=plon(ic,jdm-j-1)
         	      plat(i,jdm+j)=plat(ic,jdm-j-1)
         	      pang(i,jdm+j)=pang(ic,jdm-j-1)


	  enddo !j

        enddo !i

 	      coast(-2:0,:)=coast(idm-halo+1:idm,:)
 	      depth(-2:0,:)=depth(idm-halo+1:idm,:)
 	      uscx(-2:0,:)=uscx(idm-halo+1:idm,:)
 	      vscx(-2:0,:)=vscx(idm-halo+1:idm,:)
 	      pscx(-2:0,:)=pscx(idm-halo+1:idm,:)
 	      uscy(-2:0,:)=uscy(idm-halo+1:idm,:)
 	      vscy(-2:0,:)=vscy(idm-halo+1:idm,:)
 	      pscy(-2:0,:)=pscy(idm-halo+1:idm,:)
 	      plon(-2:0,:)=plon(idm-halo+1:idm,:)
 	      plat(-2:0,:)=plat(idm-halo+1:idm,:)
 	      pang(-2:0,:)=pang(idm-halo+1:idm,:)

 	      coast(idm+1:idm+halo,:)=coast(1:3,:)
 	      depth(idm+1:idm+halo,:)=depth(1:3,:)
 	      uscx(idm+1:idm+halo,:)=uscx(1:3,:)
 	      vscx(idm+1:idm+halo,:)=vscx(1:3,:)
 	      pscx(idm+1:idm+halo,:)=pscx(1:3,:)
 	      uscy(idm+1:idm+halo,:)=uscy(1:3,:)
 	      vscy(idm+1:idm+halo,:)=vscy(1:3,:)
 	      pscy(idm+1:idm+halo,:)=pscy(1:3,:)
 	      plon(idm+1:idm+halo,:)=plon(1:3,:)
 	      plat(idm+1:idm+halo,:)=plat(1:3,:)
 	      pang(idm+1:idm+halo,:)=pang(1:3,:)

c what are the consequences of setting these values to zero .....
 	      coast(:,1:3)=0
 	      depth(:,1:3)=0.0
 	      uscx(:,1:3)=0.0
 	      vscx(:,1:3)=0.0
 	      pscx(:,1:3)=0.0
 	      uscy(:,1:3)=0.0
 	      vscy(:,1:3)=0.0
 	      pscy(:,1:3)=0.0
 	      plon(:,1:3)=0.0
 	      plat(:,1:3)=0.0
 	      pang(:,1:3)=0.0

        do jblk=jblks,jblke !1,17
	
         do iblk=iblks,iblke !1,5

c save 2D grid variables -------------------------------------------------------------	
c coast pang depth uscx vscx pscx uscy vscy pscy plon plat 


c coast_blk_
	 outfil='griddata/coast_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) coast(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)	

c depth_blk_
	 outfil='griddata/depth_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) depth(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)	

c pang_blk_
	 outfil='griddata/pang_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) pang(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)	

c uscx_blk_
	 outfil='griddata/uscx_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) uscx(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)	

c vscx_blk_
	 outfil='griddata/vscx_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) vscx(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)

c pscx_blk_
	 outfil='griddata/pscx_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) pscx(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)

c uscy_blk_
	 outfil='griddata/uscy_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) uscy(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)	

c vscy_blk_
	 outfil='griddata/vscy_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) vscy(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)

c pscy_blk_
	 outfil='griddata/pscy_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted',
     & 		      access='append')
	   write(99) pscy(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 
           call flush(6)

c plon_blk_
	 outfil='griddata/plon_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted')
	   write(99) plon(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 

c           write(6,*) plon(1+blki*(iblk-1)-halo:1+blki*(iblk-1)-halo+3,
c     & 			  1+blkj*(jblk-1)-halo:1+blkj*(jblk-1)-halo+3)

           call flush(6)

c plat_blk_
	 outfil='griddata/plat_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &'.BinF'
	   open(unit=99,file=outfil,status='replace',form='unformatted')
	   write(99) plat(1+blki*(iblk-1)-halo:blki*iblk+halo,
     & 			  1+blkj*(jblk-1)-halo:blkj*jblk+halo)
	   close(99)
           write(6,'(2a)') 'saved to: ',  outfil
           write(6,'(a)') blank 

c           write(6,*) plat(1+blki*(iblk-1)-halo:1+blki*(iblk-1)-halo+3,
c     & 			  1+blkj*(jblk-1)-halo:1+blkj*(jblk-1)-halo+3)

           call flush(6)


        enddo

        enddo
	
      
      stop
      end
