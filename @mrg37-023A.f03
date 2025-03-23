!*** Late Version: 2025/3/23 ******************** Fortran 2003/2008 ****
!*                                                                     *
!*   ## Macro-particle Kinetic Simulation for Solar EM Simulation ##   *
!*     << Fully-implicit scheme with kinetic ions and electrons >>     *
!*                                                                     *
!*    Refs.: 1) M.Tanaka, J.Comput.Phys., vol. 79, 206 (1988).         *
!*           2) M.Tanaka, J.Comput.Phys., vol.107, 124 (1993).         *
!*           3) M.Tanaka, Comput.Phys.Comm., vol.87, 117 (1995).       *
!*           4) M.Tanaka, Comput.Phys.Comm., vol.241, 56 (2019).       *
!*           5) H.Shimazu, M.Tanaka, S.Machida, J.Geophys.Res.,        *
!*               101, 27565 (1996).                                    *
!*           6) M.Tanaka, Bulletin of Chubu University (Mar,2022).     *
!*                                                                     *
!*    Simulation files                                                 *
!*    1. @mrg37_080A.f03: simulation code, job serial number '80'      *
!*    2. param_080A.h   : parameter file                               *
!*    3. rec_3d80A      : Simulation time, box size, parameters of     *
!*                  ions and electrons, decentering parameter, etc,    *
!*                                                                     *
!*  * For kinetic ions and electrons, the time step of dt=1.2/wpe      *
!*    may be used. One should read the reference Ref.2 JCP (1993).     *
!*                                                                     *
!*  * Author: Motohiko Tanaka, Ph.D., Chikusa, Nagoya 464, Japan.      *
!*                                                                     *
!*  https://github.com/Mtanaka77/                                      *
!*  Macro-Particle_Simulation_of_Magnetic_Reconnection                 *
!*                                                                     *
!**** First Version: 7/31/1996 ************************* 09/12/2000 ****
!
!    @mrg3-A023.f03: Non-periodic in y direction (two-points) and 
!    periodic (three-points) in x and z directions.
!
!      The full-implicit plasma simulation code was created at 
!    /cfpsol/ subroutines. It was successfully applied in 2D by 1995 
!    and 3D in 1996. The present code is Fortran 2003/2008 rewritten 
!    in MPICH Ver.3 and Ver.4, whose parallel version is completed 
!    by using mpi_sendrecv with PE's mx*myA*mz/npc overlaps.
!
!    Change in 2022
!    1) The mpi routines isend and irecv are used, i.e., call isend(sendbuf),
!       call irecv(recvbuf), /wwstbm/ -> /sendrev1/,/semdrev2/.
!
!    2) The nearest integer: ip= hxi*rxl(l) +0.5 -> 0.001 - 0.999 -> 0
!      The sharing is, xx = hxi*rxl(l) -ip  L.1200. 
!    3) Graphical plots are used in real*4
!
!    4) The MPI routine mpi_allgather is with parallel processors.
!    5) Arrays are:
!     real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: ex,ey,ez,...
!     real(C_DOUBLE),dimension(3,mx,my+1,mz) :: xx in cfpsol routines.
!
!    6) The first term of an array starts at 1, np1(1)=1, nz1(1)=1
!     in the /cfpsol/. The code starts with 0, as np1(1)-1 for particles.
!
!    7) /fulmov/, /cfpsol/, /emcoef/, /wwstbi/j/k/m are in parallel         
!      execution, /poissn/, /emcof3/ are non-parallel.
!   
!    8) The files is created by /fplot3/ and /cplot3/
!       by this program. They are converted by 'pspdf fortr.77(.ps)' 
!       on Linux, and is then output as fortr.77.pdf on Windows 10.
!                                                         May 8, 2022
!
!    9) They become true if iwrt(it,5).eq.0 or mod(it,5).eq.0
! 
!     if(iwrt(it,nplot).eq.0) then
!     if(mod(it,5).eq.0) then
!          it           iwrt         mod
!           1            1            1
!           2            1            2
!           3            1            3
!           4            1            4
!           5            0            0
!
!   10) This implicit particle code is free from the Courant condition.
!      Rather it is bound by not dt*wce >> 0. It is different from
!      explicit codes, which must satisfy Delta_x/Delta_t > 1.
!
!   > Generic function names are used as sqrt <- dsqrt, max <- amax1
!
!-----------------------------------------------------------------------
!*                                                                     *
!*     /main/ ------ /trans/                                           *
!*                        /fulmov/                                     *
!*                              ---> partbc                            *
!*                              ---> srimp1-srimp2                     *
!*                        /cfpsol/,/poissn/,/diag1/                    *
!*                              ---> fplot3, cplot3                    *
!*                  /init/ ------ /loadpt/, /readpt/ initally          *
!*                                                                     *
!*      sequential:                                                    *
!*        /restrt/ ------- arrays (write /read)                        *
!*                                                                     *
!-----------------------------------------------------------------------
!* Fortran 2003 by F77/F90 is made by continuity and lower cases, 
!*  :s%/^C/!/
!*  :s%/^c/!/
!*  tr 'A-Z' 'a-z' <@mrg3.f >@mrg37.f03
!*
!* $ mpif90 -mcmodel=medium -fPIC @mrg37-080A.f03 &> log
!* $ mpiexec -n number_of_cpu's a.out &
!-----------------------------------------------------------------------
!*  Fortrtan 2003/Fortran 2008 by direct write outputs
!*               write(11,'(" arrayx,arrayy(i),arrayz=",3i6)')... 
!-----------------------------------------------------------------------
!
      program macro_particles_code
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h'
!
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
      integer(C_INT) rank,size,ipar,ierror,cl_first
      integer(C_INT) npr,kstart
!
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,vxi,vyi,vzi
      real(C_DOUBLE),dimension(np0) :: xe,ye,ze,vxe,vye,vze,mue,vpe,vhe
!-----------------------------------------------------------------------
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      real(C_DOUBLE)  walltime0,walltime1,walltime2,walltime3
      character(len=10) :: date_now
      character(len=8)  :: time_now
!*
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  & 
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_float)  plodx
      integer(C_INT) kploy,kploz
      common/plotiv/ plodx,kploy,kploz
!*
      character(len=8)  label(8)
!
!-----------------------------------------------------------------------
!*    datum0 ..... job control (not /common/ parameters).
!*    datum1 ..... parameters  (change not allowed at restart).
!*    datum2 ..... plot parameters.
!-----------------------------------------------------------------------
!  Flow:  Start: &datum0 kstart= 0 -> trans: iresrt= 1 -> call restrt
!         Restart: &datum0 kstart= 1 -> call restrt
!
      real(C_DOUBLE) tfinal,cptot,cpu1
      integer(C_INT) iresrt,i,j,k,l,istop,iwrt,nframe
      common/irest/ iresrt
!
      real(C_DOUBLE) arb,zcent,ycent1,ycent2,Ez00
      common/profl/  arb,zcent,ycent1,ycent2,Ez00
!
      integer(C_INT),dimension(mxyz) :: arrayx,arrayy,arrayz
      common/array1d/ arrayx,arrayy,arrayz
!
      integer(C_INT) ifilxs,ifilys,ifilzs
      common/damper/ ifilxs,ifilys,ifilzs
!
!*  kstart..... used in /init/, /trans/, /fulmov/, in param_080A.h.
!   Ez00   .... Ez00 x Ba
!
      namelist/datum0/  kstart,tfinal,cptot,istop
      namelist/datum1/  dt,xmax,ymax,zmax,wspec,qspec,   &
                        nspec,vdr,vbeam,Ez00,            &
                        wce_by_wpe,te_by_ti,veth,thb,    &
                        aimpl,rwd,epsln1,pi,             &
                        ifilx,ifily,ifilz,               &
                        ifilxs,ifilys,ifilzs,            &
                        itermx,iterfx,itersx,            &
                        iloadp,arb,nha
      namelist/datum2/                                   &
               vlima,vlimb,bmin,emin,npwrt,nfwrt,nplot,nhist
!
      character(len=8) :: fortr51,fortr61
      common/fortr50/ fortr51(8),fortr61(8)
!
!*----------------------------------------------------------------------
!***********************************************************************
!*    Read-in control parameters.                                      *
!***********************************************************************
!
!     open (unit=05,file=praefixs,form='formatted')  !<- L.357
!     open (unit=12,file=praefixc//'.12'//suffix1,form='unformatted') ! L.12155
      open (unit=15,file=praefixc//'.15'//suffix2,form='unformatted')
      open (unit=18,file=praefixc//'.18'//suffix2,form='unformatted')
      close(15)
      close(18)
!
!************************
!*  Initial MPI setup.  *
!************************
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      ipar = 1 + rank     !! pe #= 1,2,3...
!
      io_pe = 0
      if(ipar.eq.1) io_pe = 1
!
! Check results if it is necessary
      fortr51(1)='fortr.51'
      fortr51(2)='fortr.52'
      fortr51(3)='fortr.53'
      fortr51(4)='fortr.54'
      fortr51(5)='fortr.55'
      fortr51(6)='fortr.56'
      fortr51(7)='fortr.57'
      fortr51(8)='fortr.58'
!*
      fortr61(1)='fortr.61'
      fortr61(2)='fortr.62'
      fortr61(3)='fortr.63'
      fortr61(4)='fortr.64'
      fortr61(5)='fortr.65'
      fortr61(6)='fortr.66'
      fortr61(7)='fortr.67'
      fortr61(8)='fortr.68'
!
!     open (unit=50+ipar,file=fortr51(ipar),form='formatted')
!     open (unit=60+ipar,file=fortr61(ipar),form='formatted')
!     close(50+ipar)
!     close(60+ipar)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,form='formatted')
!
        write(11,'("### @mrg37.f03 (3-D periodic/bound system) ###")')
        write(11,'(" Macro EM particle simulation code (implicit)")')
        write(11,*)
        write(11,*) ' mx,myA,mz=',mx,myA,mz
        write(11,*) ' 3*mxyz(the total) =',3*mx*myA*mz
        write(11,*) ' 3*mxyz/kd(divided)=',3*mx*myA*kd
!
        write(11,*)
      end if
!
      label(1)='Collisio'
      label(2)='nless ma'
      label(3)='gnetic r'
      label(4)='econnect'
      label(5)='ion in 3'
      label(6)='-D space'
      label(7)=' ions/el'
      label(8)='ectrons '
!
      call date_and_time_7 (date_now,time_now)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'(8a8,/)') (label(i),i=1,8)
        write(11,'("*date: ",a10,2x,"*time: ",a8,/)') &
                                      date_now,time_now  
        close(11)
      end if
!
      call lbltop (date_now,label)
!
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!   cfpsol: np1(1) = 1, np2(1)= 3*mx*(my+1)*kd
!           nz1(1) = 1, nz2(1)= kd
!   escorr: nps1(1)= 1, nps2(1)= mx*(my+1)*kd
!   poissn:
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!     kd= mz/npc is defined in param_080A.h
!
!             number_of_cpu's
      do k= 1,npc         
      np1(k)= (k-1)*3*mx*(my+1)*kd +1    !<-- np1(1)=1
      np2(k)=     k*3*mx*(my+1)*kd       !<-- np2(2)=3*mx*myA*kd
!
      nz1(k) = (k-1)*kd +1               !<-- nz1(1)=1,  kd+1
      nz2(k) =     k*kd                  !<-- nz2(1)=kd, 2*kd
      end do
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) ' np1(k),np2(k),nz1(k),nz2(k)...'
        write(11,'(2i10,2x,2i10)') (np1(k),np2(k),nz1(k),nz2(k),k=1,npc)
        write(11,*)
!
        close(11)
      end if
!
!
!* Tables to arrays in (arrayx, arrayy, arrayz)
      npr= np0
! 
      l= 0
!*
      do k= 0,mz-1  !<<-  k=0 (start from 0)
      do j= 0,my    !     j=0
      do i= 0,mx-1  !     i=0
      l= l+1
!
      arrayx(l)= i
      arrayy(l)= j
      arrayz(l)= k
      end do
      end do
      end do
!*
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'One-dimensional array is... ',l
        write(11,'(" arrayx,arrayy(i),arrayz=",3i6)') &
                    (arrayx(i),arrayy(i),arrayz(i),i=1,30)
        close(11)
      end if
!
!********************************************************
!*  The macro-particle fully-implicit simulation code.  *
!********************************************************
!
!  All nodes must execute at the beginning
!*                       ++++++++
      open (unit=05,file=praefixs,form='formatted')
!
      read(05,datum0)
      read(05,datum1)
      read(05,datum2)
      close(05)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) '*Initial parameters... &datum0, &datum1 are:'
        write(11,datum0)
        write(11,datum1)
        write(11,datum2)
!
        close(11)
      end if
!*
!    +++++++++++++++++++++++++++
      ifilx = 1   !! &datum1 in rec_3d
      ifily = 1   !!
      ifilz = 1   !!
!
      ifilxs= 3   ! saved in /damper/
      ifilys= 3 
      ifilzs= 3 
!    +++++++++++++++++++++++++++
!
!-----------------------------
!*   Restart procedure.  
!-----------------------------
!
      if(kstart.ne.0) then
        iresrt= 1
!
        call restrt (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                     xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                     npr,ipar,size,iresrt) 
      end if
!
!**********************************************
!*  Initial loading of particles and fields.  *
!**********************************************
!
      call init (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                 xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                 npr,kstart)
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      cl_first= 1   ! measure timing
      call clocks (walltime0,size,cl_first)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'(" ** initial walltime ",f7.3," secs")') walltime0
        close(11)
      end if
!
!*************************
!*    Main part starts   *
!*************************
!
      istop= 0
!
!* Plot results forta.77.ps are initilized here....
      if(io_pe.eq.1) then
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps', &
              form='formatted')
!
        nframe= 4
        call gopen (nframe)
        close(77)
!
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) '# Plot results are output as '// &
                     praefixc//'.77'//suffix2//'.ps'
        close(11)
      end if
!
!***
      cl_first= 2
      call clocks (walltime1,size,cl_first)
!
      if(io_pe.eq.1) then
        if(iwrt(it,nplot).eq.0) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("**** it=",i5," *****")') it
        close(11)
        end if
      end if
!
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++
      call trans (xi,yi,zi,vxi,vyi,vzi,               &
                  xe,ye,ze,vxe,vye,vze,mue,vpe,vhe,   &
                  np1,np2,nz1,nz2,npr,kstart,         &
                  tfinal,cptot,ipar,size)
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!************************
!*   Diagnostic plots   *
!************************
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) '* Plot hist (final)  it=',it
        close(11)
!
!   The history plot is made.
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!   ---------------------------          ++++++
        call histry
!   ---------------------------
        close(77)
      end if
!
      cl_first= 2
      call clocks (walltime2,size,cl_first)
!
      cpu1= walltime2 - walltime0 
      call date_and_time_7 (date_now,time_now)
!
!****************************************
!*   Final procedure to close this job. *                              *
!****************************************
!
      istop= 1
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("** time used (run)",f9.2," secs.,  it=",i5, &
                                             " ****",/)') cpu1,it
        close(11) 
      end if
!
!***
!  Final procedures for future runs
!
      iresrt= 2
      call restrt (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                   xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                   npr,ipar,size,iresrt) 
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
        write(11,*) 'Write L.500 (Final procedure)'
        close(11)
      end if
!       
!***
!
      cl_first= 2
      call clocks (walltime3,size,cl_first)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("** All time (run + restrt)",f9.2," secs.")') &
                                               walltime3 -walltime0
!
        call date_and_time_7 (date_now,time_now)
        write(11,'("*Final date: ",a10,2x,"*time: ",a8,/)') &
                                      date_now,time_now  
        close(11)
!
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
        call plote
        close(77)
      end if
!
      call mpi_finalize (ierror)
!
      stop
      end program macro_particles_code
!
!
!-----------------------------------------------------------------------
      block data
!-----------------------------------------------------------------------
!*    the seed of random numbers must be odd integers.
!
!     use, intrinsic :: iso_c_binding
!
      implicit none
!
      real*4  plodx
      integer*4   kploy,kploz,ir1,ir2
      common/plotiv/ plodx,kploy,kploz
!
      common/ranfa/  ir1  !! integer for ir1
      common/ranfb/  ir2
!***
      data   ir1/3021/,ir2/7331/,                                 &
             plodx/3.99d0/,kploy/2/,kploz/3/
      end block data
!
!
!-----------------------------------------------------------------------
      integer function iwrt (it,interv) 
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
      integer(C_INT) it,interv
!
!  Hit is true if mod(it,..)= 0 
      if(mod(it,interv).eq.0) then
        iwrt= 0
      else
        iwrt= 1
      end if
!
!  Special cases
!     if(it.eq.1) then
!       iwrt= 0
!     end if
!
!     if(it.ge.77777) then
!       iwrt= 1 
!     end if
!
      return
      end function iwrt
!
!
!-----------------------------------------------------------------------
      subroutine trans (xi,yi,zi,vxi,vyi,vzi,               &
                        xe,ye,ze,vxe,vye,vze,mue,vpe,vhe,   &
                        np1,np2,nz1,nz2,npr,kstart,         &
                        tfinal,cptot,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
! 
      include 'mpif.h'
      include 'param_080A.h'
!
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
      integer(C_INT) npr,kstart,ipar,size
      real(C_DOUBLE) tfinal,cptot
!
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,vxi,vyi,vzi,           &
                                       xe,ye,ze,vxe,vye,vze,mue,vpe,vhe
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                         ex,ey,ez,bx,by,bz,        &
                                         ex0,ey0,ez0,bx0,by0,bz0,  &
                                         qix,qiy,qiz,qex,qey,qez,  &
                                         emx,emy,emz,qi,qe,amu,avhh,&
                                         pot
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/sclpot/ pot
!---------------------------------------------------------
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!* 
      real(C_DOUBLE) tdec(3000)  !<-- DOUBLE
      common/ehist/  tdec
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) dtsav,adtsav,hdtsav,gnu0,cs
      real(C_DOUBLE) ase,asb,asl,we,wb,wl
      integer(C_INT) iterm,iterf,iters
      common/emiter/ ase,asb,asl,we,wb,wl,iterm,iterf,iters
!***
      integer(C_INT) iresrt,ipc,iwrt,npl,iaverg,cl_first,i,j,k
      common/irest/ iresrt
      real(C_DOUBLE) walltime1,walltime2,walltime3,walltime4,walltime5
!
!*----------------------------------------------------------------------
!***********************************************************************
!*  1. Initial state only for it= 1                                    *
!*     the ion drift is given by vd0 along the z-direction.            *
!***********************************************************************
!
      cl_first= 1
      call clocki (walltime1,size,cl_first)
!
      if(kstart.eq.0) then  !! only at it= 0
!        +++++++++++
!
        it= 0
        t = 0
!
        iaver= 0
        ldec= 1  ! at the start time
        tdec(ldec)= t
!
        dtsav=  dt
        adtsav= adt
        hdtsav= hdt
!
        dt=  0
        adt= 0
        hdt= 0
!
!  Accumulate moments
!
        ipc= 1  !!      
        call fulmov (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                     npr,ipc,1,ipar,size)
!
        call fulmov (xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                     npr,ipc,2,ipar,size)
!
          cl_first= 2
          call clocki (walltime2,size,cl_first)
!
        call emfld0
!
          cl_first= 2
          call clocki (walltime3,size,cl_first)
          call clocki (walltime4,size,cl_first)
!
        dt=  dtsav
        adt= adtsav
        hdt= hdtsav
!*                --- for it=0, use /diag1/ ---
        kstart= 1
        go to 300
      end if
!
!-----------------------------------------------------------------------
!*  The new time-cycle starts here.                                    *
!-----------------------------------------------------------------------
!***********************************************************************
!*  2. Accumulate the source moments and solve the fully-implicit      *
!*     coupled field-particle equations.                               *
!***********************************************************************
!-----------------------------------------------------------------------
!*    Correction to the longitudinal part is made separately outside
!*    the field iteration loop (a good approximation).
!-------------------------------------------------------- date: 6/89 ---
!*    epsln1 : tolerance of the convergence (emfild, escorr).
!
 1000 cl_first= 2
      call clocki (walltime1,size,cl_first)
!
      if(mod(it,5).eq.1 .and. t.ge.tfinal) then
!
        if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
        write(11,*) '# t.ge.tfinal'
        close(11)
        end if
!
        return
      end if
!
      if((walltime1/60.d0).gt.cptot) then
!
        if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
        write(11,*) '# walltime1/60.d0).gt.cptot'
        close(11)
        end if
!
        return  !<-- in minutes
      end if
!     +++++++++++++++++++++++++++++++++++++++++++
!
      it= it +1 
      t = t +dt
!        +++++++++++++++++
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
        ldec= ldec +1  ! only at this time
        tdec(ldec)= t
      end if
!
!   Accumulate the moments and do emfild    
!
      call prefld
!
      ipc= 1  
      call fulmov (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                   npr,ipc,1,ipar,size)
!
      call fulmov (xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                   npr,ipc,2,ipar,size)
!
        cl_first= 2
        call clocki (walltime2,size,cl_first)
!
      call emfild (np1,np2,nz1,nz2,ipar)
!
!
!***********************************************************************
!*  3.  Update particles : x(n) to x(n+1)  (ipc=0).                    *
!***********************************************************************
!
        cl_first= 2
        call clocki (walltime3,size,cl_first)
!
!  For move
      ipc= 0 
      call fulmov (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                   npr,ipc,1,ipar,size)
!
      call fulmov (xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                   npr,ipc,2,ipar,size)
!
!***********************************************************************
!*  4.  Renewal : ex0 <--- ex  (after particle-move).                  *
!***********************************************************************
!  Rewrite ex() to be one-last exa() to advance forward
!
  300 continue
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ex0(i,j,k)= ex(i,j,k)
      ey0(i,j,k)= ey(i,j,k)
      ez0(i,j,k)= ez(i,j,k)
      bx0(i,j,k)= bx(i,j,k)  !<-- wce_by_wpe is not included
      by0(i,j,k)= by(i,j,k)
      bz0(i,j,k)= bz(i,j,k)
      end do
      end do
      end do
!
        cl_first= 2
        call clocki (walltime4,size,cl_first)
!
!     if(io_pe.eq.1) then
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
!
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("# it,t=",i6,f10.3,"  total, this step(sec)=", &
                   2f10.3,/,"  ful(1),em,ful(0)=",3f10.3,/)')     &
                      it,t,walltime4,walltime4-walltime1,        &
                      walltime2-walltime1,walltime3-walltime2,   &
                      walltime4-walltime3
        close(11)
      end if
!
!***
      if(iwrt(it,nhist).eq.0 .and. it.gt.1) then
!                +++++
        iresrt= 2
        call restrt (xi,yi,zi,vxi,vyi,vzi,qspec(1),wspec(1),  &
                     xe,ye,ze,vxe,vye,vze,qspec(2),wspec(2),  &
                     npr,ipar,size,iresrt) 
!
        if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.11'//suffix2,             & 
                status='unknown',position='append',form='formatted')
          write(11,*) 'Write L.840 it=',it
          close(11)
        end if
!
      end if
!       
!***
!
      if(mod(it,nha).ne.0 .or. it.eq.0) go to 1000
!               +++  ++ + 
!
!***********************************************************************
!*  5. Diagnostic routines diag1                                       *
!***********************************************************************
!   Only mod(it,nha)= 0
!
      iaver= iaver +1
      call faverg (iaver,npl)
!
!* after faverg ! 
!
      if(.false.) then
        call diag1 (xi,yi,zi,vxi,vyi,vzi,qspec(1),npr,1)
!
        call diag1 (xe,ye,ze,vxe,vye,vze,qspec(2),npr,2)
!
        call fvplot
      end if
!
!     gnu0= 1./200.
      gnu0= 0.d0
      cs=  veth/sqrt(2.d0*wspec(1))
      call dragco (gnu0,cs,iaverg,npl)
!
      if(iwrt(it,nha).eq.0) then
         iaver= 0
         call freset
      end if
!
      go to 1000
!
      return
      end subroutine trans
!
!
!-----------------------------------------------------------------------
      subroutine faverg (iaver,nav)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      integer(C_INT) iaver,nav
      real(C_DOUBLE),dimension(mxyzA) :: ex,ey,ez,bx,by,bz,        &
                                         ex0,ey0,ez0,bx0,by0,bz0
      real(C_DOUBLE),dimension(mxyzA) :: qix,qiy,qiz,qex,qey,qez,  &
                                         emx,emy,emz,qi,qe,amu,avhh  
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
!
!*    --- averaged ---        ++++
      real(C_float),dimension(mxyzA) :: cix,ciy,ciz,cex,cey,cez,       &
                                        avex,avey,avez,avbx,avby,avbz, &
                                        avqi,avqe
      common/plotav/ cix,ciy,ciz,cex,cey,cez,       &
                     avex,avey,avez,avbx,avby,avbz, &
                     avqi,avqe
      real(C_float)  rnav
      integer(C_INT) ijk
!
      if(iaver.gt.nav) return
!
      do ijk= 0,mxyzA   !<-- real*4 for plots
      cix(ijk)= cix(ijk) +qix(ijk)
      ciy(ijk)= ciy(ijk) +qiy(ijk)
      ciz(ijk)= ciz(ijk) +qiz(ijk)
      cex(ijk)= cex(ijk) +qex(ijk)
      cey(ijk)= cey(ijk) +qey(ijk)
      cez(ijk)= cez(ijk) +qez(ijk)
!
      avqi(ijk)= avqi(ijk) +qi(ijk)
      avqe(ijk)= avqe(ijk) +qe(ijk)
!
      avex(ijk)= avex(ijk) +ex(ijk)
      avey(ijk)= avey(ijk) +ey(ijk)
      avez(ijk)= avez(ijk) +ez(ijk)
      avbx(ijk)= avbx(ijk) +bx(ijk)
      avby(ijk)= avby(ijk) +by(ijk)
      avbz(ijk)= avbz(ijk) +bz(ijk)
      end do
!***
!
      if(iaver.ne.nav) return
      rnav= 1.d0/nav
!
      do ijk= 1,mxyzA  !<<-- real*4
      cix(ijk)= rnav*cix(ijk)
      ciy(ijk)= rnav*ciy(ijk)
      ciz(ijk)= rnav*ciz(ijk)
      cex(ijk)= rnav*cex(ijk)
      cey(ijk)= rnav*cey(ijk)
      cez(ijk)= rnav*cez(ijk)
!
      avqi(ijk)= rnav*avqi(ijk)
      avqe(ijk)= rnav*avqe(ijk)
!
      avex(ijk)= rnav*avex(ijk)
      avey(ijk)= rnav*avey(ijk)
      avez(ijk)= rnav*avez(ijk)
      avbx(ijk)= rnav*avbx(ijk)
      avby(ijk)= rnav*avby(ijk)
      avbz(ijk)= rnav*avbz(ijk)
      end do
!
      return
      end subroutine faverg
!
!
!-----------------------------------------------------------------------
      subroutine freset
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!          +++++++            +++++ 
      real(C_float),dimension(mxyzA) :: cix,ciy,ciz,cex,cey,cez,       &
                                        avex,avey,avez,avbx,avby,avbz, &
                                        avqi,avqe 
      common/plotav/ cix,ciy,ciz,cex,cey,cez,       &
                     avex,avey,avez,avbx,avby,avbz, &
                     avqi,avqe
      integer(C_INT) ijk
!
      do ijk= 1,mxyzA
      cix(ijk)= 0
      ciy(ijk)= 0
      ciz(ijk)= 0
      cex(ijk)= 0
      cey(ijk)= 0
      cez(ijk)= 0
!
      avqi(ijk)= 0
      avqe(ijk)= 0
!
      avex(ijk)= 0
      avey(ijk)= 0
      avez(ijk)= 0
      avbx(ijk)= 0
      avby(ijk)= 0
      avbz(ijk)= 0
      end do
!
      return
      end subroutine freset
!
!
!-----------------------------------------------------------------------
      subroutine dragco (gnu0,cs,iaver,nav)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: gnu
      common/dragcf/ gnu
      real(C_DOUBLE) gnu0,cs
      integer(C_INT) iaver,nav
!          +++++++            +++++++++
      real(C_float),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                            cix,ciy,ciz,cex,cey,cez,  &
                                            avex,avey,avez,avbx,avby, &
                                            avbz,avqi,avqe
      common/plotav/ cix,ciy,ciz,cex,cey,cez,       &
                     avex,avey,avez,avbx,avby,avbz, &
                     avqi,avqe
!
      real(C_DOUBLE) vye
      integer(C_INT) i,j,k
!
      if(iaver.eq.nav) then
!
         do k= 0,mz-1
         do j= 0,my 
         do i= 0,mx-1
         gnu(i,j,k)= 0
         end do
         end do
         end do
!
         do k= 0,mz-1
         do j= 0,my 
         do i= 0,mx-1
         vye= abs(cey(i,j,k)/avqe(i,j,k))
         if(vye.ge.cs) then
            gnu(i,j,k)= gnu0
         end if
         end do
         end do
         end do
      end if
!
      return
      end subroutine dragco
!
!
!-----------------------------------------------------------------------
      subroutine fulmov (x,y,z,vx,vy,vz,qmult,wmult,npr,ipc,ksp,ipar,size)
!-----------------------------------------------------------------------
!********************************************************** 6/02/89 ****
!*    Particle mover with /analytical/ formula for v*b rotation.       *
!*    (v(n+1/2) of the intermediate time level does not appear.)       *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      real(C_DOUBLE),dimension(np0) :: rxl,ryl,rzl,vxj,vyj,vzj
      real(C_DOUBLE) qmult,wmult
      integer(C_INT) npr,ipc,ksp,ipar,size,MPIerror 
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                             ex,ey,ez,bx,by,bz,        &
                                             ex0,ey0,ez0,bx0,by0,bz0,  &
                                             qix,qiy,qiz,qex,qey,qez,  &
                                             emx,emy,emz,qi,qe,amu,avhh
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                             exa,eya,eza,bxa,bya,bza
!---------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!----------------------------------------------------------------------
! 
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      real(C_DOUBLE) wkix,wkih,wkex,wkeh,wxsq,whsq
      common/wkinel/ wkix,wkih,wkex,wkeh
!
      real(C_DOUBLE) arb,zcent,ycent1,ycent2,Ez00
      common/profl/  arb,zcent,ycent1,ycent2,Ez00
!
      real(C_DOUBLE) hh,ht,ht2,                               &
                     xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr,   &
                     exi,eyi,ezi,bxi,byi,bzi,                 &
                     bsqi,acx,acy,acz,ach,dvx,dvy,dvz,vy0,ranfp
      integer(C_INT) l,i,il,ir,ip,j,jl,jr,jp,k,kl,kr,kp,      &
                     iwrt,syme,symb
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!***********************************************************************
!* 0. Define the e,b of time level (n+aimpl) when ipc=0 or >1.         *
!*    but, for ipc=1, exa= ex0.                                        *
!***********************************************************************
!
      do k= 0,mz-1
      do j= 0,my  
      do i= 0,mx-1
      exa(i,j,k)= aimpl*ex(i,j,k) +(1.d0-aimpl)*ex0(i,j,k)
      eya(i,j,k)= aimpl*ey(i,j,k) +(1.d0-aimpl)*ey0(i,j,k)
      eza(i,j,k)= aimpl*ez(i,j,k) +(1.d0-aimpl)*ez0(i,j,k)
!
      bxa(i,j,k)= aimpl*bx(i,j,k) +(1.d0-aimpl)*bx0(i,j,k) +bxc
      bya(i,j,k)= aimpl*by(i,j,k) +(1.d0-aimpl)*by0(i,j,k) +byc
      bza(i,j,k)= aimpl*bz(i,j,k) +(1.d0-aimpl)*bz0(i,j,k) +bzc
      end do                                   ! for particles only
      end do
      end do
!                                  !<-- extended regions
      call outmesh3 (exa,eya,eza)  !<-- Make -2,-1,...,mx,mx+1 meshes
      call outmesh3 (bxa,bya,bza)
!
      syme= -1
      call filt3e (exa,eya,eza,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,syme)
!
      symb= +1
      call filt3e (bxa,bya,bza,bxc,byc,bzc,ifilx,ifily,ifilz,symb)
!
      hh=  dt*qmult/wmult
      ht=  0.5d0*hh
      ht2= ht**2
!
!***********************************************************************
!* 1. Advance particles:                                               *
!***********************************************************************
!  << Needs the left-side grid for y, the central grids for x,z >>
!
!   To improve orbit tracking, the magnetic field in the Lorentz term
!   is evaluated at the mid-orbit positions (also in do 300).
!
      do l= ipar,npr,size
      rxl(l) = x(l) +hdt*vx(l)  !! ipc= 0,1
      ryl(l) = y(l) +hdt*vy(l) 
      rzl(l) = z(l) +hdt*vz(l)
      end do
!
      call partbcEST (rxl,ryl,rzl,npr,ipar,size)  !<-- estimated
!
!
      wkix= 0
      wkih= 0
!
      do l= ipar,npr,size 
      ip= int(hxi*rxl(l) +0.500000001d0)
      jp= int(hyi*ryl(l) +0.000000001d0)
      kp= int(hzi*rzl(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      fyl= hyi*ryl(l) -jp
      fyr= 1.d0 -fyl
!
!# These are outside of the target regions ###
!
      if(jp.ge.my) then
        jr= my+1  ! to be safe within labels
        jl= my    !
        fyr= 0
        fyl= 1.d0
      else if(jp.lt.0) then 
        jr=  0
        jl= -1
        fyr= 1.d0
        fyl= 0
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
      xx = hxi*rxl(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      zz = hzi*rzl(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      exi = fyr*  &
           ((exa(ir,jr,kr)*fxr+exa(i,jr,kr)*fxc+exa(il,jr,kr)*fxl)*fzr &
           +(exa(ir,jr,k )*fxr+exa(i,jr,k )*fxc+exa(il,jr,k )*fxl)*fzc &
           +(exa(ir,jr,kl)*fxr+exa(i,jr,kl)*fxc+exa(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((exa(ir,jl,kr)*fxr+exa(i,jl,kr)*fxc+exa(il,jl,kr)*fxl)*fzr &
           +(exa(ir,jl,k )*fxr+exa(i,jl,k )*fxc+exa(il,jl,k )*fxl)*fzc &
           +(exa(ir,jl,kl)*fxr+exa(i,jl,kl)*fxc+exa(il,jl,kl)*fxl)*fzl)
!
      eyi = fyr*  &
           ((eya(ir,jr,kr)*fxr+eya(i,jr,kr)*fxc+eya(il,jr,kr)*fxl)*fzr &
           +(eya(ir,jr,k )*fxr+eya(i,jr,k )*fxc+eya(il,jr,k )*fxl)*fzc &
           +(eya(ir,jr,kl)*fxr+eya(i,jr,kl)*fxc+eya(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((eya(ir,jl,kr)*fxr+eya(i,jl,kr)*fxc+eya(il,jl,kr)*fxl)*fzr &
           +(eya(ir,jl,k )*fxr+eya(i,jl,k )*fxc+eya(il,jl,k )*fxl)*fzc &
           +(eya(ir,jl,kl)*fxr+eya(i,jl,kl)*fxc+eya(il,jl,kl)*fxl)*fzl)
!
      ezi = fyr*  &
           ((eza(ir,jr,kr)*fxr+eza(i,jr,kr)*fxc+eza(il,jr,kr)*fxl)*fzr &
           +(eza(ir,jr,k )*fxr+eza(i,jr,k )*fxc+eza(il,jr,k )*fxl)*fzc &
           +(eza(ir,jr,kl)*fxr+eza(i,jr,kl)*fxc+eza(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((eza(ir,jl,kr)*fxr+eza(i,jl,kr)*fxc+eza(il,jl,kr)*fxl)*fzr &
           +(eza(ir,jl,k )*fxr+eza(i,jl,k )*fxc+eza(il,jl,k )*fxl)*fzc &
           +(eza(ir,jl,kl)*fxr+eza(i,jl,kl)*fxc+eza(il,jl,kl)*fxl)*fzl)
!
!  The DC magnetic term bxc-bzc is added here !
      bxi = fyr*  &
           ((bxa(ir,jr,kr)*fxr+bxa(i,jr,kr)*fxc+bxa(il,jr,kr)*fxl)*fzr & 
           +(bxa(ir,jr,k )*fxr+bxa(i,jr,k )*fxc+bxa(il,jr,k )*fxl)*fzc &
           +(bxa(ir,jr,kl)*fxr+bxa(i,jr,kl)*fxc+bxa(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((bxa(ir,jl,kr)*fxr+bxa(i,jl,kr)*fxc+bxa(il,jl,kr)*fxl)*fzr &
           +(bxa(ir,jl,k )*fxr+bxa(i,jl,k )*fxc+bxa(il,jl,k )*fxl)*fzc &
           +(bxa(ir,jl,kl)*fxr+bxa(i,jl,kl)*fxc+bxa(il,jl,kl)*fxl)*fzl)
!
      byi = fyr*  &
           ((bya(ir,jr,kr)*fxr+bya(i,jr,kr)*fxc+bya(il,jr,kr)*fxl)*fzr &
           +(bya(ir,jr,k )*fxr+bya(i,jr,k )*fxc+bya(il,jr,k )*fxl)*fzc &
           +(bya(ir,jr,kl)*fxr+bya(i,jr,kl)*fxc+bya(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((bya(ir,jl,kr)*fxr+bya(i,jl,kr)*fxc+bya(il,jl,kr)*fxl)*fzr &
           +(bya(ir,jl,k )*fxr+bya(i,jl,k )*fxc+bya(il,jl,k )*fxl)*fzc &
           +(bya(ir,jl,kl)*fxr+bya(i,jl,kl)*fxc+bya(il,jl,kl)*fxl)*fzl)
!
      bzi = fyr*  &
           ((bza(ir,jr,kr)*fxr+bza(i,jr,kr)*fxc+bza(il,jr,kr)*fxl)*fzr &
           +(bza(ir,jr,k )*fxr+bza(i,jr,k )*fxc+bza(il,jr,k )*fxl)*fzc &
           +(bza(ir,jr,kl)*fxr+bza(i,jr,kl)*fxc+bza(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((bza(ir,jl,kr)*fxr+bza(i,jl,kr)*fxc+bza(il,jl,kr)*fxl)*fzr &
           +(bza(ir,jl,k )*fxr+bza(i,jl,k )*fxc+bza(il,jl,k )*fxl)*fzc &
           +(bza(ir,jl,kl)*fxr+bza(i,jl,kl)*fxc+bza(il,jl,kl)*fxl)*fzl)
!
      bsqi= bxi**2 +byi**2 +bzi**2
      acx= exi +vy(l)*bzi -vz(l)*byi
      acy= eyi +vz(l)*bxi -vx(l)*bzi
      acz= ezi +vx(l)*byi -vy(l)*bxi
!
      ach= exi*bxi +eyi*byi +ezi*bzi
      dvx= ( acx +ht2*ach*bxi +ht*(acy*bzi -acz*byi) )/(1.d0+ht2*bsqi)
      dvy= ( acy +ht2*ach*byi +ht*(acz*bxi -acx*bzi) )/(1.d0+ht2*bsqi)
      dvz= ( acz +ht2*ach*bzi +ht*(acx*byi -acy*bxi) )/(1.d0+ht2*bsqi)
!
      wkix= wkix +0.5d0*(acx**2 +acy**2 +acz**2)
      wkih= wkih +0.5d0* ach**2
!
!  For update
!   ++++++++++++++++++++++++++++++++++++++++++++
      if(ipc.eq.0) then
!         hh= dt*qmult/wmult, x and vx are at the time t^n
      x(l) = x(l) +dt*(vx(l) +0.5d0*hh*dvx) 
      y(l) = y(l) +dt*(vy(l) +0.5d0*hh*dvy) 
      z(l) = z(l) +dt*(vz(l) +0.5d0*hh*dvz)
! 
      vx(l)= vx(l) +hh*dvx
      vy(l)= vy(l) +hh*dvy
      vz(l)= vz(l) +hh*dvz
!
!  For ipc=1, accumulate srimp 1-4. use vx(n+1/2)   
      else if(ipc.ge.1) then 
!         the decentering mechanism is here         
      vxj(l)= vx(l) +aimpl*hh*dvx  
      vyj(l)= vy(l) +aimpl*hh*dvy
      vzj(l)= vz(l) +aimpl*hh*dvz
!
      rxl(l)= x(l) +adt*( vx(l) +0.5d0*hh*dvx )
      ryl(l)= y(l) +adt*( vy(l) +0.5d0*hh*dvy )
      rzl(l)= z(l) +adt*( vz(l) +0.5d0*hh*dvz )
      end if
!   ++++++++++++++++++++++++++++++++++++++++++++
      end do
!
! Synchronize
      call mpi_allreduce (wkix,wxsq,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (wkih,whsq,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,MPIerror)
      wkix= wxsq
      wkih= whsq
!
!        +++++++++++++++++
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
      if(ksp.eq.1) then
        edec(ldec,5)= wkix
        edec(ldec,6)= wkih
      else if(ksp.eq.2) then
        edec(ldec,7)= wkix
        edec(ldec,8)= wkih
      end if
      end if
!
!***********************************************************************
!* 2. Accumulate moments for the e.m. field solver  (srimp1 - 2).      *
!***********************************************************************
!
!  For update: exact values
      if(ipc.eq.0) then
!
        call partbc (x,y,z,vx,vy,vz,npr,ipar,size)
!
!  Drive E x B  1.d-2/wce= 5.d-2 ?
!        Ez00=0.25d-2, wce=0.2
!
        do l= ipar,npr,size 
        if((abs(z(l)-zcent ).lt.0.15d0*zmax) .and.  &
          ((abs(y(l)-ycent2).lt.0.025d0*ymax) .or.  &
           (abs(y(l)-ycent1).lt.0.025d0*ymax)) ) then
!**
          ip= int(hxi*x(l) +0.500000001d0)
          jp= int(hyi*y(l) +0.000000001d0)
          kp= int(hzi*z(l) +0.500000001d0)
!**
!           if(ranfp(0.d0).gt.0.99d0) then
!         if(ranfp(0.d0).gt.0.995d0) then
          if(ranfp(0.d0).gt.0.999d0) then
          vy0= Ez00/bxa(ip,jp,kp)
!
          if(abs(y(l)-ycent2).lt.0.05d0*ymax) then
            vy(l)= vy(l) - vy0
!
          else if(abs(y(l)-ycent1).lt.0.05d0*ymax) then
            vy(l)= vy(l) + vy0
          end if
          end if
        end if
        end do
      end if
!
!
!---------------------------
!  For accumulating
!---------------------------
!
      if(ipc.ge.1) then
!                    <-- rxl-rzl and vyj
!
        call partbc (rxl,ryl,rzl,vxj,vyj,vzj,npr,ipar,size) 
!
        if(ksp.eq.1) then
          call srimp1 (rxl,ryl,rzl,vxj,vyj,vzj,qmult,qix,qiy,qiz,    &
                                            npr,ipar,size)
          call srimp2 (rxl,ryl,rzl,qmult,qi,npr,ipar,size)
!
        else if(ksp.eq.2) then
          call srimp1 (rxl,ryl,rzl,vxj,vyj,vzj,qmult,qex,qey,qez,    &
                                            npr,ipar,size)
          call srimp2 (rxl,ryl,rzl,qmult,qe,npr,ipar,size)
        end if
      end if
!
      return
      end subroutine fulmov
!
!
!-----------------------------------------------------------------------
      subroutine fulmv2 (x,y,z,vx,vy,vz,qmult,wmult,npr,ksp,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      real(C_DOUBLE),dimension(np0) :: rxl,ryl,rzl
      real(C_DOUBLE) qmult,wmult
      integer(C_INT) npr,ksp,ipar,size
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                             ex,ey,ez,bx,by,bz,       &
                                             ex0,ey0,ez0,bx0,by0,bz0, &
                                             qix,qiy,qiz,qex,qey,qez, &
                                             emx,emy,emz,qi,qe
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                             exa,eya,eza,bxa,bya,bza
!---------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!----------------------------------------------------------------------
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
                     exi,eyi,ezi,bxi,byi,bzi,               &
                     bsqi,acx,acy,acz,ach,dvx,dvy,dvz,      &
                     hh,ht,ht2
      integer(C_INT) l,i,il,ir,ip,j,jl,jr,jp,k,kl,kr,kp,    &
                     syme,symb
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!***********************************************************************
!* 1. Move particles:                                                  *
!***********************************************************************
!-----------------------------------------------------------------------
!*    Define e,b of time level (n+aimpl).
!-----------------------------------------------------------------------
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      exa(i,j,k)= aimpl*ex(i,j,k) +(1.d0-aimpl)*ex0(i,j,k)
      eya(i,j,k)= aimpl*ey(i,j,k) +(1.d0-aimpl)*ey0(i,j,k)
      eza(i,j,k)= aimpl*ez(i,j,k) +(1.d0-aimpl)*ez0(i,j,k)
!
      bxa(i,j,k)= aimpl*bx(i,j,k) +(1.d0-aimpl)*bx0(i,j,k) +bxc
      bya(i,j,k)= aimpl*by(i,j,k) +(1.d0-aimpl)*by0(i,j,k) +byc
      bza(i,j,k)= aimpl*bz(i,j,k) +(1.d0-aimpl)*bz0(i,j,k) +bzc
      end do                                 !!! for particles only
      end do
      end do
!                                  !<-- extend regions
      call outmesh3 (exa,eya,eza)
      call outmesh3 (bxa,bya,bza)
!
      syme= -1
      call filt3e (exa,eya,eza,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,syme)
!
      symb= +1
      call filt3e (bxa,bya,bza,bxc,byc,bzc,ifilx,ifily,ifilz,symb)
!
      hh= dt*qmult/wmult 
      ht= 0.5d0*hh
      ht2= ht**2
!
!
      do l= ipar,npr,size
      rxl(l) = x(l) +hdt*vx(l)
      ryl(l) = y(l) +hdt*vy(l)
      rzl(l) = z(l) +hdt*vz(l)
      end do
!
      call partbcEST (rxl,ryl,rzl,npr,ipar,size)  !<-- estimated
!
!
      do l= ipar,npr,size
      ip= int(hxi*rxl(l) +0.500000001d0)
      jp= int(hyi*ryl(l) +0.000000001d0)
      kp= int(hzi*rzl(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      fyl= hyi*ryl(l) -jp
      fyr= 1.d0 -fyl
!
      if(jp.ge.my) then
        jr= my+1  ! to be safe within labels
        jl= my    !
        fyr= 0
        fyl= 1.d0
      else if(jp.lt.0) then
        jr=  0
        jl= -1
        fyr= 1.d0
        fyl= 0
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
      xx = hxi*rxl(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      zz = hzi*rzl(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      exi = fyr* &
           ((exa(ir,jr,kr)*fxr+exa(i,jr,kr)*fxc+exa(il,jr,kr)*fxl)*fzr  &
           +(exa(ir,jr,k )*fxr+exa(i,jr,k )*fxc+exa(il,jr,k )*fxl)*fzc  &
           +(exa(ir,jr,kl)*fxr+exa(i,jr,kl)*fxc+exa(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((exa(ir,jl,kr)*fxr+exa(i,jl,kr)*fxc+exa(il,jl,kr)*fxl)*fzr  &
           +(exa(ir,jl,k )*fxr+exa(i,jl,k )*fxc+exa(il,jl,k )*fxl)*fzc  &
           +(exa(ir,jl,kl)*fxr+exa(i,jl,kl)*fxc+exa(il,jl,kl)*fxl)*fzl)
!
      eyi = fyr* &
           ((eya(ir,jr,kr)*fxr+eya(i,jr,kr)*fxc+eya(il,jr,kr)*fxl)*fzr  &
           +(eya(ir,jr,k )*fxr+eya(i,jr,k )*fxc+eya(il,jr,k )*fxl)*fzc  &
           +(eya(ir,jr,kl)*fxr+eya(i,jr,kl)*fxc+eya(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((eya(ir,jl,kr)*fxr+eya(i,jl,kr)*fxc+eya(il,jl,kr)*fxl)*fzr  &
           +(eya(ir,jl,k )*fxr+eya(i,jl,k )*fxc+eya(il,jl,k )*fxl)*fzc  &
           +(eya(ir,jl,kl)*fxr+eya(i,jl,kl)*fxc+eya(il,jl,kl)*fxl)*fzl) 
!
      ezi = fyr* &
           ((eza(ir,jr,kr)*fxr+eza(i,jr,kr)*fxc+eza(il,jr,kr)*fxl)*fzr  &
           +(eza(ir,jr,k )*fxr+eza(i,jr,k )*fxc+eza(il,jr,k )*fxl)*fzc  &
           +(eza(ir,jr,kl)*fxr+eza(i,jr,kl)*fxc+eza(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((eza(ir,jl,kr)*fxr+eza(i,jl,kr)*fxc+eza(il,jl,kr)*fxl)*fzr  &
           +(eza(ir,jl,k )*fxr+eza(i,jl,k )*fxc+eza(il,jl,k )*fxl)*fzc  &
           +(eza(ir,jl,kl)*fxr+eza(i,jl,kl)*fxc+eza(il,jl,kl)*fxl)*fzl)
!
      bxi = fyr* &
           ((bxa(ir,jr,kr)*fxr+bxa(i,jr,kr)*fxc+bxa(il,jr,kr)*fxl)*fzr  &
           +(bxa(ir,jr,k )*fxr+bxa(i,jr,k )*fxc+bxa(il,jr,k )*fxl)*fzc  &
           +(bxa(ir,jr,kl)*fxr+bxa(i,jr,kl)*fxc+bxa(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((bxa(ir,jl,kr)*fxr+bxa(i,jl,kr)*fxc+bxa(il,jl,kr)*fxl)*fzr  &
           +(bxa(ir,jl,k )*fxr+bxa(i,jl,k )*fxc+bxa(il,jl,k )*fxl)*fzc  &
           +(bxa(ir,jl,kl)*fxr+bxa(i,jl,kl)*fxc+bxa(il,jl,kl)*fxl)*fzl)
!
      byi = fyr* &
           ((bya(ir,jr,kr)*fxr+bya(i,jr,kr)*fxc+bya(il,jr,kr)*fxl)*fzr  &
           +(bya(ir,jr,k )*fxr+bya(i,jr,k )*fxc+bya(il,jr,k )*fxl)*fzc  &
           +(bya(ir,jr,kl)*fxr+bya(i,jr,kl)*fxc+bya(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((bya(ir,jl,kr)*fxr+bya(i,jl,kr)*fxc+bya(il,jl,kr)*fxl)*fzr  &
           +(bya(ir,jl,k )*fxr+bya(i,jl,k )*fxc+bya(il,jl,k )*fxl)*fzc  &
           +(bya(ir,jl,kl)*fxr+bya(i,jl,kl)*fxc+bya(il,jl,kl)*fxl)*fzl)
!
      bzi = fyr* &
           ((bza(ir,jr,kr)*fxr+bza(i,jr,kr)*fxc+bza(il,jr,kr)*fxl)*fzr  &
           +(bza(ir,jr,k )*fxr+bza(i,jr,k )*fxc+bza(il,jr,k )*fxl)*fzc  &
           +(bza(ir,jr,kl)*fxr+bza(i,jr,kl)*fxc+bza(il,jr,kl)*fxl)*fzl) &
          + fyl*                                                        &
           ((bza(ir,jl,kr)*fxr+bza(i,jl,kr)*fxc+bza(il,jl,kr)*fxl)*fzr  &
           +(bza(ir,jl,k )*fxr+bza(i,jl,k )*fxc+bza(il,jl,k )*fxl)*fzc  &
           +(bza(ir,jl,kl)*fxr+bza(i,jl,kl)*fxc+bza(il,jl,kl)*fxl)*fzl) 
      bsqi= bxi**2 +byi**2 +bzi**2
!
      acx= exi +vy(l)*bzi -vz(l)*byi
      acy= eyi +vz(l)*bxi -vx(l)*bzi
      acz= ezi +vx(l)*byi -vy(l)*bxi
!
      ach= exi*bxi +eyi*byi +ezi*bzi
      dvx= ( acx +ht2*ach*bxi +ht*(acy*bzi -acz*byi) )/(1.d0+ht2*bsqi)
      dvy= ( acy +ht2*ach*byi +ht*(acz*bxi -acx*bzi) )/(1.d0+ht2*bsqi)
      dvz= ( acz +ht2*ach*bzi +ht*(acx*byi -acy*bxi) )/(1.d0+ht2*bsqi)
!
      rxl(l)= x(l) +dt*(vx(l) +0.5d0*hh*dvx)  !<-- decentering
      ryl(l)= y(l) +dt*(vy(l) +0.5d0*hh*dvy)
      rzl(l)= z(l) +dt*(vz(l) +0.5d0*hh*dvz)
      end do
!
      call partbcEST (rxl,ryl,rzl,npr,ipar,size)  !<-- estimated
!
      if(ksp.eq.1) then
        call srimp2 (rxl,ryl,rzl,qmult,qi,npr,ipar,size)
!
      else if(ksp.eq.2) then
        call srimp2 (rxl,ryl,rzl,qmult,qe,npr,ipar,size)
      end if
!
      return
      end subroutine fulmv2
!
!
!-----------------------------------------------------------------------
      subroutine initve (x,y,z,mue,vpe,vhe,vxe,vye,vze,wmult,        &
                         npr,ipar,size)
!-----------------------------------------------------------------------
!***********************************************************************
!*    convert (vxe,vye,vze) into (mue,vp**2,vhe).                      *
!*    mue(l) = 0.5d0*vperp2/bs0, magnetic moment/mass                  *
!*    vpe(l)= vperp2                                                   *
!*    vhe(l) = vpara, parallel velocity (with sign)                    *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,mue,vpe,vhe,vxe,vye,vze
      real(C_DOUBLE) wmult 
      integer(C_INT) npr,ipar,size
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                            ex,ey,ez,bx,by,bz,       &
                                            ex0,ey0,ez0,bx0,by0,bz0, &
                                            bx1,by1,bz1
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
!----------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) wkix,wkih,wkex,wkeh,wxsq,whsq
      common/wkinel/ wkix,wkih,wkex,wkeh
! 
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr,  &
                     bxi,byi,bzi,bs0,vpara,vperp2 
      integer(C_INT) l,ip,jp,kp,i,il,ir,j,jl,jr,k,kl,kr,MPIerror
!
!              The DC magnetic field is added to bx1() !
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      bx1(i,j,k)= bx0(i,j,k) +bxc
      by1(i,j,k)= by0(i,j,k) +byc
      bz1(i,j,k)= bz0(i,j,k) +bzc
      end do
      end do
      end do
!                                  !<-- expand regions
      call outmesh3 (bx1,by1,bz1)
!
      wkex= 0
      wkeh= 0
!***
      do l= ipar,npr,size
      ip= int(hxi*x(l) +0.500000001d0)
      jp= int(hyi*y(l) +0.000000001d0)
      kp= int(hzi*z(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1
!
      fyl= hyi*y(l) -jp
      fyr= 1.d0 -fyl
!
      if(jp.ge.my) then
        jr= my+1  ! to be safe within labels
        jl= my    ! 11/06
        fyr= 0
        fyl= 1.d0
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
        fyr= 1.d0
        fyl= 0
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
      xx = hxi*x(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      zz = hzi*z(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      bxi = fyr*                                                       &
           ((bx1(ir,jr,kr)*fxr+bx1(i,jr,kr)*fxc+bx1(il,jr,kr)*fxl)*fzr &
           +(bx1(ir,jr,k )*fxr+bx1(i,jr,k )*fxc+bx1(il,jr,k )*fxl)*fzc &
           +(bx1(ir,jr,kl)*fxr+bx1(i,jr,kl)*fxc+bx1(il,jr,kl)*fxl)*fzl)&
          + fyl*                                                       &
           ((bx1(ir,jl,kr)*fxr+bx1(i,jl,kr)*fxc+bx1(il,jl,kr)*fxl)*fzr &
           +(bx1(ir,jl,k )*fxr+bx1(i,jl,k )*fxc+bx1(il,jl,k )*fxl)*fzc &
           +(bx1(ir,jl,kl)*fxr+bx1(i,jl,kl)*fxc+bx1(il,jl,kl)*fxl)*fzl)
!
      byi = fyr*                                                       &
           ((by1(ir,jr,kr)*fxr+by1(i,jr,kr)*fxc+by1(il,jr,kr)*fxl)*fzr &
           +(by1(ir,jr,k )*fxr+by1(i,jr,k )*fxc+by1(il,jr,k )*fxl)*fzc &
           +(by1(ir,jr,kl)*fxr+by1(i,jr,kl)*fxc+by1(il,jr,kl)*fxl)*fzl)&
          + fyl*                                                       &
           ((by1(ir,jl,kr)*fxr+by1(i,jl,kr)*fxc+by1(il,jl,kr)*fxl)*fzr &
           +(by1(ir,jl,k )*fxr+by1(i,jl,k )*fxc+by1(il,jl,k )*fxl)*fzc &
           +(by1(ir,jl,kl)*fxr+by1(i,jl,kl)*fxc+by1(il,jl,kl)*fxl)*fzl)
!
      bzi = fyr*                                                       &
           ((bz1(ir,jr,kr)*fxr+bz1(i,jr,kr)*fxc+bz1(il,jr,kr)*fxl)*fzr &
           +(bz1(ir,jr,k )*fxr+bz1(i,jr,k )*fxc+bz1(il,jr,k )*fxl)*fzc &
           +(bz1(ir,jr,kl)*fxr+bz1(i,jr,kl)*fxc+bz1(il,jr,kl)*fxl)*fzl)&
          + fyl*                                                       &
           ((bz1(ir,jl,kr)*fxr+bz1(i,jl,kr)*fxc+bz1(il,jl,kr)*fxl)*fzr &
           +(bz1(ir,jl,k )*fxr+bz1(i,jl,k )*fxc+bz1(il,jl,k )*fxl)*fzc &
           +(bz1(ir,jl,kl)*fxr+bz1(i,jl,kl)*fxc+bz1(il,jl,kl)*fxl)*fzl)
      bs0=  sqrt(bxi**2 +byi**2 +bzi**2)
!
!-----------------------------------------------------------------------
!*  # compute parallel velocity and magnetic moment ###
!-----------------------------------------------------------------------
!
      vpara=  (vxe(l)*bxi +vye(l)*byi +vze(l)*bzi)/bs0
      vperp2= (vxe(l) -vpara*bxi/bs0)**2 +(vye(l) -vpara*byi/bs0)**2  &
                                         +(vze(l) -vpara*bzi/bs0)**2
!
!*    vxe(l) is used always for ---> mue/mass= (1/2)v**2/b.
!     vye(i) for perpendicular**2
!     vze(i) for parallel
!
      mue(l) = 0.5d0*vperp2/bs0
      vpe(l) = vperp2
      vhe(l) = vpara
!
      wkex= wkex +wmult*0.5d0*vperp2
      wkeh= wkeh +wmult*0.5d0*vpara**2
      end do
!
! Synchronize
      call mpi_allreduce (wkex,wxsq,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (wkeh,whsq,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,MPIerror)
      wkex= wxsq
      wkeh= whsq
!
      return
      end subroutine initve
!
!
!-----------------------------------------------------------------------
      subroutine partbc (x,y,z,vx,vy,vz,npr,ipar,size) 
!-----------------------------------------------------------------------
!***********************************************************************
!*    Boundary conditions.                                             *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      integer(C_INT) npr,ipar,size
!------------------------------------------------------
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) l
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
!* The xyz coordinates 
!   the most leftside region of a three-point sharing is i=(-1,0,1),
!  that covers the region x=(-0.5*xmax,1.5*xmax).
!
      dx= hx/2  
      dz= hz/2
!
      do l= ipar,npr,size
      if(x(l).ge.xmax-dx) then
        x(l)= x(l) -xmaxe       ! a bit smaller than xmax
      else if(x(l).le.-dx) then
        x(l)= x(l) +xmaxe
      end if
!*
      if(y(l).ge.ymax) then
        y(l)= 2.d0*ymax -y(l)   !<--fold the exact point
        vy(l)= -vy(l)
      else if(y(l).le.0.d0) then
        y(l)= -y(l) 
        vy(l)= -vy(l)
      end if
!*
      if(z(l).ge.zmax-dz) then
        z(l)= z(l) -zmaxe
      else if(z(l).le.-dz) then
        z(l)= z(l) +zmaxe
      end if
      end do
!* 
      return
      end subroutine partbc
!
!
!-----------------------------------------------------------------------
      subroutine partbcEST (x,y,z,npr,ipar,size) 
!-----------------------------------------------------------------------
!***********************************************************************
!*    Boundary conditions, estimated positions                         *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z
      integer(C_INT) npr,ipar,size
!------------------------------------------------------
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) l
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
!* The xyz coordinates 
      dx= hx/2 
      dz= hz/2
!
      do l= ipar,npr,size
      if(x(l).ge.xmax-dx) then
        x(l)= x(l) -xmaxe  
      else if(x(l).le.-dx) then
        x(l)= x(l) +xmaxe
      end if
!*
      if(y(l).ge.ymax) then
        y(l)= 2.d0*ymax -y(l) 
      else if(y(l).le.0.d0) then
        y(l)= -y(l) 
      end if
!*
      if(z(l).ge.zmax-dz) then
        z(l)= z(l) -zmaxe
      else if(z(l).le.-dz) then
        z(l)= z(l) +zmaxe
      end if
      end do
!* 
      return
      end subroutine partbcEST
!
!
!-----------------------------------------------------------------------
      subroutine partdk (x,y,z,mue,vpe,vhe,npr,ipar,size) 
!-----------------------------------------------------------------------
!*********************************************** updated: 4/23/92 ******
!*    Boundary conditions - drift kinetic electrons                    *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,mue,vpe,vhe
      integer(C_INT) npr,ipar,size
!------------------------------------------------------
!     real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
!                                             ex,ey,ez,bx,by,bz,       &
!                                             ex0,ey0,ez0,bx0,by0,bz0, &
!                                             bxa,bya,bza
!     common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      integer(C_INT) l
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
!* The xyz coordinates 
      dx= hx/2  
      dz= hz/2
!
!* ksp= 2
      do l= ipar,npr,size
      if(x(l).ge.xmax-dx) then
        x(l)= x(l) -xmaxe
!
      else if(x(l).le.-dx) then
        x(l)= x(l) +xmaxe
      end if
!*
!     epara= (aex*abx +aey*aby +aez*abz)/bsa
!     ach=  qw*epara -mue(l)*grha
!     fric= gnu(i,j,k)*hdt
!     vhh= ((1.d0-0.5d0*fric)*vhe(l) +hdt*ach)/(1.d0+0.5d0*fric)
!     vxe(l)= vhe(l)*bx1(i,j,k)/bsa !+( ey0(i,j,k)*bz1(i,j,k) &
                                      !  -ez0(i,j,k)*by1(i,j,k))/bsq2
      if(y(l).ge.ymax) then
        y(l)= 2.d0*ymax -y(l)   
        vhe(l)= -vhe(l)
!
      else if(y(l).le.0.d0) then
        y(l)= -y(l) 
        vhe(l)= -vhe(l)
      end if
!*
      if(z(l).ge.zmax-dz) then
        z(l)= z(l) -zmaxe
!
      else if(z(l).le.-dz) then
        z(l)= z(l) +zmaxe
      end if
      end do
! 
      return
      end subroutine partdk
!
!
!-----------------------------------------------------------------------
      subroutine partbcT (x,y,z,vx,vy,vz,npr) 
!-----------------------------------------------------------------------
!*********************************************** updated: 4/23/92 ******
!*    <<particle>> boundary conditions.                                *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      integer(C_INT) npr
!------------------------------------------------------
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) l
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
!* The xyz coordinates 
      dx= hx/2  
      dz= hz/2
!
      do l= 1,npr
      if(x(l).ge.xmax-dx) then
        x(l)= x(l) -xmaxe
!
      else if(x(l).le.-dx) then
        x(l)= x(l) +xmaxe
      end if
!*
      if(y(l).ge.ymax) then
        y(l)= 2.d0*ymax -y(l) 
        vy(l)= -vy(l)
!
      else if(y(l).le.0.d0) then
        y(l)= -y(l)
        vy(l)= -vy(l)
      end if
!*
      if(z(l).ge.zmax-dz) then
        z(l)= z(l) -zmaxe
!
      else if(z(l).le.-dz) then
        z(l)= z(l) +zmaxe
      end if
      end do
!
      return
      end subroutine partbcT
!
!
!-----------------------------------------------------------------------
      subroutine partdkT (x,y,z,mue,vpe,vhe,npr) 
!-----------------------------------------------------------------------
!*********************************************** updated: 4/23/92 ******
!*    <<particle>> boundary conditions.                                *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,mue,vpe,vhe
      integer(C_INT) npr
!------------------------------------------------------
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      integer(C_INT) l
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
!* The xyz coordinates 
      dx= hx/2
      dz= hz/2
!
      do l= 1,npr
      if(x(l).ge.xmax-dx) then
        x(l)= x(l) -xmaxe
!
      else if(x(l).le.-dx) then
        x(l)= x(l) +xmaxe
      end if
!*
      if(y(l).ge.ymax) then
        y(l)= 2.d0*ymax -y(l)   
        vhe(l)= -vhe(l)
!
      else if(y(l).le.0.d0) then
        y(l)= -y(l) 
        vhe(l)= -vhe(l)
      end if
!*
      if(z(l).ge.zmax-dz) then
        z(l)= z(l) -zmaxe
!
      else if(z(l).le.-dz) then
        z(l)= z(l) +zmaxe
      end if
!*
      end do
! 
      return
      end subroutine partdkT
!
!
!********** << subroutines 2 >> ****************************************
!*    ## 3-d macroscale E.M. particle code ##                          *
!*         << full-implicit scheme >>                                  *
!********************************************** update : 09/18/1998 ****
!-----------------------------------------------------------------------
      subroutine srimp1 (rxm,rym,rzm,vxj,vyj,vzj,qmult,               & 
                                 qjx,qjy,qjz,npr,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: rxm,rym,rzm,vxj,vyj,vzj
      real(C_DOUBLE) qmult
!                            3-point -> left,middle,right
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: qjx,qjy,qjz
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: qqjx,qqjy,qqjz
!
      integer(C_INT) npr,ipar,size
!--------------------------------------------
!
      integer(C_INT) ierror
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) i,j,k,ip,jp,kp,l,il,ir,jl,jr,kl,kr
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr,  &
                     qq,qgamx,qgamy,qgamz
!
      do k= -2,mz+1
      do j= -1,my+1   !<-- extend regions: (-1,my+1)
      do i= -2,mx+1
      qjx(i,j,k)= 0
      qjy(i,j,k)= 0
      qjz(i,j,k)= 0
!
      qqjx(i,j,k)= 0
      qqjy(i,j,k)= 0
      qqjz(i,j,k)= 0
      end do
      end do
      end do
!
!***
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,hxi,hyi,hzi,rxm,rym,rzm,    &
!$OMP        npr,hx,hy,hz,pxc,pxl,pxr,pzc,pzl,pzr, &
!$OMP        qmult,vxj,vyj,vzj)            ,       &
!$OMP PRIVATE(l,ip,jp,kp,i,il,ir,jl,jr,k,kl,kr,    &
!$OMP         xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
!$OMP         qq,qgamx,qgamy,qgamz)                &
!$OMP REDUCTION(+:qjx,qjy,qjz)
!$OMP DO SCHEDULE(STATIC,1)
!
      do l= ipar,npr,size
      ip= int(hxi*rxm(l) +0.500000001d0)
      jp= int(hyi*rym(l) +0.000000001d0)
      kp= int(hzi*rzm(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      if(jp.ge.my) then
        jr= my+1
        jl= my
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
      xx = hxi*rxm(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      fyl= hyi*rym(l) -jp
      fyr= 1.d0 -fyl
!
      zz = hzi*rzm(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      qq = qmult 
      qgamx = qq*vxj(l)
      qgamy = qq*vyj(l)
      qgamz = qq*vzj(l)
!
      qjx(il,jl,kl)= qjx(il,jl,kl) + qgamx*fxl*fyl*fzl
      qjx(il,jl,k )= qjx(il,jl,k ) + qgamx*fxl*fyl*fzc
      qjx(il,jl,kr)= qjx(il,jl,kr) + qgamx*fxl*fyl*fzr
      qjx(i ,jl,kl)= qjx(i ,jl,kl) + qgamx*fxc*fyl*fzl
      qjx(i ,jl,k )= qjx(i ,jl,k ) + qgamx*fxc*fyl*fzc
      qjx(i ,jl,kr)= qjx(i ,jl,kr) + qgamx*fxc*fyl*fzr
      qjx(ir,jl,kl)= qjx(ir,jl,kl) + qgamx*fxr*fyl*fzl
      qjx(ir,jl,k )= qjx(ir,jl,k ) + qgamx*fxr*fyl*fzc
      qjx(ir,jl,kr)= qjx(ir,jl,kr) + qgamx*fxr*fyl*fzr
!
      qjx(il,jr,kl)= qjx(il,jr,kl) + qgamx*fxl*fyr*fzl
      qjx(il,jr,k )= qjx(il,jr,k ) + qgamx*fxl*fyr*fzc
      qjx(il,jr,kr)= qjx(il,jr,kr) + qgamx*fxl*fyr*fzr
      qjx(i ,jr,kl)= qjx(i ,jr,kl) + qgamx*fxc*fyr*fzl
      qjx(i ,jr,k )= qjx(i ,jr,k ) + qgamx*fxc*fyr*fzc
      qjx(i ,jr,kr)= qjx(i ,jr,kr) + qgamx*fxc*fyr*fzr
      qjx(ir,jr,kl)= qjx(ir,jr,kl) + qgamx*fxr*fyr*fzl
      qjx(ir,jr,k )= qjx(ir,jr,k ) + qgamx*fxr*fyr*fzc
      qjx(ir,jr,kr)= qjx(ir,jr,kr) + qgamx*fxr*fyr*fzr
!
      qjy(il,jl,kl)= qjy(il,jl,kl) + qgamy*fxl*fyl*fzl
      qjy(il,jl,k )= qjy(il,jl,k ) + qgamy*fxl*fyl*fzc
      qjy(il,jl,kr)= qjy(il,jl,kr) + qgamy*fxl*fyl*fzr
      qjy(i ,jl,kl)= qjy(i ,jl,kl) + qgamy*fxc*fyl*fzl
      qjy(i ,jl,k )= qjy(i ,jl,k ) + qgamy*fxc*fyl*fzc
      qjy(i ,jl,kr)= qjy(i ,jl,kr) + qgamy*fxc*fyl*fzr
      qjy(ir,jl,kl)= qjy(ir,jl,kl) + qgamy*fxr*fyl*fzl
      qjy(ir,jl,k )= qjy(ir,jl,k ) + qgamy*fxr*fyl*fzc
      qjy(ir,jl,kr)= qjy(ir,jl,kr) + qgamy*fxr*fyl*fzr
!
      qjy(il,jr,kl)= qjy(il,jr,kl) + qgamy*fxl*fyr*fzl
      qjy(il,jr,k )= qjy(il,jr,k ) + qgamy*fxl*fyr*fzc
      qjy(il,jr,kr)= qjy(il,jr,kr) + qgamy*fxl*fyr*fzr
      qjy(i ,jr,kl)= qjy(i ,jr,kl) + qgamy*fxc*fyr*fzl
      qjy(i ,jr,k )= qjy(i ,jr,k ) + qgamy*fxc*fyr*fzc
      qjy(i ,jr,kr)= qjy(i ,jr,kr) + qgamy*fxc*fyr*fzr
      qjy(ir,jr,kl)= qjy(ir,jr,kl) + qgamy*fxr*fyr*fzl
      qjy(ir,jr,k )= qjy(ir,jr,k ) + qgamy*fxr*fyr*fzc
      qjy(ir,jr,kr)= qjy(ir,jr,kr) + qgamy*fxr*fyr*fzr
!
      qjz(il,jl,kl)= qjz(il,jl,kl) + qgamz*fxl*fyl*fzl
      qjz(il,jl,k )= qjz(il,jl,k ) + qgamz*fxl*fyl*fzc
      qjz(il,jl,kr)= qjz(il,jl,kr) + qgamz*fxl*fyl*fzr
      qjz(i ,jl,kl)= qjz(i ,jl,kl) + qgamz*fxc*fyl*fzl
      qjz(i ,jl,k )= qjz(i ,jl,k ) + qgamz*fxc*fyl*fzc
      qjz(i ,jl,kr)= qjz(i ,jl,kr) + qgamz*fxc*fyl*fzr
      qjz(ir,jl,kl)= qjz(ir,jl,kl) + qgamz*fxr*fyl*fzl
      qjz(ir,jl,k )= qjz(ir,jl,k ) + qgamz*fxr*fyl*fzc
      qjz(ir,jl,kr)= qjz(ir,jl,kr) + qgamz*fxr*fyl*fzr
!
      qjz(il,jr,kl)= qjz(il,jr,kl) + qgamz*fxl*fyr*fzl
      qjz(il,jr,k )= qjz(il,jr,k ) + qgamz*fxl*fyr*fzc
      qjz(il,jr,kr)= qjz(il,jr,kr) + qgamz*fxl*fyr*fzr
      qjz(i ,jr,kl)= qjz(i ,jr,kl) + qgamz*fxc*fyr*fzl
      qjz(i ,jr,k )= qjz(i ,jr,k ) + qgamz*fxc*fyr*fzc
      qjz(i ,jr,kr)= qjz(i ,jr,kr) + qgamz*fxc*fyr*fzr
      qjz(ir,jr,kl)= qjz(ir,jr,kl) + qgamz*fxr*fyr*fzl
      qjz(ir,jr,k )= qjz(ir,jr,k ) + qgamz*fxr*fyr*fzc
      qjz(ir,jr,kr)= qjz(ir,jr,kr) + qgamz*fxr*fyr*fzr
      end do
!$OMP END DO
!$OMP END PARALLEL
!
!  mxyzA= (mx+4)*(myA+2)*(mz+4)
      call mpi_allreduce (qjx(-2,-1,-2),qqjx(-2,-1,-2),mxyzA,mpi_real8, &
                          mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce (qjy(-2,-1,-2),qqjy(-2,-1,-2),mxyzA,mpi_real8, &
                          mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce (qjz(-2,-1,-2),qqjz(-2,-1,-2),mxyzA,mpi_real8, &
                          mpi_sum,mpi_comm_world,ierror)
!                                                +++++
!  All data are constracted of qjx(i,j,k) - qjz(i,j,k)
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
      qjx(i,j,k)= qqjx(i,j,k)
      qjy(i,j,k)= qqjy(i,j,k)
      qjz(i,j,k)= qqjz(i,j,k)
      end do
      end do
      end do
!
!***
      call vmesh3 (qjx,qjy,qjz)  !! 3-point mesh
!
      return
      end subroutine srimp1
!
!
!-----------------------------------------------------------------------
      subroutine srimp2 (rxm,rym,rzm,qmult,q,npr,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
!---------------------------------------------------------------
      real(C_DOUBLE),dimension(np0) :: rxm,rym,rzm
      real(C_DOUBLE) qmult
!                              ++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: q
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: qq
!
      integer(C_INT) npr,ipar,size
!---------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) i,j,k,ip,jp,kp,il,ir,jl,jr,kl,kr,l,ierror
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
                     qgam
!
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
       q(i,j,k)= 0
      qq(i,j,k)= 0
      end do
      end do
      end do
!
!***
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,hxi,hyi,hzi,rxm,rym,rzm,    &
!$OMP        npr,hx,hy,hz,pxc,pxl,pxr,pzc,pzl,pzr, &
!$OMP        pyc,qmult)                            &
!$OMP PRIVATE(l,ip,jp,kp,i,il,ir,jl,jr,k,kl,kr,    &
!$OMP        xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
!$OMP        qgam)                                &
!$OMP REDUCTION(+:q)
!$OMP DO SCHEDULE(STATIC,1)
!
      do l= ipar,npr,size
      ip= int(hxi*rxm(l) +0.500000001d0)
      jp= int(hyi*rym(l) +0.000000001d0)
      kp= int(hzi*rzm(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1
!
      if(jp.ge.my) then
        jr= my+1
        jl= my
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
!
      xx = hxi*rxm(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      fyl= hyi*rym(l) -jp
      fyr= 1.d0 -fyl
!
      zz = hzi*rzm(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      qgam = qmult  !<-- 
      q(il,jl,kl)= q(il,jl,kl) + qgam*fxl*fyl*fzl
      q(il,jl,k )= q(il,jl,k ) + qgam*fxl*fyl*fzc
      q(il,jl,kr)= q(il,jl,kr) + qgam*fxl*fyl*fzr
      q(i ,jl,kl)= q(i ,jl,kl) + qgam*fxc*fyl*fzl
      q(i ,jl,k )= q(i ,jl,k ) + qgam*fxc*fyl*fzc
      q(i ,jl,kr)= q(i ,jl,kr) + qgam*fxc*fyl*fzr
      q(ir,jl,kl)= q(ir,jl,kl) + qgam*fxr*fyl*fzl
      q(ir,jl,k )= q(ir,jl,k ) + qgam*fxr*fyl*fzc
      q(ir,jl,kr)= q(ir,jl,kr) + qgam*fxr*fyl*fzr
!
      q(il,jr,kl)= q(il,jr,kl) + qgam*fxl*fyr*fzl
      q(il,jr,k )= q(il,jr,k ) + qgam*fxl*fyr*fzc
      q(il,jr,kr)= q(il,jr,kr) + qgam*fxl*fyr*fzr
      q(i ,jr,kl)= q(i ,jr,kl) + qgam*fxc*fyr*fzl
      q(i ,jr,k )= q(i ,jr,k ) + qgam*fxc*fyr*fzc
      q(i ,jr,kr)= q(i ,jr,kr) + qgam*fxc*fyr*fzr
      q(ir,jr,kl)= q(ir,jr,kl) + qgam*fxr*fyr*fzl
      q(ir,jr,k )= q(ir,jr,k ) + qgam*fxr*fyr*fzc
      q(ir,jr,kr)= q(ir,jr,kr) + qgam*fxr*fyr*fzr
      end do
!$OMP END DO
!$OMP END PARALLEL
!
      call mpi_allreduce (q(-2,-1,-2),qq(-2,-1,-2),mxyzA,mpi_real8, &
                          mpi_sum,mpi_comm_world,ierror)
!                                                +++++
      do k= -2,mz+1
      do j= -1,my+1  !<-- -1,my+1
      do i= -2,mx+1
      q(i,j,k)= qq(i,j,k)
      end do
      end do
      end do
!
      call vmesh1 (q)
!
      return
      end subroutine srimp2
!
!
!-----------------------------------------------------------------------
      subroutine srimp3 (rxm,rym,rzm,mue,amu,npr,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
!--------------------------------------------------------------
      real(C_DOUBLE),dimension(np0) :: rxm,rym,rzm,mue
!                              ++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: amu
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: aamu
!
      integer(C_INT) npr,ipar,size
!--------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) i,j,k,ip,jp,kp,il,ir,jl,jr,kl,kr,l,ierror
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr,fmue
!
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
       amu(i,j,k)= 0
      aamu(i,j,k)= 0
      end do
      end do
      end do
!
!***
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,hxi,hyi,hzi,rxm,rym,rzm,    &
!$OMP        npr,hx,hy,hz,pxc,pxl,pxr,pzc,pzl,pzr, &
!$OMP        pyc,gx,gy,gz,mue)                     &
!$OMP PRIVATE(l,ip,jp,kp,i,il,ir,jl,jr,k,kl,kr,    &
!$OMP         xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
!$OMP         fmue)                                &
!$OMP REDUCTION(+:amu)
!$OMP DO SCHEDULE(STATIC,1)
!
      do l= ipar,npr,size
      ip= int(hxi*rxm(l) +0.500000001d0)
      jp= int(hyi*rym(l) +0.000000001d0)
      kp= int(hzi*rzm(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      if(jp.ge.my) then
        jr= my+1
        jl= my
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
!
      xx = hxi*rxm(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      fyl= hyi*rym(l) -jp
      fyr= 1.d0 -fyl
!
      zz = hzi*rzm(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      fmue = mue(l) 
      amu(il,jl,kl)= amu(il,jl,kl) + fmue*fxl*fyl*fzl
      amu(il,jl,k )= amu(il,jl,k ) + fmue*fxl*fyl*fzc
      amu(il,jl,kr)= amu(il,jl,kr) + fmue*fxl*fyl*fzr
      amu(i ,jl,kl)= amu(i ,jl,kl) + fmue*fxc*fyl*fzl
      amu(i ,jl,k )= amu(i ,jl,k ) + fmue*fxc*fyl*fzc
      amu(i ,jl,kr)= amu(i ,jl,kr) + fmue*fxc*fyl*fzr
      amu(ir,jl,kl)= amu(ir,jl,kl) + fmue*fxr*fyl*fzl
      amu(ir,jl,k )= amu(ir,jl,k ) + fmue*fxr*fyl*fzc
      amu(ir,jl,kr)= amu(ir,jl,kr) + fmue*fxr*fyl*fzr
!
      amu(il,jr,kl)= amu(il,jr,kl) + fmue*fxl*fyr*fzl
      amu(il,jr,k )= amu(il,jr,k ) + fmue*fxl*fyr*fzc
      amu(il,jr,kr)= amu(il,jr,kr) + fmue*fxl*fyr*fzr
      amu(i ,jr,kl)= amu(i ,jr,kl) + fmue*fxc*fyr*fzl
      amu(i ,jr,k )= amu(i ,jr,k ) + fmue*fxc*fyr*fzc
      amu(i ,jr,kr)= amu(i ,jr,kr) + fmue*fxc*fyr*fzr
      amu(ir,jr,kl)= amu(ir,jr,kl) + fmue*fxr*fyr*fzl
      amu(ir,jr,k )= amu(ir,jr,k ) + fmue*fxr*fyr*fzc
      amu(ir,jr,kr)= amu(ir,jr,kr) + fmue*fxr*fyr*fzr
      end do
!$OMP END DO
!$OMP END PARALLEL
!
      call mpi_allreduce (amu(-2,-1,-2),aamu(-2,-1,-2),mxyzA,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
!
      do k= -2,mz+1
      do j= -1,my+1  !<-- -1,my+1
      do i= -2,mx+1
      amu(i,j,k)= aamu(i,j,k)
      end do
      end do
      end do
!
      call vmesh1 (amu)
!
      return
      end subroutine srimp3
!
!
!-----------------------------------------------------------------------
      subroutine srimp3A (rxm,rym,rzm,vhe,avhh,npr,ipar,size)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
!--------------------------------------------------------------
      real(C_DOUBLE),dimension(np0) :: rxm,rym,rzm,vhe
!                              ++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: avhh
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: aavhh
!
      integer(C_INT) npr,ipar,size
!--------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) i,j,k,ip,jp,kp,il,ir,jl,jr,kl,kr,l,ierror
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr,fvhh
!
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
       avhh(i,j,k)= 0
      aavhh(i,j,k)= 0
      end do
      end do
      end do
!
!***
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,hxi,hyi,hzi,rxm,rym,rzm,    &
!$OMP        npr,hx,hy,hz,pxc,pxl,pxr,pzc,pzl,pzr, &
!$OMP        pyc,gx,gy,gz,vhe)                     &
!$OMP PRIVATE(l,ip,jp,kp,i,il,ir,jl,jr,k,kl,kr,    &
!$OMP         xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
!$OMP         fvhh)                                &
!$OMP REDUCTION(+:avhh)
!$OMP DO SCHEDULE(STATIC,1)
!
      do l= ipar,npr,size
      ip= int(hxi*rxm(l) +0.500000001d0)
      jp= int(hyi*rym(l) +0.000000001d0)
      kp= int(hzi*rzm(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      if(jp.ge.my) then
        jr= my+1
        jl= my
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
!
      xx = hxi*rxm(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      fyl= hyi*rym(l) -jp
      fyr= 1.d0 -fyl
!
      zz = hzi*rzm(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      fvhh = vhe(l)
      avhh(il,jl,kl)= avhh(il,jl,kl) + fvhh*fxl*fyl*fzl
      avhh(il,jl,k )= avhh(il,jl,k ) + fvhh*fxl*fyl*fzc
      avhh(il,jl,kr)= avhh(il,jl,kr) + fvhh*fxl*fyl*fzr
      avhh(i ,jl,kl)= avhh(i ,jl,kl) + fvhh*fxc*fyl*fzl
      avhh(i ,jl,k )= avhh(i ,jl,k ) + fvhh*fxc*fyl*fzc
      avhh(i ,jl,kr)= avhh(i ,jl,kr) + fvhh*fxc*fyl*fzr
      avhh(ir,jl,kl)= avhh(ir,jl,kl) + fvhh*fxr*fyl*fzl
      avhh(ir,jl,k )= avhh(ir,jl,k ) + fvhh*fxr*fyl*fzc
      avhh(ir,jl,kr)= avhh(ir,jl,kr) + fvhh*fxr*fyl*fzr
!
      avhh(il,jr,kl)= avhh(il,jr,kl) + fvhh*fxl*fyr*fzl
      avhh(il,jr,k )= avhh(il,jr,k ) + fvhh*fxl*fyr*fzc
      avhh(il,jr,kr)= avhh(il,jr,kr) + fvhh*fxl*fyr*fzr
      avhh(i ,jr,kl)= avhh(i ,jr,kl) + fvhh*fxc*fyr*fzl
      avhh(i ,jr,k )= avhh(i ,jr,k ) + fvhh*fxc*fyr*fzc
      avhh(i ,jr,kr)= avhh(i ,jr,kr) + fvhh*fxc*fyr*fzr
      avhh(ir,jr,kl)= avhh(ir,jr,kl) + fvhh*fxr*fyr*fzl
      avhh(ir,jr,k )= avhh(ir,jr,k ) + fvhh*fxr*fyr*fzc
      avhh(ir,jr,kr)= avhh(ir,jr,kr) + fvhh*fxr*fyr*fzr
      end do
!$OMP END DO
!$OMP END PARALLEL
!
      call mpi_allreduce (avhh(-2,-1,-2),aavhh(-2,-1,-2),mxyzA,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
!
      do k= -2,mz+1
      do j= -1,my+1  !<-- -1,my+1
      do i= -2,mx+1
      avhh(i,j,k)= aavhh(i,j,k)
      end do
      end do
      end do
!
      call vmesh1 (avhh)
!
      return
      end subroutine srimp3A
!
!
!-----------------------------------------------------------------------
      subroutine srimp4 (rxm,rym,rzm,mue,emx,emy,emz,npr,ipar,size) 
!-----------------------------------------------------------------------
!***********************************************************************
!*    calculate  sum mue*unit<b>  for the electrons.                   *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: rxm,rym,rzm,mue
!                              ++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: emx,emy,emz
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: eemx,eemy,eemz
!
      integer(C_INT) npr,ipar,size
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                                   ex,ey,ez,bx,by,bz,     &
                                                   ex0,ey0,ez0,bx0,by0,bz0
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
!--------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) i,j,k,ip,jp,kp,il,ir,jl,jr,kl,kr,l,ierror
      real(C_DOUBLE) xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
                     bxs,bys,bzs,rbs,wmue,emuex,emuey,emuez
!
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
      emx(i,j,k)= 0
      emy(i,j,k)= 0
      emz(i,j,k)= 0
!
      eemx(i,j,k)= 0
      eemy(i,j,k)= 0
      eemz(i,j,k)= 0
      end do
      end do
      end do
!                                  !<-- expand regions
      call outmesh3 (bx0,by0,bz0)
!
!***
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,hxi,hyi,hzi,rxm,rym,rzm,    &
!$OMP        npr,hx,hy,hz,pxc,pxl,pxr,pzc,pzl,pzr, &
!$OMP        pyc,gx,gy,gz,mue,bx0,by0,bz0)         &
!$OMP PRIVATE(l,ip,jp,kp,i,il,ir,jl,jr,k,kl,kr,    &
!$OMP         xx,fxl,fxc,fxr,fyr,fyl,zz,fzl,fzc,fzr, &
!$OMP         bxs,bys,bzs,rbs,wmue,emuex,emuey,emuez)&
!$OMP REDUCTION(+:emx,emy,emz)
!$OMP DO SCHEDULE(STATIC,1)
!
      do l= ipar,npr,size
      ip= int(hxi*rxm(l) +0.500000001d0)
      jp= int(hyi*rym(l) +0.000000001d0)
      kp= int(hzi*rzm(l) +0.500000001d0)
!
      il= ip-1 
      i = ip
      ir= ip+1
!
      jl= jp
      jr= jp+1 
!
      if(jp.ge.my) then
        jr= my+1
        jl= my
      else if(jp.lt.0) then  ! on 11/08
        jr=  0
        jl= -1
      end if
!
      kl= kp-1 
      k = kp
      kr= kp+1
!
!
      xx = hxi*rxm(l) -ip
      fxl= 0.5d0*(0.5d0-xx)*(0.5d0-xx)
      fxc= 0.75d0-xx*xx
      fxr= 0.5d0*(0.5d0+xx)*(0.5d0+xx)
!
      fyl= hyi*rym(l) -jp
      fyr= 1.d0 -fyl
!
      zz = hzi*rzm(l) -kp
      fzl= 0.5d0*(0.5d0-zz)*(0.5d0-zz)
      fzc= 0.75d0-zz*zz
      fzr= 0.5d0*(0.5d0+zz)*(0.5d0+zz)
!
      bxs = fyr*  &
           ((bx0(ir,jr,kr)*fxr+bx0(i,jr,kr)*fxc+bx0(il,jr,kr)*fxl)*fzr &
           +(bx0(ir,jr,k )*fxr+bx0(i,jr,k )*fxc+bx0(il,jr,k )*fxl)*fzc &
           +(bx0(ir,jr,kl)*fxr+bx0(i,jr,kl)*fxc+bx0(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((bx0(ir,jl,kr)*fxr+bx0(i,jl,kr)*fxc+bx0(il,jl,kr)*fxl)*fzr &
           +(bx0(ir,jl,k )*fxr+bx0(i,jl,k )*fxc+bx0(il,jl,k )*fxl)*fzc &
           +(bx0(ir,jl,kl)*fxr+bx0(i,jl,kl)*fxc+bx0(il,jl,kl)*fxl)*fzl)
!
      bys = fyr*  &
           ((by0(ir,jr,kr)*fxr+by0(i,jr,kr)*fxc+by0(il,jr,kr)*fxl)*fzr &
           +(by0(ir,jr,k )*fxr+by0(i,jr,k )*fxc+by0(il,jr,k )*fxl)*fzc &
           +(by0(ir,jr,kl)*fxr+by0(i,jr,kl)*fxc+by0(il,jr,kl)*fxl)*fzl)&
          + fyl* &
           ((by0(ir,jl,kr)*fxr+by0(i,jl,kr)*fxc+by0(il,jl,kr)*fxl)*fzr &
           +(by0(ir,jl,k )*fxr+by0(i,jl,k )*fxc+by0(il,jl,k )*fxl)*fzc &
           +(by0(ir,jl,kl)*fxr+by0(i,jl,kl)*fxc+by0(il,jl,kl)*fxl)*fzl)
!
      bzs = fyr*  &
           ((bz0(ir,jr,kr)*fxr+bz0(i,jr,kr)*fxc+bz0(il,jr,kr)*fxl)*fzr &
           +(bz0(ir,jr,k )*fxr+bz0(i,jr,k )*fxc+bz0(il,jr,k )*fxl)*fzc &
           +(bz0(ir,jr,kl)*fxr+bz0(i,jr,kl)*fxc+bz0(il,jr,kl)*fxl)*fzl)&
          + fyl*  &
           ((bz0(ir,jl,kr)*fxr+bz0(i,jl,kr)*fxc+bz0(il,jl,kr)*fxl)*fzr &
           +(bz0(ir,jl,k )*fxr+bz0(i,jl,k )*fxc+bz0(il,jl,k )*fxl)*fzc &
           +(bz0(ir,jl,kl)*fxr+bz0(i,jl,kl)*fxc+bz0(il,jl,kl)*fxl)*fzl)
      rbs= 1.d0/sqrt(bxs**2 +bys**2 +bzs**2)
!
      wmue = mue(l)  !<-- it was charm(l)= 1. in @ftpw8.f
      emuex = wmue*bxs*rbs
      emuey = wmue*bys*rbs
      emuez = wmue*bzs*rbs
!
      emx(il,jl,kl)= emx(il,jl,kl) + emuex*fxl*fyl*fzl
      emx(il,jl,k )= emx(il,jl,k ) + emuex*fxl*fyl*fzc
      emx(il,jl,kr)= emx(il,jl,kr) + emuex*fxl*fyl*fzr
      emx(i ,jl,kl)= emx(i ,jl,kl) + emuex*fxc*fyl*fzl
      emx(i ,jl,k )= emx(i ,jl,k ) + emuex*fxc*fyl*fzc
      emx(i ,jl,kr)= emx(i ,jl,kr) + emuex*fxc*fyl*fzr
      emx(ir,jl,kl)= emx(ir,jl,kl) + emuex*fxr*fyl*fzl
      emx(ir,jl,k )= emx(ir,jl,k ) + emuex*fxr*fyl*fzc
      emx(ir,jl,kr)= emx(ir,jl,kr) + emuex*fxr*fyl*fzr
!
      emx(il,jr,kl)= emx(il,jr,kl) + emuex*fxl*fyr*fzl
      emx(il,jr,k )= emx(il,jr,k ) + emuex*fxl*fyr*fzc
      emx(il,jr,kr)= emx(il,jr,kr) + emuex*fxl*fyr*fzr
      emx(i ,jr,kl)= emx(i ,jr,kl) + emuex*fxc*fyr*fzl
      emx(i ,jr,k )= emx(i ,jr,k ) + emuex*fxc*fyr*fzc
      emx(i ,jr,kr)= emx(i ,jr,kr) + emuex*fxc*fyr*fzr
      emx(ir,jr,kl)= emx(ir,jr,kl) + emuex*fxr*fyr*fzl
      emx(ir,jr,k )= emx(ir,jr,k ) + emuex*fxr*fyr*fzc
      emx(ir,jr,kr)= emx(ir,jr,kr) + emuex*fxr*fyr*fzr
!
      emy(il,jl,kl)= emy(il,jl,kl) + emuey*fxl*fyl*fzl
      emy(il,jl,k )= emy(il,jl,k ) + emuey*fxl*fyl*fzc
      emy(il,jl,kr)= emy(il,jl,kr) + emuey*fxl*fyl*fzr
      emy(i ,jl,kl)= emy(i ,jl,kl) + emuey*fxc*fyl*fzl
      emy(i ,jl,k )= emy(i ,jl,k ) + emuey*fxc*fyl*fzc
      emy(i ,jl,kr)= emy(i ,jl,kr) + emuey*fxc*fyl*fzr
      emy(ir,jl,kl)= emy(ir,jl,kl) + emuey*fxr*fyl*fzl
      emy(ir,jl,k )= emy(ir,jl,k ) + emuey*fxr*fyl*fzc
      emy(ir,jl,kr)= emy(ir,jl,kr) + emuey*fxr*fyl*fzr
!
      emy(il,jr,kl)= emy(il,jr,kl) + emuey*fxl*fyr*fzl
      emy(il,jr,k )= emy(il,jr,k ) + emuey*fxl*fyr*fzc
      emy(il,jr,kr)= emy(il,jr,kr) + emuey*fxl*fyr*fzr
      emy(i ,jr,kl)= emy(i ,jr,kl) + emuey*fxc*fyr*fzl
      emy(i ,jr,k )= emy(i ,jr,k ) + emuey*fxc*fyr*fzc
      emy(i ,jr,kr)= emy(i ,jr,kr) + emuey*fxc*fyr*fzr
      emy(ir,jr,kl)= emy(ir,jr,kl) + emuey*fxr*fyr*fzl
      emy(ir,jr,k )= emy(ir,jr,k ) + emuey*fxr*fyr*fzc
      emy(ir,jr,kr)= emy(ir,jr,kr) + emuey*fxr*fyr*fzr
!
      emz(il,jl,kl)= emz(il,jl,kl) + emuez*fxl*fyl*fzl
      emz(il,jl,k )= emz(il,jl,k ) + emuez*fxl*fyl*fzc
      emz(il,jl,kr)= emz(il,jl,kr) + emuez*fxl*fyl*fzr
      emz(i ,jl,kl)= emz(i ,jl,kl) + emuez*fxc*fyl*fzl
      emz(i ,jl,k )= emz(i ,jl,k ) + emuez*fxc*fyl*fzc
      emz(i ,jl,kr)= emz(i ,jl,kr) + emuez*fxc*fyl*fzr
      emz(ir,jl,kl)= emz(ir,jl,kl) + emuez*fxr*fyl*fzl
      emz(ir,jl,k )= emz(ir,jl,k ) + emuez*fxr*fyl*fzc
      emz(ir,jl,kr)= emz(ir,jl,kr) + emuez*fxr*fyl*fzr
!
      emz(il,jr,kl)= emz(il,jr,kl) + emuez*fxl*fyr*fzl
      emz(il,jr,k )= emz(il,jr,k ) + emuez*fxl*fyr*fzc
      emz(il,jr,kr)= emz(il,jr,kr) + emuez*fxl*fyr*fzr
      emz(i ,jr,kl)= emz(i ,jr,kl) + emuez*fxc*fyr*fzl
      emz(i ,jr,k )= emz(i ,jr,k ) + emuez*fxc*fyr*fzc
      emz(i ,jr,kr)= emz(i ,jr,kr) + emuez*fxc*fyr*fzr
      emz(ir,jr,kl)= emz(ir,jr,kl) + emuez*fxr*fyr*fzl
      emz(ir,jr,k )= emz(ir,jr,k ) + emuez*fxr*fyr*fzc
      emz(ir,jr,kr)= emz(ir,jr,kr) + emuez*fxr*fyr*fzr
      end do
!$OMP END DO
!$OMP END PARALLEL
!
      call mpi_allreduce (emx(-2,-1,-2),eemx(-2,-1,-2),mxyzA,mpi_real8,mpi_sum, & 
                          mpi_comm_world,ierror)
      call mpi_allreduce (emy(-2,-1,-2),eemy(-2,-1,-2),mxyzA,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      call mpi_allreduce (emz(-2,-1,-2),eemz(-2,-1,-2),mxyzA,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
!                                              ! +++++
      do k= -2,mz+1
      do j= -1,my+1 
      do i= -2,mx+1
      emx(i,j,k)= eemx(i,j,k)
      emy(i,j,k)= eemy(i,j,k)
      emz(i,j,k)= eemz(i,j,k)
      end do
      end do
      end do
!
      call vmesh3 (emx,emy,emz)  !! 3-point mesh
!
      return
      end subroutine srimp4
!
!
!-----------------------------------------------------------------------
      subroutine outmesh3 (ax,ay,az) 
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!                              +++++++++++++++++++++++ important !!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: ax,ay,az
      integer(C_INT) i,j,k
!***
!  Make (0:mx-1,1:my+1,0:mz-1) into expanded (-2:mx+1,-1:my+1,-2:mz+1)
!  1) x-direction: -1,0,...,mx+1,mx+2
!
      do i= -2,-1
      do j= 0,my   
      do k= 0,mz-1    ! Fold back position  
      ax(i,j,k)= ax(i+mx,j,k)  ! ax(-2)->ax(mx-2), ax(-1)->ax(mx-1)
      ay(i,j,k)= ay(i+mx,j,k)
      az(i,j,k)= az(i+mx,j,k)
      end do
      end do
      end do
!
      do i= mx,mx+1
      do j= 0,my
      do k= 0,mz-1
      ax(i,j,k)= ax(i-mx,j,k)  ! ax(mx)->ax(-1), ax(mx+1)->ax(0)
      ay(i,j,k)= ay(i-mx,j,k)
      az(i,j,k)= az(i-mx,j,k)
      end do
      end do
      end do
!
!  2) y-direction
!              -1,0,1,...,my,my+1
      j= -1
      do k= 0,mz-1
      do i= -2,mx+1
      ax(i,j,k)= 0  !<-- bound
      ay(i,j,k)= 0
      az(i,j,k)= 0
      end do
      end do
!
      j= my+1
      do k= 0,mz-1
      do i= -2,mx+1
      ax(i,j,k)= 0
      ay(i,j,k)= 0
      az(i,j,k)= 0
      end do
      end do
!
!  3) z-direction
!
      do k= -2,-1
      do j= -1,my+1
      do i= -2,mx+1
      ax(i,j,k)= ax(i,j,k+mz)  ! ax(-2)->ax(mz-2), ax(-1)->ax(mz-1)
      ay(i,j,k)= ay(i,j,k+mz)
      az(i,j,k)= az(i,j,k+mz)
      end do
      end do
      end do
!
      do k= mz,mz+1
      do j= -1,my+1
      do i= -2,mx+1
      ax(i,j,k)= ax(i,j,k-mz)  ! ax(mz), ax(mz+1)
      ay(i,j,k)= ay(i,j,k-mz)
      az(i,j,k)= az(i,j,k-mz)
      end do
      end do
      end do
!
      return
      end subroutine outmesh3
!
!
!-----------------------------------------------------------------------
      subroutine outmesh1 (ax)
!-----------------------------------------------------------------------
!***********************************************************************
!*    Adjust the charge density to include "uneven mesh" effect.       *
!***********************************************************************
!  Copy outmesh3 to outmesh1
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!                              +++++++++++++++++++++++ important !!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: ax
      integer(C_INT) i,j,k
!***
!  Make (0:mx-1,0:my,0:mz-1) into expanded (-2:mx+1,-1:my+1,-2:mz+1)
!  1) x-direction: -1,0,...,mx+1,mx+2
!
      do i= -2,-1
      do j= 0,my    !<--
      do k= 0,mz-1
      ax(i,j,k)= ax(i+mx,j,k)  ! ax(-2,,), ax(-1,,)
      end do
      end do
      end do
!
      do i= mx,mx+1
      do j= 0,my
      do k= 0,mz-1
      ax(i,j,k)= ax(i-mx,j,k)  ! ax(mx), ax(mx+1)
      end do
      end do
      end do
!
!  2) y-direction
!              -1,0,1,...,my,my+1
      j= -1
      do k= 0,mz-1
      do i= -2,mx+1
      ax(i,j,k)= 0  !<-- bound
      end do
      end do
!
      j= my+1
      do k= 0,mz-1
      do i= -2,mx+1
      ax(i,j,k)= 0
      end do
      end do
!
!  3) z-direction
!
      do k= -2,-1
      do j= -1,my+1
      do i= -2,mx+1
      ax(i,j,k)= ax(i,j,k+mz)  ! ax(-2,,), ax(-1,,)
      end do
      end do
      end do
!
      do k= mz,mz+1
      do j= -1,my+1
      do i= -2,mx+1
      ax(i,j,k)= ax(i,j,k-mz)  ! ax(mz), ax(mz+1)
      end do
      end do
      end do
!
      return
      end subroutine outmesh1
!
!
!-----------------------------------------------------------------------
      subroutine vmesh3 (ax,ay,az) 
!-----------------------------------------------------------------------
!  Move two meshes (either outside left or right) into inside 
!  either mesh= 0 or mx-1. Ori kaesi (x,z), ori tatami (y)
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!                              +++++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: ax,ay,az
      integer(C_INT) i,j,k
!***
!  Collapse the initial array ax(-2:mx+1,-1:my+1,-2:mz+1)  
!  into the interior array: ax(0:mx-1,0:my,0:mz-1)
!  1) x-direction:
!
      do i= -2,-1
      do j= -1,my+1 
      do k= -2,mz+1            !<-- Fold back
      ax(mx+i,j,k)= ax(i,j,k)  !<-- ax(-1) -> ax(mx-1)
      ay(mx+i,j,k)= ay(i,j,k)  !   to ax(i+mx)
      az(mx+i,j,k)= az(i,j,k)
      end do
      end do
      end do
!
      do i= mx,mx+1
      do j= -1,my+1
      do k= -2,mz+1
      ax(i-mx,j,k)= ax(i,j,k)  ! ax(mx) -> ax(0)
      ay(i-mx,j,k)= ay(i,j,k)
      az(i-mx,j,k)= az(i,j,k)
      end do
      end do
      end do
!
!  2) y-direction
!    Bound - fold back
!
      j= -1 
      do k= -2,mz+1
      do i= 0,mx-1  
      ax(i,0,k)= ax(i,0,k) +ax(i,-1,k) ! ax(0) +ax(-1) -> ax(0)
      ay(i,0,k)= ay(i,0,k) +ay(i,-1,k)
      az(i,0,k)= az(i,0,k) +az(i,-1,k)
      end do
      end do
!
      j= my+1
      do k= -2,mz+1
      do i= 0,mx-1
      ax(i,my,k)= ax(i,my,k) +ax(i,my+1,k) ! ax(my) +ax(my+1)
      ay(i,my,k)= ay(i,my,k) +ay(i,my+1,k) !  = ax(my)
      az(i,my,k)= az(i,my,k) +az(i,my+1,k)
      end do
      end do
!
!  3) z-direction
!   Periodic
!
      do k= -2,-1
      do j= 0,my   
      do i= 0,mx-1 
      ax(i,j,mz+k)= ax(i,j,k)  ! Fold back
      ay(i,j,mz+k)= ay(i,j,k)
      az(i,j,mz+k)= az(i,j,k)
      end do
      end do
      end do
!
      do k= mz,mz+1
      do j= 0,my
      do i= 0,mx-1
      ax(i,j,k-mz)= ax(i,j,k)  ! ax(mz) -> ax(0)
      ay(i,j,k-mz)= ay(i,j,k)
      az(i,j,k-mz)= az(i,j,k)
      end do
      end do
      end do
!
      return
      end subroutine vmesh3
!
!
!-----------------------------------------------------------------------
      subroutine vmesh1 (ax)
!-----------------------------------------------------------------------
!  Copy vmesh3 -> vmesh1
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!                              +++++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: ax
      integer(C_INT) i,j,k
!***
!  Collapse the initial array ax(-2:mx+1,-1:my+1,-2:mz+1)  
!  into the end array: ax(0:mx-1,0:my,0:mz-1)
!  1) x-direction
!
      do i= -2,-1
      do j= -1,my+1 
      do k= -2,mz+1            !<-- Fold back
      ax(mx+i,j,k)= ax(i,j,k)  !<-- ax(-1) -> ax(mx-1)
      end do
      end do
      end do
!
      do i= mx,mx+1
      do j= -1,my+1
      do k= -2,mz+1
      ax(i-mx,j,k)= ax(i,j,k)  ! ax(mx) -> ax(0)
      end do
      end do
      end do
!
!  2) y-direction
!    Bound - fold back
!
      j= -1 
      do k= -2,mz+1
      do i= 0,mx-1  
      ax(i,0,k)= ax(i,0,k) +ax(i,-1,k) ! ax(0) +ax(-1) -> ax(0)
      end do
      end do
!
      j= my+1
      do k= -2,mz+1
      do i= 0,mx-1
      ax(i,my,k)= ax(i,my,k) +ax(i,my+1,k) ! ax(my) +ax(my+1)
      end do                               !  = ax(my)
      end do
!
!  3) z-direction
!   Periodic
!
      do k= -2,-1
      do j= 0,my   
      do i= 0,mx-1 
      ax(i,j,mz+k)= ax(i,j,k)  ! Fold back
      end do
      end do
      end do
!
      do k= mz,mz+1
      do j= 0,my
      do i= 0,mx-1
      ax(i,j,k-mz)= ax(i,j,k)  ! ax(mz) -> ax(0)
      end do
      end do
      end do
!
      return
      end subroutine vmesh1
!
!
!-----------------------------------------------------------------------
      subroutine emfld0 
!-----------------------------------------------------------------------
!************************************************ version: 12/18/93 ****
!*    All ex - bz0 fields are first defined in this subroutine.        *
!*    Find a static b field by ampere's law.                           *
!***********************************************************************
!  This is used only at t= 0
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!*
      integer(C_INT) nobx
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                            ex,ey,ez,bx,by,bz,       &
                                            ex0,ey0,ez0,bx0,by0,bz0, &
                                            qix,qiy,qiz,qex,qey,qez, &
                                            qi,qe,emx,emy,emz,       &
                                            amu,avhh,pot
!                              ++++++++
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/sclpot/ pot
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                            cjx,cjy,cjz,sxi,syi,szi,qq
! --------------------------------------------------------------------
!***
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,   &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),            &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) vrg1,rwd2,anbar
      common/vring/  vrg1

      integer(C_INT) i,j,k,ir,il,kr,kl,ndim,iterp,syme,iwrt
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      character(len=8) poiss_sol
      common/poiss/ poiss_sol
!
!*---------------------------------------------------------------------
!***********************************************************************
!*    Calculate Laplacian_B = -curl(cj)/q0                             *
!***********************************************************************
!*  <<note>>  the current includes only explicit one at t= 0.
!*            no staggered grids.
!     ++++++++++
      nobx= nob3  !!<--used in /emcoef3/, /cresmd/ 
!     ++++++++++
!
!   Outside meshes calculated by outmesh3 () is then arranged by 
!   vmesh3 () in srimp1-srimp4 are used. So, (1,mx) are used.
! 
      do k= 0,mz-1
      do j= 0,my    !<-- inner region: 0,my
      do i= 0,mx-1
      cjx(i,j,k)= qix(i,j,k) +qex(i,j,k)
      cjy(i,j,k)= qiy(i,j,k) +qey(i,j,k)
      cjz(i,j,k)= qiz(i,j,k) +qez(i,j,k)
      end do
      end do
      end do
!
!-----------------------------------------------------------------------
!   Magnetic field by curl*B= J -> curl.curl B= curl.J
!* ## filter the sources on the field grids. ###
!-----------------------------------------------------------------------
      syme= -1
      call outmesh3 (cjx,cjy,cjz)  !<- vmesh3 already used in srimp1 !
      call filt3e (cjx,cjy,cjz,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,syme)
!
!
!*  Laplacian_B = -curl(cj)/q0 
!
      do k= 0,mz-1    !!<-- outmesh3 is used for cjx-cjz
      do j= 1,my-1    !!<-- meshes (0,my) are used hereafter !! 
      do i= 0,mx-1
!
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!
      sxi(i,j,k)= ( (cjy(i,j,kr)  -cjy(i,j,kl))/hz2      &
                   -(cjz(i,j+1,k) -cjz(i,j-1,k))/hy2)/q0
      syi(i,j,k)= ( (cjz(ir,j,k)  -cjz(il,j,k))/hx2      &
                   -(cjx(i,j,kr)  -cjx(i,j,kl))/hz2 )/q0
      szi(i,j,k)= ( (cjx(i,j+1,k) -cjx(i,j-1,k))/hy2     &
                   -(cjy(ir,j,k)  -cjy(il,j,k))/hx2 )/q0
      end do
      end do
      end do
!
      do k= 0,mz-1
      do i= 0,mx-1
      sxi(i,0,k)= ( (cjy(i,0,kr) -cjy(i,0,kl))/hz2 )/q0
      syi(i,0,k)= ( (cjz(ir,0,k) -cjz(il,0,k))/hx2      &
                   -(cjx(i,0,kr) -cjx(i,0,kl))/hz2 )/q0
      szi(i,0,k)= (-(cjy(ir,0,k) -cjy(il,0,k))/hx2 )/q0
!
      sxi(i,my,k)= ( (cjy(i,my,kr) -cjy(i,my,kl))/hz2 )/q0
      syi(i,my,k)= ( (cjz(ir,my,k) -cjz(il,my,k))/hx2      &
                    -(cjx(i,my,kr) -cjx(i,my,kl))/hz2 )/q0
      szi(i,my,k)= (-(cjy(ir,my,k) -cjy(il,my,k))/hx2 )/q0
      end do
      end do
!
      call outmesh3 (sxi,syi,szi)
!
!
      do k= -2,mz+1 
      do j= -1,my+1  
      do i= -2,mx+1
      bx(i,j,k)=  0
      by(i,j,k)=  0 
      bz(i,j,k)=  0 
      end do
      end do
      end do
!
!--------------------
      ndim= 3
!--------------------
      poiss_sol= 'poiss-bx'
      call poissn (bx,sxi,nobx,ndim,itersx,iterp)
!
      poiss_sol= 'poiss-by'
      call poissn (by,syi,nobx,ndim,itersx,iterp)
!
      poiss_sol= 'poiss-bz'
      call poissn (bz,szi,nobx,ndim,itersx,iterp)
!
!*  Boundary field fixing.
!
      do k= 0,mz-1
      do i= 0,mx-1 
      bx(i,0,k)= bx(i,1,k) 
      by(i,0,k)= 0        
      bz(i,0,k)= bz(i,1,k)
!
      bx(i,my,k)= bx(i,my-1,k)
      by(i,my,k)= 0
      bz(i,my,k)= bz(i,my-1,k)
      end do
      end do
!
!     if(io_pe.eq.1) then
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        i= mx/2
        k= mz/2
        write(11,*) 'Magnetic field at k=',k
!
        do j= 0,my
        if(mod(j,5).eq.1) then
        write(11,'("j,k=",2i4,1p3d12.3)') &
                                j,k,bx(i,j,k),by(i,j,k),bz(i,j,k)
        end if
        end do
!
        i= mx/2
        k= 10 +mz/2
        write(11,*) 'Magnetic field at k=',k
!
        do j= 0,my
        if(mod(j,5).eq.1) then
        write(11,'("j,k=",2i4,1p3d12.3)') &
                                j,k,bx(i,j,k),by(i,j,k),bz(i,j,k)
        end if
        end do
!
        close(11)
      end if 
!
! if a periodic case
!     do j= 0,my
!     do k= 0,mz-1
!     bx(0,j,k)= 0
!     by(0,j,k)= by(mx-1,j,k)
!     bz(0,j,k)= bz(mx-1,j,k)
!
!     bx(mx-1,j,k)= 0
!     by(mx-1,j,k)= by(mx-2,j,k)
!     bz(mx-1,j,k)= bz(mx-2,j,k)
!     end do
!     end do
!
!     do i= 0,mx-1
!     do j= 0,my
!     bx(i,j,0)= bx(i,j,mz-1)
!     by(i,j,0)= by(i,j,mz-1)
!     bz(i,j,0)= 0
!
!     bx(i,j,mz-1)= bx(i,j,0)
!     by(i,j,mz-1)= by(i,j,0)
!     bz(i,j,mz-1)= 0
!     end do
!     end do
!
!*  Adding the dc component.
!
!     do 260 k= 1,mz
!     do 260 i= 1,mx
!     bx(i,k)=  bx(i,k) +bxc
!     by(i,k)=  by(i,k) +byc
!     bz(i,k)=  bz(i,k) +bzc
! 260 continue
!
!***********************************************************************
!*  2(a):  Electric field (transverse).                                *
!***********************************************************************
!     << need to take the transverse part here >>
!
      do k= 0,mz-1
      do j= 0,my 
      do i= 0,mx-1
      ex(i,j,k)= 0
      ey(i,j,k)= 0
      ez(i,j,k)= 0
      end do
      end do
      end do
!
!***********************************************************************
!*  2(b):  Electric field (longitudinal).                              *
!***********************************************************************
!   /poissn/ is used at it= 0
!
      do k= 0,mz-1
      do j= 0,my 
      do i= 0,mx-1
      qq(i,j,k)= -( qi(i,j,k) +qe(i,j,k) )/q0
      pot(i,j,k)= 0
      end do
      end do
      end do
!
      call outmesh1 (qq)  !<-- expand to (-1:mx+2,...)
!
!--------------------
      ndim= 3
!--------------------
      poiss_sol= 'full-pot'
      call poissn (pot,qq,nobx,ndim,itersx,iterp)
!
      do k= 0,mz-1
      do j= 1,my-1 
      do i= 0,mx-1
      ex(i,j,k)= ex(i,j,k) -(pot(pxr(i),j,k) -pot(pxl(i),j,k))/hx2
      ey(i,j,k)= ey(i,j,k) -(pot(i,j+1,k)    -pot(i,j-1,k))/hy2
      ez(i,j,k)= ez(i,j,k) -(pot(i,j,pzr(k)) -pot(i,j,pzl(k)))/hz2
      end do
      end do
      end do
!
      do k= 0,mz-1
      do i= 0,mx-1
      ex(i,0,k)= 0
      ey(i,0,k)= -ey(i,1,k)
      ez(i,0,k)= 0
!
      ex(i,my,k)= 0
      ey(i,my,k)= -ey(i,my-1,k)
      ez(i,my,k)= 0
      end do
      end do
!
!***********************************************************************
!*  Initial fields on the particle grids.                            *
!***********************************************************************
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ex0(i,j,k)= ex(i,j,k)
      ey0(i,j,k)= ey(i,j,k)
      ez0(i,j,k)= ez(i,j,k)
      bx0(i,j,k)= bx(i,j,k)
      by0(i,j,k)= by(i,j,k)
      bz0(i,j,k)= bz(i,j,k)
      end do
      end do
      end do
!
      return
      end subroutine emfld0
!
!
!-----------------------------------------------------------------------
      subroutine emfild (np1,np2,nz1,nz2,ipar)
!-----------------------------------------------------------------------
!*   Full implicit ions and drift-kinetic electrons
!**** 11/10/89 ******************************************* 4/18/92 *****
!*                                                                     *
!*   This field-solver is designed for inhomogeneous and magnetized    *
!*   plasmas.  since the field is solved in the real space, almost     *
!*   any deformation of the plasmas is allowed.                        *
!*                                                                     *
!*   (ref.: M.Tanaka, J.Comput.Phys. 1992.)                            *
!********************************************************* 3/09/87 *****
!*                                                                     *
!*   The e*b, grad-b and curvature drift currents of the electrons     *
!*   are incorporated in the "coupled field-particle equation".        *
!*   this equation is free from the courant condition by virtue of     *
!*   its full-implicitness.                                            *
!*                                                                     *
!*   (Ref.:  J.Comput.Phys., vol.79, p.209 (1988))                     *
!***********************************************************************
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
      integer(C_INT) ipar
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                           ex,ey,ez,bx,by,bz,       &
                                           ex0,ey0,ez0,bx0,by0,bz0, &
                                           qix,qiy,qiz,qex,qey,qez, &
                                           emx,emy,emz,qi,qe,amu,avhh, &
                                           pot,gnu
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/sclpot/ pot
      common/dragcf/ gnu
!
      real(C_DOUBLE),dimension(0:my) :: filt 
      common/bpfilt/ filt
!                              ++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                           exa,eya,eza,bxa,bya,bza, &
                                           bss,                     &
                                           cjx,cjy,cjz,crx,cry,crz, &
                                           grh,grx,gry,grz,         &
                                           rhsx,rhsy,rhsz
!          +++++++
      real(C_float),dimension(0:mx-1,0:my,0:mz-1) :: &
                                           eex,eey,eez,bbx,bby,bbz
!------------------------------------------------------------------
      real(C_DOUBLE) tdec(3000)  !<-- DOUBLE
      common/ehist/  tdec
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,   &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),            &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) wkix,wkih,wkex,wkeh
      common/wkinel/ wkix,wkih,wkex,wkeh
!
      real(C_DOUBLE) ase,asb,asl,we,wb,wl,sbp2,sep2,wbp2,wep2
      integer(C_INT) iterm,iterf,iters
      common/emiter/ ase,asb,asl,we,wb,wl,iterm,iterf,iters
!
      integer(C_INT) i,j,k,kk,ir,il,kr,kl,syme,symb,            &
                     itag,iwrt,iperio 
      real(C_DOUBLE) hxhy4,hxhz4,hyhz4,dtic,dtic2,dtice,dtice2, &
                     rax,ray,raz,aex,aey,aez,grhx,grhy,grhz,    &
                     grbx,grby,grbz,bsq2,bsa1,                  &
                     ehh,ebx,eby,ebz,drag,mue1,vhh2,            &
                     cjy1,cjy2,sb2,se2,aldt2
!
      integer(C_INT),dimension(mxyz) :: arrayx,arrayy,arrayz
      common/array1d/ arrayx,arrayy,arrayz
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!*  itermx= 1 
! 
      go to 1000
!
!--------------------------
      entry prefld
!--------------------------
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      exa(i,j,k) = aimpl*ex(i,j,k) +(1.d0-aimpl)*ex0(i,j,k)
      eya(i,j,k) = aimpl*ey(i,j,k) +(1.d0-aimpl)*ey0(i,j,k)
      eza(i,j,k) = aimpl*ez(i,j,k) +(1.d0-aimpl)*ez0(i,j,k)
      end do
      end do
      end do
!
!
      do k= 0,mz-1
      do j= 1,my-1
      do i= 0,mx-1
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!
      bx(i,j,k)= bx0(i,j,k) +dt*( (eya(i,j,kr)  -eya(i,j,kl))/hz2  &
                                 -(eza(i,j+1,k) -eza(i,j-1,k))/hy2)
!
      by(i,j,k)= by0(i,j,k) +dt*( (eza(ir,j,k) -eza(il,j,k))/hx2   &
                                 -(exa(i,j,kr) -exa(i,j,kl))/hz2 )
!
      bz(i,j,k)= bz0(i,j,k) +dt*( (exa(i,j+1,k) -exa(i,j-1,k))/hy2 &
                                 -(eya(ir,j,k) -eya(il,j,k))/hx2 )
      end do
      end do
      end do
!
!
      do k= 0,mz-1
      do i= 0,mx-1
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!*
      bx(i,0,k)= bx0(i,0,k) +dt*( (eya(i,0,kr) -eya(i,0,kl))/hz2 &
                                -2*eza(i,1,k)/hy2) 
      by(i,0,k)= 0
      bz(i,0,k)= bz0(i,0,k) +dt*(-(eya(ir,0,k) -eya(il,0,k))/hx2 &
                                +2*exa(i,1,k)/hy2) 
!
      bx(i,my,k)= bx0(i,my,k) +dt*( (eya(i,my,kr) -eya(i,my,kl))/hz2 &
                                 +2*eza(i,my,k)/hy2)
      by(i,my,k)= 0
      bz(i,my,k)= bz0(i,my,k) +dt*(-(eya(ir,my,k) -eya(il,my,k))/hx2 &
                                 -2*exa(i,my,k)/hy2) 
      end do
      end do
!
!*  Anchor field (smooth).
!
!     syme= -1
!     call filt3e (ex,ey,ez,0.,0.,0.,ifilx,ifily,ifilz,syme)
!-------------------
      return
!-------------------
!
!***********************************************************************
!* 1. Accumulate the source moments in fulmov (ipc).                   *
!***********************************************************************
!-----------------------------------------------------------------------
!*  Define rhsx-rhsz for points of the inner domain.  
!-----------------------------------------------------------------------
!  Need rhsx-rhsz: expanded later
 1000 continue
!
      hxhy4= hx2*hy2  !<<- aldt2 description is hxhy4
      hxhz4= hx2*hz2
      hyhz4= hy2*hz2 
!
      aldt2= -aimpl*(1.d0-aimpl)*dt**2
!
      do k= 0,mz-1
      do j= 1,my-1    !<-- j=1,my-1 
      do i= 0,mx-1
! 
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!
      rhsx(i,j,k)= ex0(i,j,k)  &
             +aldt2*( -(ex0(i,j+1,k) +ex0(i,j-1,k) -2.d0*ex0(i,j,k))/hysq &
                      -(ex0(i,j,kr)  +ex0(i,j,kl)  -2.d0*ex0(i,j,k))/hzsq &
                    +(  ey0(ir,j+1,k) -ey0(ir,j-1,k)                    &
                      -(ey0(il,j+1,k) -ey0(il,j-1,k)) )/hxhy4           &
                    +(  ez0(ir,j,kr) -ez0(ir,j,kl)                      &
                      -(ez0(il,j,kr) -ez0(il,j,kl)) )/hxhz4 )           &
             +dt*( (bz0(i,j+1,k) -bz0(i,j-1,k))/hy2                     &
                  -(by0(i,j,kr)  -by0(i,j,kl))/hz2 )
!
      rhsy(i,j,k)= ey0(i,j,k)  &
             +aldt2*( ( ex0(ir,j+1,k) -ex0(ir,j-1,k)                    &
                      -(ex0(il,j+1,k) -ex0(il,j-1,k)) )/hxhy4           &
                      -(ey0(ir,j,k) +ey0(il,j,k) -2.d0*ey0(i,j,k))/hxsq &
                      -(ey0(i,j,kr) +ey0(i,j,kl) -2.d0*ey0(i,j,k))/hzsq &
                     +( ez0(i,j+1,kr) -ez0(i,j-1,kr)                    &
                       -ez0(i,j+1,kl) +ez0(i,j-1,kl) )/hyhz4 )          &
             +dt*( (bx0(i,j,kr) -bx0(i,j,kl))/hz2                       &
                  -(bz0(ir,j,k) -bz0(il,j,k))/hx2 )
!
      rhsz(i,j,k)= ez0(i,j,k)  &
             +aldt2*( ( ex0(ir,j,kr) -ex0(ir,j,kl)                      &
                      -(ex0(il,j,kr) -ex0(il,j,kl)) )/hxhz4             &
                     +( ey0(i,j+1,kr) -ey0(i,j-1,kr)                    &
                       -ey0(i,j+1,kl) +ey0(i,j-1,kl) )/hyhz4            &
                      -(ez0(ir,j,k)  +ez0(il,j,k)  -2.d0*ez0(i,j,k))/hxsq  &
                      -(ez0(i,j+1,k) +ez0(i,j-1,k) -2.d0*ez0(i,j,k))/hysq) &
            +dt*( (by0(ir,j,k)  -by0(il,j,k))/hx2                       &
                 -(bx0(i,j+1,k) -bx0(i,j-1,k))/hy2 )
      end do
      end do
      end do
!
!
      do k= 0,mz-1
      do i= 0,mx-1
      rhsx(i,0,k)= 0
      rhsy(i,0,k)= -rhsy(i,1,k)
      rhsz(i,0,k)= 0
!
      rhsx(i,my,k)= 0
      rhsy(i,my,k)= -rhsy(i,my-1,k)
      rhsz(i,my,k)= 0
      end do
      end do
!
!
!***********************************************************************
!* 2. The current density: jx(t+adt) +djx(b(n+aimpl)) -djx(b(n)).      *
!***********************************************************************
!-----------------------------------------------------------------------
!*  the field  ex= ex(n) for ipc=1, but ex= ex(n+1) for ipc> 1.
!-----------------------------------------------------------------------
!  Magnetic field
!   the future E field ex= ex0 (filtered at y boundaries).
!
      call outmesh3 (ex,ey,ez)
      call outmesh3 (bx,by,bz)
      call outmesh3 (bx0,by0,bz0)
!
!***
      itermx= 1
      do kk= 1,itermx
!
      do k= 0,mz-1
      do j= 0,my 
      do i= 0,mx-1
!
      bxa(i,j,k)= aimpl*bx(i,j,k) +(1.d0-aimpl)*bx0(i,j,k) +bxc
      bya(i,j,k)= aimpl*by(i,j,k) +(1.d0-aimpl)*by0(i,j,k) +byc
      bza(i,j,k)= aimpl*bz(i,j,k) +(1.d0-aimpl)*bz0(i,j,k) +bzc
!
      bsq2= bxa(i,j,k)**2 +bya(i,j,k)**2 +bza(i,j,k)**2
      bsa1=  sqrt(bsq2)
!
      bss(i,j,k)= bsa1
!                                           wce is included here
      rax= bxa(i,j,k)/bsa1
      ray= bya(i,j,k)/bsa1
      raz= bza(i,j,k)/bsa1
!
      dtic  = hdt*qwi*bsa1
      dtic2 = dtic**2
!
      dtice = hdt*qwe*bsa1
      dtice2= dtice**2
!
!
      aex= (1.d0-aimpl)*ex0(i,j,k) 
      aey= (1.d0-aimpl)*ey0(i,j,k)
      aez= (1.d0-aimpl)*ez0(i,j,k)
!
      ehh=  aex*rax +aey*ray +aez*raz
      ebx=  aey*raz -aez*ray
      eby=  aez*rax -aex*raz
      ebz=  aex*ray -aey*rax
!
!-----------------------------------------------------------------------
!*    qix, qex include the 0-th orbit current predicted with
!     e(n+aimpl) and b(n+aimpl) fields.
!-----------------------------------------------------------------------
!
      drag=  1.d0 +aimpl*gnu(i,j,k)*hdt
!
      if(j.eq.0 .or. j.eq.my) then
        cjx(i,j,k)= 0
        cjy(i,j,k)= 0
        cjz(i,j,k)= 0
      else
!
!     +++++++++++++++++
      cjx(i,j,k)= qix(i,j,k) +qex(i,j,k) &
            +qi(i,j,k)*adt*qwi*(aex +dtic2*ehh*rax  +dtic*ebx)  &
                                                 /(1.d0+dtic2)  &
            +qe(i,j,k)*adt*qwe*(aex +dtice2*ehh*rax +dtice*ebx) &
                                                 /(1.d0+dtice2)
      cjy(i,j,k)= qiy(i,j,k) +qey(i,j,k) &
            +qi(i,j,k)*adt*qwi*(aey +dtic2*ehh*ray  +dtic*eby)  &
                                                 /(1.d0+dtic2)  &
            +qe(i,j,k)*adt*qwe*(aey +dtice2*ehh*ray +dtice*eby) &
                                                 /(1.d0+dtice2)
      cjz(i,j,k)= qiz(i,j,k) +qez(i,j,k) &
            +qi(i,j,k)*adt*qwi*(aez +dtic2*ehh*raz  +dtic*ebz)  &
                                                 /(1.d0+dtic2)  &
            +qe(i,j,k)*adt*qwe*(aez +dtice2*ehh*raz +dtice*ebz) &
                                                 /(1.d0+dtice2)
      end if
      end do 
      end do
      end do
!
!
!-----------------------------------------------------------------------
!*  The magnetic terms and magnetization current (pseudo drift current).
!-----------------------------------------------------------------------
!
      do k= 0,mz-1
      do j= 1,my-1 
      do i= 0,mx-1
!
      bsq2= bxa(i,j,k)**2 +bya(i,j,k)**2 +bza(i,j,k)**2
      bsa1=  sqrt(bsq2)
!
      rax= bxa(i,j,k)/bsa1
      ray= bya(i,j,k)/bsa1
      raz= bza(i,j,k)/bsa1
!
      dtic  = hdt*qwi*bsa1
      dtic2 = dtic**2
!
      dtice = hdt*qwe*bsa1
      dtice2= dtice**2
!
      if(j.lt.0 .or. j.ge.my) then  ! 11/08
        grbx= 0
        grby= 0
        grbz= 0
      else
        grbx= (bss(ir,j,k) -bss(il,j,k))/hx2
        grby= (bss(i,j+1,k) -bss(i,j,k))/hy  ! j=0    -> bss.L 0
        grbz= (bss(i,j,kr) -bss(i,j,kl))/hz2 ! j=my-1 -> bss.R m
      end if 
!
!  mue/(qspec(2)*B) (b x dB/dx) ...grx(i,j,k)
      grx(i,j,k)= ray*grbz -raz*grby
      gry(i,j,k)= raz*grbx -rax*grbz
      grz(i,j,k)= rax*grby -ray*grbx
!
!  vhh**2*wspec(2)/(qspec(2)*B) curl_(b x dB/dx)
!  ach=  qw*epara -(mue(l)/wspec(2))*grha  !<<- parallel of dB/dx
      grh(i,j,k)=  rax*(bxa(ir,j,k) -bxa(il,j,k))/hx2  &
                 + ray*(bya(i,j+1,k)-bya(i,j,k)) /hy   &
                 + raz*(bza(i,j,kr) -bza(i,j,kl))/hz2
!
      grhx= grh(i,j,k)*rax
      grhy= grh(i,j,k)*ray
      grhz= grh(i,j,k)*raz
!
      crx(i,j,k)= ray*grhz -raz*grhy
      cry(i,j,k)= raz*grhx -rax*grhz
      crz(i,j,k)= rax*grhy -ray*grhx
!
      mue1= amu(i,j,k)/qspec(2)
      vhh2= avhh(i,j,k)**2*wspec(2)/qspec(2)
!
!     +++++++++++++++++
      cjx(i,j,k)= cjx(i,j,k) &
            +adt*qwi*(qiy(i,j,k)*bza(i,j,k) -qiz(i,j,k)*bya(i,j,k)) &
                                                 /(1.d0+dtic2) &
            +adt*qwi*dtic*( &
                  bxa(i,j,k)*(qiy(i,j,k)*ray +qiz(i,j,k)*raz)  &
                 -bya(i,j,k)* qix(i,j,k)*ray -bza(i,j,k)*qix(i,j,k)*raz) &
                                                 /(1.d0+dtic2) &
!
            +adt*qwe*(qey(i,j,k)*bza(i,j,k) -qez(i,j,k)*bya(i,j,k)) &
                                                 /(1.d0+dtice2) &
            +adt*qwe*dtice*( &
                  bxa(i,j,k)*(qey(i,j,k)*ray +qez(i,j,k)*raz)  &
                 -bya(i,j,k)* qex(i,j,k)*ray -bza(i,j,k)*qex(i,j,k)*raz) &
                                                 /(1.d0+dtice2) 
!
      cjy(i,j,k)= cjy(i,j,k) &
            +adt*qwi*(qiz(i,j,k)*bxa(i,j,k) -qix(i,j,k)*bza(i,j,k)) &
                                                 /(1.d0+dtic2) &
            +adt*qwi*dtic*( &
                  -bxa(i,j,k)*qiy(i,j,k)*rax -bza(i,j,k)*qiy(i,j,k)*raz &
                  +bya(i,j,k)*(qix(i,j,k)*rax +qiz(i,j,k)*raz)) &
                                                 /(1.d0+dtic2)  &
            +adt*qwe*(qez(i,j,k)*bxa(i,j,k) -qex(i,j,k)*bza(i,j,k)) &
                                                 /(1.d0+dtice2) &
            +adt*qwe*dtice*( &
                  -bxa(i,j,k)*qey(i,j,k)*rax -bza(i,j,k)*qey(i,j,k)*raz &
                  +bya(i,j,k)*(qex(i,j,k)*rax +qez(i,j,k)*raz)) &
                                                 /(1.d0+dtice2) 
!
      cjz(i,j,k)= cjz(i,j,k) &
            +adt*qwi*(qix(i,j,k)*bya(i,j,k) -qiy(i,j,k)*bxa(i,j,k)) &
                                                 /(1.d0+dtic2) &
            +adt*qwi*dtic*(  &
                 -bxa(i,j,k)*qiz(i,j,k)*rax -bya(i,j,k)*qiz(i,j,k)*ray &
                 +bza(i,j,k)*(qiy(i,j,k)*ray +qix(i,j,k)*rax)) &
                                                 /(1.d0+dtic2) &
!
            +adt*qwe*(qex(i,j,k)*bya(i,j,k) -qey(i,j,k)*bxa(i,j,k)) &
                                                 /(1.d0+dtice2) &
            +adt*qwe*dtice*(  &
                 -bxa(i,j,k)*qez(i,j,k)*rax -bya(i,j,k)*qez(i,j,k)*ray &
                 +bza(i,j,k)*(qey(i,j,k)*ray +qex(i,j,k)*rax)) &
                                                 /(1.d0+dtice2) 
      end do
      end do
      end do
!
!***********************************************************************
!* 3. The RHS of the field-particle coupled equation.                  *
!***********************************************************************
!-----------------------------------------------------------------------
!*  Filter the source terms.
!-----------------------------------------------------------------------
!*  preparation before filtering.
!*  boundary values used in cfpsol: j= 0 and my
!                                   +++++++++++
!***
      do k= 0,mz-1 
      do i= 0,mx-1
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!
      cjx(i,0,k)= 0
      cjx(i,1,k)= 0
      cjx(i,2,k)= 0
!
      cjy1= ( cjy(ir,0,k) +cjy(ir,1,k) +cjy(ir,2,k)   &
            + cjy(il,0,k) +cjy(il,1,k) +cjy(il,2,k)   &
            + cjy(i, 0,k) +cjy(i, 1,k) +cjy(i, 2,k)   &
            + cjy(i,0,kr) +cjy(i,1,kr) +cjy(i,2,kr)   &
            + cjy(i,0,kl) +cjy(i,1,kl) +cjy(i,2,kl) )/15.d0
      cjy(i,0,k)=  0
      cjy(i,1,k)=  cjy1
      cjy(i,2,k)= (cjy1 +cjy(i,3,k))/2.d0
!
      cjz(i,0,k)= 0
      cjz(i,1,k)= 0
      cjz(i,2,k)= 0
!
      cjx(i,my,  k)= 0
      cjx(i,my-1,k)= 0
      cjx(i,my-2,k)= 0
!
      cjy2= ( cjy(ir,my-2,k) +cjy(ir,my-1,k)  +cjy(ir,my,k) &
            + cjy(il,my-2,k) +cjy(il,my-1,k)  +cjy(il,my,k) &
            + cjy(i, my-2,k) +cjy(i, my-1,k)  +cjy(i, my,k) &
            + cjy(i, my-2,kr) +cjy(i,my-1,kr) +cjy(i,my,kr) &
            + cjy(i, my-2,kl) +cjy(i,my-1,kl) +cjy(i,my,kl) )/15.d0
      cjy(i,my,  k)= 0
      cjy(i,my-1,k)=  cjy2
      cjy(i,my-2,k)= (cjy2 +cjy(i,my-3,k))/2.d0
!
      cjz(i,my,  k)= 0
      cjz(i,my-1,k)= 0
      cjz(i,my-2,k)= 0
      end do
      end do
!***
!
      do k= 0,mz-1
      do j= 0,my 
      do i= 0,mx-1
      rhsx(i,j,k)= rhsx(i,j,k) -dt*cjx(i,j,k)/q0
      rhsy(i,j,k)= rhsy(i,j,k) -dt*cjy(i,j,k)/q0
      rhsz(i,j,k)= rhsz(i,j,k) -dt*cjz(i,j,k)/q0
      end do
      end do
      end do
!
      syme= -1
      call outmesh3 (rhsx,rhsy,rhsz)
      call filt3e (rhsx,rhsy,rhsz,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,syme)
!
!***********************************************************************
!* 4. Solve the field-particle equation for the "given" rhs.           *
!***********************************************************************
!  This fields are for 0,mx-1, 0,my, 0,mz-1
! 
      iterf= 150  !<<- iterf must be >= 1
!
!    +++++++++++++ R ++++++ Source ++++++++++++++++++++++++++++++++
      call cfpsol (ex,ey,ez,rhsx,rhsy,rhsz,np1,np2,nz1,nz2,ipar, &
                   iterm,iterf) 
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!*
      if(mod(it,5).eq.1) then
        syme= -1 
        call outmesh3 (ex,ey,ez) 
        call filt3e (ex,ey,ez,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,syme)
      end if
!                                  <-- ex,ey,ez are filtered in /emfild/ 
!                                  <-- ex0 were filtered in the last step
!
!***********************************************************************
!* 5. The magnetic field is calculated.                                *
!***********************************************************************
!   here, the future e field ex is filtered at y boundaries.
!  Ampere's law, with filters and aimpl factor
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      exa(i,j,k)= aimpl*ex(i,j,k) +(1.d0-aimpl)*ex0(i,j,k)
      eya(i,j,k)= aimpl*ey(i,j,k) +(1.d0-aimpl)*ey0(i,j,k)
      eza(i,j,k)= aimpl*ez(i,j,k) +(1.d0-aimpl)*ez0(i,j,k)
      end do            
      end do
      end do
!
      call outmesh3 (exa,eya,eza) 
!
      do k= 0,mz-1
      do j= 1,my-1 
      do i= 0,mx-1
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!
!   new        old 
      bx(i,j,k)= bx0(i,j,k)  &
                        +dt*( (eya(i,j,kr)  -eya(i,j,kl))/hz2   &
                             -(eza(i,j+1,k) -eza(i,j-1,k))/hy2 )
      by(i,j,k)= by0(i,j,k)  &
                        +dt*( (eza(ir,j,k) -eza(il,j,k))/hx2    &
                             -(exa(i,j,kr) -exa(i,j,kl))/hz2 )
      bz(i,j,k)= bz0(i,j,k)  & 
                        +dt*( (exa(i,j+1,k) -exa(i,j-1,k))/hy2  &
                             -(eya(ir,j,k)  -eya(il,j,k))/hx2)
      end do
      end do
      end do
!
!
      do k= 0,mz-1
      do i= 0,mx-1
      ir= pxr(i)
      il= pxl(i)
!
      kr= pzr(k)
      kl= pzl(k)
!*
      bx(i,0,k)= bx0(i,0,k) +dt*( (eya(i,0,kr) -eya(i,0,kl))/hz2 &
                                -2*eza(i,1,k)/hy2) 
      by(i,0,k)= 0
      bz(i,0,k)= bz0(i,0,k) +dt*(-(eya(ir,0,k) -eya(il,0,k))/hx2 &
                                +2*exa(i,1,k)/hy2) 
!
      bx(i,my,k)= bx0(i,my,k) +dt*( (eya(i,my,kr) -eya(i,my,kl))/hz2 &
                                 +2*eza(i,my,k)/hy2)
      by(i,my,k)= 0
      bz(i,my,k)= bz0(i,my,k) +dt*(-(eya(ir,my,k) -eya(il,my,k))/hx2 &
                                 -2*exa(i,my,k)/hy2) 
      end do
      end do
!
!  Entended mesh points for particles 
      if(mod(it,5).eq.1) then
        symb= +1 
        call outmesh3 (bx,by,bz)
        call filt3e (bx,by,bz,0.d0,0.d0,0.d0,ifilx,ifily,ifilz,symb)
      end if
!                                  <-- bx,by,bz are filtered in /emfild/ 
!                                  <-- bx0 were filtered in the last step
!
!***********************************************************************
!* 6. Accuracy check of the magnetic iteration  (iterm loop).          *
!***********************************************************************
!
      sb2= 0
      se2= 0
      sbp2= 0
      sep2= 0
!
      do k= 0,mz-1
      do j= 0,my 
      do i= 0,mx-1
      sb2= sb2 +bx(i,j,k)**2 +by(i,j,k)**2 +bz(i,j,k)**2
      se2= se2 +ex(i,j,k)**2 +ey(i,j,k)**2 +ez(i,j,k)**2
!
      sbp2= sbp2 +by(i,j,k)**2 +bz(i,j,k)**2
      sep2= sep2 +ey(i,j,k)**2 +ez(i,j,k)**2
      end do
      end do
      end do
!
!-----------------------------------------------------------------------
!*    wb and we: convergence index.
!-----------------------------------------------------------------------
!
      we=    se2/float(3*mxyz)
      wb=    sb2/float(3*mxyz)
      wep2= sep2/float(2*mxyz)
      wbp2= sbp2/float(2*mxyz)
      asb= 0 
      ase= 0
!
!        +++++++++++++++++
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
        edec(ldec,1)=  wb
        edec(ldec,2)=  we
        edec(ldec,3)=  wbp2
        edec(ldec,4)=  wep2
!
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
        write(11,*) 'we, wb=',we,wb 
        close(11)
      end if
!
!***
      end do
!
!*  End of itermx loop
!***
!
      if(io_pe.eq.1) then
        if(iwrt(it,nplot).eq.0) then
!                  +++++
!       if(mod(it,10).eq.1) then
!
        call lblbot (t)
!
        do k= 0,mz-1
        do j= 0,my   
        do i= 0,mx-1
        eex(i,j,k)= ex(i,j,k)
        eey(i,j,k)= ey(i,j,k)
        eez(i,j,k)= ez(i,j,k)
        end do
        end do
        end do
!
!   (ex,ey,ez) real*4
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        iperio= 0
        itag= 4
        call fplot3 (eex,eey,eez,0.,0.,0.,real(xmax),real(ymax),real(zmax), &
                     iperio,itag,'E cfpsol',8)
!
!
        do k= 0,mz-1
        do j= 0,my  
        do i= 0,mx-1
!       bbx(i,j,k)= qix(i,j,k) +qex(i,j,k)
!
        bbx(i,j,k)= bx(i,j,k)
        bby(i,j,k)= by(i,j,k)
        bbz(i,j,k)= bz(i,j,k)
        end do
        end do
        end do
!
!   (bx,by,bz) real*4
        iperio= 0
        itag= 5
        call fplot3 (bbx,bby,bbz,0.,0.,0.,real(xmax),real(ymax),real(zmax), &
                     iperio,itag,'bx by bz',8)
        close(77)
!
        end if
      end if
!
      if(.false.) then
!     if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        i= mx/2
        j= my/2
        do k= 0,mz-1
        write(11,'(3i4,1p6d12.3)') i,j,k,ex(i,j,k),ey(i,j,k),ez(i,j,k), &
                      bx(i,j,k),by(i,j,k),bz(i,j,k)
        end do
!
        close(11)
      end if
!
      return
      end subroutine emfild
!
!
!********** << subroutines >> ******************************************
!*                                                                     *
!*     cfpsol ..... for the electromagnetic field.                     *
!*       /cfpsol/ -> /bcgstb/ -> /wwstbi/ etc.                         *
!*                                                                     *
!*** Original: 4/23/1992 ******************* Fortran 2003: 12/21/2021 **
!-----------------------------------------------------------------------
      subroutine cfpsol (ex,ey,ez,rhsx,rhsy,rhsz,np1,np2,nz1,nz2,ipar, &
                         iterm,iterf)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
      include 'mpif.h'
!
!     parameter  (nob=15,iblk=3)
!----------------------------------------------------------------------
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) ::            &
                                               ex,ey,ez,rhsx,rhsy,rhsz
!
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
      integer(C_INT) ipar,iterm,iterf,ierr
!
      integer(C_INT) ipr(10),ijk
      real(C_DOUBLE) rpr(10),rsdl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  BCG subroutines
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa         
      integer(C_INT),dimension(1:mxyz3,nob),save :: ja,na
!
      real(C_DOUBLE),dimension(1:mxyz3) :: ss,xx
      real(C_DOUBLE),dimension(1:mxyz3) :: zz  
!
      integer(C_INT),dimension(mxyz) :: arrayx,arrayy,arrayz
      common/array1d/ arrayx,arrayy,arrayz
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t8,dt,aimpl,adt,hdt,ahdt2,adtsq,     & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t8,dt,aimpl,adt,hdt,ahdt2,adtsq,     &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      integer(C_INT) i,j,k,itrmx,iwrt,ierror
      real(C_DOUBLE) wsq
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!*---------------------------------------------------------------------
!  So far single calculations in /emcoef/, /bcgstb/ 
!
      call emcoef (aa,na)
!                     ++
!
      do ijk= 1,mxyz  
      i= arrayx(ijk)
      j= arrayy(ijk)
      k= arrayz(ijk)
! 
      xx(3*ijk-2) =  ex(i,j,k)       !! (xx(3*ijk-2),,)= 
      xx(3*ijk-1) =  ey(i,j,k)       !!        (ex(i,j,k),,)
      xx(3*ijk  ) =  ez(i,j,k)
!
      ss(3*ijk-2) = rhsx(i,j,k)
      ss(3*ijk-1) = rhsy(i,j,k)
      ss(3*ijk  ) = rhsz(i,j,k)
      end do
!
!***
!  Bound in y: cfpsol, escorr, poissn 
!
      do ijk= 1,mxyz 
      j= arrayy(ijk)
! 
      if(j.eq.0 .or. j.eq.my) then 
        ss(3*ijk-2)= 0
        ss(3*ijk-1)= 0
        ss(3*ijk  )= 0
      end if
      end do
!
!-----------------------------------------------------------------------
!  The solution of the 3D CFP equations.
!-----------------------------------------------------------------------
!* 1)  choose /bcgstb/ when the core 3*3 matrix is off-diagonal dominant
!     (if no convergence with /bcgsts/.)              date: 9/15/92.
!  2)  the field boundary values are treated in the matrix calculation.
!
!***
!                                    set the parameters
      ipr(1) = iterfx  !<-- in /cfpsol/ 
      ipr(3) = 0       !<-- iprc=0, no preconditioning
      ipr(4) = 3       !<-- iblk
      rpr(1) = 1.d-5   !<-- eps
!***
!  zz() is used in mpi_allgather
!                           s  r 
      call bcgstb (aa,ja,na,ss,xx,np1,np2,nz1,nz2,ipar, &
                   ipr,rpr,ierr)
!
!     integer(C_INT) la,ma,n
!     parameter (la=mxyz3,ma=nob,n=mxyz3)
!     call bcgstb_sgl (aa,la,n,ma,jm,ja,na,ss,xx,  &
!                      ipr,rpr,iw,w,wp,ierr)
!
!  After bcgstb
!     rpr(2)  rsdl
!     ipr(2)  iterf 
!----------------------------------------------------
!*  write back to ex, ey and ez.       3/07/1999
!----------------------------------------------------
!   First address of (ex,ey,ez) : 1+(k-1)*mxyz3/npc
!    We only need mx*myA*mz/npc, but the full index is kept !
!
!                            +++++++++
      call mpi_allgather (xx(np1(ipar)),mxyz3/npc,mpi_real8,  & !! xx()
                          zz           ,mxyz3/npc,mpi_real8,  & !! zz(1) 
                          mpi_comm_world,ierror) 
!     -----------------------------------------------------------
!
      do ijk= 1,mxyz 
      i= arrayx(ijk)
      j= arrayy(ijk)
      k= arrayz(ijk)
!
      if(j.eq.0 .or. j.eq.my) then
        ex(i,j,k)= 0
        ey(i,j,k)= 0
        ez(i,j,k)= 0
      else                        !!--> parallel, zz( )
        ex(i,j,k)= zz(3*ijk-2)        ! zz(3*ijk-2)>= 1
        ey(i,j,k)= zz(3*ijk-1)        !    3*ijk-1 >= 2
        ez(i,j,k)= zz(3*ijk  )        !    3*ijk   >= 3
      end if
      end do
!*
!
      wsq= 0
!
      do ijk= 1,mxyz 
      i= arrayx(ijk)
      j= arrayy(ijk)
      k= arrayz(ijk)
!
      wsq= wsq +ex(i,j,k)**2 +ey(i,j,k)**2 +ez(i,j,k)**2
      end do
!
      wsq= sqrt(wsq/float(3*mxyz))
!
!     if(ierr.ne.0 .and. io_pe.eq.1) then
      if(.false.) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
!                     iterf  cond  sec   rsdl
        iterf= ipr(2)
        rsdl = rpr(2)
        write(11,'("#(bcgstb-em) it=",i5," iterf,ierr=",i3,i5,    &
               " (",f6.3," sec); rsdl=",1pd13.5,"; <e>=",d13.5)') &
                                     it,iterf,ierr,rpr(3),rsdl,wsq
        close(11)
      end if
!
      return
      end subroutine cfpsol
!
!
!----------------------------------------------------------------------
      subroutine emcoef (aa,na) 
!----------------------------------------------------------------------
!*   /emcoef/ in /cfpsol/
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
!***********************************************************************
!*    Define the whole CFP matrix.                                     *
!***********************************************************************
!*    aa (component, grid number, band number)
!          (x,y,z)    1 - mxyz     1 - nob
!
!     parameter  (nob=15)
! ---------------------------------------------------------------------
!  3 arguments
      real(C_DOUBLE),dimension(3,1:mxyz,nob) :: aa 
      integer(C_INT),dimension(3,1:mxyz,nob) :: na 
! 
      real(C_DOUBLE),dimension(3,nob) :: ca
      integer(C_INT),dimension(3,nob) :: lai,laj,lak,loc
      integer(C_INT) ipar,lxyz
! ---------------------------------------------------------------------
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                          ex,ey,ez,bx,by,bz,       &
                                          ex0,ey0,ez0,bx0,by0,bz0, &
                                          qix,qiy,qiz,qex,qey,qez, &
                                          emx,emy,emz,qi,qe,amu,avhh, &
                                          gnu
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/dragcf/ gnu
!------------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      real(C_DOUBLE),dimension(0:my) :: filt 
      common/bpfilt/ filt
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      integer(C_INT) i,j,k,ii,jj,kk,m
      real(C_DOUBLE) hxhy4,hxhz4,hyhz4,bxa,bya,bza,bss,     &
                     dtic,dtic2,fmain,fpare,fexbe,          &
                     dtice,dtice2,fmaine,rax,ray,raz,       &
                     a11,a22,a33,a12,a13,a21,a23,a31,a32
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!-----------------------------------------------------------------------
!***********************************************************************
!*  (1) The curl*curl terms.                                           *
!***********************************************************************
!*  The coefficients must be given so that na becomes an acending
!*    order in terms of index m  (error in /wwstbc/ of /bcgstb/).
!
!!    if(first) then
!       first= .true.  !! .false.
!
!   The boundary conditions are given in subroutines emcoef, escsol, emcof3
!***********************************************
!   Start with 1 (one) for xx(), ss() and zz()  <--emcoef
!
      do k= 0,mz-1
      do j= 0,my 
!
! ------- boundary points ------------
! (a) j= 0
      if(j.eq.0) then 
!***
        do i= 0,mx-1 
!* side 1:                 for ex(i,0,k): even to the left
!          l m 
        ca(1,1) =  1.d0
        ca(1,2) = -1.d0   ! even for ex(i,0,k) -ex(i,1,k)= 0
!*  
        ca(2,1) =  1.d0
        ca(2,2) =  1.d0   ! odd for ey(i,0,k) +ey(i,1,k)= 0
!*  
        ca(3,1) =  1.d0
        ca(3,2) = -1.d0   ! even for ez(i,0,k) -ez(i,1,k)= 0
! 
!         if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
!         open (unit=11,file=praefixc//'.11'//suffix2,             & 
!               status='unknown',position='append',form='formatted')
!         if(i.eq.0 .and. j.eq.0 .and. k.eq.0) then
!         write(11,941) i,j,k,lxyz+1,aa(1,lxyz+1,1),na(1,lxyz+1,1)
! 941     format('941: 1,1) i,j,k,lxyz+1,aa,na=',3i3,i8,d12.3,i10)
!         end if
!*
!  The number of LHS terms [continuous, (mxyz3,nob)]
!  The order is: (1,lxyz+1,1), (2,lxyz+2,1), (3,lxyz+3,1),
!                (1,lxyz+1,2), (2,lxyz+2,2), (3,lxyz+3,2),... 
! 
        lxyz = i +mx*(0 +myA*k)
! 
        do m= 1,nob
        if(m.eq.1) then
        aa(1,lxyz+1,1) = ca(1,1)  ! aa( ,m=1) 
        na(1,lxyz+1,1) = 3*(i +mx*(0 +myA*k)) +1
!
        else if(m.eq.2) then
        aa(1,lxyz+1,2) = ca(1,2)  ! aa( ,m=2)
        na(1,lxyz+1,2) = 3*(i +mx*(1 +myA*k)) +1 
!
        else
          aa(1,lxyz+1,m) = 0
          na(1,lxyz+1,m) = 3*(i +mx*(0 +myA*k)) +1 
        end if
        end do
!
        do m= 1,nob
        if(m.eq.1) then
        aa(2,lxyz+1,1) = ca(2,1) 
        na(2,lxyz+1,1) = 3*(i +mx*(0 +myA*k)) +2
!
        else if(m.eq.2) then
        aa(2,lxyz+1,2) = ca(2,2) 
        na(2,lxyz+1,2) = 3*(i +mx*(1 +myA*k)) +2
!
        else
          aa(2,lxyz+1,m) = 0
          na(2,lxyz+1,m) = 3*(i +mx*(0 +myA*k)) +2
        end if
        end do
!
        do m= 1,nob
        if(m.eq.1) then
        aa(3,lxyz+1,1) = ca(3,1) 
        na(3,lxyz+1,1) = 3*(i +mx*(0 +myA*k)) +3
!
        else if(m.eq.2) then
        aa(3,lxyz+1,2) = ca(3,2) 
        na(3,lxyz+1,2) = 3*(i +mx*(1 +myA*k)) +3
!
        else
          aa(3,lxyz+1,m)= 0
          na(3,lxyz+1,m)= 3*(i +mx*(0 +myA*k)) +3
        end if
        end do
        end do 
!
!***
! (b) 0 < j < my
!
      else if(j.gt.0 .and. j.lt.my) then
!
      hxhy4= hx2*hy2  !<-- hy2
      hxhz4= hx2*hz2
      hyhz4= hy2*hz2 
!
      do i= 0,mx-1  
!
!*  lai,laj,lak and loc,ca; there are 15 non-zero terms 
!
!-- sx(i,k) --
!                                for (i-1,j,k-1)
      lai(1,1) = i-1
      laj(1,1) = j
      lak(1,1) = k-1
      loc(1,1) = 3
       ca(1,1) = adtsq/hxhz4
!                      ..... 1
!                                for (i,j,k-1)
      lai(1,2) = i
      laj(1,2) = j
      lak(1,2) = k-1
      loc(1,2) = 1
       ca(1,2) = -adtsq/hzsq
!                       .... 2
!                         
      lai(1,3) = i+1
      laj(1,3) = j
      lak(1,3) = k-1
      loc(1,3) = 3
       ca(1,3) = -adtsq/hxhz4
!                       ..... 3
!                                for (i-1,j-1,k)
      lai(1,4) = i-1
      laj(1,4) = j-1
      lak(1,4) = k
      loc(1,4) = 2
       ca(1,4) =  adtsq/hxhy4
!                       ..... 4
!                                for (i,j-1,k)
      lai(1,5) = i
      laj(1,5) = j-1
      lak(1,5) = k
      loc(1,5) = 1
       ca(1,5) = -adtsq/hysq 
!                       .... 6
!                     
!                                for (i+1,j-1,k)
      lai(1,6) = i+1
      laj(1,6) = j-1
      lak(1,6) = k
      loc(1,6) = 2
       ca(1,6) = -adtsq/hxhy4
!                       ..... 11
!*------------------------------ for (i,j,k)... diagonal -------------
      lai(1,7) = i
      laj(1,7) = j
      lak(1,7) = k
      loc(1,7) = 1
       ca(1,7) =  1.d0 + adtsq*(2.d0/hysq +2.d0/hzsq)
!                  .... 7,8,9 
      lai(1,8) = i
      laj(1,8) = j
      lak(1,8) = k
      loc(1,8) = 2
       ca(1,8) = 0.d0
!
      lai(1,9) = i
      laj(1,9) = j
      lak(1,9) = k
      loc(1,9) = 3
       ca(1,9) = 0.d0
!*--------------------------------------------------------------------
!                                for (i-1,j+1,k)
      lai(1,10) = i-1
      laj(1,10) = j+1
      lak(1,10) = k
      loc(1,10) = 2
       ca(1,10) = -adtsq/hxhy4
!                        ..... 5
!                                for (i,j+1,k)
      lai(1,11) = i
      laj(1,11) = j+1
      lak(1,11) = k
      loc(1,11) = 1
       ca(1,11) = -adtsq/hysq
!                        .... 10
!                                for (i+1,j+1,k)
      lai(1,12) = i+1
      laj(1,12) = j+1
      lak(1,12) = k
      loc(1,12) = 2
       ca(1,12) = adtsq/hxhy4
!                       ..... 12
!                                for (i-1,j,k+1)
      lai(1,13) = i-1
      laj(1,13) = j
      lak(1,13) = k+1
      loc(1,13) = 3
       ca(1,13) = -adtsq/hxhz4
!                        ..... 13
!                                for (i,j,k+1)
      lai(1,14) = i
      laj(1,14) = j
      lak(1,14) = k+1
      loc(1,14) = 1
       ca(1,14) = -adtsq/hzsq
!                        .... 14
!                                for (i+1,j,k+1)
      lai(1,15) = i+1
      laj(1,15) = j
      lak(1,15) = k+1
      loc(1,15) = 3
       ca(1,15) = adtsq/hxhz4
!                       ..... 15
!-- sy(k) --
!                             
!                                for (i,j-1,k-1)
      lai(2,1) = i
      laj(2,1) = j-1 
      lak(2,1) = k-1
      loc(2,1) = 3
       ca(2,1) = adtsq/hyhz4
!                      ..... 1
!                                for (i,j,k-1)
      lai(2,2) = i
      laj(2,2) = j
      lak(2,2) = k-1
      loc(2,2) = 2
       ca(2,2) = -adtsq/hzsq
!                       .... 2
!                                for (i,j+1,k-1)
      lai(2,3) = i
      laj(2,3) = j+1
      lak(2,3) = k-1
      loc(2,3) = 3
       ca(2,3) = -adtsq/hyhz4
!                       ..... 3
!                                for (i-1,j-1,k)
      lai(2,4) = i-1
      laj(2,4) = j-1
      lak(2,4) = k
      loc(2,4) = 1
       ca(2,4) = adtsq/hxhy4
!                      ..... 4
!                                for (i+1,j-1,k)
      lai(2,5) = i+1
      laj(2,5) = j-1
      lak(2,5) = k
      loc(2,5) = 1
       ca(2,5) = -adtsq/hxhy4
!                       ..... 10
!                                for (i-1,j,k)
      lai(2,6) = i-1
      laj(2,6) = j
      lak(2,6) = k
      loc(2,6) = 2
       ca(2,6) = -adtsq/hxsq
!                       .... 5
!*------------------------------ for (j,k)... diagonal -------------
      lai(2,7) = i
      laj(2,7) = j
      lak(2,7) = k
      loc(2,7) = 1
       ca(2,7) = 0.d0
!
      lai(2,8) = i
      laj(2,8) = j
      lak(2,8) = k
      loc(2,8) = 2
       ca(2,8) = 1.d0 +adtsq*(2.d0/hxsq +2.d0/hzsq)
!                 ..... 7,8,9
      lai(2,9) = i
      laj(2,9) = j
      lak(2,9) = k
      loc(2,9) = 3
       ca(2,9) = 0.d0
!*--------------------------------------------------------------------
!                                for (i+1,j,k)
      lai(2,10) = i+1
      laj(2,10) = j 
      lak(2,10) = k
      loc(2,10) = 2
       ca(2,10) = -adtsq/hxsq
!                        .... 11
!                                for (i-1,j+1,k)
      lai(2,11) = i-1
      laj(2,11) = j+1
      lak(2,11) = k
      loc(2,11) = 1
       ca(2,11) = -adtsq/hxhy4
!                        ,,,,, 6
!                                for (i+1,j+1,k)
      lai(2,12) = i+1
      laj(2,12) = j+1
      lak(2,12) = k
      loc(2,12) = 1
       ca(2,12) = adtsq/hxhy4
!                       ..... 12
!                                for (i,j-1,k+1)
      lai(2,13) = i
      laj(2,13) = j-1
      lak(2,13) = k+1
      loc(2,13) = 3
       ca(2,13) = -adtsq/hyhz4
!                        ..... 13
!                                for (i,j,k+1)
      lai(2,14) = i
      laj(2,14) = j
      lak(2,14) = k+1
      loc(2,14) = 2
       ca(2,14) = -adtsq/hzsq
!                        .... 14
!                                for (i,j+1,k+1)
      lai(2,15) = i
      laj(2,15) = j+1
      lak(2,15) = k+1
      loc(2,15) = 3
       ca(2,15) = adtsq/hyhz4
!                       ..... 15
!
!-- sz(k) --
!
!                                for (i,j-1,k-1)
      lai(3,1) = i
      laj(3,1) = j-1
      lak(3,1) = k-1
      loc(3,1) = 2
       ca(3,1) = adtsq/hyhz4
!                      ..... 2
!                                for (i-1,j,k-1)
      lai(3,2) = i-1
      laj(3,2) = j
      lak(3,2) = k-1
      loc(3,2) = 1
       ca(3,2) = adtsq/hxhz4
!                      ..... 1
!                                for (i+1,j,k-1)
      lai(3,3) = i+1
      laj(3,3) = j 
      lak(3,3) = k-1
      loc(3,3) = 1
       ca(3,3) = -adtsq/hxhz4
!                       ..... 4
!                                for (i,j+1,k-1)
      lai(3,4) = i
      laj(3,4) = j+1
      lak(3,4) = k-1
      loc(3,4) = 2
       ca(3,4) = -adtsq/hyhz4
!                       ..... 3
!                                for (i,j-1,k)
      lai(3,5) = i
      laj(3,5) = j-1
      lak(3,5) = k
      loc(3,5) = 3
       ca(3,5) = -adtsq/hysq
!                       .... 6
!                                for (i-1,j,k)
      lai(3,6) = i-1
      laj(3,6) = j
      lak(3,6) = k
      loc(3,6) = 3
       ca(3,6) = -adtsq/hxsq
!                       .... 5
!*------------------------------ for (i,j,k)... diagonal -------------
      lai(3,7) = i
      laj(3,7) = j
      lak(3,7) = k
      loc(3,7) = 1
       ca(3,7) = 0.d0
!
      lai(3,8) = i
      laj(3,8) = j
      lak(3,8) = k
      loc(3,8) = 2
       ca(3,8) = 0.d0
!
      lai(3,9) = i
      laj(3,9) = j
      lak(3,9) = k
      loc(3,9) = 3
       ca(3,9) = 1.d0 + adtsq*(2.d0/hxsq +2.d0/hysq)
!                 .... 7,8,9
!*--------------------------------------------------------------------
!                                for (i+1,j,k)
      lai(3,10) = i+1
      laj(3,10) = j 
      lak(3,10) = k
      loc(3,10) = 3
       ca(3,10) = -adtsq/hxsq
!                        .... 11
!                                for (i,j+1,k)
      lai(3,11) = i
      laj(3,11) = j+1
      lak(3,11) = k
      loc(3,11) = 3
       ca(3,11) = -adtsq/hysq
!                        .... 10
!                                for (i-1,j,k+1) 
      lai(3,12) = i-1
      laj(3,12) = j
      lak(3,12) = k+1
      loc(3,12) = 1
       ca(3,12) = -adtsq/hxhz4
!                        ..... 12
!                                for (i,j,k+1)
      lai(3,13) = i
      laj(3,13) = j-1 
      lak(3,13) = k+1
      loc(3,13) = 2 
       ca(3,13) = - adtsq/hyhz4
!                         ..... 13
!                                for (i+1,j,k+1)
      lai(3,14) = i+1
      laj(3,14) = j 
      lak(3,14) = k+1
      loc(3,14) = 1
       ca(3,14) = adtsq/hxhz4
!                       ..... 15       
!                                for (i,j+1,k+1)
      lai(3,15) = i 
      laj(3,15) = j+1
      lak(3,15) = k+1
      loc(3,15) = 2
       ca(3,15) = adtsq/hyhz4
!                       ..... 14
!
!***********************************************************************
!*  (2) Add plasma response to the diagonal terms (i,j,k).             *
!***********************************************************************
!
      bxa = aimpl*bx(i,j,k) +(1.d0-aimpl)*bx0(i,j,k) +bxc
      bya = aimpl*by(i,j,k) +(1.d0-aimpl)*by0(i,j,k) +byc
      bza = aimpl*bz(i,j,k) +(1.d0-aimpl)*bz0(i,j,k) +bzc
      bss = sqrt(bxa**2 +bya**2 +bza**2)       !!<-- wc is added 
!
      dtic  = hdt*qwi*bss
      dtic2 = dtic**2
!
      dtice = hdt*qwe*bss
      dtice2= dtice**2
!
!  adt*ex()/q0
      fmain = adtsq*qwi*qi(i,j,k)/((1.d0 +dtic2)*q0) 
      fmaine= adtsq*qwe*qe(i,j,k)/((1.d0 +dtice2)*q0) 
!
      rax= bxa/bss
      ray= bya/bss
      raz= bza/bss
!
!     +++++++++++++++++
      a11=  fmain*(1.d0 +dtic2*rax**2) +fmaine*(1.d0 +dtice2*rax**2)
      a22=  fmain*(1.d0 +dtic2*ray**2) +fmaine*(1.d0 +dtice2*ray**2)
      a33=  fmain*(1.d0 +dtic2*raz**2) +fmaine*(1.d0 +dtice2*raz**2)
!
      a12=  fmain*(dtic2*rax*ray +dtic*raz) +fmaine*(dtice2*rax*ray +dtice*raz)
      a13=  fmain*(dtic2*rax*raz -dtic*ray) +fmaine*(dtice2*rax*raz -dtice*ray)
!
      a21=  fmain*(dtic2*ray*rax -dtic*raz) +fmaine*(dtice2*ray*rax -dtice*raz)
      a23=  fmain*(dtic2*ray*raz +dtic*rax) +fmaine*(dtice2*ray*raz +dtice*rax)
!
      a31=  fmain*(dtic2*raz*rax +dtic*ray) +fmaine*(dtice2*raz*rax +dtice*ray)
      a32=  fmain*(dtic2*raz*ray -dtic*rax) +fmaine*(dtice2*raz*ray -dtice*rax)
!
!**********************************************************
!*  Core matrix: l= 1-3 and m= 9-11 are added.            *
!*    diagonal terms are --- (1,9) (2,10) (3,11)          *
!*    ca(l,m)   :  vacuum (second derivative) terms.      *
!*    a11 - a33 :  plasma response terms.                 *
!**********************************************************
!
      ca(1,7) = ca(1,7) + a11
      ca(1,8) = ca(1,8) + a12
      ca(1,9) = ca(1,9) + a13
!
      ca(2,7) = ca(2,7) + a21
      ca(2,8) = ca(2,8) + a22
      ca(2,9) = ca(2,9) + a23
!
      ca(3,7) = ca(3,7) + a31
      ca(3,8) = ca(3,8) + a32
      ca(3,9) = ca(3,9) + a33
!
!  Outside the own region is also calculated !!
!  For lai()= i-1 and i=0 on the ipar=1 node, it becomes 
!   pxc(lai(i,l,m))= pxc(-1)= mx-1, which is the ipar=npc 
!   node due to periodicity.
!
      lxyz = i +mx*(j +myA*k)
!
      do m= 1,nob
      ii= pxc(lai(1,m))
      jj=     laj(1,m)          !!<-- bound, laj()= 0-my
      kk= pzc(lak(1,m))  !lak(1,m)  ! pzc(lak(1,m))
! 
      aa(1,lxyz+1,m) = ca(1,m)
      na(1,lxyz+1,m) = 3*(ii +mx*(jj +myA*kk)) +loc(1,m) 
      end do
!                     
      do m= 1,nob
      ii= pxc(lai(2,m))
      jj=     laj(2,m) 
      kk= pzc(lak(2,m))  !lak(2,m)  ! pzc(lak(2,m))
!
      aa(2,lxyz+1,m) = ca(2,m)
      na(2,lxyz+1,m) = 3*(ii +mx*(jj +myA*kk)) +loc(2,m)
      end do
!
      do m= 1,nob
      ii= pxc(lai(3,m))
      jj=     laj(3,m) 
      kk= pzc(lak(3,m))  !lak(3,m)  ! pzc(lak(3,m))
!
      aa(3,lxyz+1,m) = ca(3,m)
      na(3,lxyz+1,m) = 3*(ii +mx*(jj +myA*kk)) +loc(3,m)
      end do
      end do 
!
!***
! (c) j= my
      else if(j.eq.my) then
!
!* side 2:                 for ex(my,k): to the right
        do i= 0,mx-1
        ca(1,1) =  1.d0
        ca(1,2) = -1.d0
!*                         for ey(my,k)
        ca(2,1) =  1.d0
        ca(2,2) =  1.d0
!*                         for ez(my,k)
        ca(3,1) =  1.d0
        ca(3,2) = -1.d0
!*
        lxyz = i +mx*(my +myA*k)
!
        do m= 1,nob
        if(m.eq.1) then
        aa(1,lxyz+1,1) = ca(1,1)  ! aa( ,m=1)
        na(1,lxyz+1,1) = 3*(i +mx*(my-1 +myA*k)) +1
!
        else if(m.eq.2) then
        aa(1,lxyz+1,2) = ca(1,2)  ! aa( ,m=2)
        na(1,lxyz+1,2) = 3*(i +mx*(my   +myA*k)) +1
!
        else
          aa(1,lxyz+1,m) = 0
          na(1,lxyz+1,m) = 3*(i +mx*(my   +myA*k)) +1
        end if
        end do
!
        do m= 1,nob
        if(m.eq.1) then
        aa(2,lxyz+1,1) = ca(2,1)  ! aa( ,m=1)
        na(2,lxyz+1,1) = 3*(i +mx*(my-1 +myA*k)) +2
!                                 +             ++
        else if(m.eq.2) then
        aa(2,lxyz+1,2) = ca(2,2)  ! aa( ,m=2)
        na(2,lxyz+1,2) = 3*(i +mx*(my   +myA*k)) +2
!
        else
          aa(2,lxyz+1,m) = 0
          na(2,lxyz+1,m) = 3*(i +mx*(my   +myA*k)) +2
        end if
        end do
!
        do m= 1,nob
        if(m.eq.1) then
        aa(3,lxyz+1,1) = ca(3,1)  ! aa( ,m=1)
        na(3,lxyz+1,1) = 3*(i +mx*(my-1 +myA*k)) +3
!                                 +             ++
        else if(m.eq.2) then
        aa(3,lxyz+1,2) = ca(3,2)  ! aa( ,m=2)
        na(3,lxyz+1,2) = 3*(i +mx*(my   +myA*k)) +3
!
        else
          aa(3,lxyz+1,m) = 0
          na(3,lxyz+1,m) = 3*(i +mx*(my   +myA*k)) +3
        end if
        end do 
        end do  
      end if
!***
      end do
      end do 
!
      return
      end subroutine emcoef
!
!
!***********************************************************************
!*                                                                     *
!*    Bi-conjugate gradient method: block type.                        *
!*                                                                     *
!***********************************************************************
!-----------------------------------------------------------------------
      subroutine bcgstb (aa,ja,na,b,u,np1,np2,nz1,nz2,ipar,  &
                         ipr,rpr,ierr)
!-----------------------------------------------------------------------
! Purpose :
!    This subroutine solves non symmetric linear systems
!    by bi-cgstab method with block precondtioning.
!
! Format :
!                 -- -- - -- -- -   --- ---
!     call bcgstb(aa,la,n,ma,ja,b,u,ipr,rpr,iw,w,wp,ierr)
!                                 - --- ---         ----
! Arguments :
!
!  aa(la,ma) (real*8) ------------ coefficient matrix
!  ja(la,ma) (integer*4) --------- column no. table
!                                  ja(i,jj) is column no.of aa(i,jj)
!  n (integer*4) ----------------- vector size
!  b(n) (real*8) ----------------- right-hand side's vector
!  u(n) (real*8) ----------------- solution vector
!  ipr(10) (integer*4) ----------- integer parameter
!                 ipr(1) (in) .... maximum number of iteration
!                 ipr(2) (out) ... final iteration number
!                 ipr(3) (in) .... type of preconditioning
!                 ipr(4) (in) .... size of block matrix
!  rpr(10) (real*8) -------------- real parameter
!                 rpr(1) (in) .... admissible residual norm
!                 rpr(2) (out) ... final residual norm
!                 rpr(3) (out) ... cpu time
!  iw(n,2) (integer*4) ----------- integer work area
!  w(n,6+ipr(4)) (real*8) -------- real work area
!                 w(*,1) ......... residual vector
!  ierr (integer*4) -------------- return code
!                 .eq. 0 ......... normal termination
!                 .ge. 1000 ...... warning: useing default values
!                 .ge. 3000 ...... fatal error: abnormal termination
!
!               Copyright : All rights are reserved by NEC (1993).
!----------------------------------------------------------------------
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa 
      integer(C_INT),dimension(1:mxyz3,nob) :: ja,na
!
      real(C_DOUBLE),dimension(1:mxyz3) :: b,u
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
!
      integer(C_INT),dimension(10) :: ipr
      real(C_DOUBLE),dimension(10) :: rpr
      integer(C_INT),save :: istart=0
      integer(C_INT) :: jm,ipar
!                       ++++++++
!
      integer(C_INT),dimension(1:mxyz3,2) :: iw
      real(C_DOUBLE),dimension(1:mxyz3,3) :: wp
!
      integer(C_INT) ierr,itrm,iflg,i,iwrt
      real(C_DOUBLE) eps,bb
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
      character(len=8) :: fortr51,fortr61
      common/fortr50/ fortr51(8),fortr61(8)
!
      real :: t1,t2,t3,t4,t5
      call cpu_time (t1)
!
!
      ierr = 0
!                                      parameter check
      call wwstbp (ipr,rpr,ierr)
!
!                                      parameter set
      itrm = ipr(1)  !<-- iterf in /cfpsol/
!     iprc = ipr(3)  !<-- iprc= 0
!     iblk = ipr(4)  !<-- iblk= 3
      eps  = rpr(1)
!
!                                      matrix check
      if(istart.eq.0) then
        istart= 1
!
!  Define iw() array only for /wwstbc,wwstbs,wwstbk/
!
        call wwstbc (aa,na,iw,jm,iflg,ierr)
!
        if ( ierr.ge.3000 ) goto 4000
!
!                                      matrix sort if necessary
        if ( iflg.eq.1 ) then
!
!  All data are defined in /wwstbs/. jm is used for it=0
!
          call wwstbs (aa,ja,na,jm,iw,ierr)
! 
          ierr = max(ierr,1200)
        endif
        if ( ierr.ge.3000 ) goto 4000
!                    ++++
      else
!       istart=1
!       Parallel version in /wwstbs3/ and /wwstbm/
!
        call wwstbs3 (aa,ja,na,iw,np1,np2,ipar,ierr)
      end if                         ! iblk in param_080A.h
      
      call cpu_time (t2)
!                                      preconditioning
!
!  Array wp() is calculated here,and used in /wwstbj/ and /wwstbk/
!
!                           ++
      call wwstbi (aa,ja,wp,iw,np1,np2,ipar,ierr) 
!                        ++ 
!
      call cpu_time (t3) 
! 
!                                      make initial value
!  Here b and u come in with wp. 
!
      if ( ierr.eq.1000 .or. ierr.eq.1010 ) then
         do i= np1(ipar),np2(ipar) 
         u(i) = b(i) 
         end do 
!
!        if(io_pe.eq.1) then
         if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
!
         open (unit=11,file=praefixc//'.11'//suffix2,             & 
               status='unknown',position='append',form='formatted')
         write(11,*) '+ After wwstbi in ierr(1000 or 1010)=',ierr
!
         close(11)
         end if
!
      else   
!  
!*  u = m^(-1)*b 
!                     + +   
         call wwstbj (b,u,wp,np1,np2,ipar) 
      end if 
! 
!                                       bi-cgstab's iteration
      call cpu_time (t4) 
!   
!  The arrays b,u are valid for i=np1(ipar)-1,np2(ipar)+1
      call wwstbk (aa,ja,b,u,np1,np2,nz1,nz2,ipar,itrm,eps, &
                   iw,wp,ierr)
!
!                                       ending
 4000 continue
      call cpu_time (t5)
!
      ipr(2) = itrm
      rpr(2) = eps
      rpr(3) = t5 - t1
!
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
!     if(.false.) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
!       bb= 0
!       do i= np1(ipar),np2(ipar)
!       bb= bb +u(i)**2
!       end do
!
!       bb= sqrt(bb/(np2(ipar)-np1(ipar)))
!       write(11,700) bb
! 700   format('bcgstb: bb=',1pd12.3)
!
!
        write(11,'("bcgstb: cpu time, bs,bi,bj,bk=",1p4d12.3)') &
                                          t2-t1,t3-t2,t4-t3,t5-t4
        close(11)
      end if
!
      return
      end subroutine bcgstb
!
!                                                         iprc=0 
!----------------------------------------------------------------------
      subroutine wwstbk (aa,ja,b,u,np1,np2,nz1,nz2,ipar,itrm,eps, &
                         jd,wp,ierr)
!----------------------------------------------------------------------
!     Purpose : solve the linear equation by bi-cgstab method.
!----------------------------------------------------------------------
!  Arrays are used in np1(ipar) and np2(ipar) (void otherwise), 
!  while they are copied by mpi_allreduce
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
      include 'mpif.h'
!
!     parameter  (nob=15)
!*---------------------------------------------------------
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa 
      integer(C_INT),dimension(1:mxyz3,nob) :: ja
!
      real(C_DOUBLE),dimension(1:mxyz3) :: b,u,r,p,q,s,t,v
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
!
      integer(C_INT),dimension(1:mxyz3,2) :: jd
      real(C_DOUBLE),dimension(1:mxyz3,3) :: wp
!
      real(C_DOUBLE) eps
      integer(C_INT) ipar,itrm,ierr
!
      real(C_DOUBLE) bb,qrofp,qrofn,qv,alpha,beta,omega, &
                     bnrm,tr,tt,rr,uu,rsdl,bbb
      real(C_DOUBLE) ttt(2),ttt2(2),qqrofn,qv1
      integer(C_INT) itr,nstb,ierror,i,iwrt
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!*---------------------------------------------------------
!
      do i= 1,mxyz3
      r(i)= 0 
      p(i)= 0 
      q(i)= 0 
      s(i)= 0 
      t(i)= 0 
      v(i)= 0 
      end do 
!
!
      itr= 0
      bb = 0
!
      do i= np1(ipar),np2(ipar) 
      bb= bb +b(i)**2
      end do
!
      call mpi_allreduce (bb,bbb,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,ierror)
      bb= bbb
      bnrm = sqrt(bb)
!
!*  Source check.
!
      if(bb.lt.1.d-72) then
        ierr= 1300
!
        do i= np1(ipar),np2(ipar)
        u(i)= 0
        end do
        go to 2000
!
      else if(bb.gt.1.d+72) then
         ierr= 3500
         go to 2000
      endif
!
!*  r = a*u 
!             
      call wwstbm (u,r,aa,ja,np1,np2,nz1,nz2,ipar)
!                  + r
!
      do i= np1(ipar),np2(ipar) 
      r(i)= b(i) -r(i)  !<-- minimize b(i)-r(i)
      end do
!
!
      do i= np1(ipar),np2(ipar) 
      q(i)= r(i)
      p(i)= 0
      v(i)= 0
      end do
!
      alpha= 1.d0
      omega= 1.d0
      qrofn= 1.d0
!                                       iteration
!*-----------------------------------------***--------------------------
 1000 itr = itr + 1
!
      qrofp= qrofn
      qrofn= 0
!
      do i= np1(ipar),np2(ipar) 
      qrofn= qrofn +q(i)*r(i)
      end do
! 
      call mpi_allreduce (qrofn,qqrofn,1,mpi_real8,mpi_sum,  &
                          mpi_comm_world,ierror)
      qrofn= qqrofn
!
!
      if(abs(qrofn).lt.1.d-72) then
         ierr= 3600
      endif
!                                            beta = (qrofn*alpa)/
!                                                  /(qrofp*omeg)
!*  To avoid zero-divide.  11/24/95
!
!     if(abs(qrofp).gt.1.d-10) then
         beta = (qrofn*alpha)/(qrofp*omega)
!     else
!        beta = alpha/omega  !<-- bit differ
!     end if
!
!
      do i= np1(ipar),np2(ipar) 
      v(i)= r(i) -beta*omega*v(i)
      end do
!
!*  s = wp^(-1)*v           
!                  + +
      call wwstbj (v,s,wp,np1,np2,ipar) 
!                          wp() is calculated outside 
!
      do i= np1(ipar),np2(ipar) 
      p(i)= s(i) +beta*p(i)
      end do
!
!*  v = a*p
!
      call wwstbm (p,v,aa,ja,np1,np2,nz1,nz2,ipar)
!                  + r
!
      qv= 0
!
      do i= np1(ipar),np2(ipar) 
      qv= qv +q(i)*v(i)
      end do
! 
      call mpi_allreduce (qv,qv1,1,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      qv= qv1
!
!
      if(abs(qv).lt.1.d-72) then
        ierr= 3700
        go to 2000
      endif
!
      alpha = qrofn/qv
!
! 
      do i= np1(ipar),np2(ipar) 
      r(i)= r(i) -alpha*v(i)
      end do
!
!*  s = wp^(-1)*r                ! wp() from outside
!
      call wwstbj (r,s,wp,np1,np2,ipar)
!                  + r 
!*  t = a*s
! 
      call wwstbm (s,t,aa,ja,np1,np2,nz1,nz2,ipar)
!                  + r
!
      tr= 0
      tt= 0
!
      do i= np1(ipar),np2(ipar) 
      tr= tr +t(i)*r(i)
      tt= tt +t(i)*t(i)
      end do
! 
      ttt(1)= tr
      ttt(2)= tt
      call mpi_allreduce (ttt(1),ttt2(1),2,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      tr= ttt2(1)
      tt= ttt2(2)
!
!
      if(tt.lt.1.d-72) then
        ierr = 3700
        go to 2000
      endif
!
      omega = tr/tt
!
      do i= np1(ipar),np2(ipar) 
      u(i)= u(i) + alpha*p(i) + omega*s(i) 
      r(i)= r(i) - omega*t(i)
      end do
!
!
      rr = 0
      uu = 0
!
      do i= np1(ipar),np2(ipar) 
      rr= rr + r(i)**2
      uu= uu + u(i)**2
      end do
! 
      ttt(1)= rr
      ttt(2)= uu
      call mpi_allreduce (ttt(1),ttt2(1),2,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      rr= ttt2(1)
      uu= ttt2(2)
!
!
      if(rr.gt.1.d+72) then
         ierr= 3800
         go to 2000
      endif
!
      rsdl= sqrt(rr)/bnrm
      uu  = sqrt(uu/float(mxyz3)) 
!
!     ***********************************
      if(rsdl.lt.eps) go to 2000
!
      if(itr.le.itrm) then
        go to 1000     !!<--iteration continues
      else
        ierr= 2000
        nstb= nstb +1  !!<--pass this point
      endif
!     ***********************************
!*-----------------------------------------***--------------------------
!                                       end of iteration
 2000 itrm = itr
      eps  = rsdl
!
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("++wwstbk: itr=",i4,"  uu,rsdl=",1p2d13.4)') &
                                                      itr,uu,rsdl
        close(11)
      end if
!     end if
!
      return
      end subroutine wwstbk
!
!
!----------------------------------------------------------------------
      subroutine wwstbp (ipr,rpr,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
      integer(C_INT) ipr(10),ierr,n
      real(C_DOUBLE) rpr(10)
!
!     if ( mxyz3.le.0 ) then
!        ierr = 3000
!        return
!     endif
!
!     if ( mxyz3.lt.mxyz3 ) then
!        ierr = 3100
!        return
!     endif
!
      if ( ipr(1).le.0 ) then
         ipr(1) = 5*int(sqrt(real(mxyz3)))
         ierr   = 1100
      endif
!
      if ( ipr(3).lt.0 .or. ipr(3).gt.1 ) then
         ipr(3) = 0
         ierr   = 1110
      endif
!
      if ( ipr(4).le.1 .or. ipr(4).ge.4 ) then
         ierr   = 3400
      endif
!
      if ( rpr(1).lt.1.d-72 ) then
         rpr(1) = 1.d-12
         ierr   = 1150
      endif
!
      return
      end subroutine wwstbp
!
!
!----------------------------------------------------------------------
      subroutine wwstbc (aa,na,jd,jm,iflg,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
      integer(C_INT),dimension(npc) :: np1,np2
      integer(C_INT) ipar
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa 
      integer(C_INT),dimension(1:mxyz3,nob) :: na
      integer(C_INT),dimension(1:mxyz3,2) :: jd
!
      integer(C_INT) iflg,jm,i,j,jpst,j0,j1,jj,ierr
!
      iflg = 0
      jm = 0
!
      do i= 1,mxyz3
         jpst = 0
         j0 = ((i-1)/iblk)*iblk + 1
         j1 = ((i-1)/iblk)*iblk + iblk
         jd(i,1) = 0
!
         do jj= 1,nob
            j = na(i,jj)
            if ( j.lt.1 .or. j.gt.mxyz3 ) then
               aa(i,jj) = 0.d0
               na(i,jj) = i
               goto 110
            endif
            if ( j.eq.i .and. abs(aa(i,jj)).lt.1.d-72 ) goto 110
!
            jm = max(jm,jj)  !<<-- this defines max jm value
            if ( iflg.eq.1 ) goto 110
!
            if ( j.gt.jpst ) then
               jpst = j
            else
               iflg = 1
               goto 110
            endif
!
            if ( jd(i,1).eq.0 .and. j.ge.j0 ) then
               jd(i,1) = jj
            endif
            if ( jd(i,1).ne.0 .and. j.le.j1 ) then
               jd(i,2) = jj
            endif
  110    end do     
!
         if ( jd(i,2)-jd(i,1).ge.iblk ) then
            iflg = 1
         endif
!
      end do
!
      return
      end subroutine wwstbc
!
!
!*vocl total, scalar
!----------------------------------------------------------------------
      subroutine wwstbs (aa,ja,na,jm,jd,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'  !<-- mx,myA,mz,iblk in param_080A.h
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa,ax 
      integer(C_INT),dimension(1:mxyz3,nob) :: ja,na
      integer(C_INT),dimension(1:mxyz3,2) :: jd
      integer(C_INT) jm,ierr  
!
      real(C_DOUBLE) a0
      integer(C_INT) i,j,j0,j1,jc,ii,jj
      logical :: first= .true.
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
! !NEC$ NOVECTOR
      do i= 1,mxyz3 
!
      do j= 1,nob   ! ma
      ja(i,j)= na(i,j)  !<-- ja(i,) from original na(i,)'s 
        ax(i,j)= aa(i,j)  !!<-- memorize
      end do
!     ++++++++++
!
      do j0=  1,jm
      j  = ja(i,j0)
      a0 = aa(i,j0)
!
        if ( j.eq.i .and. abs(a0).lt.1.d-72 ) then
          ja(i,j0) = mxyz3 + 1
          goto 110
        endif
!
      do j1= j0-1,1,-1
      if ( ja(i,j1).gt.j ) then
        ja(i,j1+1) = ja(i,j1)
        aa(i,j1+1) = aa(i,j1)
      else
        ja(i,j1+1) = j
        aa(i,j1+1) = a0  !<-- hit matches
        goto 110
      endif
      end do
!
      ja(i,1) = j
      aa(i,1) = a0
  110 continue
      end do
!     ++++++
!
      j0 = ((i-1)/iblk)*iblk + 1
      j1 = ((i-1)/iblk)*iblk + iblk
      jd(i,1) = 0
!
!     ++++++++++
      do jj= 1,jm
      j = ja(i,jj)
!
       if ( j.eq.mxyz3+1 ) then  !<-- i-th line is given
         ja(i,jj) = i            !  y-end
         goto 130
       endif
!                         jj= 1,2,...
!                           j0 =< ja(i,jj) 
      if ( jd(i,1).eq.0 .and. j.ge.j0 ) then
        jd(i,1) = jj
      endif
!                           ja(i,jj) <= j1    
      if ( jd(i,1).ne.0 .and. j.le.j1 ) then
        jd(i,2) = jj
      endif
!
  130 continue
      end do
!     ++++++++++
!
        if ( jd(i,2)-jd(i,1).ge.iblk ) then
           ierr = 3200
           return
        endif
      end do     
!***
!
      if(.false.) then
!     if(io_pe.eq.1 .and. first) then
!!        first= .false.
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        do ii= 1,3
        i= ii
        do j= 1,nob
        jc= na(i,j)
        jj= ja(i,j)
        write(11,'("# i,j,jc,ax,jj,aa=",2i6,2x,i10,d12.3,2x,i10,d12.3)') &
                                                i,j,jc,ax(i,j),jj,aa(i,j)
        end do
        end do
!
        do ii= 1,3
        i= 3*mx +ii
        do j= 1,nob
        jc= na(i,j)
        jj= ja(i,j)
        write(11,'("# i,j,jc,ax,jj,aa=",2i6,2x,i10,d12.3,2x,i10,d12.3)') &
                                                i,j,jc,ax(i,j),jj,aa(i,j)
        end do
        end do
!
        close(11)
      end if
!
      return
      end subroutine wwstbs
!
! 
!----------------------------------------------------------------------
      subroutine wwstbs3 (aa,ja,na,jd,np1,np2,ipar,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'  !<-- mx,myA,mz,iblk in param_080A.h
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa,ax 
      integer(C_INT),dimension(1:mxyz3,nob) :: ja,na
      integer(C_INT),dimension(1:mxyz3,2) :: jd
!
      integer(C_INT),dimension(npc) :: np1,np2
      integer(C_INT) ipar,ierr  ! jm
!
      real(C_DOUBLE) a0
      integer(C_INT) i,j,j0,j1,jc,ii,jj
      logical :: first= .true.
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      do i= np1(ipar),np2(ipar)
!
      do j= 1,nob 
      ja(i,j)= na(i,j)  !<-- ja(i,) copied from original na(i,)
        ax(i,j)= aa(i,j)  !!<-- output
      end do
!     ++++++++++
!
      do j0=  1,nob  ! jm
      j  = ja(i,j0)
      a0 = aa(i,j0)
!
        if ( j.eq.i .and. abs(a0).lt.1.d-72 ) then
          ja(i,j0) = mxyz3 + 1
          goto 110
        endif
!
      do j1= j0-1,1,-1
      if ( ja(i,j1).gt.j ) then
        ja(i,j1+1) = ja(i,j1)
        aa(i,j1+1) = aa(i,j1)
      else
        ja(i,j1+1) = j
        aa(i,j1+1) = a0 
        goto 110
      endif
      end do
!
      ja(i,1) = j
      aa(i,1) = a0
  110 continue
      end do
!     ++++++
!
      j0 = ((i-1)/iblk)*iblk + 1
      j1 = ((i-1)/iblk)*iblk + iblk
      jd(i,1) = 0
!
!     ++++++++++
      do jj= 1,nob  ! jm
      j = ja(i,jj)
!
      if ( j.eq.mxyz3+1 ) then  !<-- above 110
        ja(i,jj) = i   
        goto 130
      endif
!                           j0 =< ja(i,jj) 
      if ( jd(i,1).eq.0 .and. j.ge.j0 ) then
        jd(i,1) = jj
      endif
!                           ja(i,jj) <= j1    
      if ( jd(i,1).ne.0 .and. j.le.j1 ) then
        jd(i,2) = jj
      endif
!
  130 continue
      end do
!     ++++++++++
!
        if ( jd(i,2)-jd(i,1).ge.iblk ) then
           ierr = 3200
           return
        endif
      end do     
!***
!
      if(.false.) then
!     if(io_pe.eq.1 .and. first) then
        first= .false.
!
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        do ii= 1,3
        i= ii
        do j= 1,nob
        jc= na(i,j)
        jj= ja(i,j)
        write(11,'("# i,j,jc,ax,jj,aa=",2i6,2x,i10,d12.3,2x,i10,d12.3)') &
                                              i,j,jc,ax(i,j),jj,aa(i,j)
        end do
        end do
!
        do ii= 1,3
        i= 3*mx +ii
        do j= 1,nob
        jc= na(i,j)
        jj= ja(i,j)
        write(11,'("# i,j,jc,ax,jj,aa=",2i6,2x,i10,d12.3,2x,i10,d12.3)') &
                                              i,j,jc,ax(i,j),jj,aa(i,j)
        end do
        end do
!
        close(11)
      end if
!
      return
      end subroutine wwstbs3
!
!
!----------------------------------------------------------------------
      subroutine wwstbi (aa,ja,wp,jd,np1,np2,ipar,ierr)
!----------------------------- **-- wp is defined outside -------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!                                     
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa
      integer(C_INT),dimension(1:mxyz3,nob) :: ja
!
      real(C_DOUBLE),dimension(1:mxyz3,3) :: wp
      integer(C_INT),dimension(1:mxyz3,2) :: jd
!
      integer(C_INT),dimension(npc) :: np1,np2
      integer(C_INT) ipar,ierr,i,j
!
      real(C_DOUBLE) w0(3*3),w1(3,3),det
      equivalence   (w0(1),w1(1,1))
!
      integer(C_INT) i0,ii,ij,jj
!
!
!     if ( iprc.eq.1 ) goto 1000
!
      do i0= np1(ipar)-1,np2(ipar)-1,iblk    ! 0,n-1,iblk
!*
      do ij= 1,iblk*iblk
      w0(ij) = 0.d0
      end do
!
      do ii= 1,iblk
      i = i0 + ii
!
!  i0=0 ; j0=((i-1)/3)*3+1=1; j1=0+3=3; 
!    jj=1,3
!     ja(1,1)-(1-1)-1=ja(1,m)-1
!     ja(2, )
!     ja(3, )
      do jj= jd(i,1),jd(i,2)
      j = ja(i,jj)
      w0(ii+(j-i0-1)*iblk) = aa(i,jj)
      end do
      end do
!
!   call wwstbb (w0,wp,np1,np2,ipar,i0,ierr)
!
      if ( iblk.eq.2 ) then
!
        det = w1(1,1)*w1(2,2) - w1(1,2)*w1(2,1)
!
        if ( abs(det).lt.1.d-72 ) then
          ierr = 3300
          return
        else
          det = 1.d0/det
        endif
!
        wp(i0+1,1) =  w1(2,2)*det
        wp(i0+1,2) = -w1(1,2)*det
        wp(i0+2,1) = -w1(2,1)*det
        wp(i0+2,2) =  w1(1,1)*det
!
      elseif ( iblk.eq.3 ) then
!
        det = w1(1,1)*w1(2,2)*w1(3,3) + w1(1,2)*w1(2,3)*w1(3,1)  & 
            + w1(1,3)*w1(3,2)*w1(2,1) - w1(1,3)*w1(2,2)*w1(3,1)  &
            - w1(1,2)*w1(2,1)*w1(3,3) - w1(1,1)*w1(3,2)*w1(2,3)
!
        if ( abs(det).lt.1.d-72 ) then
          ierr = 3300
          return
        else
          det = 1.d0/det
        endif
!
        wp(i0+1,1) =  (w1(2,2)*w1(3,3) -w1(2,3)*w1(3,2))*det
        wp(i0+1,2) = -(w1(1,2)*w1(3,3) -w1(1,3)*w1(3,2))*det
        wp(i0+1,3) =  (w1(1,2)*w1(2,3) -w1(1,3)*w1(2,2))*det
!
        wp(i0+2,1) = -(w1(2,1)*w1(3,3) -w1(2,3)*w1(3,1))*det
        wp(i0+2,2) =  (w1(1,1)*w1(3,3) -w1(1,3)*w1(3,1))*det
        wp(i0+2,3) = -(w1(1,1)*w1(2,3) -w1(1,3)*w1(2,1))*det
!
        wp(i0+3,1) =  (w1(2,1)*w1(3,2) -w1(2,2)*w1(3,1))*det
        wp(i0+3,2) = -(w1(1,1)*w1(3,2) -w1(1,2)*w1(3,1))*det
        wp(i0+3,3) =  (w1(1,1)*w1(2,2) -w1(1,2)*w1(2,1))*det
      end if
      end do
!*             wp( ) are copied to outside 
      return
!
!1000 continue
!     write(11,*) 'Not incremented. RETURN !'
!
      end subroutine wwstbi
!
!                                              iprc=0
!-------------------------------------------------------------------------
      subroutine wwstbj (v,w,wp,np1,np2,ipar)
!------------------------*-*----------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(1:mxyz3) :: v,w
!
      real(C_DOUBLE),dimension(1:mxyz3,3) ::wp
      integer(C_INT),dimension(npc) :: np1,np2
      integer(C_INT) ipar,i
!
      if (iblk.eq.2) then
!
         do i= np1(ipar),np2(ipar),2   ! 1,n,2
         w(i)   = wp(i,1)  *v(i) + wp(i,2)  *v(i+1)
         w(i+1) = wp(i+1,1)*v(i) + wp(i+1,2)*v(i+1)
         end do
!
      elseif (iblk.eq.3) then
!
         do i= np1(ipar),np2(ipar),3   ! 1,n,3
         w(i)   = wp(i,  1)*v(i) + wp(i,  2)*v(i+1) + wp(i,  3)*v(i+2)
         w(i+1) = wp(i+1,1)*v(i) + wp(i+1,2)*v(i+1) + wp(i+1,3)*v(i+2)
         w(i+2) = wp(i+2,1)*v(i) + wp(i+2,2)*v(i+1) + wp(i+2,3)*v(i+2)
         end do
!
!     elseif ( iprc.eq.1 ) then
!        write(11,*) 'Not incremented. RETURN !'
!
      end if
!
      return
      end subroutine wwstbj
!
!
!----------------------------------------------------------------------
      subroutine wwstbm (v,w,aa,ja,np1,np2,nz1,nz2,ipar) 
!------------------------s-r-------------------------------------------
!**********************************************************************
!     Calculate:  w(i) = sum(j)| a(i,j)*v(na(i,j)).                   *
!**********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
!
!     parameter  (nob=15)
!----------------------------------------------------------------------
      real(C_DOUBLE),dimension(1:mxyz3) :: v,w,vl
!
      real(C_DOUBLE),dimension(1:mxyz3,nob) :: aa
      integer(C_INT),dimension(1:mxyz3,nob) :: ja
! 
      integer(C_INT),dimension(npc) :: np1,np2,nz1,nz2
      integer(C_INT) ipar,rank,i,j,k,ij,ijk 
!
      real(C_DOUBLE),dimension(mxyz/npc) :: vgx,vgy,vgz,vggx,vggy,vggz
!
      integer(C_INT) io_pe,ijm
      common/iope66/ io_pe
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      character(len=8) :: fortr51,fortr61
      common/fortr50/ fortr51(8),fortr61(8)
!*---------------------------------------------------------------------
!* Copy an array to overlap the original.
!
      do i= np1(ipar),np2(ipar)
      vl(i)= v(i)
      end do
!
!  to leftward side
!       call mpi_irecv (recvbuf,mxy,mpi_real8,rank-1,mpi_any_tag, &
!                       mpi_comm_world,recv_request,ierror)
!       call mpi_isend (sendbuf,mxy,mpi_real8,rank+1,0,           &
!                       mpi_comm_world,send_request,ierror)
!  to rightward side 
!       call mpi_irecv (recvbuf,mxy,mpi_real8,rank+1,mpi_any_tag, &
!                       mpi_comm_world,recv_request,ierror)
!       call mpi_isend (sendbuf,mxy,mpi_real8,rank-1,0,           &
!                       mpi_comm_world,send_request,ierror)
!
      rank= ipar -1
!     +++++++++++++
!*
!  To the left boundary.
!
      ij= 0
      do k= nz1(ipar),nz2(ipar)
      do j= 1,my+1
      do i= 1,mx 
!
      ij= ij +1
      ijk= 3*(i-1 +mx*(j-1 +myA*(k-1)))
! 
      vgx(ij)= v(1+ijk) 
      vgy(ij)= v(2+ijk)
      vgz(ij)= v(3+ijk)
      end do
      end do
      end do
!
      call sendrev1 (vgx,vggx,rank,mxyz/npc) 
      call sendrev1 (vgy,vggy,rank,mxyz/npc) 
      call sendrev1 (vgz,vggz,rank,mxyz/npc) 
!
      if(ipar.eq.1) then
        ij= 0
        do k= nz1(npc),nz2(npc)
        do j= 1,my+1
        do i= 1,mx
!
        ij= ij +1
        ijk= 3*(i-1 +mx*(j-1 +myA*(k-1)))
!
        vl(ijk+1)= vggx(ij) 
        vl(ijk+2)= vggy(ij) 
        vl(ijk+3)= vggz(ij) 
        end do
        end do
        end do
!
      else
        ij= 0
        do k= nz1(ipar-1),nz2(ipar-1)
        do j= 1,my+1
        do i= 1,mx
!
        ij= ij +1
        ijk= 3*(i-1 +mx*(j-1 +myA*(k-1)))
!
        vl(ijk+1)= vggx(ij) 
        vl(ijk+2)= vggy(ij) 
        vl(ijk+3)= vggz(ij) 
        end do
        end do
        end do
      end if
!
!*
!  To the right boundary 
!
      ij= 0
      do k= nz1(ipar),nz2(ipar) 
      do j= 1,my+1
      do i= 1,mx
!
      ij= ij +1
      ijk= 3*(i-1 +mx*(j-1 +myA*(k-1)))
! 
      vgx(ij)= v(1+ijk) 
      vgy(ij)= v(2+ijk)
      vgz(ij)= v(3+ijk)
      end do 
      end do
      end do
!
      call sendrev2 (vgx,vggx,rank,mxyz/npc) 
      call sendrev2 (vgy,vggy,rank,mxyz/npc) 
      call sendrev2 (vgz,vggz,rank,mxyz/npc) 
!
      if(ipar.eq.npc) then
        ij= 0
        do k= nz1(1),nz2(1) 
        do j= 1,my+1
        do i= 1,mx
!
        ij= ij +1
        ijk= 3*(i-1 +mx*(j-1 +myA*(k-1))) ! Can be left or right
! 
        vl(ijk+1)= vggx(ij) 
        vl(ijk+2)= vggy(ij) 
        vl(ijk+3)= vggz(ij)
        end do
        end do
        end do
!
      else
        ij= 0
        do k= nz1(ipar+1),nz2(ipar+1) 
        do j= 1,my+1
        do i= 1,mx
!
        ij= ij +1
        ijk= 3*(i-1 +mx*(j-1 +myA*(k-1))) ! Can be left or right
! 
        vl(ijk+1)= vggx(ij) 
        vl(ijk+2)= vggy(ij) 
        vl(ijk+3)= vggz(ij)
        end do
        end do
        end do
      end if
!
!
!***********************************************************************
!* (2) Calculate a*v on each pe's.                                     *
!***********************************************************************
!  In the expanded pe's nodes including left and right pe's.
!
      do i= np1(ipar),np2(ipar)
      w(i)= 0.d0
      end do
!
      do i= np1(ipar),np2(ipar)
      do j= 1,nob 
      w(i)= w(i) + aa(i,j)*vl(ja(i,j))
      end do  
      end do
!
      return
      end subroutine wwstbm
!
!
!---------------------------------------------------------------
      subroutine sendrev1 (sendbuf,recvbuf,rank,mxy)
!---------------------------------------------------------------
! Isend/Irecv: periodic case
!   https://wiki.calculquebec.ca/w/Isend-Irecv_non_bloquant
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'mpif.h'
      include   'param_080A.h'
!
      real(C_DOUBLE),dimension(mxy) :: sendbuf,recvbuf
      integer(C_INT) rank,mxy,ierror
!
      integer(kind=4),dimension(MPI_STATUS_SIZE)               &
                                   :: send_status, recv_status
      integer(kind=4),dimension(1) :: send_request, recv_request
!+++
! Transfer to upward: periodic case
!   sendbuff: 1 (fourth element) 
!   recvbuff: npc-1, the same size is assumed
!
      if( rank.eq.0 ) then
        call MPI_Irecv (recvbuf,mxy,mpi_real8,npc-1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,1,0,               &
                        mpi_comm_world,send_request,ierror)  
!
      else if( rank.eq.(npc-1) ) then
        call MPI_Irecv (recvbuf,mxy,mpi_real8,npc-2,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,0,0,               &
                        mpi_comm_world,send_request,ierror)
!
      else
        call MPI_Irecv (recvbuf,mxy,mpi_real8,rank-1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,rank+1,0,           &
                        mpi_comm_world,send_request,ierror)
      end if
!
      call mpi_wait (send_request,send_status,ierror)
      call mpi_wait (recv_request,recv_status,ierror)
!
      return
      end subroutine sendrev1
!
!
!---------------------------------------------------------------
      subroutine sendrev2 (sendbuf,recvbuf,rank,mxy)
!---------------------------------------------------------------
! Isend/Irecv: periodic case
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'mpif.h'
      include   'param_080A.h'
!
      real(C_DOUBLE),dimension(mxy) :: sendbuf,recvbuf
      integer(C_INT) rank,mxy,ierror
!
      integer(kind=4),dimension(MPI_STATUS_SIZE)               &
                                   :: send_status, recv_status
      integer(kind=4),dimension(1) :: send_request, recv_request
!+++
! Transfer to downward: 
!
      if ( rank.eq.0 ) then
        call MPI_Irecv (recvbuf,mxy,mpi_real8,1,mpi_any_tag,  &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,npc-1,0,        &
                        mpi_comm_world,send_request,ierror)  
!
      else if( rank.eq.(npc-1) ) then
        call MPI_Irecv (recvbuf,mxy,mpi_real8,0,mpi_any_tag,  &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,npc-2,0,        &
                        mpi_comm_world,send_request,ierror)
!
      else
        call MPI_Irecv (recvbuf,mxy,mpi_real8,rank+1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,mxy,mpi_real8,rank-1,0,           &
                        mpi_comm_world,send_request,ierror)
      end if
!
      call mpi_wait (send_request,send_status,ierror)
      call mpi_wait (recv_request,recv_status,ierror)
!
      return
      end subroutine sendrev2
!
!
!----------------------------------------------------------------------
      subroutine wwstbb (w0,wp,np1,np2,ipar,i0,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
      integer(C_INT),dimension(npc) :: np1,np2
      integer(C_INT) ipar,i0,ierr
!
      real(C_DOUBLE) w0(iblk,iblk),det
      real(C_DOUBLE) wp(1:mxyz3,3)
!
      if ( iblk.eq.2 ) then
         det = w0(1,1)*w0(2,2) - w0(1,2)*w0(2,1)
         if ( abs(det).lt.1.d-72 ) then
            ierr = 3300
            return
         else
            det = 1.d0/det
         endif
!
         wp(i0+1,1) =  w0(2,2)*det
         wp(i0+1,2) = -w0(1,2)*det
         wp(i0+2,1) = -w0(2,1)*det
         wp(i0+2,2) =  w0(1,1)*det
!
      elseif ( iblk.eq.3 ) then
         det = w0(1,1)*w0(2,2)*w0(3,3) + w0(1,2)*w0(2,3)*w0(3,1)  & 
             + w0(1,3)*w0(3,2)*w0(2,1) - w0(1,3)*w0(2,2)*w0(3,1)  &
             - w0(1,2)*w0(2,1)*w0(3,3) - w0(1,1)*w0(3,2)*w0(2,3)
!
         if ( abs(det).lt.1.d-72 ) then
            ierr = 3300
            return
         else
            det = 1.d0/det
         endif
!
         wp(i0+1,1) =  (w0(2,2)*w0(3,3) -w0(2,3)*w0(3,2))*det
         wp(i0+1,2) = -(w0(1,2)*w0(3,3) -w0(1,3)*w0(3,2))*det
         wp(i0+1,3) =  (w0(1,2)*w0(2,3) -w0(1,3)*w0(2,2))*det
         wp(i0+2,1) = -(w0(2,1)*w0(3,3) -w0(2,3)*w0(3,1))*det
         wp(i0+2,2) =  (w0(1,1)*w0(3,3) -w0(1,3)*w0(3,1))*det
         wp(i0+2,3) = -(w0(1,1)*w0(2,3) -w0(1,3)*w0(2,1))*det
         wp(i0+3,1) =  (w0(2,1)*w0(3,2) -w0(2,2)*w0(3,1))*det
         wp(i0+3,2) = -(w0(1,1)*w0(3,2) -w0(1,2)*w0(3,1))*det
         wp(i0+3,3) =  (w0(1,1)*w0(2,2) -w0(1,2)*w0(2,1))*det
      endif
!
      return
      end subroutine wwstbb
!
!
!----------------------------------------------------------------------
      subroutine wwstba (w0,w1,w2)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit real(C_DOUBLE) (a-h,o-z)
!
      include 'param_080A.h'
      dimension w0(iblk,iblk),w1(iblk,iblk),w2(iblk,iblk)
!*                         <-- papam_A13X.h
!*vdir novector
!*vocl loop, scalar
!
      do i=1,iblk
!*
!*vdir novector
!*vocl loop, scalar
!
         do j=1,iblk
         w2(i,j) = 0.d0
!*
!*vdir novector
!*vocl loop, scalar
!
            do k=1,iblk
            w2(i,j) = w2(i,j) + w0(i,k)*w1(k,j)
            end do
         end do
      end do
!
      return
      end subroutine wwstba
!
!
!***********************************************************************
!*    The Poisson solver which uses /bcgsts/ matrix solver.            *
!***********************************************************************
!-----------------------------------------------------------------------
      subroutine poissn (phi,rho,nobx,ndim,itrmax,iterp) 
!-----------------------------------------------------------------------
!***********************************************************************
!*    << Poisson solver in 3-d : cresmd method >>                      *
!*                                                                     *
!*        Laplacian*phi(x,y,z) = rho(x,y,z).                           *
!*                                                                     *
!*   rho(i,j,k)....... source term.             <----- input           *
!*                                                                     *
!*   phi(i,j,k)... solution in (x,y,z) space.                          *
!*            phi comes with a guess.           <----- input,output    *
!*                                                                     *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
!
!     parameter  (nob3=7)
!*----------------------------------------------------------------------
!                              +++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: phi,rho
      real(C_DOUBLE),dimension(mxyz) :: ss,xx
!
      integer(C_INT) ipar,size,ndim,itrmax,iterp
!***                                               poisson routine
      real(C_DOUBLE),dimension(mxyz,nob3) :: aa       !!<-- first defined
      integer(C_INT),dimension(mxyz,nob3) :: na 
!
      integer(C_INT) nobx,ipr(10)
      real(C_DOUBLE) rpr(10)
!
      real(C_DOUBLE) eps,rsdl
!     common/cresm1/ aa  !!<-- poissn, nobx= 7
!!    common/cresm1/ aa  !!<-- escsol, nobx= 19 
!*----------------------------------------------------------------------
!
      character(len=8) poiss_sol
      common/poiss/ poiss_sol
!
      integer(C_INT) itrm,ijk,i,j,k,ierr
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!***
      do k= 0,mz-1 
      do j= 0,my
      do i= 0,mx-1
      ijk= i +mx*(j +myA*k)        !<-- ijk=0
!
      ss(ijk+1)= rho(i,j,k)
      xx(ijk+1)= 0
      end do
      end do
      end do
!
!***
!  Bound by y: poisson
!
      do k= 0,mz-1 
      do j= 0,my
      do i= 0,mx-1
      if(j.eq.0 .or. j.eq.my) then
        ijk= i +mx*(j +myA*k) 
!
        ss(ijk+1)= 0
      end if
      end do
      end do
      end do
!***
!                                       define the cfp matrix
      call emcof3 (aa,na,ndim)
!     ++++++++++++++++++++++++  
!                                       squeeze zero elements.
!                                         aa + na --> aa' + na'; la
      itrm = 300  !itrmax
      eps  = 1.0d-5
!
      ipr(1) = itrm
      rpr(1) = eps
!
!                                       solve the equation
!                  ++ R             +++++++
      call cresmd (ss,xx,aa,na,nobx,ipr,rpr,ierr)
!     +++++++++++++++++++++++++++++++++++++++++++
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ijk= i +mx*(j +myA*k) 
!
      if(j.eq.0 .or. j.eq.my) then
        phi(i,j,k)= 0
      else
        phi(i,j,k)= xx(ijk+1) 
      end if
      end do
      end do
      end do
!
      iterp = ipr(2)
      rsdl  = rpr(2)
!
!     if(ierr.ne.0) then
!        write (6 ,*) '#(bcg-ps)  iter, rsdl= ',ipr(2),rpr(2)
!     end if
!
      return
      end subroutine poissn
!
!
!-----------------------------------------------------------------------
      subroutine emcof3 (aa,na,ndim) 
!-----------------------------------*******-----------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
!     parameter  (nob3=7)
!*----------------------------------------------------------------
      real(C_DOUBLE),dimension(mxyz,nob3) :: aa
      integer(C_INT),dimension(mxyz,nob3) :: na
      integer(C_INT) ndim
!                                   ******
      real(C_DOUBLE),dimension(nob3) :: ca
      integer(C_INT),dimension(nob3) :: lai,laj,lak
!
      integer(C_INT) i,j,k,ii,jj,kk,m,lxyz
!     real(C_DOUBLE) ca2,ca4
!*---------------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      integer(C_INT) io_pe
      common/iope66/ io_pe
!----------------------------------------------------------------------
!
!     if(ndim.eq.3) then
!        ca2 =  1.d0/hysq 
!        ca4 =  2.d0/hysq +2.d0/hzsq 
!
!     else if(ndim.eq.2) then
!        ca2 =  0
!        ca4 =  2.d0/hzsq
!     end if
!
!*                               !!<--emcof3
      do k= 0,mz-1         ! All k !!
      do j= 0,my           ! 200 j=0:my
!***
      if(j.eq.0) then
!**
        do i= 0,mx-1
!* (i,j-1,k)
        lai(1) =  i      ! number i=0: lai(1,0)=0
        laj(1) =  0      !  na(,)=1    laj(1,0)=0
        lak(1) =  k      !             lak(1,0)=0
         ca(1) =  1.d0   !            1:nob2 0:mx-1
!
!* (i,j+1,k)
        lai(2)=   i
        laj(2)=   1
        lak(2)=   k
         ca(2)=   1.d0 
!*
        lxyz= i +mx*(0 +myA*k)
!
        aa(lxyz+1,1) = ca(1)  ! aa(1,m=1)
        na(lxyz+1,1) = i +mx*(0 +myA*k) +1  ! na(1,m=1)=1
!            l    m           +          +
        aa(lxyz+1,2) = ca(2)  ! aa(1,m=2) <- lxyz+1 !!
        na(lxyz+1,2) = i +mx*(1 +myA*k) +1
! 
        do m= 3,nob3
        aa(lxyz+1,m)= 0
        na(lxyz+1,m)= lxyz +1 
        end do
        end do 
!***
! 
      else if(j.gt.0 .and. j.lt.my) then
!
      do i= 0,mx-1    ! 300
!* (i,j,k-1)
      lai(1)=   i
      laj(1)=   j
      lak(1)=   k-1
       ca(1)=   1.d0/hzsq
!
!* (j-1,k)
      lai(2)=   i
      laj(2)=   j-1
      lak(2)=   k
       ca(2)=   1.d0/hysq
!
!* (i-1,j,k)
      lai(3)=   i-1
      laj(3)=   j
      lak(3)=   k
       ca(3)=   1.d0/hxsq
!
!* (i,j,k)
      lai(4)=   i
      laj(4)=   j
      lak(4)=   k
       ca(4)=   -2.d0*(1/hxsq +1/hysq +1/hzsq)
!
!* (i+1,j,k)
      lai(5)=   i+1
      laj(5)=   j
      lak(5)=   k
       ca(5)=   1.d0/hxsq
!
!* (i,j+1,k)
      lai(6)=   i
      laj(6)=   j+1
      lak(6)=   k
       ca(6)=   1.d0/hysq
!
!* (i,j,k+1)
      lai(7)=   i
      laj(7)=   j
      lak(7)=   k+1
       ca(7)=   1.d0/hzsq
!
!
      lxyz= i +mx*(j +myA*k)
!
      do m= 1,nob3
      ii= pxc(lai(m))            !     periodic
      jj=     laj(m)             !!<-- bound, laj()=0-my
      kk= pzc(lak(m))            !     periodic
!
      aa(lxyz+1,m) = ca(m) 
      na(lxyz+1,m) = ii +mx*(jj +myA*kk) +1
      end do                            ! +
      end do 
!
!***
!
      else if(j.eq.my) then
!
        do i= 0,mx-1
!* (i,j-1,k)
        lai(1)=   i
        laj(1)=   my-1  !<--
        lak(1)=   k
         ca(1)=   1.d0
!
!* (j+1,k)
        lai(2)=   i
        laj(2)=   my
        lak(2)=   k
         ca(2)=   1.d0
!
!
        lxyz= i +mx*(my +myA*k)  !<-- (i,my,k) !!
!
        aa(lxyz+1,1) = ca(1)  ! aa( ,m=1)
        na(lxyz+1,1) = i +mx*(my-1 +myA*k) +1
! 
        aa(lxyz+1,2) = ca(2)  ! aa( ,m=2)
        na(lxyz+1,2) = i +mx*(my   +myA*k) +1  !<-- lxyz+1 !!
! 
        do m= 3,nob3
        aa(lxyz+1,m)= 0
        na(lxyz+1,m)= lxyz +1 
        end do
        end do  
!***
      end if
      end do 
      end do 
!
      return
      end subroutine emcof3
!
!
!******************************************************* 7/15/93 *******
!*    Conjugate residual method  (for non-symmetric matrix).           *
!*      na( ,nobx) is 19 or 7                                          *
!**************************************** parallel version: 11/4/94 ****
!---------------------------------------***-***-------------------------
      subroutine cresmd (b,x,aa,na,nobx,ipr,rpr,ierr)
!------------------------*-R-**********---------------------------------
!   The order of two arguments ss, xx is the same as those of @ftp8.f. 
!   One component is used to solve: A*x = b.
!
!     call cresmd (ierr)
!     common/cresm1/ aa(mxyz,nob2)
!     common/cresma/ na(mxyz,nob2)
!     common/cresmb/ nobx
!     common/cresm2/ b(mxyz),x(mxyz)
!     common/cresm3/ ipr(10),rpr(10)
!
!    arguments :
!
!    aa(mxyz,nob2) (real*4, 8) --------- coefficient matrix
!    na(mxyz,nob2) (integer*4) --------- column no. table
!            na(i,jj) gives a location of the base vector for a(i,jj).
!    note: nob2= 9; but, any nobx .le. nob2 is acceptable.
!
!    b(mxyz)     (real*4, 8) -------------- source vector
!    x(mxyz)     (real*4, 8) -------------- solution (unknown)
!
!    ipr(10) (integer*4) ----------- integer parameter
!                   ipr(1) (in) .... maximum number of iteration
!                   ipr(2) (out) ... final iteration number
!    rpr(10) (real*8) -------------- real parameter
!                   rpr(1) (in) .... tolerance of convergence
!                   rpr(2) (out) ... residual norm
!                   rpr(3) (out) ... cpu time
!    ierr (integer*4) -------------- return code
!                     = 0    ......... normal termination
!                     = 1000 ......... source vector = 0.
!                     = 2000 ......... no convergence within itrm
!                     = 3000 ......... blow up (x > 1.d+50)
!
!                                   Copyright : M.Tanaka  7/15/1993.
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
!     parameter  (nob2=19) or (nob3=7)
! -----------------------------------------------------------------
      real(C_DOUBLE),dimension(mxyz) :: b,x 
      integer(C_INT) ierr
!
      real(C_DOUBLE),dimension(mxyz,nobx) :: aa
      integer(C_INT),dimension(mxyz,nobx) :: na
      integer(C_INT) nobx
!
      integer(C_INT),dimension(10) :: ipr
      real(C_DOUBLE),dimension(10) :: rpr 
!
      real(C_DOUBLE) eps,sq0,sq1,ss,bb
      integer(C_INT) cl_first,ncrs,itrm,syme,i
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!***
!*---------------------------------------------------------------------
!*   Restart procedure is taken in /cressl/ iteration if convergence  *
!*  slows down (ierr= 2000).                                          *
!*                             (also in /wwstbk/)  ** 04/22/1998 **   *
!*---------------------------------------------------------------------
!  Start timing
      real(C_float) t1,t2
!
      cl_first= 1
      call cpu_time (t1)
!
      ncrs= 1
!
 2000 continue
      itrm = ipr(1)
      eps  = rpr(1)
!
!                  + R             
      call cresin (b,x,sq0,sq1,aa,na,nobx,ierr)
!     +++++++++++++++++++++++++++++++++++++++++
!
!  One pass only
      if(ierr.ne.0) go to 1000  !!<-- ierr=1000 for x() the 0 case
!
!                  + R                     
      call cressl (b,x,sq0,sq1,eps,aa,na,nobx,ncrs,itrm,ierr)
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  
!* this can be offline 
!     if(ierr.eq.2000) then
!       if(ncrs.le.2) then
!          syme= -1
!          call filt1l (x,syme,ipar)
!          ierr= 1700 +ncrs
!
!          write(11,*) ' reset and restart (cressl), ncrs=',ncrs
!          go to 2000
!       end if
!     end if
!
 1000 cl_first= 2
      call cpu_time (t2)
!
      ipr(2) = itrm
      rpr(2) = eps
      rpr(3) = t2 -t1
!
      bb= 0
      ss= 0
!
      do i= 1,mxyz 
      bb= bb +b(i)**2
      ss= ss +x(i)**2
      end do
!
      bb= bb/mxyz
      ss= ss/mxyz
!
!     if(iwrt(it,nplot).eq.0 .and. io_pe.eq.1) then
      if(.false.) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("#(cresmd) <bb>,<ss>=",2d12.3,"itrm,eps,ierr=", &
                i5,d12.3,i5)') bb,ss,itrm,eps,ierr
        close(11)
      end if
!
      return
      end subroutine cresmd
!
!
!------------------------*-R-------------------------------------------
      subroutine cresin (b,x,sq0,sq1,aa,na,nobx,ierr)
!----------------------------------------------------------------------
!**********************************************************************
!*    initialization.                                                 *
!**********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
      include 'mpif.h'
!
!     parameter  (nob2=19) or (nob3=7)
!*--------------------------------------------------------------------
      real(C_DOUBLE),dimension(mxyz) :: b,x
      integer(C_INT) ierr
!
      real(C_DOUBLE),dimension(mxyz,nobx) :: aa
      integer(C_INT),dimension(mxyz,nobx) :: na
      integer(C_INT) nobx
!
      real(C_DOUBLE),dimension(mxyz) :: p,q,r,s,p1,p2,q1,q2
      common/cresm4/ p,q,r,s,p1,p2,q1,q2
!
      real(C_DOUBLE) sq0,sq1,swq,sswq 
      integer(C_INT) i
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!*--------------------------------------------------------------------
!
!*  s = a*x  
      call avmult (x,s,aa,na,nobx)
!                  + R        
!
      do i= 1,mxyz 
      p1(i)= 0
      q1(i)= 0
!
      r(i)= b(i) -s(i)
      p(i)= r(i)
      end do
!
!* q = a*r
! 
      call avmult (p,q,aa,na,nobx)
!                  + R     
!
      swq= 0.d0
! 
      do i= 1,mxyz 
      swq= swq + q(i)**2
      end do
!
      if(swq.gt.1.d-30) then
        sq0= swq
        ierr= 0
      else
        sq0= 0
        ierr= 1000
      end if
!
      sq1= 1.d0
!
      return
      end subroutine cresin
!
!
!------------------------*-R-------------------------------------------
      subroutine cressl (b,x,sq0,sq1,eps,aa,na,nobx,ncrs,itrm,ierr)
!----------------------------------------------------------------------
!**********************************************************************
!*    Iteration step by conjugate residual method.                    *
!**********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
      include 'mpif.h'
!
!     parameter  (nob2=19) or (nob3=7)
!*--------------------------------------------------------------------
      real(C_DOUBLE),dimension(mxyz) :: x,b
!
      real(C_DOUBLE),dimension(mxyz,nobx) :: aa
      integer(C_INT),dimension(mxyz,nobx) :: na
      integer(C_INT) nobx
      real(C_DOUBLE) sq0,sq1,sq2,eps
      integer(C_INT) ncrs,itrm,ierr 
!
      real(C_DOUBLE),dimension(mxyz) :: p,q,r,s,p1,p2,q1,q2
      common/cresm4/ p,q,r,s,p1,p2,q1,q2
!*--------------------------------------------------------------------
      real(C_DOUBLE) alpha,beta1,beta2,wb1,wrq,wr1,wxs, &
                     sw1,sw2,wq
      real(C_DOUBLE) srhs,sres,avex,rcc,conv0,conv
      integer(C_INT) i,itr
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!---------------------------------------------------------------------
!*
      wb1= 0.d0
!
!  Start 
      do i= 1,mxyz 
      wb1= wb1 +b(i)**2
      end do
!
      srhs=  sqrt(wb1/float(mxyz))
!
!                        ** iteration starts **
      itr= 0
 1000 itr= itr +1
!
      sq2= sq1
      sq1= sq0
!
      wrq= 0.d0
!
      do i= 1,mxyz 
      p2(i)= p1(i)
      p1(i)=  p(i)
      q2(i)= q1(i)
      q1(i)=  q(i)
!
      wrq= wrq + r(i)*q1(i)
      end do
!
      alpha= wrq/sq1
!
      wr1= 0.d0
      wxs= 0.d0
!
      do i= 1,mxyz 
      x(i)= x(i) +alpha*p1(i)
      r(i)= r(i) -alpha*q1(i)
!
      wr1= wr1 +r(i)**2
      wxs= wxs +x(i)**2
      end do
!
!**********************************************************************
!*    Convergence check.  (compare residual with source)              *
!**********************************************************************
!----------------------------------------------------------------------
!
      sres= sqrt(wr1/float(mxyz))
      avex= sqrt(wxs/float(mxyz))
!
!     conv0= 1.d0 ! conv ??
      conv= sres/(srhs +1.d-15)  ! residue /rhs(b)
!
!     write(11,900) itr,sres,conv,avex
! 900 format(1h ,'#(crm)  itr=',i3,'  res=',1pd12.5,', res/rhs=',
!    *       d12.5,', <x>=',d12.5)
!
!                                 ** successfully converged **
      if(conv.lt.eps) then
         ierr= 0
         go to 7000  !!<--success
      end if
!                                 ** source too small **
      if(srhs.lt.1.d-30) then
         ierr= 1000
         go to 7000
      end if
!                                 ** blow up **
      if(avex.gt.1.d+50) then
        ierr= 3000
        go to 7000
      end if
!----------------------------------------------------------------------
!                                       s = a*r
      call avmult (r,s,aa,na,nobx)
!                  + R 
!
      sw1= 0.d0
      sw2= 0.d0
!
      do i= 1,mxyz 
      sw1= sw1 +s(i)*q1(i)
      sw2= sw2 +s(i)*q2(i)
      end do
!
!
      beta1= -sw1/sq1
      beta2= -sw2/sq2
!
!
      wq= 0.d0
!
      do i= 1,mxyz 
      p(i)= r(i) +beta1*p1(i) +beta2*p2(i)
      q(i)= s(i) +beta1*q1(i) +beta2*q2(i)
!
      wq= wq + q(i)**2
      end do
!
      sq0= wq
!
!-------------------------------------------------------------------
!*  Round-off error is dominating, thus restart iterations.
!--------------------------------------------- added: 04/22/1998 ---
!*  Conditionally: if itr >= itrm, or itr >50 and eps < 1.d-4 
!
!   conv= sres/(srhs +1.d-15) = residue / rhs(b)
!     rcc= abs(conv -conv0)/(conv0 +1.d-10)
!
      if(itr.gt.50 .and. eps.lt.1.d-3) go to 6000
!     if(itr.gt.30 .and. rcc.lt.1.d-4) go to 6000
!     +++++++++++++++++++++++++++++++++++++++++
!                 <-- itr exceeds itrm, go to 6000 
      if(itr.lt.itrm) go to 1000
!
 6000 ierr= 2000
      ncrs= ncrs +1
!
 7000 eps=  conv
      itrm= itr
!
      return
      end subroutine cressl
!
!
!----------------------------------------------------------------------
      subroutine avmult (v,w,aa,na,nobx)
!------------------------*-R------------*********-*******--------------
!**********************************************************************
!     Calculate:  w(i) = sum(j)| a(i,j)*v(na(i,j)).                   *
!**********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
!     parameter  (nob2=19)(nob3=7)
!***                           one array mxyz
      real(C_DOUBLE),dimension(mxyz) :: v,w
!
      real(C_DOUBLE),dimension(mxyz,nobx) :: aa
      integer(C_INT),dimension(mxyz,nobx) :: na
      integer(C_INT) nobx,i,j 
!***
      do i= 1,mxyz 
      w(i)= 0.d0
      end do
!
      do i= 1,mxyz 
      do j= 1,nobx
      w(i)= w(i) + aa(i,j)*v(na(i,j))
      end do
      end do
!
      return
      end subroutine avmult
!
!
!-----------------------------------------------------------------------
      subroutine filt3e (ex,ey,ez,exc,eyc,ezc,ifilx,ifily,ifilz,sym)
!-----------------------------------------------------------------------
!***********************************************************************
!*    << Filtering & boosting for wall-bound system >>                 *
!*     subroutines:  filt3e, filt3b, filt1p                            *
!*                   sym= -1  sym=1  sym= 1 or -1                      *
!*                                                                     *
!*    Successive filtering & boosting with the following weights       *
!*        filtering : ( 0.25, 0.50, 0.25)                              *
!*        boosting  : (-0.25, 0.50,-0.25)                              *
!*    result in  (-1/16, 4/16, 10/16, 4/16, -1/16) weighting.          *
!***********************************************************************
!*   Wall boundary condition in y, a short array (mxyz3/npc) in z      *
!***********************************************************************
!  The array is: ax(-2:mx+1,-1:my+1,-2:mz+1) (real array)
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                                     ex,ey,ez,ax,ay,az
!
      real(C_DOUBLE) exc,eyc,ezc
      integer(C_INT) ifilx,ifily,ifilz,sym
!-----------------------------------------------------------------------     
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) i,j,k,ntz,kr,kl,krr,kll,ntx,ir,il,irr,ill,    &
                     nty,js
!
!-----------------------------------------------------------------------
!*  Filtering & boosting steps. 5-point-average scheme.
!-----------------------------------------------------------------------
!  The array starts with zero: ax(-2:mx+1,-1:my+1,-2:mz+1) 
!  The dc component must be subtracted (restored before return).
!
      do k= 0,mz-1
      do j= 0,my     !<-- the inner points
      do i= 0,mx-1
      ex(i,j,k)= ex(i,j,k) -exc
      ey(i,j,k)= ey(i,j,k) -eyc
      ez(i,j,k)= ez(i,j,k) -ezc
      end do
      end do
      end do
!
!***********************************************************************
!*  Filtering in the z-direction.                                      *
!***********************************************************************
!
      do ntz= 1,ifilz
!*
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ax(i,j,k)= ex(i,j,k)
      ay(i,j,k)= ey(i,j,k)
      az(i,j,k)= ez(i,j,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      kr = pzr(k)
      kl = pzl(k)
      krr= pzr(kr)
      kll= pzl(kl)
!***             folding back by pzr(k) and pzl(k)
      ex(i,j,k)= -0.0625d0*ax(i,j,krr) +0.25d0*ax(i,j,kr) &
                  +0.625d0*ax(i,j,k)   +0.25d0*ax(i,j,kl) -0.0625d0*ax(i,j,kll)
      ey(i,j,k)= -0.0625d0*ay(i,j,krr) +0.25d0*ay(i,j,kr) &
                  +0.625d0*ay(i,j,k)   +0.25d0*ay(i,j,kl) -0.0625d0*ay(i,j,kll)
      ez(i,j,k)= -0.0625d0*az(i,j,krr) +0.25d0*az(i,j,kr) &
                  +0.625d0*az(i,j,k)   +0.25d0*az(i,j,kl) -0.0625d0*az(i,j,kll)
      end do
      end do
      end do
!*
      end do
!
!***********************************************************************
!*  Filtering in the x-direction.                                      *
!***********************************************************************
!
      ntx= 0
 1000 ntx= ntx +1
      if(ntx.gt.ifilx) go to 2000
!*
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ax(i,j,k)= ex(i,j,k)
      ay(i,j,k)= ey(i,j,k)
      az(i,j,k)= ez(i,j,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ir=  pxr(i)
      il=  pxl(i)
      irr= pxr(ir)
      ill= pxl(il)
!                 Folding back by pxr(i) and pxl(i)
      ex(i,j,k)= -0.0625d0*ax(ill,j,k) +0.25d0*ax(il,j,k) +0.625d0*ax(i,j,k) &
                   +0.25d0*ax(ir,j,k) -0.0625d0*ax(irr,j,k)
      ey(i,j,k)= -0.0625d0*ay(ill,j,k) +0.25d0*ay(il,j,k) +0.625d0*ay(i,j,k) &
                   +0.25d0*ay(ir,j,k) -0.0625d0*ay(irr,j,k)
      ez(i,j,k)= -0.0625d0*az(ill,j,k) +0.25d0*az(il,j,k) +0.625d0*az(i,j,k) &
                   +0.25d0*az(ir,j,k) -0.0625d0*az(irr,j,k)
      end do
      end do
      end do
!*
      go to 1000
 2000 continue
!
!  Bound condition in y
!
      nty= 0
 3000 nty= nty +1
      if(nty.gt.ifily) go to 4000
!* 
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ax(i,j,k)= ex(i,j,k)
      ay(i,j,k)= ey(i,j,k)
      az(i,j,k)= ez(i,j,k)
      end do
      end do
      end do
!
!                       sym= -1: (odd, even,odd)   -1 0 (1) 2 3
!                       sym=  1: (even,odd,even)
!***
      js= 1         ! a(-2) <-> e(2), a(-1) <-> e(1)
      do k= 0,mz-1  !  ex() to ax() of opposite sign
      do i= 0,mx-1
      ax(i,-js,k)=  sym* ex(i,js,k)
      ay(i,-js,k)= -sym* ey(i,js,k)
      az(i,-js,k)=  sym* ez(i,js,k)
      end do
      end do
!                       my-2 my-1 (my) my+1 my+2
      js= 1        ! a(my+2) <-> e(my-2), a(my+1) <-> e(my-1)
      do k= 0,mz-1 
      do i= 0,mx-1
      ax(i,my+js,k)=  sym* ex(i,my-js,k)
      ay(i,my+js,k)= -sym* ey(i,my-js,k)
      az(i,my+js,k)=  sym* ez(i,my-js,k)
      end do
      end do
!***
!
      do k= 0,mz-1
      do j= 1,my-1
      do i= 0,mx-1
!                 j-2 -> -1, j+2 -> my+3
      ex(i,j,k)= -0.0625d0*ax(i,j+2,k) +0.25d0*ax(i,j+1,k) &
                  +0.625d0*ax(i,j,k)   +0.25d0*ax(i,j-1,k) &
                 -0.0625d0*ax(i,j-2,k)
      ey(i,j,k)= -0.0625d0*ay(i,j+2,k) +0.25d0*ay(i,j+1,k) &
                  +0.625d0*ay(i,j,k)   +0.25d0*ay(i,j-1,k) &
                 -0.0625d0*ay(i,j-2,k)
      ez(i,j,k)= -0.0625d0*az(i,j+2,k) +0.25d0*az(i,j+1,k) &
                  +0.625d0*az(i,j,k)   +0.25d0*az(i,j-1,k) &
                 -0.0625d0*az(i,j-2,k)
      end do
      end do
      end do
!
      go to 3000
 4000 continue
!
!--------------------------------------
!* the dc component is restored.
!--------------------------------------
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ex(i,j,k)= ex(i,j,k) +exc
      ey(i,j,k)= ey(i,j,k) +eyc
      ez(i,j,k)= ez(i,j,k) +ezc
      end do
      end do
      end do
!
      return
      end subroutine filt3e
!
!
!-----------------------------------------------------------------------
      subroutine filt1p (q,ifilx,ifily,ifilz,sym)
!-----------------------------------------------------------------------
!***********************************************************************
!*    for duplicate local array  (e/m field).                          *
!***********************************************************************
!  The array starts with zero: ax(-2:mx+1,-1:my+1,-2:mz+1) 
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!                              +++++++++++++++++++++++
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: q,a
      integer(C_INT) ifilx,ifily,ifilz,sym
!----------------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) i,j,k,ntz,kr,kl,krr,kll,ntx, &
                     ir,il,irr,ill,nty,js
!
!***********************************************************************
!*    Filter in the z direction.                                       *
!***********************************************************************
!
      do ntz= 1,ifilz
!*
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      kr=  pzr(k)
      kl=  pzl(k)
      krr= pzr(kr)
      kll= pzl(kl)
!                folding back
      q(i,j,k)= -0.0625d0*a(i,j,krr) +0.25d0*a(i,j,kr) +0.625d0*a(i,j,k) &
                  +0.25d0*a(i,j,kl) -0.0625d0*a(i,j,kll)
      end do
      end do
      end do
!*
      end do
!
!
      ntx= 0
 1000 ntx= ntx +1
      if(ntx.gt.ifilx) go to 2000
!*
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      ir = pxr(i)
      il = pxl(i)
      irr= pxr(ir)
      ill= pxl(il)
!                folding back
      q(i,j,k)= -0.0625d0*a(ill,j,k) +0.25d0*a(il,j,k) +0.625d0*a(i,j,k) &
                  +0.25d0*a(ir,j,k) -0.0625d0*a(irr,j,k)
      end do
      end do
      end do
!*
      go to 1000
 2000 continue
!
!  Bound condition in y
!
      nty= 0
 3000 nty= nty +1
      if(nty.gt.ifily) go to 4000
!*
      do k= 0,mz-1
      do j= 0,my
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!***
      js= 1         ! a(-2) <-> e(2), a(-1) <-> e(1)
      do k= 0,mz-1  !  ex() to ax() of opposite sign
      do i= 0,mx-1
      ax(i,-js,k)=  sym* ex(i,js,k)
      end do
      end do
!                       my-2 my-1 (my) my+1 my+2
      js= 1        ! a(my+2) <-> e(my-2), a(my+1) <-> e(my-1)
      do k= 0,mz-1 
      do i= 0,mx-1
      ax(i,my+js,k)=  sym* ex(i,my-js,k)
      end do
      end do
!***
!
      do k= 0,mz-1
      do j= 1,my-1
      do i= 0,mx-1
!                j-2 -> -1, j+2 -> my+3
      q(i,j,k)= -0.0625d0*a(i,j+2,k) +0.25d0*a(i,j+1,k) +0.625d0*a(i,j,k) &
                  +0.25d0*a(i,j-1,k) -0.0625d0*a(i,j-2,k)
      end do
      end do
      end do
!*
      go to 3000
 4000 continue
!
      return
      end subroutine filt1p
!
!
!-----------------------------------------------------------------------
      subroutine diag1 (x,y,z,vx,vy,vz,qmult,npr,ksp)
!-----------------------------------------------------------------------
!  Plot x-vz and/or qix-qez
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
      include 'mpif.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      real(C_DOUBLE) qmult
      integer(C_INT) npr,ksp,MPIerror
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
! 
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                              ex,ey,ez,bx,by,bz,       &
                                              ex0,ey0,ez0,bx0,by0,bz0, &
                                              qix,qiy,qiz,qex,qey,qez, &
                                              emx,emy,emz,qi,qe,amu,avhh
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
!
      real(C_DOUBLE),dimension(np0) :: ggg
!---------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) wkix,wkih,wkex,wkeh
      common/wkinel/ wkix,wkih,wkex,wkeh
!
!     real(C_float),dimension(7) :: vxav,vyav,vzav,vxrmsd,vyrmsd,   &
!                                   vzrmsd,vsqr,anpt        !<-- real*4
!     common/diagpr/ vxav,vyav,vzav,vxrmsd,vyrmsd,vzrmsd,vsqr,anpt
!
      real(C_DOUBLE),dimension(101,7,2) :: fvx,fvy,fvz
      common/diagp2/ fvx,fvy,fvz
!
      integer(C_INT) iwrt,lskip,jj,ix,iy,iz,i,j,k,l,itag,iperio
      real(C_DOUBLE) vlimx,vlimz,delx,delz,x1,x2,radi,  &
                     xcnt,ycnt,zcnt,ddx,ddy,ddz,        &
                     anp,bs
!                             +++++++++
      real(C_float),dimension(0:mx-1,0:my,0:mz-1) :: &
                                        eex,eey,eez,qqi,qqe,qqq
!
!*  ####### diag1 #######
!
      if(iwrt(it,nha).ne.0) return
!            .and. iwrt(it,nplot).ne.0) return
!     ++++++++++++++++++++++++++++++++++++++++++++++++
!
      call lblbot (t)
!
!* [1] ******************************************
      if(iwrt(it,nplot).eq.0 .and. io_pe.eq.1) then
!             ++++++++
      do j=1,7
      do i=1,101
      fvx(i,j,ksp)= 1.d-2
      fvy(i,j,ksp)= 1.d-2
      fvz(i,j,ksp)= 1.d-2
      end do
      end do
!
      if(ksp.eq.1) then
        vlimx= vlima
        vlimz= vlima
        delx= 0.04d0*vlima
        delz= delx
!
      else
        vlimx= vlimb
        vlimz= vlimb
        delx= 0.04d0*vlimb
        delz= delx
      end if
!
      x1= 0.45d0*ymax ! ycent1 of pplot3
      x2= 0.55d0*ymax ! ycent2
      radi= rwd
      lskip=  npr/8000.
!
!-----------------------------------------------------------------------
!*  **Plasma current**
!-----------------------------------------------------------------------
!
!     if(iwrt(it,nplot).eq.0 .and. io_pe.eq.1) then
      if(.false.) then
!
      if(ksp.eq.1) then
!
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        eex(i,j,k)= qix(i,j,k)
        eey(i,j,k)= qiy(i,j,k)
        eez(i,j,k)= qiz(i,j,k)
        qqi(i,j,k)=  qi(i,j,k)
        end do
        end do
        end do
!
!  (qix,qiy,qiz) real*4 
!
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        iperio= 0
        itag= 1
!       call fplot3 (eex,eey,eez,0.,0.,0.,real(xmax),real(ymax),real(zmax), &
!                    iperio,itag,'Ji.total',8)
!       call cplot3 (qqi,real(xmax),real(ymax),real(zmax),'ion.dens',8)
!       call pplot3 (x,y,z,vx,vy,vz,real(xmax),real(ymax),real(zmax),  &
!                    real(x1),real(x2),real(radi),npr,lskip,itag,      &
!                    'ions    ',8)
        close(77)
! 
      else if(ksp.eq.2) then
!
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        eex(i,j,k)= qix(i,j,k) +qex(i,j,k)
        eey(i,j,k)= qiy(i,j,k) +qey(i,j,k)
        eez(i,j,k)= qiz(i,j,k) +qez(i,j,k)
        qqq(i,j,k)=  qi(i,j,k) +qe(i,j,k)
        end do
        end do
        end do
!
!  (qqx,qqy,qqz) real*4
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        iperio= 0
        itag= 2
!       call fplot3 (eex,eey,eez,0.,0.,0.,real(xmax),real(ymax),real(zmax), &
!                    iperio,itag,'J.total ',8)
!       call cplot3 (qqq,real(xmax),real(ymax),real(zmax),'el.dens ',8)
!       call pplot3 (x,y,z,vx,vy,vz,real(xmax),real(ymax),real(zmax),   &
!                    real(x1),real(x2),real(radi),npr,lskip,itag,       &
!                    'elec    ',8)
        close(77)
      end if
      end if
!***
      end if
!* [1] ******************************************
!
!-----------------------------------------------------------------------
!*  Temperature.
!-----------------------------------------------------------------------
!
!     if(mod(it,nha).eq.0 .and. io_pe.eq.1) then
!       open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
!             status='unknown',position='append',form='formatted')
!
!       if(ksp.eq.1) then
!         call temper (x,y,z,vx,vy,vz,qmult,npr,th,tx,1)
!
!         call cplot3 (th,xmax,ymax,zmax,11,'tpara.i ',8)
!         call cplot3 (tx,xmax,ymax,zmax,12,'tperp.i ',8)
!       else
!         call temper (x,y,z,vx,vy,vz,qmult,npr,th,tx,2)
!         call cplot3 (th,xmax,ymax,zmax,13,'tpara.el',8)
!       end if
!
!       close(77)
!     end if
!
!-----------------------------------------------------------------------
!   Velocity distribution function.
!-----------------------------------------------------------------------
!* [2] *********************************************
      if(iwrt(it,nha).eq.0) then
!                +++
!*
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        if(ksp.eq.1) then
          write(11,'("t=",f7.1," wkix,wkih=",1p2d12.3)') t,wkix,wkih
        else
          write(11,'("t=",f7.1," wkex,wkeh=",1p2d12.3)') t,wkex,wkeh
        end if
!
        close(11)
      end if
!
!  Reset ggg()=0 and reduce the size intervals to all components
!   Sums of x, y, z, vx,vy, vz
! [1]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (x(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      x(l)= ggg(l)
      end do
! [2]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (y(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      y(l)= ggg(l)
      end do
! [3]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (z(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      z(l)= ggg(l)
      end do
! [4]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (vx(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      vx(l)= ggg(l)
      end do
! [5]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (vy(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      vy(l)= ggg(l)
      end do
! [6]
      do l= 1,np0
      ggg(l)= 0
      end do
!
      call mpi_allreduce (vz(1),ggg(1),np0,mpi_real8,mpi_sum,  & 
                          mpi_comm_world,MPIerror)
      do l= 1,np0
      vz(l)= ggg(l)
      end do
!***
!
      if(io_pe.eq.1) then
!        ++++++++++
!
      ycnt= 0.5*ymax
      zcnt= 0.5*zmax
      ddy= 15.
      ddz= zmax/6. 
!
!*
      do l= 1,np0
      if(abs(y(l)-ycnt).gt.ddy) go to 332
      if(abs(z(l)-zcnt).gt.ddz) go to 332  
!
      jj= (z(l) -zcnt)/ddz +3.51           ! z slices
!        1,2,3, 4, 5,6,7
      if((jj.ge.1).and.(jj.le.7)) then
!
      ix= (vx(l) +vlimx)/delx +1.0
      iy= (vy(l) +vlimz)/delz +1.0
      iz= (vz(l) +vlimz)/delz +1.0
!
      if(iabs(ix-26).gt.50)  go to 333
        fvx(ix,jj,ksp)= fvx(ix,jj,ksp) +qmult
!
  333 if(iabs(iy-26).gt.50)  go to 334
        fvy(iy,jj,ksp)= fvy(iy,jj,ksp) +qmult
!
  334 if(iabs(iz-26).gt.50)  go to 332
        fvz(iz,jj,ksp)= fvz(iz,jj,ksp) +qmult
      end if
!
  332 continue
      end do 
!
      end if  ! if(io_pe.eq.1) then
      end if  ! if(iwrt(it,nha).eq.0) then
!
!* [2] *********************************************
!
      return
      end subroutine diag1
!
!
!-----------------------------------------------------------------------
      subroutine fvplot 
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(101,7,2) :: fvx,fvy,fvz
      common/diagp2/ fvx,fvy,fvz
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE),dimension(101) :: vsa,vsb
      real(C_float)  fimax,fimin,femax,femin
      integer(C_INT) k,ILN,iwrt
!
      if(iwrt(it,nplot).eq.0 .and. io_pe.eq.1) then
!**
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',        &
              status='unknown',position='append',form='unformatted')
!
!*  fv-plots for the zones y3 - y4  (total of 7 zones)
!
        do k= 1,101
        vsa(k)= vlima*(-1.d0 +(k-1)/50.d0) 
        vsb(k)= vlimb*(-1.d0 +(k-1)/50.d0) 
        end do
!
!     call lplot(-2,4,ldec,tdec,edec(1,1),ILN,elab(1),8, &
!     subroutine lplot (ix,iy,npts,xsc,q,IL,lab1,n1,lab2,n2, & 
!
        ILN= 1
        call lplmax3 (fvx(1,3,1),fvy(1,3,1),fvz(1,3,1),fimax,fimin,101)
        call lplmax3 (fvx(1,3,2),fvy(1,3,2),fvz(1,3,2),femax,femin,101)
!
        call hplot (2,4,101,vsa,fvx(1,3,1),0.,fimax,ILN,'io.5/8  ',8, &
                    '   vx   ',8,'        ',8)
        call hplot (2,5,101,vsa,fvy(1,3,1),0.,fimax,ILN,'io.5/8  ',8, &
                    '   vy   ',8,'        ',8)
        call hplot (2,6,101,vsa,fvz(1,3,1),0.,fimax,ILN,'io.5/8  ',8, &
                    '   vz   ',8,'        ',8)
!
        call hplot (3,4,101,vsb,fvx(1,3,2),0.,femax,ILN,'el.5/8  ',8, &
                    '   mu   ',8,'        ',8)
        call hplot (3,5,101,vsb,fvz(1,3,2),0.,femax,ILN,'el.5/8  ',8, &
                    '   vh   ',8,'        ',8)
        call chart
!
        call hplot (2,4,101,vsa,fvx(1,4,1),0.,fimax,ILN,'io.6/8  ',8, &
                    '   vx   ',8,'        ',8)
        call hplot (2,5,101,vsa,fvy(1,4,1),0.,fimax,ILN,'io.6/8  ',8, &
                    '   vy   ',8,'        ',8)
        call hplot (2,6,101,vsa,fvz(1,4,1),0.,fimax,ILN,'io.6/8  ',8, &
                    '   vz   ',8,'        ',8)
!
        call hplot (3,4,101,vsb,fvx(1,4,2),0.,femax,ILN,'el.6/8  ',8, &
                    '   mu   ',8,'        ',8)
        call hplot (3,5,101,vsb,fvz(1,4,2),0.,femax,ILN,'el.6/8  ',8, &
                    '   vh   ',8,'        ',8)
        call chart
!
        close(77)
      end if
!
!**
      if(iwrt(it,nha).eq.0 .and. io_pe.eq.1) then
!
        if(ldec.ge.3000) then
          call rehist
        end if
      end if
!
      return
      end subroutine fvplot
!
!
!-----------------------------------------------------------------------
      subroutine temper (x,y,z,vx,vy,vz,qmult,npr,th,tp,ksp)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                                    th,tp,vh,vp,an
      real(C_DOUBLE) qmult
      integer(C_INT) npr,ksp 
! 
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                              ex,ey,ez,bx,by,bz,    &
                                              ex0,ey0,ez0,bx0,by0,bz0
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
!-------------------------------------------------------------
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      real(C_DOUBLE) bxa,bya,bza,bs,vpara,vpx,ran
      integer(C_INT) ip,jp,kp,i,j,k,l
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
!
      if(io_pe.ne.1) return
!***
      if(ksp.eq.1) then
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        th(i,j,k)= 0
        tp(i,j,k)= 0
        vh(i,j,k)= 0
        vp(i,j,k)= 0
        an(i,j,k)= 0
        end do
        end do
        end do
!
        do l= 1,npr 
        ip= int(hxi*x(l) +0.500000001d0)  ! il, ip, ir
        jp= int(hyi*y(l) +0.000000001d0)  !  jp, jp+1
        kp= int(hzi*z(l) +0.500000001d0)  ! kl, kp, kr
!
        i = pxc(ip)  !<--???
        j = jp
        k = pzc(kp)
!
        bxa= bx(i,j,k) +bxc
        bya= by(i,j,k) +byc
        bza= bz(i,j,k) +bzc
        bs= sqrt(bxa**2 +bya**2 +bza**2)
!
        vpara= (vx(l)*bxa +vy(l)*bya +vz(l)*bza)/bs
        vpx= vx(l)
!
        vh(i,j,k)= vh(i,j,k) +qmult*vpara
        vp(i,j,k)= vp(i,j,k) +qmult*vpx
        th(i,j,k)= th(i,j,k) +qmult*vpara**2
        tp(i,j,k)= tp(i,j,k) +qmult*vpx**2
        an(i,j,k)= an(i,j,k) +qmult
        end do 
!
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        if(an(i,j,k).gt.0.3) then
          ran= 1.d0/an(i,j,k)
        else
          ran= 0
        end if
!
        vh(i,j,k)= ran*vh(i,j,k)
        vp(i,j,k)= ran*vp(i,j,k)
        th(i,j,k)= ran*th(i,j,k)
        tp(i,j,k)= ran*tp(i,j,k)
!
        th(i,j,k)= sqrt(th(i,j,k) -vh(i,j,k)**2)
        tp(i,j,k)= sqrt(tp(i,j,k) -vp(i,j,k)**2)
        end do
        end do
        end do
!
      else
!
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        th(i,j,k)= 0
        vh(i,j,k)= 0
        an(i,j,k)= 0
        end do
        end do
        end do
!
        do l= 1,npr 
        ip= int(hxi*x(l) +0.500000001d0)  ! il, ip, ir
        jp= int(hyi*y(l) +0.000000001d0)  !  jp, jp+1
        kp= int(hzi*z(l) +0.500000001d0)  ! kl, kp, kr
!
        i = pxc(ip)
        j = jp
        k = pzc(kp)
!
        vh(i,j,k)= vh(i,j,k) +qmult*vz(l)
        th(i,j,k)= th(i,j,k) +qmult*vz(l)**2
        an(i,j,k)= an(i,j,k) +qmult
        end do
!
        do k= 0,mz-1
        do j= 0,my
        do i= 0,mx-1
        if(an(i,j,k).gt.0.3) then
          ran= 1.d0/an(i,j,k)
        else
          ran= 0
        end if
!
        vh(i,j,k)= ran*vh(i,j,k)
        th(i,j,k)= ran*th(i,j,k)
!
        th(i,j,k)= sqrt(th(i,j,k) -vh(i,j,k)**2)
        end do
        end do
        end do
      end if
!
      return
      end subroutine temper
!
!
!-----------------------------------------------------------------------
      subroutine init (xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti,  &  
                       xe,ye,ze,vxe,vye,vze,qmulte,wmulte,  &
                       npr,kstart)
!-----------------------------------------------------------------------
!***********************************************************************
!*    At initial particle loading, the Grad-Shafranov solution may be  *
!*    plugged in via ft19f001 if iloadp= 1.                            *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,vxi,vyi,vzi, &
                                       xe,ye,ze,vxe,vye,vze
      real(C_DOUBLE) qmulti,qmulte,wmulti,wmulte
      integer(C_INT) npr,kstart
!------------------------------------------------------------------- 
      real(C_DOUBLE),dimension(mxyzA) :: ex,ey,ez,bx,by,bz,       &
                                         ex0,ey0,ez0,bx0,by0,bz0, &
                                         qix,qiy,qiz,qex,qey,qez, &
                                         emx,emy,emz,qi,qe,amu,avhh, &
                                         pot,gnu
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/sclpot/ pot
      common/dragcf/ gnu
!
      real(C_DOUBLE),dimension(0:my) :: filt
      common/bpfilt/ filt
!------------------------------------------------------------------- 
!***
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!***
      integer(C_INT) i,j,k,l,ijk
      real(C_DOUBLE) vith,wg  !thb1,sn,cs
!
!***********************************************************************
!* (1) Following index tables are used to assign grids to particels.   *
!***********************************************************************
!*  filter for boundary points.
!
      do j= 0,my 
      filt(j)= 1.d0
      end do
!
!       if(.false.) then  ! 8/01
!     filt(0)= 0
!     filt(1)= 0.25d0
!     filt(2)= 0.50d0
!     filt(3)= 0.75d0
!
!     filt(my-3)= 0.75d0
!     filt(my-2)= 0.50d0
!     filt(my-1)= 0.25d0
!     filt(my  )= 0
!       end if
!
!-----------------------------------------------------------------------
!* (2) Particle index: x-direction
!-----------------------------------------------------------------------
!* Periodic:  base is pxc(0)  3/01/2022
!
      do i= 0,mx-1  !<--
      pxl(i)= i-1
      pxc(i)= i 
      pxr(i)= i+1
      end do
!                      ! particle table - cyclic
      pxl(0)= mx-1     !<-- pxl(0)= -1 or mx-1
      pxr(mx-1)= 0     !<-- pxr(mx-1)= mx or 0
!
      do i= -mx,-1
      pxl(i)= i-1 +mx
      pxc(i)= i   +mx  !<-- pxc(-1)= 0 or mx-1
      pxr(i)= i+1 +mx
      end do
!
      pxl(-mx)= mx-1
      pxr(-1) = 0
!
      do i= mx,mx+mx-1
      pxl(i)= i-1 -mx
      pxc(i)= i   -mx
      pxr(i)= i+1 -mx
      end do
!
      pxl(mx)=      mx-1
      pxr(mx+mx-1)= 0
!
!-----------------------------------------------------------------------
!*    y-direction
!-----------------------------------------------------------------------
!* Bound boundary: pyc= 0 - my
!
!     do j= -5,-1
      pyl(-1)= -1
      pyc(-1)= -1
      pyr(-1)= 0
!     end do
!
      do j= 0,my
      pyl(j)= j-1
      pyc(j)= j
      pyr(j)= j+1
      end do
!
      pyl(0)=  0
      pyr(my)= my
! 
!     do j= my+1,my+5
      pyl(my+1)= my
      pyc(my+1)= my+1
      pyr(my+1)= my+1
!     end do
!
!-----------------------------------------------------------------------
!*    z-direction
!-----------------------------------------------------------------------
!* Periodic (so far)
!
      do k= 0,mz-1
      pzl(k)= k-1
      pzc(k)= k
      pzr(k)= k+1
      end do
!
      pzl(0)= mz-1
      pzr(mz-1)= 0
!
      do k= -mz,-1
      pzl(k)= k-1 +mz
      pzc(k)= k   +mz
      pzr(k)= k+1 +mz
      end do
!
      pzl(-mz)= mz-1
      pzr(-1) = 0
!
      do k= mz,mz+mz-1
      pzl(k)= k-1 -mz
      pzc(k)= k   -mz
      pzr(k)= k+1 -mz
      end do
!
      pzl(mz)=      mz-1
      pzr(mz+mz-1)= 0
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) ' The table for i= -1,0,1,...,mx-1,mx,mx+1 '
        write(11,*) '               j= -1,0,1,...,my-1,my,my+1'
        write(11,*) '               k= -1,0,1,...,mz-1,mz,mz+1 '
        write(11,*) ' i, pxc,  j, pyc,  k, pzc.....'
        write(11,701) (i,pxc(i),i,pyc(i),i,pzc(i),i=-1,mx+1)
        write(11,703) (i,pyc(i),i=mx+2,my+1)   ! my+1
  701   format(' pxc(',i3,')=',i4,' pyc(',i3,')=',i4,  &
               ' pzc(',i3,')=',i4)
  703   format(14x,' pyc(',i3,')=',i4)
        close(11)
      end if
!
!***********************************************************************
!* (3) Define uneven mesh quantities.                                  *
!***********************************************************************
!-----------------------------------------------------------------------
!*    Generate the grid spacing first.
!     /caution/  gx(mx+1)= gx(mx) +hx must be defined for pcx-def.
!                this is also used to get hx(mx)= gx(mx+1) -gx(mx).
!-----------------------------------------------------------------------
!  *** x-direction ***
!* Periodic
!
      hx= xmax/mx  !<--
!*
      wg= -hx/2    !<--
      gx(0)= wg
      gx(-1)= gx(0) -hx
!
      do i= 1,mx+2
      wg= wg +hx
      gx(i)= wg
      end do
!
!  *** y-direction ***
!
      hy= ymax/my  !<-- yidirection are 0,..,my
!*
      wg= 0.d0     !<-- 
      gy(0)= wg
      gy(-1)= gy(0) -hy
!
      do j= 1,my
      wg= wg +hy
      gy(j)= wg 
      end do
!
      gy(my+1)= ymax +hy
! 
!  *** z-direction ***
!* Periodic (so far)
!*       define the guard cell (mz+1).
!
      hz= zmax/mz
!*
      wg= -hz/2     !<--
      gz(0)= wg     !<-- initial point
      gz(-1)= gz(0) -hz
!
      do k= 1,mz+2
      wg= wg +hz
      gz(k)= wg      !<-- inner points: (1,mz) range: -hz/2,zmax-hz/2
      end do         !!<-- guard point: (0,mz+1) range: 0,hz*mz
!
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,290)
        write(11,*) 'inner cells: mx,my+1,mz=',mx,my+1,mz
        write(11,291) hx,hy,hz
        write(11,292)
  290   format(//)
  291   format('# hx=',1pd13.5,',  hy=',d13.5,',  hz=',d13.5,/)
  292   format(' i,     gx,      hx;     gy,    hy;     gz,    hz',/)
!
        write(11,293) (i,gx(i),gy(i),gz(i),i=-1,0) 
        write(11,293) (i,gx(i),gy(i),gz(i),i=1,mx+1) 
        write(11,297) (i,gy(i),i=mx+2,my+1)
  293   format(i5,1p3d12.3)
  297   format(i5,12x,1pd12.3)
!
        close(11)
      end if
!
!***********************************************************************
!* (4) Particle position ---> grid assignment table.                   *
!***********************************************************************
!-----------------------------------------------------------------------
!*  1) find the left-side grid.
!-----------------------------------------------------------------------
!     i0 ............... field grid number.
!     ii of ptxl(ii) .... particle grid number.
!
!     i0 = 1
!     ptxl(1)= 1
!     shx= 1.d-2*hx
!
!     do 400 ii= 2,2000
!     x0= hx*(ii-1)
!     if(x0.gt.xmax) go to 400
!
!     if(x0.le.(gx(i0+1)-shx)) then
!        ptxl(ii)= i0
!     else
!        i0= i0 +1
!        ptxl(ii)= i0
!     end if
! 400 continue
!
!     do 420 ii= 1,2000
!     if(ptxl(ii).gt.mx) ptxl(ii)= 1
! 420 continue
!
!-----------------------------------------------------------------------
!*  (5) Find the nearest grids ptx, pty and ptz for x, y and z
!-----------------------------------------------------------------------
!*   Periodic cycle.
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'(" hx=",1pd15.7,",  hy=",d15.7,",  hz=",d15.7,/)') &
                                                          hx,hy,hz
      end if
!
!***********************************************************************
!* (6) Define constants in /parm2/.                                    *
!***********************************************************************
!
      hx2 = 2.d0*hx
      hy2 = 2.d0*hy
      hz2 = 2.d0*hz
!
      hxi= 0.9999999999999d0/hx
      hyi= 0.9999999999999d0/hy
      hzi= 0.9999999999999d0/hz
!
      hxsq= hx**2
      hysq= hy**2
      hzsq= hz**2
!
      xmaxe= 0.9999999999999d0*xmax
      ymaxe= 0.9999999999999d0*ymax
      zmaxe= 0.9999999999999d0*zmax
!
      adt  = aimpl*dt
      hdt  = 0.5d0*dt
!
      adtsq= adt**2 
      ahdt2= adt*hdt
!
!-----------------------------------------------------------------------
!  thermal speeds & constant magnetic field.
!-----------------------------------------------------------------------
!
      vith= veth/sqrt(te_by_ti*wspec(1))
      vthx(1)= vith
      vthz(1)= vith
!
      vthx(2)= veth
      vthz(2)= veth
!
!     thb1= pi*thb/180.d0
!     sn = sin(thb1)
!     cs = cos(thb1)
!
!  >> bxc is included here !!
      bxc= wce_by_wpe !! wce_by_wpe is given in rec3d
      byc= 0      !! (By, Bz)
      bzc= 0
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'wce/wpe=',wce_by_wpe,' in the x direction.'
        write(11,*) '  bxc,byc,bzc=',bxc,byc,bzc
        close(11)
      end if
!
!***********************************************************************
      if(kstart.ne.0) return
!***********************************************************************
!-----------------------------------------------------------------------
!*    the initial fields  ex - bz, ex0 - bz0 are re-defined
!*    in emfld0  (called in trans).
!-----------------------------------------------------------------------
!     to avoid overflow in movers at t= 0, needs some fields.
!
      do ijk= 1,mxyzA  !<-- All points included extended ones
      ex(ijk)= 0
      ey(ijk)= 0
      ez(ijk)= 0
!
      ex0(ijk)= 0
      ey0(ijk)= 0
      ez0(ijk)= 0
!
      bx(ijk)= 0
      by(ijk)= 0
      bz(ijk)= 0
!
      bx0(ijk)= 0
      by0(ijk)= 0
      bz0(ijk)= 0
      end do
!
!***********************************************************************
!* (7) Particle loading: if iloadp= 1, then read-in prescribed data.   *
!******************************************************** 11/10/90 *****
!* use same random seed for ions and electrons.  09-21-1998
!
      if(iloadp.eq.0) then  !! particle loading is executed here
!  Arrays here are continous
!
         call rantbl
         call loadpt (xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti,npr,1)
         nspec(1)= npr
!
         call rantbl
         call loadpt (xe,ye,ze,vxe,vye,vze,qmulte,wmulte,npr,2)
         nspec(2)= npr
!                              +++++++++++ xyz
!
         call ipleql (xi,yi,zi,xe,ye,ze,npr)
         q0= abs(np0)/float(mxyz0)  !<-- density
!
      else
         call readpt (xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti,  &
                      xe,ye,ze,vxe,vye,vze,qmulte,wmulte,  &
                      npr)
      end if
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'(" >> density  q0= ",1pd12.3," <<")') q0
!
        write(11,*) 'total number of ions  npr=',npr
        do l= 1,npr
        if(l.le.10) then
          write(11,'(i8,6d12.3)') &
                          l,xi(l),yi(l),zi(l),vxi(l),vyi(l),vzi(l)
        end if 
        end do
 
        write(11,*) 'total number of elec.  npr=',npr
        do l= 1,npr
        if(l.le.10) then
          write(11,'(i8,6d12.3)') &
                          l,xe(l),ye(l),ze(l),vxe(l),vye(l),vze(l)
        end if 
        end do
 
        close(11)
      end if
!
!-----------------------------------------------------------------------
!*    constants used in the field solution.
!-----------------------------------------------------------------------
!*  <<note>>  q0  is defined in loadpt  (just above).
!
      qi0=  qspec(1)*q0   !<<- not used
      qe0=  qspec(2)*q0   !<<- not used, qe0 < 0
      qwi=  qspec(1)/wspec(1)
      qwe=  qspec(2)/wspec(2)
      qqwi= qspec(1)*qwi
      qqwe= qspec(2)*qwe  !<<-  qspec(2)**2/wspec(2) > 0
!
      aqi0= aimpl*qi0     !<<- not used
      aqe0= aimpl*qe0
      aqwi= aimpl*qwi
      aqwe= aimpl*qwe
!
!-----------------------------------------------------------------------
!*    zero-flush the field arrays  (since not all areas are defined).
!-----------------------------------------------------------------------
!
      do ijk= 1,mxyzA 
      qix(ijk)= 0
      qiy(ijk)= 0
      qiz(ijk)= 0
      qex(ijk)= 0
      qey(ijk)= 0
      qez(ijk)= 0
       qi(ijk)= 0
       qe(ijk)= 0
      amu(ijk)= 0
      gnu(ijk)= 0
      end do
!
      do ijk= 1,mxyzA 
      pot(ijk)= 0
      end do
!
      return
      end subroutine init
!
!
!-----------------------------------------------------------------------
      subroutine loadpt (x,y,z,vx,vy,vz,qmult,wmult,npr,ksp)
!-----------------------------------------------------------------------
!***********************************************************************
!*    Particle loading in Cartesian coordinate                         *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_080A.h'
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
      real(C_DOUBLE) qmult,wmult
      integer(C_INT) npr,ksp
!------------------------------------------------------------
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_DOUBLE) gx,gy,gz,hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq, &
                     dx,dy,dz
      common/ptable/ gx(-1:mx+2),gy(-1:my+1),gz(-1:mz+2),      &
                     hx,hx2,hxsq,hy,hy2,hysq,hz,hz2,hzsq,dx,dy,dz
! 
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!***
      real(C_DOUBLE) arb,zcent,ycent1,ycent2,Ez00,vrg1
      common/profl/  arb,zcent,ycent1,ycent2,Ez00
      common/vring/  vrg1
!
      real(C_DOUBLE) fv1(101),fv2(101),fdr(101)  !! real*8
      real(C_DOUBLE) radi,ddx,xx1,xx,s,sdx,             &!pi2,th,ph
                     dv,v1,dv1,vv,v2,dv2,sdv,eps,y1,y2, &
                     x2,vmag,vxout,vyout,vzout,         &
                     pi2,vdi,vde,dzcent,dzsmt,dycent,   &
                     ycnt1,ycnt2,rry,rrz,vdrift,        &
                     rwd2,anbar,                        &
                     ranfp,funr,fun2,ranf,vx1,vy1,vz1
      integer(C_INT) ns,i,j,k,l,m,k2,nn,npr1,npr2
!
!-----------------------------------------------------------------------
!   Loading particle positions and velosities:
!               funr..... to specify the density profile.
!-----------------------------------------------------------------------
!
      pi2= 2.d0*pi
!
!* funr() profile
      radi= min(ymax,zmax)  !! min() <-- use ymax, zmax
      ddx=  radi/100.d0
!
      fdr(1)= 0
      xx1= 0
      xx= xx1
!
      do j=1,100 
      s= 0
      ns= 1000
      k2= ns/2
      sdx= ddx/float(ns)
!
      do k=1,k2  ! 303
      xx= xx +2.d0*sdx
      s= s +4.d0*funr(xx-sdx) +2.d0*funr(xx)
      end do
!
      s= (s +4.d0*funr(xx+sdx) +funr(xx+2.d0*sdx))*sdx/3.d0
      fdr(j+1)= fdr(j) +s
      end do
!
      do j= 1,100
      fdr(j)= fdr(j)/fdr(101)
      end do
      fdr(101)= 1.d0
!
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Loadpt for ksp=',ksp
        write(11,*) '  vthx,vthz,vbeam,vdr(ksp)=',  &
                    vthx(ksp),vthz(ksp),vbeam(ksp),vdr(ksp)
        close(11)
      end if
!
!-----------------------------------------------------------------------
!*  fun2:  used for the parallel distribution (better shape).
!-----------------------------------------------------------------------
!
      fv1(1)= 0
      fv2(1)= 0
!
      vrg1= 0
!
      vv= 0.d0
      dv= 6.d0/100.d0
!
      v1= vv*vthx(ksp)
      dv1=dv*vthx(ksp)
!
      do j= 1,100
      s= 0
      ns= 1000
      k2= ns/2
      sdv=dv/float(ns)
!
      do k=1,k2
      vv= vv +2.d0*sdv
      s= s +4.d0*fun2(vv-sdv) +2.d0*fun2(vv)
      end do
!
      s= (s +4.d0*fun2(vv+sdv) +fun2(vv+2.d0*sdv))*sdv/3.d0
      fv1(j+1)= fv1(j) +s
      end do
!
      do j=1,101
      fv1(j)= fv1(j)/fv1(101)
      end do
!
!-----------------------------------------------------------------------
!     The perpendicular velocity (vy,vz).
!*    fun2 ..... ring distribution if vdr(ksp) > 0.
!-----------------------------------------------------------------------
!
      vrg1= vdr(ksp)/vthx(ksp)
      vv= max(-3.d0,-vrg1)  !! max()
      dv= (3.d0-vv)/100.d0
!
      v2 = vv*vthx(ksp)
      dv2= dv*vthx(ksp)
!
      do j= 1,100
      s= 0.d0
      ns= 1000
      k2= ns/2
      sdv=dv/float(ns)
!
      do k=1,k2
      vv= vv +2.d0*sdv
      s= s +4.d0*fun2(vv-sdv) +2.d0*fun2(vv)
      end do
!
      s= (s +4.d0*fun2(vv+sdv) +fun2(vv+2.d0*sdv))*sdv/3.d0
      fv2(j+1)= fv2(j) +s
      end do
!
      do j=1,101
      fv2(j)= fv2(j)/fv2(101)
      end do
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'("** Distribution fdr, fv1, fv2 (accumulated) ** ",/)')
!
!       do j= 1,101
!       write(11,662) j,fdr(j),fv1(j),fv2(j)
! 662   format(i6,3d12.3)
!       end do
        close(11)
      end if
!
!***********************************************************************
!*    Particles are loaded in pairs of 16.                             *
!***********************************************************************
!*    use equally-spaced eps for initial particle positioning
!*    (inhomogeneous case)
!-------------------------
!* ksp=1 and ksp=2 
!-------------------------
!     wg= -hx/2 gx(0)= wg
!     wg= 0.d0  gy(0)= wg
!***
      l= 0
!     ++++
      do k= 1,mz 
      do j= 1,my  !!!<-- use inside the region
      do i= 1,mx
!
      do m= 1,32  !! 32 particles per cube
      l= l +1
!
!     xout= gx(i-1) +hx/2 +hx*ranfp(0.d0)  !<-- gx(0)= -hx/2 !!
!     yout= gy(j-1)       +hy*ranfp(0.d0)  !    gx(mx-1)= hx*mx -hx/2 
!     zout= gz(k-1) +hz/2 +hz*ranfp(0.d0)  !     the rightside border
      x(l)= xmax*ranfp(0.d0) -hx/2
      y(l)= ymax*ranfp(0.d0)
      z(l)= zmax*ranfp(0.d0) -hz/2
      end do
!
      end do
      end do
      end do 
!
!     ++++++
      npr= l
!     ++++++
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'init of npr= ',npr
        write(11,*) 'init of ksp= ',ksp
        write(11,*) 'The number of particles= ',l
        write(11,*) ' qmulti=',qspec(ksp)
        write(11,*)
        close(11)
      end if
!
!***********************************************************************
!*    Velocity (vx,vy,vz) is given three-dimensionally.                *
!***********************************************************************
!
      do l= 1,npr
      eps= ranf(0.d0)
!
      do k= 1,100  ! 750
      k2= k
      if(fv2(k).gt.eps) go to 760
      end do
!
  760 y1= fv2(k2-1)
      y2= fv2(k2)
      x2= (eps-y2)/(y2-y1)+k2
!
      vmag= v2 +dv2*(x2-1.d0) +vdr(ksp)
!
!     ph=  pi2*ranf(0.d0)
!     th=   pi*ranf(0.d0)
!
      vxout= vmag*(ranf(0.d0) -0.5d0)  !sin(th)
      vyout= vmag*(ranf(0.d0) -0.5d0)  !cos(th)*cos(ph) 
      vzout= vmag*(ranf(0.d0) -0.5d0)  !cos(th)*sin(ph)  
!
!--------------------------------
!*  Give initial drifts.
!--------------------------------
!
      zcent=  0.50d0 *zmax !! one cycle
      dzcent= 0.125d0*zmax !!
      dzsmt = 0.15d0 *zmax !! larger half the region
!
      ycent1= 0.30d0*ymax  !! two cycles
      ycent2= 0.70d0*ymax  !!
!     ycent1= 0.25d0*ymax  !! two cycles
!     ycent2= 0.75d0*ymax  !!
      dycent= 0.05d0*ymax  !!
!
! The current shape in (y,z)
      if(abs(z(l) -zcent).lt.dzcent) then  ! |z-0.5*zmax| < 0.125*zcent
        ycnt1= ycent1 +dycent  ! ycnt1= (0.25+0.05)*ymax 
        ycnt2= ycent2 -dycent  ! ycnt2= (0.25-0.05)*ymax
      else 
        if(abs(z(l) -zcent).lt.dzsmt) then  ! |z-0.5*zmax| < 0.15*zmax
          ycnt1= ycent1 +dycent !*(1.d0 -abs(z(l)-zcent)/dzsmt) 
          ycnt2= ycent2 -dycent !*(1.d0 -abs(z(l)-zcent)/dzsmt) 
        else                                ! else
          ycnt1= ycent1        ! ycnt1= 0.25*ymax
          ycnt2= ycent2        ! ycnt2= 0.75*ymax
        end if
      end if
!
! The vx drift for ksp=1 and 2
      rrz= 0.25d0 *zmax 
      rry= 0.075d0*ymax 
!
      vdrift= 0
!                               ycnt1, or ycnt2
      if((abs(z(l) -zcent).le.rrz) .and. &
         (abs(y(l) -ycnt1).le.rry  .or. abs(y(l) -ycnt2).le.rry)) then
!
          vdrift= vbeam(ksp)
      end if
!
      vx(l)=  vxout + vdrift  !<-- <vx>= vdrift
      vy(l)=  vyout           !<-- <vy>= 0
      vz(l)=  vzout
      end do 
!
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        nn= 0
        vx1= 0
        vy1= 0
        vz1= 0
!
        do l= 1,npr
        nn= nn +1
        vx1= vx1 +vx(l)
        vy1= vy1 +vy(l)
        vz1= vz1 +vz(l)
        end do
!
        write(11,*) 'loadpt ksp=',ksp
        write(11,'("  the total npr=",i10,/,        &
                   "  <vx1>, <vy1>, <vz1>=",1p3d12.3)') &
                                  nn,vx1/nn,vy1/nn,vz1/nn
        close(11)
      end if
!
!*  For write out.
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,666) npr,qmult,wmult
  666   format(/,' ## npr=',i10,'  qspec(ksp)=',f7.3,',  wspec(ksp)=',f7.3)
        write(11,667) vthx(ksp),vthz(ksp),vdr(ksp),vbeam(ksp)
  667   format('    vthx,vthz(ksp)=',2f8.4,'  vdr,vbeam(ksp)=',2f8.4)
        close(11)
      end if 
!
      return
      end subroutine loadpt
!
!
!-----------------------------------------------------------------------
      subroutine readpt (xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti,  &
                         xe,ye,ze,vxe,vye,vze,qmulte,wmulte,  &
                         npr)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
      include 'mpif.h'
!
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,vxi,vyi,vzi, &
                                       xe,ye,ze,vxe,vye,vze
      real(C_DOUBLE) qmulti,qmulte,wmulti,wmulte
      integer(C_INT) npr,MPIerror
!
!----------------------------------------------------------------------
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!***
      real(C_DOUBLE) sq0,ssq0
      integer(C_INT) l
!----------------------------------------------------------------------
!  Arrays are continuos
      read(19) npr
      read(19) xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti
      read(19) xe,ye,ze,vxe,vye,vze,qmulte,wmulte
!
      sq0= 0
      do l= 1,npr
      sq0= sq0 +qmulti
      end do
!
      call mpi_allreduce (sq0,ssq0,1,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      sq0= ssq0
!
      q0= sq0/float(mxyz0)
!
!  The xyz coordinates
      call partbcT (xi,yi,zi,vxi,vyi,vzi,npr)
      call partbcT (xe,ye,ze,vxe,vye,vze,npr)
!
      return
      end subroutine readpt
!
!
!-----------------------------------------------------------------------
      subroutine ipleql (xi,yi,zi,xe,ye,ze,npr)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
!----------------------------------------------------------------------
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,xe,ye,ze
      integer(C_INT) npr,l
!----------------------------------------------------------------------
!
      do l= 1,npr
      xe(l)= xi(l)
      ye(l)= yi(l)
      ze(l)= zi(l)
      end do
!
      return
      end subroutine ipleql
!
!
!-----------------------------------------------------------------------
      function fun1(v)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE) fun1,v
!
      fun1= exp(-v**2)
!
      return
      end function fun1
!
!
!-----------------------------------------------------------------------
      function fun2(v)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      real(C_DOUBLE) fun2,v,vrg1
      common /vring/ vrg1
!
      fun2= exp(-v**2)*(v+vrg1)
!
      return
      end function fun2
!
!
!-----------------------------------------------------------------------
      function funr(r)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
!  prof= max(1.d0 -arb*(r/rwd)**2, 0.d0) 
      real(C_DOUBLE) funr,r,prof
      real(C_DOUBLE) arb,zcent,ycent1,ycent2,Ez00,vrg1
      common/profl/  arb,zcent,ycent1,ycent2,Ez00
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
!*    volume element r*dr*dtheta; r is included in funr.
!
      prof= max(1.d0 -arb*(r/rwd)**2, 0.d0)  !! generic name
      funr= r*prof
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine rantbl
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      integer(C_INT) ir,iq
      common/ranfa/ ir
      common/ranfb/ iq
!
      ir= 3021 
      iq= 7331
!
      return
      end subroutine rantbl
!
!
!-----------------------------------------------------------------------
      function ranf(x)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      integer(C_INT) ir  ! iand is a generic function
      real(C_DOUBLE) ranf,x
      common/ranfa/ ir
!
      real(C_DOUBLE) invm  ! real*8
      integer(C_INT) lambda,mask
      parameter  (lambda=48828125,mask=2**30+(2**30-1))
      parameter  (invm= 0.5d0**31)
!
      ir= iand( lambda*ir, mask) !! compare by a bit
      ranf= ir*invm              !! 1 if both are 1, otherwise 0
!
      return
      end function ranf
!
!
!-----------------------------------------------------------------------
      function ranfp(x)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h'
!
      integer(C_INT) ir
      real(C_DOUBLE) ranfp,x
      common/ranfb/ ir
!
      real(C_DOUBLE) invm  ! real*8
      integer(C_INT) lambda,mask      
      parameter  (lambda=48828125,mask=2**30+(2**30-1))
      parameter  (invm= 0.5d0**31)
!
      ir= iand( lambda*ir, mask)
      ranfp= ir*invm
!
      return
      end function ranfp
!
!
!-----------------------------------------------------------------------
      subroutine histry
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
!     parameter   (nhistm=54)
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      real(C_DOUBLE) tdec(3000)  !<-- DOUBLE
      common/ehist/  tdec
!
      character(len=54) :: elab(nhistm)
!*
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!                         +++++++
      real(C_float)  emax1a,emin1a,emax2a,emin2a,emax1,emin1, &
                     emax3a,emin3a,emax4a,emin4a,emax3,emin3, &
                     emax5a,emin5a,emax6a,emin6a,emax5,emin5, &
                     emax7a,emin7a,emax8a,emin8a,emax7,emin7
      integer(C_INT) ILN,ILG,i,k
!
!    +++++++++++++++++++++++
      if(io_pe.ne.1) return
!    +++++++++++++++++++++++
!
      call lblbot(t)
!     tdec(ldec)= t
!
!**
      open (unit=11,file=praefixc//'.11'//suffix2,             & 
            status='unknown',position='append',form='formatted')
!
      write(11,*) '*histry: ldec, tdec(ldec)=',ldec,tdec(ldec)
      write(11,*)
!
      do i= 1,ldec ! in every nha 
      write(11,'("tdec=",f7.1,8d11.3)') &
                  tdec(i),edec(i,1),edec(i,2),edec(i,3),edec(i,4), &
                  edec(i,5),edec(i,6),edec(i,7),edec(i,8)
                  
      end do
!
      close(11)
!
! 
      open (unit=77,file=praefixc//'.77'//suffix2//'.ps',        &
            status='unknown',position='append',form='formatted')
!
      ILN= 1
      ILG= 2
!
      call lplmax (edec(1,1),emax1a,emin1a,ldec)
      call lplmax (edec(1,2),emax2a,emin2a,ldec)
      emax1 = max(emax1a,emax2a)
      emin1 = 0  ! min(emin1a,emin2a)
!
      call lplmax (edec(1,3),emax3a,emin3a,ldec)
      call lplmax (edec(1,4),emax4a,emin4a,ldec)
      emax3 = max(emax3a,emax4a)
      emin3 = 0  ! min(emin3a,emin4a)
!
      call lplmax (edec(1,5),emax5a,emin5a,ldec)
      call lplmax (edec(1,6),emax6a,emin6a,ldec)
      emax5 = max(emax5a,emax5a)
      emin5 = 0  ! min(emin3a,emin4a)
!
      call lplot (2,4,ldec,tdec,edec(1,1),emax1a,0.,ILN,'B2 Histr',8,&
                 '        ',8,'        ',8)
      call lplot (2,5,ldec,tdec,edec(1,2),emax2a,0.,ILN,'E2 Histr',8,&
                 '        ',8,'        ',8)
      call lplot (2,6,ldec,tdec,edec(1,3),emax3a,0.,ILN,'B2P Hist',8,&
                 '        ',8,'        ',8)
      call lplot (3,4,ldec,tdec,edec(1,4),emax4a,0.,ILN,'E2P Hist',8,&
                 '        ',8,'        ',8)
      call lplot (3,5,ldec,tdec,edec(1,5),emax5a,0.,ILN,'ionx His',8,&
                 '        ',8,'        ',8)
      call lplot (3,6,ldec,tdec,edec(1,6),emax6a,0.,ILN,'ionh His',8,&
                 '        ',8,'        ',8)
!   ++++++++++++++
      call chart
!   ++++++++++++++
!
      call lplmax (edec(1,7),emax7a,emin7a,ldec)
      call lplmax (edec(1,8),emax8a,emin8a,ldec)
      emax7 = max(emax7a,emax8a)
      emin7 = 0  ! min(emin7a,emin8a)
!
      call lplot (2,4,ldec,tdec,edec(1,7),emax7a,0.,ILN,'elex His',8,&
                 '        ',8,'        ',8)
      call lplot (2,5,ldec,tdec,edec(1,8),emax8a,0.,ILN,'eleh His',8,&
                 '        ',8,'        ',8)
!   ++++++++++++++
      call chart
!   ++++++++++++++
!**
      close(77)
!
      return
      end subroutine histry
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      integer(C_INT) i,is
      real(C_DOUBLE) f(is)
      real(C_float)  fmax,fmin
!
      fmax= -1.e10
      fmin=  1.e10
!
      do i= 1,is
      fmax= max(fmax,real(f(i)))
      fmin= min(fmin,real(f(i)))
      end do
!
      return
      end subroutine lplmax
!
!
!------------------------------------------------------
      subroutine lplmax3 (f1,f2,f3,fmax,fmin,n)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      integer(C_INT) i,n
      real(C_DOUBLE) f1(n),f2(n),f3(n)
      real(C_float)  fmax,fmin
!
      fmax= -1.e10
      fmin=  1.e10
!
      do i= 1,n
      fmax= max(fmax,real(f1(i)),real(f2(i)),real(f3(i)))
      fmin= min(fmin,real(f1(i)),real(f2(i)),real(f3(i)))
      end do
!
      return
      end subroutine lplmax3
!
!
!-----------------------------------------------------------------------
      subroutine rehist
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
!
      integer(C_INT) j,l,l1
!
      nha = 2*nha
      ldec= ldec/2
!
      do j=1,nhistm
      do l=1,ldec
      l1=2*l-1
      edec(l,j)=edec(l1,j)
      end do
      end do
!
      return
      end subroutine rehist
!
!
!-----------------------------------------------------------------------
      subroutine restrt (xi,yi,zi,vxi,vyi,vzi,qmulti,wmulti,  &
                         xe,ye,ze,vxe,vye,vze,qmulte,wmulte,  &
                         npr,ipar,size,iresrt) 
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      implicit none
      include 'param_080A.h' 
      include 'mpif.h'
!
!--------------------------------------------
      real(C_DOUBLE),dimension(np0) :: xi,yi,zi,vxi,vyi,vzi
      real(C_DOUBLE),dimension(np0) :: xe,ye,ze,vxe,vye,vze
!
      real(C_DOUBLE),dimension(np0) :: xxi,yyi,zzi,vvxi,vvyi,vvzi
      real(C_DOUBLE),dimension(np0) :: xxe,yye,zze,vvxe,vvye,vvze
!
      real(C_DOUBLE) qmulti,qmulte,wmulti,wmulte
      integer(C_INT) npr
      integer(C_INT) rank,ipar,size,iresrt,l,MPIerror
!
      real(C_DOUBLE),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                         ex,ey,ez,bx,by,bz,       &
                                         ex0,ey0,ez0,bx0,by0,bz0, &
                                         qix,qiy,qiz,qex,qey,qez, &
                                         emx,emy,emz,qi,qe,amu,avhh, &
                                         pot,gnu
      common/fields/ ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      common/srimp7/ qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      common/srimp8/ amu,avhh
      common/sclpot/ pot
      common/dragcf/ gnu
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      integer(C_INT) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec,nfwrt,npwrt,         &
                     nha,nplot,nhist
      common/parm1/  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,         &
                     itermx,iterfx,itersx,nspec(4),nfwrt,npwrt,      &
                     nha,nplot,nhist
!
      real(C_DOUBLE) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec,wspec,veth,te_by_ti,wce_by_wpe,thb,       &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      & 
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx,vthz,vdr,vbeam,                  &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec
      common/parm2/  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,   &
                     qspec(4),wspec(4),veth,te_by_ti,wce_by_wpe,thb, &
                     rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,      &
                     q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,  &
                     qqwi,qqwe,vthx(4),vthz(4),vdr(4),vbeam(4),      &
                     efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,&
                     edec(3000,12)
      real(C_float ) ctrans
!
      real(C_DOUBLE) tdec(3000)  !<-- DOUBLE
      common/ehist/  tdec
!
!          +++++++                    +++++++
      real(C_float),dimension(-2:mx+1,-1:my+1,-2:mz+1) :: &
                                           cix,ciy,ciz,cex,cey,cez, &
                                           avex,avey,avez,avbx,avby,&
                                           avbz,avqi,avqe
      common/plotav/ cix,ciy,ciz,cex,cey,cez,       &
                     avex,avey,avez,avbx,avby,avbz, &
                     avqi,avqe
      integer(C_INT) i,j,k
!----------------------------------------------------------------------
!
      rank= ipar -1
      if(iresrt.eq.1) go to 1000
!
!***********************************************************************
!*  Create the dump data file for closing the run                      *
!***********************************************************************
!  iresrt= 2 
!    Synchronization of all data with mod(l-1,size).ne.rank.
!    An n-step jumping is not arranged like mpi_gather !
!**
      do l= 1,np0
      xxi(l) = 0
      yyi(l) = 0
      zzi(l) = 0
      vvxi(l)= 0
      vvyi(l)= 0
      vvzi(l)= 0
!
      xxe(l) = 0 
      yye(l) = 0
      zze(l) = 0
      vvxe(l)= 0
      vvye(l)= 0
      vvze(l)= 0
      end do
!
!                   l-1= 0,1,...,7; rank= ipar-1= 0,1,..,7
      do l= 1,np0    
      if(mod(l-1,size).ne.rank) then
        xi(l) = 0 
        yi(l) = 0
        zi(l) = 0
        vxi(l)= 0
        vyi(l)= 0
        vzi(l)= 0
      end if
      end do
!
      do l= 1,np0
      if(mod(l-1,size).ne.rank) then
        xe(l) = 0 
        ye(l) = 0
        ze(l) = 0
        vxe(l)= 0
        vye(l)= 0
        vze(l)= 0
      end if
      end do
!
!  All synchronize from xi in a node to xxi
      call mpi_allreduce (xi(1),xxi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (yi(1),yyi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (zi(1),zzi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vxi(1),vvxi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vyi(1),vvyi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vzi(1),vvzi(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
!
      call mpi_allreduce (xe(1),xxe(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (ye(1),yye(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (ze(1),zze(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vxe(1),vvxe(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vye(1),vvye(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
      call mpi_allreduce (vze(1),vvze(1),np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,MPIerror)
!
!   Only main node
      if(io_pe.eq.1) then
!**
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'it,t=',it,t
        write(11,*) 'ions...by allreduce'
!
        do i= 1,12
        write(11,'("i=",i8,1p6d12.3)') &
                      i,xxi(i),yyi(i),zzi(i),vvxi(i),vvyi(i),vvzi(i)
        end do
!
        write(11,*)
        write(11,*) 'elec...by allreduce'
!
        do i= 1,12
        write(11,'("i=",i8,1p6d12.3)') &
                      i,xxe(i),yye(i),zze(i),vvxe(i),vvye(i),vvze(i)
        end do
!**
! Data dump for restart
!
        open (unit=12,file=praefixc//'.12'//suffix2,           &
                            status='replace',form='unformatted')
!                            ++++++++ for write
!
        write(12)  it,ldec,iaver,ifilx,ifily,ifilz,iloadp,           &
                   itermx,iterfx,itersx,(nspec(i),i=1,4),            &
                   nfwrt,npwrt,nha,nplot,nhist
        write(12)  xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,     &
                   (qspec(i),i=1,4),(wspec(i),i=1,4),                &
                   veth,te_by_ti,wce_by_wpe,thb,         &
                   rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,        &
                   q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,    &
                   qqwi,qqwe,(vthx(i),i=1,4),(vthz(i),i=1,4),        &
                   (vdr(i),i=1,4),(vbeam(i),i=1,4),&
                   efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,  &
                   ((edec(i,j),i=1,3000),j=1,4)
        write(12)  tdec
!
        write(12)  ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
        write(12)  qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe 
        write(12)  amu,avhh
        write(12)  pot,gnu
!
        write(12)  cix,ciy,ciz,cex,cey,cez,avex,avey,avez, & ! real*4
                   avbx,avby,avbz,avqi,avqe
!
        write(12)  qmulti,wmulti,qmulte,wmulte
        write(12)  npr 
!
        write(12)  xxi,yyi,zzi,vvxi,vvyi,vvzi   !<- complete if all PE
        write(12)  xxe,yye,zze,vvxe,vvye,vvze
!
        close(12)
!
        write(11,'("** Restrt file is created FT12...",/)') 
        close(11)
      end if
!**
      return
!
!
!***********************************************************************
!*  Restart the run from now                                           *
!***********************************************************************
!* iresrt= 1
!
!   All nodes must do
 1000 continue
      open (unit=12,file=praefixc//'.12'//suffix1,form='unformatted')
!                                   ++++++++ old
!
      read(12) it,ldec,iaver,ifilx,ifily,ifilz,iloadp,           &
               itermx,iterfx,itersx,(nspec(i),i=1,4),            &
               nfwrt,npwrt,nha,nplot,nhist
      read(12) xmax,ymax,zmax,hxi,hyi,hzi,xmaxe,ymaxe,zmaxe,     &
               (qspec(i),i=1,4),(wspec(i),i=1,4),                &
               veth,te_by_ti,wce_by_wpe,thb,         &
               rwd,pi,ait,t,dt,aimpl,adt,hdt,ahdt2,adtsq,        &
               q0,qi0,qe0,aqi0,aqe0,epsln1,qwi,qwe,aqwi,aqwe,    &
               qqwi,qqwe,(vthx(i),i=1,4),(vthz(i),i=1,4),        &
               (vdr(i),i=1,4),(vbeam(i),i=1,4),&
               efe,efb,etot0,bxc,byc,bzc,vlima,vlimb,bmin,emin,  &
               ((edec(i,j),i=1,3000),j=1,4)
      read(12)  tdec
!
      read(12)  ex,ey,ez,bx,by,bz,ex0,ey0,ez0,bx0,by0,bz0
      read(12)  qix,qiy,qiz,qex,qey,qez,emx,emy,emz,qi,qe
      read(12)  amu,avhh
      read(12)  pot,gnu
!
      read(12)  cix,ciy,ciz,cex,cey,cez,avex,avey,avez, & ! real*4
                avbx,avby,avbz,avqi,avqe
!
      read(12)  qmulti,wmulti,qmulte,wmulte
      read(12)  npr
!
      read(12)  xi,yi,zi,vxi,vyi,vzi    !<- complete set
      read(12)  xe,ye,ze,vxe,vye,vze
!
      close(12)
!
!*
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) '# At restart #'
        write(11,*) 'it,t=',it,t
        write(11,*)
!
        write(11,*) 'nspec(1),qspec(1),wspec(1)=',nspec(1),qspec(1),wspec(1)
        write(11,*) 'nspec(2),qspec(2),wspec(2)=',nspec(2),qspec(2),wspec(2)
!
        write(11,*) 'xmax,ymax,zmax=',xmax,ymax,zmax
        write(11,*)
!
        write(11,*) 'ions...at restart'
!
        do i= 1,12
        write(11,'("i=",i8,1p6d12.3)') &
                         i,xi(i),yi(i),zi(i),vxi(i),vyi(i),vzi(i)
        end do
!
        write(11,*)
        write(11,*) 'elec...at restart'
!
        do i= 1,12
        write(11,'("i=",i8,1p6d12.3)') &
                         i,xe(i),ye(i),ze(i),vxe(i),vye(i),vze(i)
        end do
!
!
        write(11,*) 'emfild...'
        do k= 0,3
        do j= 0,3
        do i= 0,3
        write(11,'(3i3,1p6d12.3)') &
                         i,j,k,ex(i,j,k),ey(i,j,k),ez(i,j,k),bx(i,j,k), &
                         by(i,j,k),bz(i,j,k)
        end do
        end do
        end do
!
        close(11)
      end if
!
!  Only the ipar node must have real data, is null otherwise
!  This part is actually not touched.
!
!     do l= 1,np0
!     if(mod(l-1,size).ne.rank) then
!     xi(l) = 0 
!     yi(l) = 0
!     zi(l) = 0
!     vxi(l)= 0
!     vyi(l)= 0
!     vzi(l)= 0
!     end if
!     end do
!
!     do l= 1,np0
!     if(mod(l-1,size).ne.rank) then
!     xe(l) = 0 
!     ye(l) = 0
!     ze(l) = 0
!     vxe(l)= 0
!     vye(l)= 0
!     vze(l)= 0
!     end if
!     end do
!*
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,'(" *** Start files are loaded ***",/)') 
        close(11)
      end if  
!**
      return
      end subroutine restrt
!
!
!-----------------------------------------------------------------------
      subroutine lbltop (date_now,label)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      character(len=8)  :: label(8),label1(8)
      character(len=10) :: date_now,date_now1
      common/headr1/ label1,date_now1
!
      label1(1)= label(1)
      label1(2)= label(2)
      label1(3)= label(3)
      label1(4)= label(4)
      label1(5)= label(5)
      label1(6)= label(6)
      label1(7)= label(7)
      label1(8)= label(8)
      date_now1= date_now
!
      return
      end subroutine lbltop
!
!
!-----------------------------------------------------------------------
      subroutine lblbot (ts)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE) time1,ts
      common/headr2/ time1
!
      time1= ts
!
      return
      end subroutine lblbot
!
!
!-----------------------------------------------------------------------
      subroutine endrun(h)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      character(len=8) h
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,20) h
        write(11,20) h
   20   format(///,'## run terminated ----- ',a8)
        close(11)
      end if
!
      if(.true.) stop
!
      return
      end subroutine endrun
!
!
!------------------------------------------------------
      subroutine clocks (walltime,size,cl_first)
!------------------------------------------------------
!*  measure both cpu and elapsed times 
!   needs "walltime,size,cl_first" instead of ts
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
      integer(C_INT) size,cl_first,MPIerror
      real(C_DOUBLE) walltime,walltime0,buffer1,buffer2
      save       walltime0
!
      if(cl_first.eq.1) then
        walltime0= mpi_wtime() 
!
        buffer1= walltime0
        call mpi_allreduce (buffer1,buffer2,1,mpi_real8,mpi_sum, &
                            MPI_COMM_WORLD,MPIerror)
        walltime0 = buffer2/size
      end if
!
      walltime = mpi_wtime() - walltime0
!
      buffer1= walltime
      call mpi_allreduce (buffer1,buffer2,1,mpi_real8,mpi_sum, &
                          MPI_COMM_WORLD,MPIerror)
      walltime = buffer2/size
!
      return
      end subroutine clocks
!
!
!------------------------------------------------------
      subroutine clocki (walltime,size,cl_first)
!------------------------------------------------------
!*  measure both cpu and elapsed times 
!   needs "walltime,size,cl_first" instead of ts
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'mpif.h'
      include 'param_080A.h' 
!
      integer(C_INT) size,cl_first,MPIerror
      real(C_DOUBLE) walltime,walltime0,buffer1,buffer2
      save       walltime0
!
      if(cl_first.eq.1) then
        walltime0= mpi_wtime() 
!
        buffer1= walltime0
        call mpi_allreduce (buffer1,buffer2,1,mpi_real8,mpi_sum, &
                            MPI_COMM_WORLD,MPIerror)
        walltime0 = buffer2/size
      end if
!
      walltime = mpi_wtime() - walltime0
!
      buffer1= walltime
      call mpi_allreduce (buffer1,buffer2,1,mpi_real8,mpi_sum, &
                          MPI_COMM_WORLD,MPIerror)
      walltime = buffer2/size
!
      return
      end subroutine clocki
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding  ! <-
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time_7
!
!
!-----------------------------------------------------------------------
      subroutine lplot (ix,iy,npts,xsc,q,ymax,ymin,IL,lab1,n1,lab2,n2, & 
                        lab3,n3)
!-----------------------------------------------------------------------
!* IFR= 7 ---> -ymin < y < ymax plot.
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      real(C_DOUBLE),dimension(3000) :: xsc,q
      real(C_DOUBLE),dimension(101)  :: vsc,ff
      real(C_float), dimension(3000) :: u,v
!          +++++++
      integer(C_INT) ix,iy,npts,IFR,IL,n1,n2,n3
      integer(C_INT) i1,j1,isc,i,j,k,nfine,iplot
!
      real(C_float) ymin,ymax
      real(C_float) xcm(6),ycm(6),pl(6),pr(6),ql(6),qr(6)
!
      character(len=8) lab1,lab2,lab3
      character(len=8) label(8),date_now*10,cax*1
      common/headr1/ label,date_now
      common/headr2/ time,xp_leng
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10),  &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
!   for Fujitsu.
!     data  xcm/18.46,2*9.867,3*6.18/,
!    *      ycm/16.85,2*7.435,3*4.381/,
!    *      pl/2*2.00,15.132,2.00,8.00,18.20/,
!    *      ql/1.95,10.885,1.95,13.832,7.891,1.95/
!
!   for NEC.
      data  xcm/21.0, 2*10.00, 3*7.00/,  &
            ycm/15.0, 2*6.80, 3*3.90/,   &
            pl/2.0,  2.0,14.0, 1.0,9.0,17.0/, &
            ql/2.3, 10.5,2.3, 12.9,7.6,2.3/
      logical  lab_skip
!
      real(C_float) hh,hhs,xmax,xmin,dx,dy,x0,y0, &
                    pl1,pr1,ql1,qr1,scx,scy,time,xp_leng,   &
                    xmin1,xmax1,ymin1,ymax1,x1,x2,x3,x4,    &
                    y1,y2,y3,y4,xc,xd,xl,xu,yc,yl,yr
!***
      iplot= 1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot (ix,iy,npts,vsc,ff,ymin,ymax,IL,lab1,n1,lab2,n2,  &
                   lab3,n3)
!-----------------------------------------------------------------------
      iplot=2
!
    1 isc= 1
!
      do i=1,6
      pr(i)= pl(i) +xcm(i)
      end do
!
      do j=1,6
      qr(j)= ql(j) +ycm(j)
      end do
!
      lab_skip= .false.
      if(il.eq.7) lab_skip= .true.
!
!                 ******************************************************
!*                **  Make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh = 0.70
      hhs= 0.60
!
      i1= abs(ix)
      j1= abs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (0.1,18.0,hh,label(1),0.,8)
!     call symbol (3.1,18.0,hh,date_now, 0.,10)
      call symbol (15.9,0.1,hh,'t =',0.,3)
      call number (999.0,999.0,hh,time,0.,5)
!
   10 continue
!
      if(iplot.eq.1) then
!++
      do i= 1,npts
      u(i)= xsc(i)
      end do
      xmin= u(1)
      xmax= u(npts)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(IL.eq.1) then
        v(1)   = q(1)
        v(npts)= q(npts)
!
        do i= 2,npts-1
        v(i)=   q(i)
!       v(i)= 0.33333*(q(i-1)+q(i)+q(i+1))
        end do
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(abs(IL).eq.2) then
         do i= 1,npts
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
         end do
      end if
      end if
!                                                ****************
!                                                **  iplot= 2  **
!                                                ****************
      if(iplot.eq.2) then
        do i= 1,npts
        u(i)= vsc(i)
        v(i)=  ff(i)
        end do
      end if
!
!                                **************************************
!                                ** set a new scale and draw a frame **
!                                **************************************
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      call scalex (pl(i1),ql(j1),x0,y0,dx,dy,isc)
!
      pl1(isc)= pl(i1)
      pr1(isc)= pr(i1)
      ql1(isc)= ql(j1)
      qr1(isc)= qr(j1)
      xmin1(isc)= xmin
      xmax1(isc)= xmax
      ymax1(isc)= ymax
      ymin1(isc)= ymin
!                                                      *************
!                                                      **  frame  **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  tick marks  **
!                                                    ******************
      scx= xcm(i1)/5.0
      scy= ycm(j1)/4.0
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do k=1,4
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
      end do
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do k=1,3
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
      end do
!                                                     **************
!                                                     ** numbers. **
!                                                     **************
!
      if(.not.lab_skip) then
        call number (pl(i1)-0.5,ql(j1)-0.45,hhs,xmin,0.,101)
        call number (pr(i1)-1.5,ql(j1)-0.45,hhs,xmax,0.,101)
!
        call number (pl(i1)-2.0,ql(j1)     ,hhs,ymin,0.,101)
        call number (pl(i1)-2.0,qr(j1)-0.30,hhs,ymax,0.,101)
      end if
!
!                                                     **************
!                                                     **  labels. **
!                                                     **************
      xc= 0.5*(pl(i1)+pr(i1))
      xu= xc -1.60
      xd= xc -0.20*n2/2
!
      yr= qr(j1)+0.15
      yl= ql(j1)-0.70
!
      call symbol (xu,yr,hh,lab1,0.,n1)
      call symbol (xd,yl,hh,lab2,0.,n2)
!
      xl= pl(i1)-1.50
      yc= 0.5*(ql(j1)+qr(j1))
      call symbol (xl,yc,hh,lab3,0.,n3)
!                                     **********************************
!                                     **  no plot is made if npts < 0 **
!                                     **********************************
   70 if(npts.lt.0) return
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do i=1,npts
         call plotl (u(i),v(i),isc,2)
         end do
      else
         do i=1,npts
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
         end do
      end if
!**
      call plotl (u(npts),v(npts),isc,3)
!
      return
      end subroutine lplot
!
!
!-----------------------------------------------------------------------
      subroutine pplot (x,y,z,nr,lskip,isc)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_float),dimension(100000) :: x,y,z
!          +++++++    
      integer(C_INT) nr,lskip,isc,iplot,ipen,itag,i
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      if(io_pe.eq.1) then
!***
      open (unit=18,file=praefixc//'.18'//suffix2,                   &
                status='unknown',position='append',form='unformatted')
      iplot= 3
!
      ipen= 0
      itag= 0
      write(18) iplot,ipen,itag
      write(18) nr,lskip,isc
      write(18) (x(i),i=1,nr,lskip),(y(i),i=1,nr,lskip)
!          <-- merge to real*4 and all x,y,z ??
!
      close(18)
      end if
!
      return
      end subroutine pplot
!
!
!-----------------------------------------------------------------------
      subroutine pplot3 (x,y,z,vx,vy,vz,xmax,ymax,zmax,x1,x2,radi,  &
                         nr,lskip,iptag,char1,n1)
!-----------------------------------------------------------------------
!  real*4 except real*8 x - vz
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_DOUBLE),dimension(np0) :: x,y,z,vx,vy,vz
!          +++++++
      real(C_float) xmax,ymax,zmax,x1,x2,radi,time
      integer(C_INT) nr,lskip,iptag,n1,iplot,i
!
      character(len=8) char1
      character(len=8) label(8)
      character(len=10) date_now
      common/headr1/ label,date_now
!
      real(C_DOUBLE) time1 
      common/headr2/ time1
!
      integer(C_INT) ifcopy,ipen,iksy
      common/plots/  ifcopy,ipen,iksy
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
      time= time1
!
      if(io_pe.ne.1) return
      iplot= 6
!***
      open (unit=18,file=praefixc//'.18'//suffix2,                   &
                status='unknown',position='append',form='unformatted')
!
      write(18) iplot,ipen,iptag
      write(18) time,xmax,ymax,zmax,x1,x2,radi,nr,lskip,n1
      write(18) (x(i),y(i),z(i),vx(i),vy(i),vz(i),i=1,nr,lskip)
!         merge to all x-vz ?
      write(18) label,date_now,char1
!
      close(18)
!
      return
      end subroutine pplot3
!* data are written here 
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!
       write(77,*) 'stroke'
!
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
!
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
!
       return
       end subroutine newcolor 
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  open circle centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
!
      x0= x +a
      y0= y
      call plot (x0,y0,3)
!
      do 100 j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
  100 continue
!
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
!
      if(ic.eq.1) return
!------------------------------------
!*  filled circle centered at (x,y).
!------------------------------------
!
      write(77,*) " 3.0 setlinewidth"
!
      nc= 5
      dth= pi/(2*nc +1)
!
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      x2= x1
      y2= 2.*y -y1
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
!
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end subroutine circle 
!
!
!------------------------------------------------
      integer function iwrta (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
!
      return
      end function iwrta
!
!
!------------------------------------------------
      integer function iwrtb (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
!
      return
      end function iwrtb
!
!
!------------------------------------------------
      integer function iwrtc (t,twr)
!------------------------------------------------
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwc) then
        iwc= iw
        iwrtc= 0
      else
        iwrtc= 1
      end if
!
      return
      end function iwrtc
!
!
!------------------------------------------------
      function ranff (x)
!------------------------------------------------
!*  ranf= (0,1)
!
      use, intrinsic :: iso_c_binding
      common/ranfff/ ir,iq
!
      real(C_DOUBLE)  ranff,invm
      parameter  (mask=2**30+(2**30-1),invm= 0.5d0**31)
      parameter  (lambda=48828125)
!
      ir= iand( lambda*ir, mask)
      ranff= ir*invm
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir= amod( ambda*ir +qq, ask)
!     ranff= ir/ask
!
      return
      end function ranff 
!
!
!-----------------------------------------------------------------------
      subroutine cplt2a (q,xmax,ymax,zmax,iperio,iptag,char1,n1)
!-----------------------------------------------------------------------
!  real*4
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!*
      integer  n1
      real(C_DOUBLE)   time1
      character(len=8) char1,label(8),date_now*10
      common/headr1/ label,date_now
      common/headr2/ time1
!
      integer        pxl,pxc,pxr,pyl,pyc,pyr,pzl,pzc,pzr
      common/table/  pxl(-mx:2*mx-1),pxc(-mx:2*mx-1),pxr(-mx:2*mx-1), &
                     pyl(-1:my+1),pyc(-1:my+1),pyr(-1:my+1),          &
                     pzl(-mz:2*mz-1),pzc(-mz:2*mz-1),pzr(-mz:2*mz-1)
!
      real(C_float)  q(mx,myA,mz),qc,qc2
      real(C_float)  a(11000),b(11000),aa(11000),ww(11000),as(11000)
!
      real(C_float) time,zleng,xmin,ymin,zmin,hh,hs,      &
                    yl,yr,zl,zr,yls,yrs,zls,zrs,am1,am2,  & 
                    am12,y,z,gdy,gdz,xmax,ymax,zmax
      integer(C_INT) i1,j,k,kl,kr,jk,npy,npz,iplot,ny,nz, &
                     nyz,nzp,iperio,iptag,ncontr
!
      time= time1
      zleng= 8. 
!
      xmin= 0
      ymin= 0
      zmin= 0
!
!-----------------------------------------------------------------------
!*   use both  k=1 and k=nz data for good continuation on the sides.
!*  clip-edge plot  /or wall-to-wall plot.
!-----------------------------------------------------------------------
!* 1. plot at j= j1: subscript k first.
!-----------------------------------------------------------------------
!
      i1= mx/2
!
      qc2= 1./6.
      jk= 0
!***
      npy= 0
      do j= 2,my,2
      npy= npy +1
!
        npz= 0
        do k= 2,mz-1,2
        npz= npz +1
        jk= jk +1
!
        kr= pzr(k)
        kl= pzl(k)
!
        a(jk)= qc2*( q(i1,j+1,kr) +q(i1,j+1,k) +q(i1,j+1,kl)  &
                 +   q(i1,j  ,kr) +q(i1,j  ,k) +q(i1,j  ,kl))
        end do
        end do
!
!---------------------------------------------------------
      iplot= 5
!
!     write(18) iplot,iperio,iptag
!     write(18) time,xmax,ymax,zmax,npx,npy,npz,n1
!     write(18) a,b
!     write(18) label,date_now,char1
! 
!---------------------------------------------------------
!********************************
!*  make postscript plot file.  *
!********************************
!
      nz= npz
      ny= npy
!
      hh= 0.60
      hs= 0.40
      call symbol (0.1,18.2,hh,label(1),0.,8)
!     call symbol (0.1, 0.1,hh,date_now,0.,10)
      call symbol (15.9,0.1,hh,'t=',0.,2)
      call number (16.5,0.1,hh,time,0.,101)
!
      zl= 2. +zleng +1.
      zr= zl +zleng
!
      yl= 2.5
      yr= yl +amin1(17.,zleng*ymax/zmax)
!
      zls= zl
      zrs= zr
      yls= yr  +0.5
      yrs= yls +2.0
!
      nyz= ny*nz
!
      am1= -1.e30
      am2=  1.e30
!
      do jk= 1,nyz
      am1= max(am1,a(jk))
      am2= min(am2,a(jk))
      end do
!
      call symbol (zl+1.0,yl-1.5,hh,'max=',0.,4)
      call number (zl+2.7,yl-1.5,hh,am1,0.,101)
      call symbol (zl+1.0,yl-2.0,hh,'min=',0.,4)
      call number (zl+2.7,yl-2.0,hh,am2,0.,101)
!
      call setscl (0.,0.,zmax,ymax,zl,yl,zr,yr,gdz,gdy,  &
                   n1,char1, 8,' scalar ', hh,           &
                   1,'Z',0.4, 1,'X',0.4, 1)
!
      call number (zl-1.0,yl-0.6,hs,0.  ,0.,100)
      call number (zl-2.0,yr-0.3,hs,ymax,0.,100)
      call number (zr-1.0,yl-0.6,hs,zmax,0.,100)
!
!
      jk= 0
      do 210 j= 1,ny
      do 230 k= 1,nz
      jk= jk +1
      aa(jk)= a(k+nz*(j-1))
  230 continue
!
      if(iperio.eq.1) then
        jk= jk +1
        aa(jk)= a(1+nz*(j-1))
      end if
  210 continue
!
      nzp= nz
      if(iperio.eq.1) nzp= nz +1
!
      ncontr=  7
      call eqcntr (aa,ww,nzp,ny,zl,yl,zr,yr,am2,0.0,am1,  &
                   7700,ncontr,1) 
!
!--------------------------
!  (2) section plot.
!--------------------------
        if(.false.) then
      k= nz/2
!
      do j= 1,ny
      as(j)= a(k +nz*(j-1))
      end do
!
      am1= -1.e30
      am2=  1.e30
!
      do j= 1,ny
      am1= max(am1,as(j))
      am2= min(am2,as(j))
      end do
!
      am12 = am1 -am2
      if(am12.lt.1.e-10) am12= 999.
!
      call symbol (zls-1.9,yls-1.5,hh,'max=',0.,4)
      call number (zls-0.3,yls-1.5,hh,am1,0.,101)
      call symbol (zls-1.9,yls-2.0,hh,'min=',0.,4)
      call number (zls-0.3,yls-2.0,hh,am2,0.,101)
!
      call plot (zls,yls,3)
      call plot (zls,yrs,2)
      call plot (zrs,yrs,2)
      call plot (zrs,yls,2)
      call plot (zls,yls,2)
!
      y= yls + zleng*0.5/ny
      z= zls + 2.0*(as(1) -am2)/am12
      call plot (z,y,3)
!
      do 360 j= 2,ny
      y= yls + zleng*(j-0.5)/ny
      z= zls + 2.0*(as(j) -am2)/am12
      call plot (z,y,2)
  360 continue
        end if
!
!-----------------------
      call chart
!-----------------------
      return
      end subroutine cplt2a
!
!
!-----------------------------------------------------------------------
      subroutine fplot3 (ex,ey,ez,exc,eyc,ezc,xmax,ymax,zmax,         &
                                                 iperio,iptag,char1,n1)
!-----------------------------------------------------------------------
!   real*4  Vector plots in (x,z) plane  (wall bound in y).
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_float)  ex(mx,myA,mz),ey(mx,myA,mz),ez(mx,myA,mz),  &
                     exc,eyc,ezc,xmax,ymax,zmax
      real(C_float)  ax(11000),ay(11000),az(11000),aa(11000),ww(11000)
!*
      integer(C_INT) n1,iperio,iptag
      real(C_DOUBLE) time1
      character(len=8) label(8),date_now*10,char1
      common/headr1/ label,date_now
      common/headr2/ time1
!
      integer(C_INT) i1,iplot,npy,npz,nyz,nzp,isel,             &
                     j,jk,k,kl,kr,ncontr,ny,nz,nynz
      real(C_float)  am1,am11,am12,am2,amy,arw1,arw2,gdy,gdz,   &
                     hh,hs, qc,rarrow,shfy,shfz,time,           &
                     yl,yr,zleng,ymin,zc,zl,zr
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      time= time1
      zleng= 8. 
!
!-----------------------------------------------------------------------
!* 1. plot at i= i1: subscript k first.
!-----------------------------------------------------------------------
!
      i1= mx/2
!
      qc= 1./16.
      jk= 0
!
      npy= 0
      do j= 2,my,2  !<-- j+1= my+1
      npy= npy +1
!
        npz= 0
        do k= 2,mz-1,2  !<-- k+1= mz
        npz= npz +1
!
        kr= k+1 
        kl= k-1 
!
        jk= jk+1
        ax(jk)= qc*( ex(i1,j+1,kr) +2.*ex(i1,j+1,k)    +ex(i1,j+1,kl)  &
                 +2.*ex(i1,j  ,kr) +4.*ex(i1,j  ,k) +2.*ex(i1,j  ,kl)  &
                    +ex(i1,j-1,kr) +2.*ex(i1,j-1,k)    +ex(i1,j-1,kl)) -exc
        ay(jk)= qc*( ey(i1,j+1,kr) +2.*ey(i1,j+1,k)    +ey(i1,j+1,kl)  &
                 +2.*ey(i1,j  ,kr) +4.*ey(i1,j  ,k) +2.*ey(i1,j  ,kl)  &
                    +ey(i1,j-1,kr) +2.*ey(i1,j-1,k)    +ey(i1,j-1,kl)) -eyc
        az(jk)= qc*( ez(i1,j+1,kr) +2.*ez(i1,j+1,k)    +ez(i1,j+1,kl)  &
                 +2.*ez(i1,j  ,kr) +4.*ez(i1,j  ,k) +2.*ez(i1,j  ,kl)  &
                    +ez(i1,j-1,kr) +2.*ez(i1,j-1,k)    +ez(i1,j-1,kl)) -ezc
        end do
      end do
!
!---------------------------------------------------------
      iplot= 4
!
!     write(18) iplot,iperio,iptag
!     write(18) time,xmax,zmax,npx,npz,n1
!     write(18) ax,ay,az
!     write(18) label,date_now,char1
! 
!---------------------------------------------------------
!********************************
!*  make postscript plot file.  *
!********************************
!
      hh= 0.60
      hs= 0.40
      call symbol (0.1,18.2,hh,label(1),0.,8)
!     call symbol (0.1, 0.1,hh,date_now,0.,10)
      call symbol (15.9,0.1,hh,'t=',0.,2)
      call number (16.5,0.1,hh,time,0.,101)
!
      nz= npz
      ny= npy
!
      yl= 2.0
      yr= yl +amin1(17.,zleng*ymax/zmax)
!
      zl= 2.5
      zr= zl +zleng 
      zc= (zr+zl)/2
!
!*  (1): Vector plot of the y-z components.
!
      call setscl (0.,0.,zmax,ymax,zl,yl,zr,yr,gdz,gdy, &
                   n1,char1, 2,'YZ', hh,                &
                   1,'Z',0.4, 1,'Y',0.4, 1)
        ymin= 0.
      call number (zl-1.0,yl-0.6,hs,ymin,0.,100)
      call number (zl-2.0,yr-0.3,hs,ymax,0.,100)
      call number (zr-1.0,yl-0.6,hs,zmax,0.,100)
!
!*  **maximum of the A vectors**
!
      am1= 0.
      am2= 0.
      nyz= ny*nz
!
      do jk= 1,nyz
      amy= sqrt(ax(jk)**2 +ay(jk)**2 +az(jk)**2) 
      am1= max(am1,sqrt(abs(ay(jk)**2 +az(jk)**2))) !amy)
      am2= max(am2,abs(ax(jk)))
      end do
!
      if(am1.lt.1.e-10) am1=999.0
      if(am2.lt.1.e-10) am2=999.0
!
!*  **plot arrows**
!
      arw1  = 0.65
      rarrow= 0.50
!
      arw2= arw1
      am11= am1
      am12= am1
      shfz= 0.5
      shfy= 0.5
      isel= 2
!
!   72/2 * 144/2 = 2592
      call pltarw (az,ay,nz,ny,zl,yl,zr,yr,arw1,arw2,rarrow, &
                   am11,am12,shfz,shfy,7000,isel)
!
      call symbol (zc-3.7,0.6,hh,'Pol=',0.,4)
      call number (zc-2.0,0.6,hh,am1,0.,101)
      call symbol (zc-3.7,0.1,hh,'Tor=',0.,4)
      call number (zc-2.0,0.1,hh,am2,0.,101)
!
!
!*  (2): Contour plot of the ax-component.
!
      yl= 2.0
      yr= yl +amin1(17.,zleng*ymax/zmax)
!
!   zr= zl +zleng -> space 2 cm 
      zl= zr +2.0
      zr= zl +zleng 
!
      call setscl (0.,0.,zmax,ymax,zl,yl,zr,yr,gdz,gdy,  &
                   n1,char1, 2,'AX',hh,                  &
                   1,'Z',0.4, 1,'Y',0.4, 1)
        ymin= 0.
      call number (zl-1.0,yl-0.6,hs,ymin,0.,100)
      call number (zl-2.0,yr-0.3,hs,ymax,0.,100)
      call number (zr-1.0,yl-0.6,hs,zmax,0.,100)
!
!
      am1= -1.e+30
      am2=  1.e+30
      jk= 0
!
      do j= 1,ny
      do k= 1,nz
      jk= jk +1
      am1= max(am1,ax(jk)) !<- ok
      am2= min(am2,ax(jk))
      end do
      end do
!
      jk= 0
      do j= 1,ny
      do k= 1,nz  !<-- k first
      jk= jk +1
      aa(jk)= ax(k +nz*(j-1)) 
      end do
!
!     if(iperio.eq.1) then
!       jk= jk +1
!       aa(jk)= ax(1+nz*(j-1)) !<-- nz+1 ?
!     end if
      end do
!
      nzp= nz
!     if(iperio.eq.1) nzp= nz+1
!
      ncontr=  7
      call eqcntr (aa,ww,nzp,ny,zl,yl,zr,yr,am2,0.,am1, &
                   7000,ncontr,1) 
!
!-----------------------
      call chart
!-----------------------
      return
      end subroutine fplot3
!
!
!-----------------------------------------------------------------------
      subroutine cplot3 (q,xmax,ymax,zmax,char1,n1)
!-----------------------------------------------------------------------
!***********************************************************************
!*   Contour plots of scalar quantities.                               *
!***********************************************************************
!  All real*4
!
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_080A.h' 
!
      real(C_float) q(mx,my+1,mz),xmax,ymax,zmax
      real(C_float) a(7000),b(7000),ww(7000),cut(7000,4)
!
      integer(C_INT)  n1,ik,jk,i1,j1,k1,ii,i,j,k,npx,npy,npz, &
                      il,ir,jl,jr,kl,kr,nyz,nxz,npz2,ncontr
      real(C_float)   time,xl2,xr2,yl,yr,zl,zr,hh,qc,qc2,       &
                      am21,am41,am22,am42,ams1,ams2,ams,          &
                      gdx,gdy,gdz,wamin,wamax,xmin,ymin,zmin,zleng

      real(C_DOUBLE)  time1
      character(len=8) label(8),date_now*10,char1
      common/headr1/  label,date_now
      common/headr2/  time1
!
      integer(C_INT) io_pe
      common/iope66/ io_pe
!
      time= time1
      zleng= 8.
!
      xmin= 0
      ymin= 0
      zmin= 0
!
      i1= mx/2
      j1= my/2 
      k1= mz/2 
!
!* 1. Plot at constant y: subscript i first.
!     real(C_float) q(mx,my+1,mz),xmax,ymax,zmax
!
      qc2 = 1./16.
      jk= 0
!
      npy= 0
      do j= 2,my,2  !<-- j+1= my+1
      npy= npy +1
!
        npz= 0
        do k= 2,mz-1,2  !<-- k+1= mz
        npz= npz +1
!
        kr= k+1 
        kl= k-1
!
        jk= jk +1
        a(jk)= qc2*(    q(i1,j+1,kr) +2*q(i1,j+1,k)   +q(i1,j+1,kl)  &
                    + 2*q(i1,j  ,kr) +4*q(i1,j  ,k) +2*q(i1,j  ,kl)  &
                    +   q(i1,j-1,kr) +2*q(i1,j-1,k)   +q(i1,j-1,kl) )
        end do
      end do
!
!* 2. Plot at constant z: subscript i first.
!
      qc = 1./16.
      ik= 0
!
      npx= 0
      do i= 2,mx-1,2  !<-- ir= mx
      npx= npx +1
!
      ir= i+1 
      il= i-1 
!
        npz2= 0
        do k= 2,mz-1,2  !<-- kr= mz
        npz2= npz2 +1
!
        kr= k+1 
        kl= k-1 
!
        ik= ik +1
        b(ik)= qc*(   q(ir,j1,kr) +2.*q(ir,j1,k)    +q(ir,j1,kl)  &
                  +2.*q(i ,j1,kr) +4.*q(i ,j1,k) +2.*q(i ,j1,kl)  &
                  +   q(il,j1,kr) +2.*q(il,j1,k)    +q(il,j1,kl) )
        end do
      end do
!
!  a(jk)  y z
!  b(ik)  x z
!
      zl=  2.0        
      zr=  zl +zleng 
!
      yl=  2.0 
      yr=  yl +amin1(17.,zleng*ymax/zmax) 
!                          <--- limit elongated y-length.
      hh = 0.70
      call symbol (0.1,18.2,hh,label(1),0.,8)
!     call symbol (0.1, 0.1,hh,date_now,0.,10)
      call symbol (15.9,0.1,hh,'t=',0.,2)
      call number (16.5,0.1,hh,time,0.,101)
!
!---------------------------------------------
!*  **maximum of the vectors**
!---------------------------------------------
!
      am21= -1.d+10
      am22=  1.d+10
!
      do jk= 1,npy*npz
      am21= max(am21,a(jk))
      am22= min(am22,a(jk))
      end do
!
!
      am41= -1.e+10
      am42=  1.e+10
!
      do ik= 1,npx*npz2
      am41= max(am41,b(ik))
      am42= min(am42,b(ik))
      end do
!
      ams1= max(am21,am41)
      ams2= min(am22,am42)
      ams = max(ams1,-ams2)
!     if(ams.lt.1.d-10) ams=999.0
!
      call symbol (2.0,0.80,hh,'max=',0.,4)
      call number (5.0,0.80,hh,ams1,0.,101)
      call symbol (2.0,0.10,hh,'min=',0.,4)
      call number (5.0,0.10,hh,ams2,0.,101)
!
!---------------------------------------------
!*  (1): Contours in (x,z) plane.
!---------------------------------------------
!
      call setscl (0.,0.,zmax,ymax,zl,yl,zr,yr,gdz,gdy, &
                   n1,char1, 2,'YZ', 0.6,               &
                   1,'Z',0.6, 1,'Y',0.6, 1)
!
      call number (zl-1.3,yl-0.3, hh,zmin,0.,5)
      call number (zl-1.3,yr-0.3, hh,zmax,0.,5)
      call number (zr-1.3,yl-0.5, hh,ymax,0.,5)
!
      nyz= npy*npz
      call daisho (a,nyz,wamin,wamax)
!
      ncontr= 7 
      call eqcntr (a,ww,npz,npy,zl,yl,zr,yr,wamin,0.,wamax, &
                   7000,ncontr,1) 
!
!---------------------------------------------
!*  (2): Contours in (x,z) plane.
!---------------------------------------------
!  b(ik)  x z
!
!   zr= zl +zleng -> space 2 cm 
      zl= zr +2.0 
      zr= zl +amin1(17.,zleng*xmax/zmax)
! 
      xl2= 2. 
      xr2= xl2 +zleng 
!
      call setscl (0.,0.,zmax,xmax,zl,xl2,zr,xr2,gdz,gdx, &
                   n1,char1, 2,'XZ', 0.6,                 &
                   1,'Z',0.6, 1,'X',0.6, 1)
!                                           <-- ierr= 5
      call number (zl-0.45,xl2-0.5,hh,zmin,0.,5)
      call number (zl-1.3, xr2-0.3,hh,zmax,0.,5)
      call number (zl-1.3, xl2-0.5,hh,xmax,0.,5)
!
      nxz= npx*npz2
      call daisho (b,nxz,wamin,wamax)
!
      ncontr= 7  ! 11
      call eqcntr (b,ww,npz2,npx,zl,xl2,zr,xr2,wamin,0.0,wamax, &
                   7000,ncontr,1) 
!
!---------------------
      call chart
!---------------------
      return
      end subroutine cplot3
!
!
!-----------------------------------------------------------------------
      subroutine pltarw (x,y,nx,ny,xl,yl,xr,yr,arx,ary,rarw,   &
                         anormx,anormy,shiftx,shifty,ndim,isel)
!-----------------------------------------------------------------------
!  << pltarw >>
!     1. function
!        (1) to arrow plot
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y    (nx,ny)    (i)       plotted value of x,y direction
!        (2) xl,yl,xr,yr       (i)       absolute coordinate value
!        (3) arx,ary           (i)
!        (4) rarw              (i)       length of the hat
!                                           (0.< <1.) *(axis length).
!        (5) anorm             (i)       normalization factor
!        (6) scl               (i)
!        (7) shiftx,shifty     (i)
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
      dimension  x(ndim),y(ndim)
!
      if (isel.eq.2) isel = 1
      nxy = nx*ny
      if (xr.lt.xl) return
      if (yr.lt.yl) return
      scl = float(nx)/float(ny)
      rar= 0.4*rarw
      gdx = (xr-xl)/float(nx)
      gdy = (yr-yl)/float(ny)
!
      call plot (xl,yl,3)
      call plot (xl,yr,2)
      call plot (xr,yr,2)
      call plot (xr,yl,2)
      call plot (xl,yl,2)
      call plot (xl,yl,3)
!
      if (anormx.lt.1.e-6) return
      if (anormy.lt.1.e-6) return
      if (anormx.eq.999.0) return
      if (anormy.eq.999.0) return
!
      if (isel.eq.2) goto 2000
!*
      do 3200 ij= 1,nxy
      iy= (ij-1)/nx+1
      ix= ij-nx*(iy-1)
!
      x1= ( ix-shiftx-arx*x(ij)/anormx)*gdx + xl
      x2= ( ix-shiftx+arx*x(ij)/anormx)*gdx + xl
      y1= ( iy-shifty-ary*y(ij)/anormy)*gdy + yl
      y2= ( iy-shifty+ary*y(ij)/anormy)*gdy + yl
!
      xa= rarw*x1 +(1.0-rarw)*x2
      ya= rarw*y1 +(1.0-rarw)*y2
!
      wx=-(y2-y1)*rar*scl
      wy= (x2-x1)*rar
!
      x3= xa + wx
      x4= xa - wx
      y3= ya + wy
      y4= ya - wy
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
      call plot (x3,y3,2)
      call plot (x2,y2,3)
      call plot (x4,y4,2)
 3200 continue
      return
!*
 2000 continue
      ram = 0.
      eps = 1.e-10
      anorm = max(anormx,anormy)
      if (anorm.gt.eps) ram = 1./anorm
      do 2200 ij= 1,nxy
      iy= (ij-1)/nx+1
      ix= ij-nx*(iy-1)
!
      ax1 = x(ij)/scl
      ay1 = y(ij)*scl
      asq = ax1*ax1 + ay1*ay1
      if (asq.gt.eps) then
        anrm = ram*sqrt( (x(ij)**2 + y(ij)**2) / asq )
      else
        anrm = 0.
      end if
      x1= ( ix-shiftx-arx*ax1*anrm)*gdx + xl
      x2= ( ix-shiftx+arx*ax1*anrm)*gdx + xl
      y1= ( iy-shifty-ary*ay1*anrm)*gdy + yl
      y2= ( iy-shifty+ary*ay1*anrm)*gdy + yl
!
      xa= rarw*x1 +(1.0-rarw)*x2
      ya= rarw*y1 +(1.0-rarw)*y2
!
      wx=-(y2-y1)*rar
      wy= (x2-x1)*rar
!
      x3= xa + wx
      x4= xa - wx
      y3= ya + wy
      y4= ya - wy
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
      call plot (x3,y3,2)
      call plot (x2,y2,3)
      call plot (x4,y4,2)
 2200 continue
!
      return
      end subroutine pltarw
!
!
! ++ Plotting subroutines ++
!-----------------------------------------------------------------------
      subroutine eqcntr (u,w,nx,ny,xl,yl,xr,yr,umin,ubund,umax,       &
                         ndim,lank,iwaku)
!-----------------------------------------------------------------------
!     1. function
!        (1) to draw tokosen
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) u       nx,ny     (i)       world value
!        (2) w       nx,ny     (i)       work array (real*4)
!        (3) xl,xr,yl,yr       (i)       absolute coordinate value
!        (4) umin,umax         (i)       height of max & min
!        (5) ubund             (i)       draw dash line (u < ubund)
!        (6) lank              (i)       number of draw lines
!        (7) iwaku             (i)       =1 : draw frame
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
      dimension u(ndim),w(ndim)
!
      if (nx.lt.2)  return
      if (ny.lt.2)  return
      if (xr.lt.xl) return
      if (yr.lt.yl) return
!
!------------------------------------------------
      if(umax.le.(1.000001*umin))  return
!------------------------------------------------
!
      nxm1 = nx -1
      nym1 = ny -1
!
      dx = (xr-xl)/float(nxm1)
      dy = (yr-yl)/float(nym1)
!
!
      do 10 i = 1,nx*ny
      w(i) = u(i) - umin
   10 continue
!
      if(iwaku.eq.1) then
        call plot(xl,yl,3)
        call plot(xr,yl,2)
        call plot(xr,yr,2)
        call plot(xl,yr,2)
        call plot(xl,yl,2)
        call plot(xl,yl,3)
      end if
!
      uld = float(lank+1) /(umax -umin)
      eps = 1.0e-3
!
!*-------------------------------------------------------
!*  draw contours cell by cell.
!*-------------------------------------------------------
!
      do 9000  ij= 1,nxm1*nym1
      j = (ij-1)/nxm1 + 1
      i = ij -(j-1)*nxm1
!
!*  four vertices around (i,j).
!
      i1 = i +nx*(j -1)
      i2 = i1 +1
      i3 = i1 +1 +nx
      i4 = i1 +nx
!
!*  normalized to (0., lank).
!
      u1 =  w(i1) *uld
      u2 =  w(i2) *uld
      u3 =  w(i3) *uld
      u4 =  w(i4) *uld
!
      k1 = ifix(u1)
      k2 = ifix(u2)
      k3 = ifix(u3)
      k4 = ifix(u4)
!
!*  used to judge possible intersection by contours.
!
      j1 = abs(k2-k1)
      j2 = abs(k3-k2)
      j3 = abs(k4-k3)
!
      u21= u2 -u1
      u32= u3 -u2
      u34= u3 -u4
      u41= u4 -u1
!
      ru21= 0.
      ru32= 0.
      ru34= 0.
      ru41= 0.
      if(abs(u21).gt.eps) ru21= 1./u21
      if(abs(u32).gt.eps) ru32= 1./u32
      if(abs(u34).gt.eps) ru34= 1./u34
      if(abs(u41).gt.eps) ru41= 1./u41
!
!* (1)
      if(j1.ne.0) then
         do 1000 ll = 1,j1
!
         u0 = float(ll) +float(min0(k1,k2))
         ujouge = u0/uld +umin
         if (ujouge.lt.ubund) then
             jouge = 4
         else
             jouge = 1
         end if
         if(abs(u21).lt.eps) go to 1000
!
             x1 = xl + dx*((u0-u1)*ru21 +float(i-1))
             y1 = yl + dy*float(j-1)
!
             if(((u3-u0)*(u2-u0)).gt.0. )            go to 1100
             if(((u0-u2).gt.0.).and.((u0-u4).gt.0.)) go to 1100
             if(abs(u32).lt.eps )                    go to 1100
!
             x2 = xl + dx*float(i)
             y2 = yl + dy*((u0-u2)*ru32 +float(j-1))
             call wdash (x1,y1,x2,y2,jouge)
!
 1100    continue
         if(((u4-u0)*(u3-u0)).gt.0. )  go to 1200
         if(((u1-u0)*(u3-u0)).gt.0. )  go to 1200
         if(((u2-u0)*(u4-u0)).gt.0. )  go to 1200
         if(abs(u34).lt.eps )          go to 1200
!
             x2 = xl + dx*( (u0-u4)*ru34 +float(i-1))
             y2 = yl + dy*float(j)
             call wdash (x1,y1,x2,y2,jouge)
!
 1200    continue
         if(((u1-u0)*(u4-u0)).gt.0. )            go to 1000
         if(((u0-u1).gt.0.).and.((u0-u3).gt.0.)) go to 1000
         if(abs(u41).lt.eps )                    go to 1000
!
             x2 = xl + dx*float(i-1)
             y2 = yl + dy*( (u0-u1)*ru41 +float(j-1))
             call wdash (x1,y1,x2,y2,jouge)
!
 1000    continue
      end if
!
!* (2)
      if(j2.ne.0) then
         do 2000 ll = 1,j2
!
         u0 = float(ll) + float(min0(k2,k3))
         ujouge = u0/uld + umin
         if (ujouge.lt.ubund) then
             jouge = 4
         else
             jouge = 1
         end if
         if(abs(u32).lt.eps)  go to 2000
!
              x1 = xl + dx*float(i)
              y1 = yl + dy*( (u0-u2)*ru32 +float(j-1) )
!
              if(((u4-u0)*(u3-u0)).gt.0. )            go to 2100
              if(((u0-u1).gt.0.).and.((u0-u3).gt.0.)) go to 2100
              if(abs(u34).lt.eps )                    go to 2100
!
              x2 = xl + dx*((u0-u4)*ru34 +float(i-1))
              y2 = yl + dy*float(j)
              call wdash (x1,y1,x2,y2,jouge)
!
 2100    continue
         if(((u1-u0)*(u4-u0)).gt.0. )  go to 2000
         if(((u1-u0)*(u3-u0)).gt.0. )  go to 2000
         if(((u2-u0)*(u4-u0)).gt.0. )  go to 2000
         if(abs(u41).lt.eps )          go to 2000
!
              x2 = xl + dx*float(i-1)
              y2 = yl + dy*((u0-u1)*ru41 +float(j-1))
              call wdash (x1,y1,x2,y2,jouge)
!
 2000    continue
      end if
!
!* (3)
      if(j3.ne.0) then
         do 3000 ll = 1,j3
!
         u0 = float(ll) + float(min0(k3,k4))
         ujouge = u0/uld + umin
         if (ujouge.lt.ubund) then
             jouge = 4
         else
             jouge = 1
         end if
         if(abs(u34).lt.eps)  go to 3000
!
              x1 = xl + dx*((u0-u4)*ru34 + float(i-1))
              y1 = yl + dy*float(j)
!
              if(((u1-u0)*(u4-u0)).gt.0. )            go to 3000
              if(((u0-u2).gt.0.).and.((u0-u4).gt.0.)) go to 3000
              if(abs(u41).lt.eps )                    go to 3000
!
              x2 = xl + dx*float(i-1)
              y2 = yl + dy*((u0-u1)*ru41 +float(j-1))
              call wdash (x1,y1,x2,y2,jouge)
!
 3000    continue
      end if
!
 9000 continue
      return
      end subroutine eqcntr
!
!
!-----------------------------------------------------------------------
      subroutine wdashl (wx ,wy ,nr,wybund,                           &
                         wminx,wminy,wmaxx,wmaxy,                     &
                         xl,yl,xr,yr,kline,iwaku)
!-----------------------------------------------------------------------
!
!     1. function
!        (1) to draw poly-line by multi line type
!     2. arguments   size     (i/o)     (meaning)
!        (1) wx,wy   (nr)      (i)       world coordinate value
!        (2) xl,xr,yl,yr       (i)       absolute coordinate scale value
!        (3) wminx,wminy       (i)       draw region
!            wmaxx,wmaxy
!                                        wmin>wmax : automatic control
!        (4) kline             (i)       pen type of 'wdash'
!        (5) iwaku             (i)       draw frame  (0:off ; 1:on)
!     3. called by
!             (** nothing **)
!     4. calls
!             (** wdash  **)
!-----------------------------------------------------------------------
      dimension wx(1),wy(1)
!
      if (nr.lt.2) return
      if (xr.le.xl) goto 9999
      if (yr.le.yl) goto 9999
      if (iwaku.eq.1) then
        call plot(xl,yl,3)
        call plot(xr,yl,2)
        call plot(xr,yr,2)
        call plot(xl,yr,2)
        call plot(xl,yl,2)
      end if
      gdx = (xr-xl)/(wmaxx - wminx)
      gdy = (yr-yl)/(wmaxy - wminy)
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
      if (wminy.le.wybund .and. wmaxy.gt.wybund) then
        y000 = (wybund-wminy)*gdy + yl
        if (y000.gt.yr) y000 = yr
        if (y000.lt.yl) y000 = yl
        call wdash(xl,y000,xr,y000,6)
      end if
      ir = 1
      xir =gdx*(wx(1)-wminx) + xl
      if (xir.lt.xl) then
        x1 = xl
      else if (xir.gt.xr) then
        x1 = xr
      else
        x1 = xir
      end if
      yir =gdy*(wy(1)-wminy) + yl
      if (yir.lt.yl) then
        y1 = yl
      else if (yir.gt.yr) then
        y1 = yr
      else
        y1 = yir
      end if
      k = - 1
      call plot ( x1 , y1 , 3 )
 1000 if (ir.eq.nr) goto 9000
      ir = ir + 1
      xir =gdx*(wx(ir)-wminx)+xl
      if (xir.lt.xl) then
        x2 = xl
      else if (xir.gt.xr) then
        x2 = xr
      else
        x2 = xir
      end if
      yir =gdy*(wy(ir)-wminy) + yl
      if (yir.lt.yl) then
        y2 = yl
      else if (yir.gt.yr) then
        y2 = yr
      else
        y2 = yir
      end if
      if(kline.lt.2) then
        go to 900
      else if(kline.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (kline.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (kline.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (kline.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (kline.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
      if(kline.lt.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 900
        if(rleng.lt.hh1) goto 900
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
  200   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hhh
          d=d+hh1
          goto 200
        end if
      else if (kline.eq.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 900
        if(rleng.lt.hh1) goto 900
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
  500   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hhh
          d=d+hh1
          goto 500
        end if
      else if (kline.eq.6) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 900
        if(rleng.lt.hh1) goto 900
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
  600   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hh5
          hh5 = hh6
          hh6 = hhh
          d=d+hh1
          goto 600
        end if
      end if
  900 call plot ( x2 , y2 , ( 5 + k ) / 2 )
      if (kline.gt.1 .and. kline.lt.7) k = - k
      x1 = x2
      y1 = y2
      goto 1000
 9000 call plot ( x2 , y2 , 3)
      return
 9999 continue
      call chart
      call symbol(1.0,11.0,0.2,' procedure = wdashl ',0.,20)
      call symbol(1.0,10.0,0.2,' abnormal world coordinate call',0.,31)
      call symbol(1.0,09.0,0.2,' wmaxx =',0.,8)
      call number(999.0,999.0,0.2,wmaxx,0.,2)
      call symbol(1.0,08.5,0.2,' wminx =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,08.0,0.2,' wmaxy =',0.,8)
      call number(999.0,999.0,0.2,wmaxy,0.,2)
      call symbol(1.0,07.5,0.2,' wminy =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,06.5,0.2,' xleft =',0.,8)
      call number(999.0,999.0,0.2,xl,0.,2)
      call symbol(1.0,06.0,0.2,' yleft =',0.,8)
      call number(999.0,999.0,0.2,yl,0.,2)
      call symbol(1.0,05.5,0.2,' xright=',0.,8)
      call number(999.0,999.0,0.2,xr,0.,2)
      call symbol(1.0,05.0,0.2,' yright=',0.,8)
      call number(999.0,999.0,0.2,yr,0.,2)
      write(6,*) '**********  abnormal world coordinate ********'
      write(6,*) '     procedure = wdashl'
      write(6,*) '    wmaxx =',wmaxx,' wminx = ',wminx
      write(6,*) '    wmaxy =',wmaxy,' wminy = ',wminy
      write(6,*) '    xl,yl,xr,yr =',xl,yl,xr,yr
      if(.true.) stop
      return
      end subroutine wdashl
!
!
!-----------------------------------------------------------------------
      subroutine wdash (x1,y1,x2,y2,ipen)
!-----------------------------------------------------------------------
!  << wdash  >>                      ver 2.00   16.mar.1990
!
!     1. function
!        (1) to draw line from (x1,y1) to (x2,y2) by wdash
!                            in absolute coordinate
!     2. arguments            (i/o)     (meaning)
!        (1) x1,x2,y1,y2       (i)       absolute coordinate value
!        (2) ipen              (i)       pen type of 'wdash'
!     3. called by
!             (** eqcntr  **)
!             (** wdashl  **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
!       ipen : meaning           - : 0.05 (cm)
!        1   :       line     -------------------
!        2   :  dash line     --- --- --- --- ---
!        3   :  dash line     -- -- -- -- -- -- --
!        4   :  dash line     - - - - - - - - - -
!        5   :  1 point dash  ---- - ---- - ---- -
!        6   :  2 point dash  --2.0-- - - --2.0--
!   otherwise:  line          ---------------------
!-----------------------------------------------------------------------
!
      call plot (x1,y1,3)
      k = -1
      if(ipen.lt.2) go to 999
!
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
!
      if(ipen.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (ipen.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (ipen.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (ipen.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (ipen.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
!
      if(ipen.lt.5) then
        rleng = sqrt((x2-x1)**2 + (y2-y1)**2)
        if(rleng.lt.hh1) goto 999
!
        costh = (x2 -x1)/rleng
        sinth = (y2 -y1)/rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
!
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
!
  200   if(d.le.rleng) then
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
!
        k = - k
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
        d=d +hh1
        goto 200
        end if
!
      else if (ipen.eq.5) then
        rleng = sqrt((x2-x1)**2 + (y2-y1)**2)
        if(rleng.lt.hh1) goto 999
!
        costh = ( x2 -x1 )/rleng
        sinth = ( y2 -y1 )/rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
!
  500   if(d.le.rleng) then
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
        k = - k
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
        d=d+hh1
        goto 500
        end if
!
      else if (ipen.eq.6) then
        rleng = sqrt((x2-x1)**2 + (y2-y1)**2)
        if(rleng.lt.hh1) goto 999
!
        costh = ( x2 - x1 )/rleng
        sinth = ( y2 - y1 )/rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
!
  600   if(d.le.rleng) then
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot (x,y,(5+k)/2)
        k = - k
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
        d=d+hh1
        goto 600
        end if
      end if
!
  999 call plot (x2,y2,(5+k)/2)
      call plot (x2,y2,3)
!
      return
      end subroutine wdash
!
!
!-----------------------------------------------------------------------
      subroutine setscl (wminx,wminy,wmaxx,wmaxy, xl,yl,xr,yr,gdx,gdy, &
                 n1,char1, n2,char2, hight1,                           &
                 nnx,charx,hightx, nny,chary,highty, iwaku)
!-----------------------------------------------------------------------
!  << setscl >>                   /char1/
!                          wmaxy  +--------------------+  (xl,yl)
!                               y |             (xr,yr)|  (xr,yr) on 0
!                               r |                    |
!                               a |                    |
!    (wminx,wminy)              h |                    |
!    (wmaxx,wmaxy) on is        c |                    |
!                                 |(xl,yl)             |
!                          wminy  +--------+--+--------+
!                                 wminx  /charx/       wmaxx
!-----------------------------------------------------------------------
!
!     setscl
!
!     1. function
!        (1) to scale the graphics by calcomp specifications
!     2. arguments            (i/o)     (meaning)
!        (1) wminx,wmaxx,
!            wminy,wmaxy       (i)       world coordinate value
!        (2) xl,xr,yl,yr       (i)       absolute coordinate value
!        (3) gdx,gdy           (o)       scaling factor of coordinate
!                                        from world to absolute
!        (4) char1,charx,cahry (i)       title on graph,x-axis,y-axis
!        (5) iwaku             (i)       draw frame (0:off ; 1:on)
!                                         999 : write out only title,
!                                                    not draw otherwise
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!             (** symbol **)
!             (** number **)
!-----------------------------------------------------------------------
!  NEC style
      character*8  char1(*),char2(*),charx(*),chary(*)
!
      if (wmaxx.le.wminy) goto 9999
      if (wmaxx.le.wminy) goto 9999
      if (xr.le.xl)       goto 9999
      if (yr.le.yl)       goto 9999
!
      gdx= (xr-xl)/(wmaxx-wminx)
      gdy= (yr-yl)/(wmaxy-wminy)
!
      xc = 0.5*( xr + xl )
      yc = 0.5*( yr + yl )
!
      if (n1 .gt.0) then
        if (hight1.gt.0) then
          xs1= xc -0.5*n1*hight1
          xs2= xs1 +(n1+1)*hight1
          call symbol(xs1,yr+0.1,hight1,char1(1),0.,n1)
          call symbol(xs2,yr+0.1,hight1,char2(1),0.,n2)
        end if
      end if
!-----------------------------------------------------------------------
      if (iwaku.eq.999) return
!-----------------------------------------------------------------------
!
      if (iwaku.eq.1) then
        call plot (xl,yl,3)
        call plot (xl,yr,2)
        call plot (xr,yr,2)
        call plot (xr,yl,2)
        call plot (xl,yl,2)
        call plot (999.,999.0,3)
      end if
!
      if (nnx.gt.0) then
        if (hightx.gt.0) then
          call symbol(xc-0.5*hightx*nnx,yl-0.5,hightx,charx(1),0.,1)
          do nnx1=2,nnx
          call symbol(999.0,999.0,hightx,charx(nnx1),0.,1)
          end do
        end if
      end if
      if (nny.gt.0) then
        if (highty.gt.0) then
          call symbol(xl-0.5,yc-0.5*highty*nny,highty,chary(1),0.,1)
          do nny1=2,nny
          call symbol(999.0,999.0,highty,chary(nny1),0.,1)
          end do
        end if
      else if(nny.lt.0) then
        if (highty.gt.0) then
          call symbol(xc-0.5*highty*nny,yc,highty,chary(1),0.,1)
          do nny1=2,nny
          call symbol(999.0,999.0,highty,chary(nny1),0.,1)
          end do
        end if
      end if
!
      return
!
!-----------------------------------------------------------------------
!
 9999 continue
      write(6,*) '**********  abnormal world coordinate ********'
      write(6,*) '      '
      write(6,*) '    wmaxx =',wmaxx,' wminx = ',wminx
      write(6,*) '    wmaxy =',wmaxy,' wminy = ',wminy
      write(6,*) '    xl,yl,xr,yr =',xl,yl,xr,yr
      write(6,*) '    fctr  =',fctr
      write(6,*) '      '
      call chart
      call symbol(1.0,10.0,0.2,' abnormal world coordinate call',0.,31)
      call symbol(1.0,09.0,0.2,' wmaxx =',0.,8)
      call number(999.0,999.0,0.2,wmaxx,0.,2)
      call symbol(1.0,08.5,0.2,' wminx =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,08.0,0.2,' wmaxy =',0.,8)
      call number(999.0,999.0,0.2,wmaxy,0.,2)
      call symbol(1.0,07.5,0.2,' wminy =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,07.0,0.2,' fctr  =',0.,8)
      call number(999.0,999.0,0.2,fctr,0.,2)
      call symbol(1.0,06.5,0.2,' xleft =',0.,8)
      call number(999.0,999.0,0.2,xl,0.,2)
      call symbol(1.0,06.0,0.2,' yleft =',0.,8)
      call number(999.0,999.0,0.2,yl,0.,2)
      call symbol(1.0,05.5,0.2,' xright=',0.,8)
      call number(999.0,999.0,0.2,xr,0.,2)
      call symbol(1.0,05.0,0.2,' yright=',0.,8)
      call number(999.0,999.0,0.2,yr,0.,2)
      stop
      return
      end subroutine setscl
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      x0(isc)= x00
      y0(isc)= y00
      dxi(isc)= 1./dx
      dyi(isc)= 1./dy
!
      xl(isc)= xcm
      yl(isc)= ycm
!
      return
      end subroutine scalex
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end subroutine plotl
!
!
!-----------------------------------------------------------------------
      subroutine values (x,y,height,val,theta,ifmat)
!-----------------------------------------------------------------------
!  << values >>
!     1. function
!        (1) to draw variable
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y               (i)       absolute coordinate value
!        (2) height            (i)       draw out size on paper
!        (3) val               (i)       variable
!        (4) theta             (i)       angle
!        (5) ifmat             (i)       format type
!     3. called by
!             (** nothing **)
!     4. calls
!             (** number **)
!             (** symbol **)
!-----------------------------------------------------------------------
!        ifmat = (n100)*100 + keta
!        n100 = 0 : integer format
!        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!        n100 = 2 : e format ::
!        n100 = 3 : power of ten format
!        n100 = othewise : not write out
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      real(C_float) val
      character        chr13*13,chr12*12,chr3*3
      character(len=1) minus,zero,blank
!
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!
      if (n100.eq.0) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,keta)
      else if (n100.eq.2) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
        call symbol(999.,999.,height,chr3,theta,numsy1)
      else if (n100.eq.3) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        sint = sin(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
!
!                                             *******************
!                                             ** exponent part **
!                                             *******************
!
        h2 = height * 5./7.
        x1 = (numsym+1)* height * ratio
        y1 = height * 4./7.
        if (abs(theta).lt.1e-04) then
          x1 = x + x1
          y1 = y + y1
        else
          x2 =     x1 * cost - y1 * sint
          y1 = y + x1 * sint + y1 * cost + h2*cost
          x1 = x + x2                    - h2*sint
        end if
        call symbol(x1,y1,h2,chr3,theta,numsy1)
      end if
      return
      end subroutine values
!
!
!-----------------------------------------------------------------------
      subroutine daisho(x ,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      dimension x(1)
!
      xmax1= x(1)
      xmin1= x(1)
      do 100 i=2,nx
      xmax1= max(xmax1,x(i) )
      xmin1= min(xmin1,x(i) )
  100 continue
      return
      end subroutine daisho
!
!
!***************************************************************
!*     this program package generates a unix postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      postscript header by fortran
!        t. ogino (nagoya university) february 27, 1992
!      modified to conform gsipp commands
!        motohiko tanaka (nifs)       november 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!     this ps-adobe-2.0 header allows us full paging features in
!     the ghostview.  to scroll up the page (backward), click the 
!     page number and press two buttons of mouse simultaneously.
!
!     consult: a.saitou (kyoto u.)  the definition of /@eop  
!    needs stroke for line drawings (not in the tex header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  this is an adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!ps-adobe-2.0',/       &
              '%%pages: (atend)',/     &
              '%%pageorder: ascend',/  &
              '%%endcomments',/        &
              '%%begindocument')
!
!%%%%%%%%%%%%%%%%%%% procedure defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%boundingbox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/times-roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/             &
             '/se {setfont} bind def',/               &
             '/ro {rotate}  bind def',/               &
             '/tl {translate} bind def',/             &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/  &
             '{erasepage newpath initgraphics',/                & 
             '/saveimage save def',/                            &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/  &
             '{stroke showpage',/                     &
             ' saveimage restore',/                   &
             '} bind def')
!
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/  &
             ' /end load def')
!
      write(77,27) 
   27 format('/dir 0 def')
!
      write(77,29) 
   29 format('/s             % string s -- show the string',/  &
             '{dir 1 eq',/                                     &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/  &
             ' {show} ifelse',/   &
             '} bind def')
!
      write(77,31)
   31 format('%%enddocument',/        &
             '%%endprolog',/          &
             '%%beginsetup',/         &
             '/resolution 300 def',/  &
             '/#copies 1 def',/       &
             '%%endsetup')
!
!%%%%%%%%%%%%%%%%%%% end of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%beginpagesetup',/  &
              '%%endpagesetup',/    &
              '@bop')
!
!*  set magnifying factor (gsipp to sun coordinate).
!   rotate and translate to output on a4-l paper.
!      left corner ...... (  0.,  0.)
!      right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  if nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end subroutine gopen
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end subroutine gclose
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end subroutine plote
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%pagetrailer    % need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%beginpagesetup',/  &
                 '%%endpagesetup',/    &
                 '@bop')
!
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!-----------------------------------------------------
!      first cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   frames 2-4:
!
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
!
       return
       end subroutine chart
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end subroutine factor
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end subroutine newpen
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
       return
       end subroutine linee
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
       return
       end subroutine plot
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(C_float) x0,y0,x,y,h
       real(C_float) h0,ang
       character(*)   isymb   !!! NEC: it must hitsuyo
       integer(C_INT) n0,n,i
!
       character    ica*80,ich(80)*1
       equivalence (ica,ich(1))
!
       x= x0
       y= y0
       h= h0
       n= n0
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
!
       return
       end subroutine symbol 
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(C_float)  x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character      isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       if(abs(anu).gt.1.e1 .or.  &
          abs(anu).lt.1.e-1) then
        write(isymb,'(1pe9.2)') anu    ! e9.2
       else
        write(isymb,'(f7.2)') anu      ! f7.2
       end if
!
       if(.true.) go to 300
       if(abs(anu).lt.10000.) then  !! 5 digits
         if(abs(anu).gt.0.1) then
           write(isymb,'(f6.1)') anu        ! f6.1
         else
           if(abs(anu).gt.0.001) then 
             write(isymb,'(f6.3)') anu      ! f6.3
           else
             if(abs(anu).gt.0.001) then 
               write(isymb,'(1pe9.2)') anu  ! e9.2
             else
               write(isymb,'(f6.1)') anu    ! f6.1  0.0
             end if
           end if
         end if
!
       else
         if(abs(anu).lt.100000.) then
           write(isymb,'(f7.1)') anu       ! f7.1
         else
           write(isymb,'(1pe9.2)') anu     ! e9.2
         end if
       end if
  300  continue
!
       write(77,*) '(',isymb,') s'
!
       return
       end subroutine number
!
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(C_float)  x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character      isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
      if(n0.eq.1) write(isymb,41) anu
      if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
!
       write(77,*) '(',isymb,') s'
!
       return
       end subroutine number2 
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +abs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end subroutine sunscl
