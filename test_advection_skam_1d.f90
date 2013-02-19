program execute
  implicit none

  logical :: doposlimit, domonlimit, dosellimit, dosemilagr, dopcm, &
       dowenosplit, outputtestname
  real(kind=8) :: lambdamax, epslambda
  integer :: nmethod

  dosemilagr = .false.
  dopcm = .false.

  call test1dweno(1,30,2,4,0.d0,-1,0.5d0)
  call test1dweno(2,30,2,4,0.d0,-1,0.5d0)
  call test1dweno(3,30,2,4,0.d0,-1,0.5d0)
  call test1dweno(4,30,2,4,0.d0,-1,0.5d0)

contains

  subroutine test1dweno(ntest,nx0,nscale,nlev,stretch,noutput,maxcfl)
    implicit none

    integer, intent(in) :: ntest, nx0, nscale, nlev, noutput
    real (kind=8), intent(in) :: maxcfl, stretch

    real (kind=8), allocatable, dimension(:,:) :: q, q0, dqdt, u
    real (kind=8), allocatable, dimension(:) :: dx, qi, qiface, x, xf, upos
    real (kind=8), allocatable, dimension(:) :: lambda, monlimit
   
    real(kind=8), dimension(nlev) :: e1, e2, ei
    integer, parameter :: ny = 1
    integer :: nx, nstep, i ,j , n, p, ierr, imethod, nout, npad, n1, n2
    real (kind=8) :: dt, tfinal, time, pi
    character(len=40) :: cdf_out 
    character(len=8) :: outdir

    pi = atan2(0.d0,-1.d0)

    if(nlev.lt.1) STOP 'nlev should be at least 1 in test1dweno'

    n1 = 5 ! Selective PPM FCT
    n2 = 11 ! Selective PPM FCT

    do imethod = n1,n2,6
       dosellimit = .false.
       domonlimit = .false.
       doposlimit = .false.
       dopcm = .false.  ! use PPM by default
       dowenosplit = .false.  ! use PPM by default
       select case(imethod)
       case(1)
          write(*,*) 'PPM FCT, No limiting'
          outdir = 'pfctnon/'
          nmethod = 1
       case(2)
          write(*,*) 'PPM FCT, Positive'
          doposlimit = .true.
          outdir = 'pfctpos/'
          nmethod = 2
       case(3)
          write(*,*) 'PPM FCT, Monotonic'
          domonlimit = .true.
          outdir = 'pfctmon/'
          nmethod = 3
       case(4)
          write(*,*) 'PPM FCT, Monotonic, Positive'
          domonlimit = .true.
          doposlimit = .true.
          outdir = 'pfctpmn/'
          nmethod = 4
       case(5)
          write(*,*) 'PPM FCT, Selective'
          dosellimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'pfctsel/'
          nmethod = 5
       case(6)
          write(*,*) 'PPM FCT, Selective, Positive'
          dosellimit = .true. 
          doposlimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'pfctpse/'
          nmethod = 6
       case(7)
          write(*,*) 'PPM PMOD, No limiting'
          outdir = 'ppmdnon/'
          nmethod = 21
       case(8)
          write(*,*) 'PPM PMOD, Positive'
          doposlimit = .true.
          outdir = 'ppmdpos/'
          nmethod = 22
       case(9)
          write(*,*) 'PPM PMOD, Monotonic'
          domonlimit = .true.
          outdir = 'ppmdmon/'
          nmethod = 23
       case(10)
          write(*,*) 'PPM PMOD, Monotonic, Positive'
          domonlimit = .true.
          doposlimit = .true.
          outdir = 'ppmdpmn/'
          nmethod = 24
       case(11)
          write(*,*) 'PPM PMOD, Selective'
          dosellimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdsel/'
          nmethod = 25
       case(12)
          write(*,*) 'PPM PMOD, Selective, Positive'
          dosellimit = .true. 
          doposlimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdpse/'
          nmethod = 26
       case(29)
          write(*,*) 'PPM PMOD, Selective, Clean'
          dosellimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdsec/'
          nmethod = 81
       case(30)
          write(*,*) 'PPM PMOD, Selective, Positive, Clean'
          dosellimit = .true. 
          doposlimit = .true.
          lambdamax=20.
          epslambda=1.e-8
          outdir = 'ppmdpsc/'
          nmethod = 82
       end select

       if(dopcm.and.domonlimit) then
          npad = 6
       elseif(dosellimit.or.doposlimit) then
          npad = 4
       else
          npad = 3
       end if

       outputtestname = .true.
    do p = 1,nlev
       nx = nx0*nscale**(p-1)
       allocate(q(1-npad:nx+npad,1:ny),q0(1:nx,1:ny),dqdt(1:nx,1:ny), &
            u(3-npad:nx-3+npad,1:ny), &
            dx(1:nx),x(1:nx),xf(0:nx),lambda(0:nx),monlimit(0:nx),STAT=ierr)
       
       ! set up grid
       do i = 1,nx
          dx(i) = 1.d0 + stretch*COS(2.d0*pi*DBLE(i)/DBLE(nx))
       end do
       dx = dx/SUM(dx) ! normalize for unit length
       xf(0) = 0.d0
       do i = 1,nx
          xf(i) = xf(i-1) + dx(i)
          x(i) = 0.5d0*(xf(i-1)+xf(i))
       end do

       ! set up scalar and velocity fields
       lambda = 0.d0
       monlimit = 0.d0
       q = 0.d0
       dqdt = 0.d0
       u = 0.d0
       call init1d(ntest,nx,q0,u(0:nx,1),x,tfinal,cdf_out)
       q(1:nx,:) = q0

       if (ABS(MAXVAL(dx)-MINVAL(dx)).gt.10.d0*EPSILON(dx)) then
          cdf_out = 'stretch_' // cdf_out
       end if

       cdf_out = outdir // cdf_out

       ! set up time step
       time = 0.d0
       if(noutput.eq.-1) then
          if(p.eq.1) then
             nstep = ceiling(tfinal*MAXVAL(ABS(u))/MINVAL(dx)/maxcfl)
             nout = nstep
          else
             nstep=nout &
               *ceiling(tfinal*MAXVAL(ABS(u))/MINVAL(dx)/maxcfl/DBLE(nout))
       end if
       else
          nstep=noutput &
               *ceiling(tfinal*MAXVAL(ABS(u))/MINVAL(dx)/maxcfl/DBLE(noutput))
          nout = noutput
       end if
       if((dowenosplit).and.(noutput.ne.-1)) nstep=nstep*nscale**(2.d0*DBLE(p-1)/3.d0)
       dt = tfinal/dble(nstep)



       if (p.eq.1) call output1d(q(1:nx,1),u(0:nx,1),lambda,monlimit,nx,x,xf,tfinal,-1,cdf_out,nout)
       
       call output1d(q(1:nx,1),u(0:nx,1),lambda,monlimit,nx,x,xf,time,0,cdf_out,p)

       do n = 1,nstep
          call skamstep_1d(q,dqdt,u,lambda,monlimit,nx,ny,npad,dt,dx)
          time = time + dt
          if (mod(n,nstep/nout).eq.0) then
             call output1d(q(1:nx,1),u(0:nx,1),lambda,monlimit,nx,x,xf,time,1,cdf_out,p)
          end if
       end do

       e1(p) = SUM(ABS(q(1:nx,:)-q0))/DBLE(nx)
       e2(p) = SQRT(SUM((q(1:nx,:)-q0)**2)/DBLE(nx))
       ei(p) = MAXVAL(ABS(q(1:nx,:)-q0))
       if (p.eq.1) then
          write(*,*) &
               '  nx       E1          E2         Einf   convergence rate overshoot  undershoot'
          write(*,990) nx, e1(p), e2(p), ei(p), 0.0, 0.0, 0.0, &
               MINVAL(q0)-MINVAL(q), MAXVAL(q)-MAXVAL(q0)
       else
          write(*,990) nx, e1(p), e2(p), ei(p), &
               -log(e1(p)/e1(p-1))/log(dble(nscale)), &
               -log(e2(p)/e2(p-1))/log(dble(nscale)), &
               -log(ei(p)/ei(p-1))/log(dble(nscale)), &
               MINVAL(q0)-MINVAL(q), MAXVAL(q)-MAXVAL(q0)
       end if

       if (p.eq.nlev) call output1d(q(1:nx,1),u(0:nx,1),lambda,monlimit,nx,x,xf,time,2,cdf_out,1)

       deallocate(q,q0,dqdt,u,dx,x,xf,lambda,monlimit,STAT=ierr)
    end do
 end do

990    format(i6,3e12.4,3f5.2,2e12.4)
992    format(12f8.4)
  end subroutine test1dweno

  subroutine init1d(ntest,nx,q,u,x,tfinal,cdf_out)
    implicit none
    integer, intent(in) :: ntest,nx
    real(kind=8), dimension(1:nx), intent(out) :: q
    real(kind=8), dimension(0:nx), intent(out) :: u
    real(kind=8), dimension(1:nx), intent(in) :: x
    real(kind=8), intent(out) :: tfinal
    character(len=40), intent(out) :: cdf_out 

    integer :: i
    real(kind=8) :: pi, a, z, alpha, beta, delta
    real(kind=8), dimension(1:nx) :: xx

    pi = datan2(0.d0,-1.d0)

    u = 1.d0 ! uniform velocity
    tfinal = 1.d0

    select case(ntest)
    case (1) ! sine wave
       q = sin(2.d0*pi*x)
       cdf_out = 'weno1d_sine_adv.nc'
    case (2)
       q = 0.5d0*(sin(6.d0*pi*x)+sin(8.d0*pi*x))
       cdf_out = 'weno1d_twowave_adv.nc'
    case (12)
       u = -u
       q = 0.5d0*(sin(6.d0*pi*x)+sin(8.d0*pi*x))
       cdf_out = 'weno1d_twowave_adv.nc'
    case (3)
       q = MAX(0.d0,0.5d0*(sin(6.d0*pi*x)+sin(8.d0*pi*x)))
       cdf_out = 'weno1d_twowave_pos_adv.nc'
    case (4)
       q = 0.d0
       where ((x.gt.0.4d0).and.(x.lt.0.6d0))
          q = 1.d0
       end where
       cdf_out = 'weno1d_tophat_adv.nc'
    case (5) ! sine^4 wave
       q = sin(pi*x)**4
       cdf_out = 'weno1d_sine4_adv.nc'
    case (6) ! smooth top hat
       q = 1.d0/(1.d0 + exp(80.d0*(abs(x-0.5d0) - 0.15d0)))
       cdf_out = 'weno1d_smooth_tophat_adv.nc'
    case (7) ! sine wave -- twenty times through domain
       tfinal = 20.d0
       q = sin(2.d0*pi*x)
       cdf_out = 'weno1d_sine_20xadv.nc'
    case (8) ! top hat -- five times through domain
       tfinal = 5.d0
       q = 0.d0
       where ((x.gt.0.4d0).and.(x.lt.0.6d0))
          q = 1.d0
       end where
       cdf_out = 'weno1d_tophat_5xadv.nc'
    case (11) ! top hat -- five times through domain
       tfinal = 5.d0
       q = 1.d0
       where ((x.gt.0.4d0).and.(x.lt.0.6d0))
          q = 2.d0
       end where
       cdf_out = 'weno1d_tophat2_5xadv.nc'
    case (9) ! test profile with discontinuities in 0th, 1st and 2nd deriv
       tfinal = x(2) - x(1)
       where (abs(x-0.5d0).ge.0.375d0)
          q = 0.d0
       elsewhere ((x.gt.0.5d0).and.(x.le.0.875d0))
          q = 0.5d0
!!$          q = -1.d0 + 4.d0*x - 2.d0*x**2
       elsewhere
          q = 0.5d0*(1.d0 - sin(4*pi*x))
       end where
       cdf_out = 'weno1d_testprof_adv.nc'
    case (10) ! test profile from Jiang and Shu JCP (1996)
       ! four times through domain
       tfinal = 4.d0

       ! rescale x for their x \in [-1,1] domain
       xx = 2.d0*x - 1.d0

       ! some constants
       a = 0.5d0
       z = -0.7d0
       delta = 0.005d0
       alpha = 10.d0
       beta = log(2.d0)/36.d0/delta**2

       where (abs(xx+0.7d0).le.0.1d0)
          ! gaussian (smooth q)
          q = (1.d0/6.d0) &
               *(exp(-beta*(xx-z+delta)**2) &
               + exp(-beta*(xx-z-delta)**2) &
               + 4.d0*exp(-beta*(xx-z)**2))
       elsewhere (abs(xx+0.3d0).le.0.1d0)
          ! top hat (q discontinuous)
          q = 1.d0
       elsewhere (abs(xx-0.1d0).le.0.1d0)
          ! sharp cone (q' discontinuous)
          q = 1.d0 - 10.d0*abs(xx-0.1d0)
       elsewhere (abs(xx-0.5d0).le.0.1d0)
          ! curved cone (singularity in q')
          q = (1.d0/6.d0) &
               *(sqrt(max(0.d0,1.d0-alpha**2*(xx-a+delta)**2)) &
               + sqrt(max(0.d0,1.d0-alpha**2*(xx-a-delta)**2)) &
               + 4.d0*sqrt(max(0.d0,1.d0-alpha**2*(xx-a)**2)))
       elsewhere (abs(xx+0.3d0).le.0.1d0)
       elsewhere
          q = 0.d0
       end where
       cdf_out = 'weno1d_jiangshu_adv.nc'
    end select

    if(outputtestname) write(*,*) '================= ', TRIM(cdf_out),  ' ================='
    outputtestname = .false.

  end subroutine init1d
  
  subroutine skamstep_1d(q,dqdt,u,lambda,monlimit,nx,ny,npad,dt,dx)
    !  use forward-in-time Skamarock (2006) scheme

    implicit none

    integer, intent(in) :: nx, ny, npad
    real (kind=8), intent(in)  :: dt

    real (kind=8), dimension(1-npad:nx+npad,1:ny), intent(inout) :: q
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: dqdt ! extra tendency
    real (kind=8), dimension(3-npad:nx-3+npad,ny), intent(inout) :: u
    real (kind=8), dimension(1:nx), intent(in) :: dx

    integer, parameter :: nmaxcfl = 1

    ! local variables, one-dimensional arrays in x-direction
    real (kind=8), dimension(0:nx) :: monlimit
    real (kind=8), dimension(0:nx) :: lambda
    real (kind=8), dimension(0:nx+1) :: rhoq1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: q1dx
    real (kind=8), dimension(-2:nx+2) :: rhouh1d
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rho1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rhoprime1dx

    ! local variables
    real (kind=8), dimension(1:nx,1:ny) :: q0,q1
    real (kind=8), dimension(0:nx+1) :: rhoq
    real (kind=8), dimension(-2:nx+2) :: rhou
    real (kind=8), dimension(-1:nx+2) :: rhoprime
    real (kind=8), dimension(-1:nx+2) :: rho
    real(kind=8) :: scale
    integer, dimension(2) :: xbctype
    real (kind=8), dimension(2) :: fxbc ! specified bndry flux 
    real (kind=8), dimension(0:nx  ,ny) :: xflx

    integer i,j, npad2, npadq, npadrp, npadrho

    ! set up boundary types, fluxes, limiting flags
!!$    dosellimit = .false. ! selective monotonic limiting
!!$    domonlimit = .true. ! monotonic limiting everywhere
!!$    doposlimit = .false. ! positive definite limiting
    xbctype = 2           ! periodic boundary conditions
    fxbc = 0.d0          ! fixed fluxes unused, but initialize anyways.

    scale = 1.

    q1dx = 0.d0
    rhoprime1dx = 0.d0
    rho1dx = 0.d0
    rhouh1d = 0.d0

    do j = 1,ny
          q1dx(1:nx)     = q(1:nx,j)
          rhoq1dx(1:nx)  = q(1:nx,j)*dx(:)

          rhoprime1dx(1:nx) = dx(:)
          rho1dx(1:nx)      = dx(:)

          ! fix periodic bcs on rhoq
          rhoq1dx(0) = rhoq1dx(nx)
          rhoq1dx(nx+1) = rhoq1dx(1)

          ! fix periodic bcs on q
          npadq = npad+nmaxcfl-1
          q1dx(1-npadq:0) = q1dx(nx+1-npadq:nx)
          q1dx(nx+1:nx+npadq) = q1dx(1:npadq)

          ! fix periodic bcs on rhoprime
          npadrp = nmaxcfl+2
          rhoprime1dx(1-npadrp:0) = rhoprime1dx(nx+1-npadrp:nx)
          rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(1:npadrp)

          ! fix periodic bcs on rho
          npadrho = npad-1 + (nmaxcfl+1)/2
          rho1dx(1-npadrho:0) = rho1dx(nx+1-npadrho:nx)
          rho1dx(nx+1:nx+npadrho) = rho1dx(1:npadrho)

       ! 
       rhouh1d(0:nx) = u(0:nx,j)
       rhouh1d(-2:-1) = rhouh1d(nx-1:nx)
       rhouh1d(nx+1:nx+2) = rhouh1d(1:2) 

       call ppmwrap(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx,xflx(0,j),dt, &
            nx,npad,nmaxcfl,xbctype,fxbc, &
            dosemilagr,dosellimit,domonlimit,doposlimit, &
            dopcm, dowenosplit, &
            scale,nmethod,lambdamax,epslambda, &
            lambda,monlimit)

!!$       npad2 = npad + nmaxcfl
!!$       if(dosellimit) then
!!$          call ppmsweep_select(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx, &
!!$               xflx(0,j), dt,nx,npad2,xbctype,fxbc,doposlimit,scale, &
!!$               lambda,monlimit)
!!$       elseif(dowenosplit) then
!!$          call wenosweep(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx, &
!!$               xflx(0,j), dt,nx,npad2,xbctype,fxbc,doposlimit,scale, &
!!$               lambda,monlimit)
!!$       else
!!$          call ppmsweep(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx, &
!!$               xflx(0,j), dt,nx,npad2,xbctype,fxbc,domonlimit,doposlimit, &
!!$               lambda,monlimit)
!!$       end if

       ! update solution
       q(1:nx,j) = rhoq1dx(1:nx)/rhoprime1dx(1:nx)

    end do

  end subroutine skamstep_1d

  subroutine output1d(q,u,lam,mnl,nx,x,xf,time,status,cdf_out,ilevel)
    implicit none

    !     VMS include statement (on UW Vax machines)
    !     include 'local_root:[include]netcdf.inc'

    integer, intent(in) :: nx, status, ilevel
    real (kind=8), intent(in) :: time
    real (kind=8), dimension(1:nx), intent(in) :: q, x
    real (kind=8), dimension(0:nx), intent(in) :: u, xf, lam, mnl
    character(len=40), intent(in) :: cdf_out ! Name of the netCDF file

    integer msize,n1d,i,j,ierr,idq,idu,idv,idt,idl,idm

    real (kind=4) :: time4
    real (kind=8) :: tfcn
    real (kind=4), dimension(nx)   :: temp
    real (kind=4), dimension(0:nx) :: temp1
    real (kind=4), dimension(:), allocatable :: temp2

    !  netCDF file declaration
    integer cdfid               ! ID for the netCDF file to be created
    character(len=8):: xname,xfname,tname,qname,uname,nxname,nxfname,ntname, &
         lname, mname

    !  netCDF variables declaration

    integer :: &
    &    start(4) &             ! netCDF stuff;
    &   ,count(4) &             ! netCDF stuff;
    &   ,vdim(4)  &             ! netCDF stuff;
    &   ,ndim     &             ! netCDF stuff;
    &   ,idx,idy,idxf,idyf,idnx,idny,idnxf,idnyf,idnt
    data       &             
    &    start /1, 1, 1, 1/, count /1, 1, 1, 1/

    integer nxp1
    save cdfid, idq, idu, idl, idm, idt, idnt

    !     UNIX include statement (on UW Unix machines)
    include 'netcdf.inc'

   !  Begin

    if (status.eq.-1) then
       !                  Create netCDF file
       ierr = NF_CREATE(cdf_out, NF_CLOBBER, cdfid)

       ! set up time dimension
       tname = 'time'
       ntname = 'nt'
       ierr = NF_REDEF(cdfid)
       ierr = NF_DEF_DIM(cdfid,TRIM(ntname),ilevel+1,idnt)
       ierr = NF_DEF_VAR(cdfid,TRIM(tname), NF_FLOAT, 1, idnt, idt)
       ierr = NF_ENDDEF(cdfid)

       ALLOCATE(temp2(1:ilevel+1),STAT=ierr)
       temp2(1) = 0.d0
       do i = 1,ilevel
          temp2(i+1) = DBLE(i)*time/DBLE(ilevel)
       end do
       ierr = NF_PUT_VAR_REAL(cdfid, idt, temp2)
       DEALLOCATE(temp2)
       
       return

    elseif (status.eq.0) then

       write(xname,'(a1,i1)') 'x', ilevel
       write(xfname,'(a2,i1)') 'xf', ilevel
       write(nxname,'(a2,i1)') 'nx', ilevel
       write(nxfname,'(a3,i1)') 'nxf', ilevel
       write(qname,'(a1,i1)') 'Q', ilevel
       write(uname,'(a1,i1)') 'U', ilevel
       write(lname,'(a6,i1)') 'LAMBDA', ilevel
       write(mname,'(a6,i1)') 'MONLMT', ilevel

       nxp1=nx+1

       ierr = NF_REDEF(cdfid)

       ierr = NF_DEF_DIM(cdfid,TRIM(nxname),nx,idnx)
       ierr = NF_DEF_DIM(cdfid,TRIM(nxfname),nxp1,idnxf)

       ierr = NF_DEF_VAR(cdfid,TRIM(xname), NF_FLOAT, 1, idnx, idx)
       ierr = NF_DEF_VAR(cdfid,TRIM(xfname), NF_FLOAT, 1, idnxf, idxf)

       ndim = 2
       vdim(1) = idnx
       vdim(2) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(qname),NF_FLOAT,ndim,vdim,idq)

       ndim = 2
       vdim(1) = idnxf
       vdim(2) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(uname),NF_FLOAT,ndim,vdim,idu)

       ndim = 2
       vdim(1) = idnxf
       vdim(2) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(lname),NF_FLOAT,ndim,vdim,idl)

       ndim = 2
       vdim(1) = idnxf
       vdim(2) = idnt
       ierr = NF_DEF_VAR(cdfid,TRIM(mname),NF_FLOAT,ndim,vdim,idm)

       ierr = NF_ENDDEF(cdfid)

       temp = x
       ierr = NF_PUT_VAR_REAL(cdfid, idx, temp)

       temp1 = xf
       ierr = NF_PUT_VAR_REAL(cdfid, idxf, temp1)

       start(2) = 1

       !                 Open netCDF file:

!!$       idq=ncvid(cdfid,TRIM(qname), ierr)
!!$       idu=ncvid(cdfid,TRIM(uname), ierr)
!!$       idt=ncvid(cdfid,TRIM(tname) ,ierr)

    else if(status.eq.2) then
       !         Close netCDF
       call ncclos(cdfid, ierr)
       return
    end if

    !               Output fields

    !              concentration field
    count(1)=nx
    temp = q(:)
    call ncvpt(cdfid, idq , start, count, temp, ierr)

    !              u field
    count(1)=nx+1
    temp1 = u(:)
    call ncvpt(cdfid, idu , start, count, temp1, ierr)

    count(1)=nx+1
    temp1 = lam(:)
    call ncvpt(cdfid, idl , start, count, temp1, ierr)

    count(1)=nx+1
    temp1 = mnl(:)
    call ncvpt(cdfid, idm , start, count, temp1, ierr)

!!$    if(ilevel.eq.1) then
!!$       write(*,*) time, start(2)
!!$       time4=time
!!$       call ncvpt1(cdfid, idt, start(2), time4, ierr)
!!$    end if
    start(2) = start(2) + 1

  end subroutine output1d

end program execute
