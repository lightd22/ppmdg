subroutine skamstep_2d(q,dqdt,u0,v0,u2,v2,rho_in,rhoq,rhoprime,nx,ny,npad,dt,jcbn,&
            xlambda,xmonlimit,ylambda,ymonlimit)
    !  use forward-in-time Skamarock (2006) scheme

    implicit none

    integer, intent(in) :: nx, ny, npad
    real (kind=8), intent(in)  :: dt

    real (kind=8), dimension(1-npad:nx+npad,1-npad:ny+npad), intent(inout) :: q
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: dqdt ! extra tendency
    real (kind=8), dimension(0:nx,1:ny), intent(in) :: u0, u2
    real (kind=8), dimension(1:nx,0:ny), intent(in) :: v0, v2
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: jcbn
    real (kind=8), dimension(0:nx,1:ny), intent(out) :: xlambda,xmonlimit
    real (kind=8), dimension(1:nx,0:ny), intent(out) :: ylambda,ymonlimit

    ! local variables, two-dimensional arrays
    real (kind=8), dimension(1:nx,1:ny), intent(in) :: rho_in
    real (kind=8), dimension(1:nx,1:ny), intent(inout) :: rhoq,rhoprime
    real (kind=8), dimension(0:nx,1:ny) :: u, uh
    real (kind=8), dimension(1:nx,0:ny) :: v, vh

!bloss    integer, parameter :: nmaxcfl = 20

    
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rho
    real (kind=8), dimension(0-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhouh
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl, &
                             0-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhovh

    ! local variables, one-dimensional arrays in x-direction
    real (kind=8), dimension(0:nx+1) :: rhoq1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: q1dx
    real (kind=8), dimension(-2:nx+2) :: rhouh1d
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rho1dx
    real (kind=8), dimension(1-npad-nmaxcfl:nx+npad+nmaxcfl) :: rhoprime1dx

    ! local variables, one-dimensional arrays in y-direction
    real (kind=8), dimension(0:ny+1) :: rhoq1dy
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: q1dy
    real (kind=8), dimension(-2:ny+2) :: rhovh1d
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rho1dy
    real (kind=8), dimension(1-npad-nmaxcfl:ny+npad+nmaxcfl) :: rhoprime1dy

    real(kind=8) :: scale
    integer, dimension(2) :: xbctype, ybctype
    real (kind=8), dimension(2) :: fxbc, fybc ! specified bndry flux 
    real (kind=8), dimension(0:nx  ,ny) :: xflx
    real (kind=8), dimension(nx  ,0:ny) :: yflx

    integer i,j, ii, npadq, npadu, npadrho, npadrp, npad2

    real(kind=8), external :: tfcn
    real(kind=8) :: t_temp, tmpflx(0:ny), tmplam(0:ny), tmpmon(0:ny)

    ! set up boundary types, fluxes, limiting flags
    if(nofluxew) then
       xbctype = 1 ! try open bc
    else
       xbctype = 2           ! periodic boundary conditions
    end if
    if(nofluxns) then
       ybctype = 1 ! try open bc
    else
       ybctype = 2           ! periodic boundary conditions
    end if

    fxbc = 0.d0          ! initialize in case fixed flux is used.
    fybc = 0.d0          ! initialize in case fixed flux is used.

!!$    if(nofluxew) xbctype = 3
!!$    if(nofluxns) xbctype = 3


    scale = 1.

    t_temp = time+0.5d0*dt
    uh = u0
    vh = v0
    if(transient) then
       uh = uh*tfcn(t_temp)
       vh = vh*tfcn(t_temp)
    end if

    if(transient2) then
       uh = u0*tfcn(t_temp)**2 + u2*(1.d0-tfcn(t_temp)**2)
       vh = v0*tfcn(t_temp)**2 + v2*(1.d0-tfcn(t_temp)**2)
    end if

    if(transient3) then
       t_temp = 2.d0*t_temp
       uh = u0 + u2*tfcn(t_temp)
       vh = v0 + v2*tfcn(t_temp)
    end if

    ! set padding size for each variable
    npadrho = npad-1 + (nmaxcfl+1)/2
    npadu = 2 + (nmaxcfl+1)/2

    npadq = npad+nmaxcfl-1
    npadrp = nmaxcfl+2

    ! set up density array for back trajectory computation
    rho(:,:) = 0.d0 ! zero out
    rho(1:nx,1:ny) = rho_in(:,:)

    ! x-mass flux at time t + dt/2
    rhouh(:,:) = 0.d0 ! zero out
    rhouh(0:nx,1:ny) = uh(:,:)

    ! y-mass flux at time t + dt/2
    rhovh(:,:) = 0.d0 ! zero out
    rhovh(1:nx,0:ny) = vh(:,:)
    
    if(nofluxew) then ! no scalar flux at x boundaries
       do i = 1,npadrho
          rho(1-i,1:ny) = rho(1,1:ny)
          rho(nx+i,1:ny) = rho(nx,1:ny)
       end do
       do i = 1,npadu
          rhouh(-i,1:ny) = rhouh(0,1:ny)
          rhouh(nx+i,1:ny) = rhouh(nx,1:ny)
       end do
    else  ! periodic in x
       rho(1-npadrho:0,1:ny) = rho(nx+1-npadrho:nx,1:ny)
       rho(nx+1:nx+npadrho,1:ny) = rho(1:npadrho,1:ny)
       rhouh(-npadu:-1,1:ny) = rhouh(nx+1-npadu:nx,1:ny)
       rhouh(nx+1:nx+npadu,1:ny) = rhouh(1:npadu,1:ny)
    end if
    
    if(nofluxns) then ! no scalar flux at y boundaries
       do i = 1,npadrho
          rho(:,1-i) = rho(:,1)
          rho(:,ny+i) = rho(:,ny)
       end do
       do i = 1,npadu
          rhovh(:,-i) = rhovh(:,0)
          rhovh(:,ny+i) = rhovh(:,ny)
       end do
    else ! periodic in y
       rho(:,1-npadrho:0) = rho(:,ny+1-npadrho:ny)
       rho(:,ny+1:ny+npadrho) = rho(:,1:npadrho)
       rhovh(:,-npadu:-1) = rhovh(:,ny+1-npadu:ny)
       rhovh(:,ny+1:ny+npadu) = rhovh(:,1:npadu)
    end if

    q1dx = 0.d0
    rhoprime1dx = 0.d0
    rho1dx = 0.d0
    rhouh1d = 0.d0

    q1dy = 0.d0
    rhoprime1dy = 0.d0
    rho1dy = 0.d0
    rhovh1d = 0.d0

    if(oddstep) then

       ! perform sweeps in x-direction
       do j = 1,ny
          q1dx(1:nx)     = q(1:nx,j)
          rhoq1dx(1:nx)  = rhoq(1:nx,j)

          rhoprime1dx(1:nx) = rhoprime(1:nx,j)

          if(nofluxew) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dx(1-ii) = MAX(0.d0,q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2)))
                   q1dx(nx+ii) = MAX(0.d0,q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1)))
                end do
             else
                do ii = 1,npadq
                   q1dx(1-ii) = q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2))
                   q1dx(nx+ii) = q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dx(1-npadrp:0) = rhoprime1dx(1)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(nx)
             
             ! fill in rhoq values 
             rhoq1dx(0) = rhoprime1dx(0)*q1dx(0)
             rhoq1dx(nx+1) = rhoprime1dx(nx+1)*q1dx(nx+1)

          else ! periodic in x
             rhoq1dx(0) = rhoq1dx(nx)
             rhoq1dx(nx+1) = rhoq1dx(1)

             q1dx(1-npadq:0) = q1dx(nx+1-npadq:nx)
             q1dx(nx+1:nx+npadq) = q1dx(1:npadq)

             rhoprime1dx(1-npadrp:0) = rhoprime1dx(nx+1-npadrp:nx)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dx(1-npadrho:nx+npadrho) = rho(1-npadrho:nx+npadrho,j)
          rhouh1d(-2:nx+2) = rhouh(-2:nx+2,j)

          call ppmwrap(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx,xflx(0,j),dt, &
                       nx,npad,nmaxcfl,xbctype,fxbc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda, &
                       xlambda(0,j),xmonlimit(0,j))

          ! update solution
          q(1:nx,j) = rhoq1dx(1:nx)/rhoprime1dx(1:nx)
          rhoprime(1:nx,j) = rhoprime1dx(1:nx)
          rhoq(1:nx,j) = rhoq1dx(1:nx)

       end do

       ! perform sweeps in y-direction
       do i = 1,nx
          q1dy(1:ny)     = q(i,1:ny)
          rhoq1dy(1:ny)  = rhoq(i,1:ny)

          rhoprime1dy(1:ny) = rhoprime(i,1:ny)

          if(nofluxns) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dy(1-ii) = MAX(0.d0,q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2)))
                   q1dy(ny+ii) = MAX(0.d0,q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1)))
                end do
             else
                do ii = 1,npadq
                   q1dy(1-ii) = q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2))
                   q1dy(ny+ii) = q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dy(1-npadrp:0) = rhoprime1dy(1)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(ny)
             
             ! fill in rhoq values 
             rhoq1dy(0) = rhoprime1dy(0)*q1dy(0)
             rhoq1dy(ny+1) = rhoprime1dy(ny+1)*q1dy(ny+1)

          else ! periodic in y
             rhoq1dy(0) = rhoq1dy(ny)
             rhoq1dy(ny+1) = rhoq1dy(1)

             q1dy(1-npadq:0) = q1dy(ny+1-npadq:ny)
             q1dy(ny+1:ny+npadq) = q1dy(1:npadq)

             rhoprime1dy(1-npadrp:0) = rhoprime1dy(ny+1-npadrp:ny)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dy(1-npadrho:ny+npadrho) = rho(i,1-npadrho:ny+npadrho)
          rhovh1d(-2:ny+2) = rhovh(i,-2:ny+2)

          call ppmwrap(rhoq1dy,q1dy,rhovh1d,rho1dy,rhoprime1dy,tmpflx,dt, &
                       ny,npad,nmaxcfl,ybctype,fybc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda,tmplam, tmpmon)

          yflx(i,:) = tmpflx(:)
          ylambda(i,:) = tmplam(:)
          ymonlimit(i,:) = tmpmon(:)

          ! update solution
          q(i,1:ny) = rhoq1dy(1:ny)/rhoprime1dy(1:ny)
          rhoprime(i,1:ny) = rhoprime1dy(1:ny)
          rhoq(i,1:ny) = rhoq1dy(1:ny)

       end do

    else

       ! perform sweeps in y-direction
       do i = 1,nx
          q1dy(1:ny)     = q(i,1:ny)
          rhoq1dy(1:ny)  = rhoq(i,1:ny)

          rhoprime1dy(1:ny) = rhoprime(i,1:ny)

          if(nofluxns) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dy(1-ii) = MAX(0.d0,q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2)))
                   q1dy(ny+ii) = MAX(0.d0,q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1)))
                end do
             else
                do ii = 1,npadq
                   q1dy(1-ii) = q1dy(1) + DBLE(ii)*(q1dy(1)-q1dy(2))
                   q1dy(ny+ii) = q1dy(ny) + DBLE(ii)*(q1dy(ny)-q1dy(ny-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dy(1-npadrp:0) = rhoprime1dy(1)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(ny)
             
             ! fill in rhoq values 
             rhoq1dy(0) = rhoprime1dy(0)*q1dy(0)
             rhoq1dy(ny+1) = rhoprime1dy(ny+1)*q1dy(ny+1)

          else ! periodic in y
             rhoq1dy(0) = rhoq1dy(ny)
             rhoq1dy(ny+1) = rhoq1dy(1)

             q1dy(1-npadq:0) = q1dy(ny+1-npadq:ny)
             q1dy(ny+1:ny+npadq) = q1dy(1:npadq)

             rhoprime1dy(1-npadrp:0) = rhoprime1dy(ny+1-npadrp:ny)
             rhoprime1dy(ny+1:ny+npadrp) = rhoprime1dy(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dy(1-npadrho:ny+npadrho) = rho(i,1-npadrho:ny+npadrho)
          rhovh1d(-2:ny+2) = rhovh(i,-2:ny+2)

          call ppmwrap(rhoq1dy,q1dy,rhovh1d,rho1dy,rhoprime1dy,tmpflx,dt, &
                       ny,npad,nmaxcfl,ybctype,fybc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda,tmplam, tmpmon)

          yflx(i,:) = tmpflx(:)
          ylambda(i,:) = tmplam(:)
          ymonlimit(i,:) = tmpmon(:)

          ! update solution
          q(i,1:ny) = rhoq1dy(1:ny)/rhoprime1dy(1:ny)
          rhoprime(i,1:ny) = rhoprime1dy(1:ny)
          rhoq(i,1:ny) = rhoq1dy(1:ny)

       end do

       ! perform sweeps in x-direction
       do j = 1,ny
          q1dx(1:nx)     = q(1:nx,j)
          rhoq1dx(1:nx)  = rhoq(1:nx,j)

          rhoprime1dx(1:nx) = rhoprime(1:nx,j)

          if(nofluxew) then
             ! linear extrapolation
             if(doposlimit) then
                do ii = 1,npadq
                   q1dx(1-ii) = MAX(0.d0,q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2)))
                   q1dx(nx+ii) = MAX(0.d0,q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1)))
                end do
             else
                do ii = 1,npadq
                   q1dx(1-ii) = q1dx(1) + DBLE(ii)*(q1dx(1)-q1dx(2))
                   q1dx(nx+ii) = q1dx(nx) + DBLE(ii)*(q1dx(nx)-q1dx(nx-1))
                end do
             end if

             ! constant extrapolation of rhoprime
             rhoprime1dx(1-npadrp:0) = rhoprime1dx(1)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(nx)
             
             ! fill in rhoq values 
             rhoq1dx(0) = rhoprime1dx(0)*q1dx(0)
             rhoq1dx(nx+1) = rhoprime1dx(nx+1)*q1dx(nx+1)

          else ! periodic in x
             rhoq1dx(0) = rhoq1dx(nx)
             rhoq1dx(nx+1) = rhoq1dx(1)

             q1dx(1-npadq:0) = q1dx(nx+1-npadq:nx)
             q1dx(nx+1:nx+npadq) = q1dx(1:npadq)

             rhoprime1dx(1-npadrp:0) = rhoprime1dx(nx+1-npadrp:nx)
             rhoprime1dx(nx+1:nx+npadrp) = rhoprime1dx(1:npadrp)
          end if

          ! get rho and rhou from 2d arrays padded above
          rho1dx(1-npadrho:nx+npadrho) = rho(1-npadrho:nx+npadrho,j)
          rhouh1d(-2:nx+2) = rhouh(-2:nx+2,j)

          call ppmwrap(rhoq1dx,q1dx,rhouh1d,rho1dx,rhoprime1dx,xflx(0,j),dt, &
                       nx,npad,nmaxcfl,xbctype,fxbc, &
                       dosemilagr,dosellimit,domonlimit,doposlimit, &
                       dopcm, dowenosplit, &
                       scale,nmethod,lambdamax,epslambda, &
                       xlambda(0,j),xmonlimit(0,j))

          ! update solution
          q(1:nx,j) = rhoq1dx(1:nx)/rhoprime1dx(1:nx)
          rhoprime(1:nx,j) = rhoprime1dx(1:nx)
          rhoq(1:nx,j) = rhoq1dx(1:nx)

       end do

    end if

  end subroutine skamstep_2d
