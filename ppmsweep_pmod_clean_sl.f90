subroutine ppmsweep_pmod_select_clean(q,rhou,rhoprime,fx,dt,  &
     N,npad,bctype,fbc,doposlimit,scale)
  !     
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(0:N)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale

  !     # in/outputs
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(0:N+1)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 

  !     # local variables
  integer i
  real(kind=8) :: tmpeps, tmpx, a0, a1, a2
  real(kind=8) :: qface(-1:N+1), gamma(-1:N+2), qavg(0:N)
  real(kind=8) :: dql, dqr, qup, cr
  integer :: iup, il, ir

  !     # local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0
  real(kind=8), parameter :: maxlambda = 20., epshybrid = 1.e-8
  !     
  !     FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bctype(1).eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bctype(2).eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)
  !     
  !     # compute ppm estimate for flux
  !===============================
  !     
  !     # re-scale epshybrid by the square of scale.  This will prevent 
  !     #   epshybrid from overwhelming the smoothness indicators when 
  !     #   q is small.
  tmpeps = epshybrid*scale*scale

  do i = -1,N+2
    gamma(i) = (q(i-1)-q(i))**2 + (q(i)-q(i+1))**2
  end do

  if(doposlimit) then
    do i = -1,N+1
      qface(i) = MAX(0.d0, &
	           fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
    end do
  else
    do i = -1,N+1
      qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
    end do
  end if

  do i = 0,N

    if(rhou(i).ge.0.d0) then
      !     # courant number
      iup = i
      il = i+1
      ir = i-1
      dql = qface(i) - q(i)
      dqr = qface(i-1) - q(i)
    else
      iup = i+1
      il = i
      ir = i+2
      dql = qface(i) - q(i+1)
      dqr = qface(i+1) - q(i+1)
    end if

    cr = ABS(rhou(i))*dt/rhoprime(iup)
    qup = q(iup)
    a0 = qface(i)
    a1 = -4.d0*dql - 2.d0*dqr
    a2 =  3.d0*dql + 3.d0*dqr

    !     check whether lambda > maxlambda
    if(MAX(gamma(il),gamma(iup),gamma(ir)) &
	        .gt.maxlambda*(tmpeps &
	        + MIN(gamma(il),gamma(iup),gamma(ir)) &
	        )) then
      ! selective monotonicity preservation has been triggered
      ! because lambda > maxlambda

      ! make sure the edge values are monotonic
      dql = MAX(MIN(q(il)-q(iup),0.d0), &
	           MIN(MAX(q(il)-q(iup),0.d0),dql))
      dqr = MAX(MIN(q(ir)-q(iup),0.d0), &
	           MIN(MAX(q(ir)-q(iup),0.d0),dqr))

      ! check whether (qbar-ql)*(qbar-qr)>0
      if(dql*dqr.gt.0.d0) then
        ! if so, use piecewise constant reconstruction
        a0 = qup
        a1 = 0.d0
        a2 = 0.d0
      else
        ! next, reconstruct the coefficients (in case dql, dqr changed)
        a0 = qup + dql
        a1 = -4.d0*dql - 2.d0*dqr
        a2 =  3.d0*dql + 3.d0*dqr
        ! check whether the polynomial has an extrema within the interval
        if(abs(a1+a2).lt.abs(a2)) then            
          ! remove the extrema by modifying either ql or qr
          if(ABS(dql).le.ABS(dqr)) then
            !     extrema is on the ql side of the interval from ql to qr, 
            !     modify qr to remove extrema.  Note that a0 doesn't change.
            a1 = 0.d0
            a2 = -3.d0*dql
          else
            !     extrema is on the qr side of the interval from ql to qr, 
            !     modify ql to remove extrema
            a0 = -2.d0*dqr + qup 
            a1 =  6.d0*dqr
            a2 = -3.d0*dqr
          end if
        end if
      end if

      !     limit parabola to prevent negative values
    elseif(doposlimit.AND.(abs(a1+a2).lt.abs(a2))) then            
      !     # extrema is within interval, locate it and
      ! test for sign of parabola at extrema
      tmpx = -a1/(2.d0*a2)
      if(a0 + tmpx*(a1 + tmpx*a2).lt.0.d0) then
        ! if extrema is negative, use piecewise-constant reconstruction
        a0 = qup
        a1 = 0.d0
        a2 = 0.d0
      end if
    end if

    !     compute outgoing fluxes at right and left edges of cell
    qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
    fx(i) = rhou(i)*qavg(i)
  end do

  !===============================================================
  !========DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !==============LEFT BOUNDARY ==============
  if((bctype(1).eq.1).and.(rhou(0).gt.0.d0)) then
    !     inflow at left boundary, impose no scalar gradient across cell next to 
    !     boundary.  Use fbc(1) to account for possible mean gradient of scalar
    fx(0) = rhou(0)*(qavg(1) - fbc(1))
  elseif(bctype(1).eq.3) then
    !     specify fixed flux at left boundary by modifying flux correction.
    fx(0) = fbc(1)
  end if

  !==============RIGHT BOUNDARY ==============
  if((bctype(2).eq.1).and.(rhou(N).lt.0.d0)) then
    !     inflow at right boundary, impose no scalar gradient across cell next to 
    !     boundary.  Use fbc(2) to account for possible mean gradient of scalar
    fx(N) = rhou(N)*(qavg(N-1) + fbc(2))
  elseif(bctype(2).eq.3) then
    !     specify fixed flux at right boundary
    fx(N) = fbc(2)
  end if

  !     # update solution (positive-definite) along with rhoprime
  !==============================================================
  !     update solution (q now holds rho*q)
  if(doposlimit) then
    do i = 1,N
      q(i) = MAX(0.d0,rhoprime(i)*q(i) + dt*(fx(i-1) - fx(i)))
    end do
  else
    do i = 1,N
      q(i) = rhoprime(i)*q(i) + dt*(fx(i-1) - fx(i))
    end do
  end if

  !     update rhoprime
  do i = 1,N
    rhoprime(i) = rhoprime(i) + dt*(rhou(i-1)-rhou(i))
  end do

  !     compute new value for q (dividing through by rhoprime)
  do i = 1,N
    q(i) = q(i)/rhoprime(i)
  end do

end subroutine ppmsweep_pmod_select_clean

subroutine ppmsweep_pmod_select_clean_sl(q,rhou,rhoprime,fx,dt,  &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale)
  ! flux form semi-Lagrangian version
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: nmaxcfl
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(0:N)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale

  !     # in/outputs
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 

  !     # local variables
  integer i
  real(kind=8) :: tmpeps, tmpx, a0, a1, a2
  real(kind=8) :: qface(2-npad:N+npad-2)
  real(kind=8) :: gamma(2-npad:N+npad-1)
  real(kind=8) :: qavg(0:N)
  real(kind=8) :: dql, dqr, qup, cr, rhoufrac
  integer :: iup, il, ir, k

  !     # local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0
  real(kind=8), parameter :: maxlambda = 20., epshybrid = 1.e-8
  !     
  !     FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bctype(1).eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bctype(2).eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)
  !     
  !     # compute ppm estimate for flux
  !===============================
  !     
  !     # re-scale epshybrid by the square of scale.  This will prevent 
  !     #   epshybrid from overwhelming the smoothness indicators when 
  !     #   q is small.
  tmpeps = epshybrid*scale*scale

  do i = 2-npad,N+npad-1
    gamma(i) = (q(i-1)-q(i))**2 + (q(i)-q(i+1))**2
  end do

  if(doposlimit) then
    do i = 2-npad,N+npad-2
      qface(i) = MAX(0.d0, &
	           fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
    end do
  else
    do i = 2-npad,N+npad-2
      qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
    end do
  end if

  do i = 0,N

    fx(i) = 0.d0
    rhoufrac = rhou(i)*dt ! mass fluxed through face i over the time step

    if(rhou(i).ge.0.d0) then
      ! use flux-form semi-lagrangian method.
      ! accumulate flux first by accumulating rhoprime*q in cells
      ! traversed by the integer part of the cfl 
      iup = i
      do k = 1,nmaxcfl+1
        if(rhoufrac.lt.rhoprime(iup)) EXIT
        fx(i) = fx(i) + rhoprime(iup)*q(iup) ! integer cfl part
        rhoufrac = rhoufrac - rhoprime(iup) ! integer cfl mass flux
        iup = iup - 1 ! change index of cell being added to flux
      end do
      !     # courant number
      iup = iup
      il = iup+1
      ir = iup-1
      dql = qface(iup) - q(iup)
      dqr = qface(ir) - q(iup)
    else
      iup = i+1
      do k = 1,nmaxcfl+1
        if(-rhoufrac.lt.rhoprime(iup)) EXIT
        fx(i) = fx(i) - rhoprime(iup)*q(iup) ! integer cfl part
        rhoufrac = rhoufrac + rhoprime(iup) ! integer cfl mass flux
        iup = iup + 1 ! change index of cell being added to flux
      end do
      !     
      !# compute upwind and PCM correction flux
      !     
      il = iup-1
      ir = iup+1
      dql = qface(il) - q(iup)
      dqr = qface(iup) - q(iup)
    end if

    cr = ABS(rhoufrac)/rhoprime(iup) ! note: factor of dt already in rhoufrac
    qup = q(iup)
    a0 = q(iup) + dql
    a1 = -4.d0*dql - 2.d0*dqr
    a2 =  3.d0*dql + 3.d0*dqr

    !     check whether lambda > maxlambda
    if(MAX(gamma(il),gamma(iup),gamma(ir)) &
	        .gt.maxlambda*(tmpeps &
	        + MIN(gamma(il),gamma(iup),gamma(ir)) &
	        )) then
      ! selective monotonicity preservation has been triggered
      ! because lambda > maxlambda

      ! make sure the edge values are monotonic
      dql = MAX(MIN(q(il)-q(iup),0.d0), &
	           MIN(MAX(q(il)-q(iup),0.d0),dql))
      dqr = MAX(MIN(q(ir)-q(iup),0.d0), &
	           MIN(MAX(q(ir)-q(iup),0.d0),dqr))

      ! check whether (qbar-ql)*(qbar-qr)>0
      if(dql*dqr.gt.0.d0) then
        ! if so, use piecewise constant reconstruction
        a0 = qup
        a1 = 0.d0
        a2 = 0.d0
      else
        ! next, reconstruct the coefficients (in case dql, dqr changed)
        a0 = qup + dql
        a1 = -4.d0*dql - 2.d0*dqr
        a2 =  3.d0*dql + 3.d0*dqr
        ! check whether the polynomial has an extrema within the interval
        if(abs(a1+a2).lt.abs(a2)) then            
          ! remove the extrema by modifying either ql or qr
          if(ABS(dql).le.ABS(dqr)) then
            !     extrema is on the ql side of the interval from ql to qr, 
            !     modify qr to remove extrema.  Note that a0 doesn't change.
            a1 = 0.d0
            a2 = -3.d0*dql
          else
            !     extrema is on the qr side of the interval from ql to qr, 
            !     modify ql to remove extrema
            a0 = -2.d0*dqr + qup 
            a1 =  6.d0*dqr
            a2 = -3.d0*dqr
          end if
        end if
      end if

      !     limit parabola to prevent negative values
    elseif(doposlimit.AND.(abs(a1+a2).lt.abs(a2))) then            
      !     # extrema is within interval, locate it and
      ! test for sign of parabola at extrema
      tmpx = -a1/(2.d0*a2)
      if(a0 + tmpx*(a1 + tmpx*a2).lt.0.d0) then
        ! if extrema is negative, use piecewise-constant reconstruction
        a0 = qup
        a1 = 0.d0
        a2 = 0.d0
      end if
    end if

    !     compute outgoing fluxes at right and left edges of cell
    qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
    fx(i) = fx(i) + rhoufrac*qavg(i) ! includes factor of dt.
  end do

  !===============================================================
  !========DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !==============LEFT BOUNDARY ==============
  if((bctype(1).eq.1).and.(rhou(0).gt.0.d0)) then
    !     inflow at left boundary, impose no scalar gradient across cell next to 
    !     boundary.  Use fbc(1) to account for possible mean gradient of scalar
    if(rhou(1).eq.0.) then
      fx(0) = dt*rhou(0)*(qavg(1) - fbc(1))
    else
      fx(0) = rhou(0)*fx(1)/rhou(1)
    end if
  elseif(bctype(1).eq.3) then
    !     specify fixed flux at left boundary by modifying flux correction.
    fx(0) = dt*fbc(1)
  end if

  !==============RIGHT BOUNDARY ==============
  if((bctype(2).eq.1).and.(rhou(N).lt.0.d0)) then
    !     inflow at right boundary, impose no scalar gradient across cell next to 
    !     boundary.  Use fbc(2) to account for possible mean gradient of scalar
    if(rhou(N-1).eq.0.) then
      fx(N) = dt*rhou(N)*(qavg(N-1) + fbc(2))
    else
      fx(N) = rhou(N)*fx(N-1)/rhou(N-1)
    end if
  elseif(bctype(2).eq.3) then
    !     specify fixed flux at right boundary
    fx(N) = dt*fbc(2)
  end if

  !     # update solution (positive-definite) along with rhoprime
  !==============================================================
  !     update solution (q now holds rho*q)
  if(doposlimit) then
    do i = 1,N
      q(i) = MAX(0.d0,rhoprime(i)*q(i) + fx(i-1) - fx(i))
    end do
  else
    do i = 1,N
      q(i) = rhoprime(i)*q(i) + fx(i-1) - fx(i)
    end do
  end if

  !     update rhoprime
  do i = 1,N
    rhoprime(i) = rhoprime(i) + dt*(rhou(i-1)-rhou(i))
  end do

  !     compute new value for q (dividing through by rhoprime)
  do i = 1,N
    q(i) = q(i)/rhoprime(i)
  end do

  fx(0:N) = fx(0:N)/dt

end subroutine ppmsweep_pmod_select_clean_sl

