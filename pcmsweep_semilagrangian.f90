! ============================================================
subroutine pcmsweep_fct_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !# inputs
  logical, intent(in) :: doposlimit
  integer, intent(in) :: N, npad, nmaxcfl, bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: dt, scale
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !# in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhoprime(1-npad-nmaxcfl:N+npad+nmaxcfl)

  !# outputs
  real(kind=8) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, imin, imax, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) ::  crface, crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl), qtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-2:N+2), rhoprime_new(-1:N+2)
  real(kind=8) :: beta1, beta2, beta3, tmpeps, lambda(-2:N+2), rhoqtd(-1:N+2)
  real(kind=8) :: secnd(-2-nmaxcfl:N+3+nmaxcfl), tmp_fxcorr(-2:N+2)
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: qup(-2:N+2), qcorr(-2:N+2)
  integer :: monlimit(-3:N+3), itmp, ib, ic
  real(kind=8) :: sdir
  integer :: idir, iup, idn, ii, iu, iu2, ishift(-2:N+2)
  real(kind=8) :: cr, dq, ql, qr, qbar, a0, a1, a2, a3
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16, &
       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  monlimit(:) = 0
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad-nmaxcfl:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad+nmaxcfl) = q(N)
  !     
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  !# compute ppm estimate for flux
  !===============================
  !     
  do i = -1-nmaxcfl,-2+nmaxcfl
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do
  
  do i = -2,N+2
     !# qface holds interpolant of q at faces.
     qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+1+nmaxcfl)) &
          + fac2*(q(i-1+nmaxcfl) + q(i+2+nmaxcfl))

     if(rhou(i).ge.0.d0) then
        fxup(i) = 0.d0
        rhoufrac = rhou(i)*dt
        iup = i
        ii = iup
        do k = 1,nmaxcfl
           if(rhoufrac.lt.rhoprime(ii)) EXIT
           fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
           rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
           ii = ii - 1 ! change index of cell being added to flux
        end do
        !
        !# compute upwind and PCM correction flux
        !     
        ishift(i) = iup - ii                ! size of integer shift at face i
        fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
        cr = rhoufrac/rhoprime(ii)
        dq = fac3*(q(ii-1)-q(ii+1)) + fac4*(q(ii-2)-q(ii+2)) ! reverse sign
        ql = qface(ii) ! note: ql and qr flipped here.
        qr = qface(ii-1)
        qbar = q(ii)
     else
        fxup(i) = 0.d0
        rhoufrac = rhou(i)*dt
        iup = i+1
        ii = iup
        do k = 1,nmaxcfl
           if(-rhoufrac.lt.rhoprime(ii)) EXIT
           fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
           rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
           ii = ii + 1 ! change index of cell being added to flux
        end do
        !     
        !# compute upwind and PCM correction flux
        !     
        ishift(i) = iup - ii                ! size of integer shift at face i
        fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
        cr = ABS(rhoufrac/rhoprime(ii))
        dq = fac3*(q(ii+1)-q(ii-1)) + fac4*(q(ii+2)-q(ii-2))
        ql = qface(ii-1)
        qr = qface(ii)
        qbar = q(ii)
     end if
     !     
     !# compute ratio of largest to smallest smoothness indicator
     beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)-q(ii-1))**2
     beta2 = (q(ii)-q(ii-1))**2 + (q(ii)-q(ii+1))**2
     beta3 = (q(ii+1)-q(ii))**2 + (q(ii+2)-q(ii+1))**2
     !     
     !# compute ratio of largest to smallest smoothness indicator.
     !     
     lambda(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
     !
     !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
     ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
     a0 =       ql                                  ! constant term
     a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
     a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
     a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

     !# compute PCM flux by integrating across reconstruction.
     qup(i) = qbar
     qcorr(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0)) - qbar
     fxcorr(i) = rhoufrac*qcorr(i)
     tmp_fxcorr(i) = fxcorr(i)
  end do
  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  imin = 0
  imax = N

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).or.(bcleft.eq.3)) then
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(1) to account for possible mean gradient of scalar
        fxup(0) = dt*rhou(0)*(qup(1) - fbc(1))
        fxcorr(0) = dt*rhou(0)*qcorr(1)

        lambda(0) = lambda(1)
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fxcorr(0) = dt*fxleft - fxup(0) 
        imin = 1 ! no limiting of left boundary flux
     end if

     ! if open or fixed flux bcs at the left end, specify q and fxup
     !   so that cells -1 and 0 match cell 1 (essentially constraining 
     !   any flux correction to be applied only using values from the
     !   domain interior
     fxup(-2:-1) = fxup(0)
     fxcorr(-2:-1) = fxcorr(0)
     tmp_fxcorr(-2:0) = fxcorr(-2:0)
     q(-1:0) = q(1)
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).or.(bcright.eq.3)) then
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(2) to account for possible mean gradient of scalar
        fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
        fxcorr(N) = dt*rhou(N)*qcorr(N-1)

        lambda(N) = lambda(N-1)
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fxcorr(N) = dt*fxright - fxup(N)
        imax = N-1 ! no limiting of right boundary flux
     end if

     ! if open or fixed flux bcs at the left end, specify q and fxup
     !   so that cells -1 and 0 match cell 1 (essentially constraining 
     !   any flux correction to be applied only using values from the
     !   domain interior)
     fxup(N+1:N+2) = fxup(N)
     fxcorr(N+1:N+2) = fxcorr(N)
     tmp_fxcorr(N:N+2) = fxcorr(N:N+2)
     q(N+1:N+2) = q(N)
  end if


  !
  !# update rhoprime -- kept to ensure mass consistency
  rhoprime_new(-1:N+2) = rhoprime(-1:N+2) + dt*(rhou(-2:N+1) - rhou(-1:N+2))
  rhoqtd(-1:N+2) = rhoprime(-1:N+2)*q(-1:N+2) + fxup(-2:N+1) - fxup(-1:N+2)
  qtd(-1:N+2) = rhoqtd(-1:N+2)/rhoprime_new(-1:N+2)

  tmpeps = eps2*scale ! a small parameter, used to prevent division by zero

  !============= selective, positive flux correction =============
  do i = imin,imax

     ! perform monotonic flux correction if lambda>lambda_max at this face.
     if(lambda(i).gt.maxlambda) then
        monlimit(i) = 1

        if(tmp_fxcorr(i).ge.0.d0) then

           ii = i - ishift(i)

           !check whether flux out of cell i will induce new minimum in q
           qmn = min(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)

           !check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
           fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)

        else

           ii = i - ishift(i)

           !check whether flux into cell i will induce new maximum in q 
           qmx = max(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
           fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
           ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)

           !check whether flux out of cell i+1 will induce new min in q
           qmn = min(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
           fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
           ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)

        end if

        !# limit flux correction 
        fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)

     end if

  end do

  if(doposlimit) then

     do i = imin,imax
        ! perform flux correction for non-negativity in all locations
        if(tmp_fxcorr(i).ge.0.d0) then
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoqtd(i)
        else
           fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
           ratiom = rhoqtd(i+1)
        end if
        if(fluxout.gt.ratiom) then
           fxcorr(i) = tmp_fxcorr(i)*max(0.d0,min(1.d0, &
                ratiom/(tmpeps+fluxout), &      ! ratiom for positivity
                fxcorr(i)/(tmp_fxcorr(i)+tmpeps))) ! ratiom for selective monotonicity
        end if
     end do

     !# update solution, eliminating any remaining negative values
     !============================================================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = rhoqtd(i) + fxcorr(i-1) - fxcorr(i) !max(0.d0,rhoq(i) + fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do

  else

     !# update solution
     !=================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do
  end if ! doposlimit

  fx(:) = fx(:)/dt
  limitpts(:) = monlimit(0:N)
  lambdavec(:) = lambda(0:N)

end subroutine pcmsweep_fct_select_semilagrangian
!
! ============================================================
! ============================================================
subroutine pcmsweep_pmod_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !# inputs
  logical, intent(in) :: doposlimit
  integer, intent(in) :: N, npad, nmaxcfl, bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: dt, scale
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !# in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhoprime(1-npad-nmaxcfl:N+npad+nmaxcfl)

  !# outputs
  real(kind=8) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, imin, imax, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: limit(-1:N+1), lambda(-1:N+1)
  real(kind=8) :: tmpeps, qmn, qmx, qavg(0:N)
  real(kind=8) :: ql, qr, qbar, a0, a1, a2, a3, cr, dq
  real(kind=8) :: discr, x1, x2, x1dis, x2dis
  real(kind=8) ::  crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: beta1, beta2, beta3
  real(kind=8) :: sdir, rhoqtd, fluxout, tmp_fxcorr
  real(kind=8) :: fxup(-1:N+1), fxcorr(-1:N+1)
  real(kind=8) :: qup(-1:N+1), qcorr(-1:N+1)
  integer :: idir, iup, idn, ii, iu, iu2, ishift(-2:N+2)
  integer :: ilft, ibar, irgt
  logical :: poscorr
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16, &
       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
!!$  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, &
!!$!!$       eps2=1.d-16, epshybrid=1.d-4, maxlambda=20.d0
!!$       eps2=1.d-16, epshybrid=1.d-8, maxlambda=100.d0
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     q(1-npad-nmaxcfl:0) = q(1)
     rhoq(0) = rhoq(1)
  end if
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     q(N+1:N+npad+nmaxcfl) = q(N)
     rhoq(N+1) = rhoq(N)
  end if
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  limit(:) = 0
  !     
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  if(doposlimit) then
     
     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -2-nmaxcfl,nmaxcfl
        qface(i) = MAX(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do

     do i = -1,N+1
        ! initialize qface and fdbeta to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+nmaxcfl) = MAX(0.d0,fac1*(q(i+nmaxcfl) + q(i+nmaxcfl+1)) &
                            + fac2*(q(i+nmaxcfl-1) + q(i+nmaxcfl+2)))

        if(rhou(i).ge.0.d0) then
           fxup(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i
           do k = 1,nmaxcfl
              if(rhoufrac.lt.rhoprime(ii)) EXIT
              fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
              ii = ii - 1 ! change index of cell being added to flux
           end do
           !
           !# compute upwind and PCM correction flux
           !     
           cr = rhoufrac/rhoprime(ii)
           dq = fac3*(q(ii-1)-q(ii+1)) + fac4*(q(ii-2)-q(ii+2)) ! reverse sign
           ql = qface(ii) ! note: ql and qr flipped here.
           qr = qface(ii-1)
           qbar = q(ii)
           ilft = ii+1
           ibar = ii
           irgt = ii-1
        else
           fxup(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i+1
           do k = 1,nmaxcfl
              if(-rhoufrac.lt.rhoprime(ii)) EXIT
              fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
              ii = ii + 1 ! change index of cell being added to flux
           end do
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhoufrac/rhoprime(ii))
           dq = fac3*(q(ii+1)-q(ii-1)) + fac4*(q(ii+2)-q(ii-2))
           ql = qface(ii-1)
           qr = qface(ii)
           qbar = q(ii)
           ilft = ii-1
           ibar = ii
           irgt = ii+1
        end if
        !     
        !# compute ratio of largest to smallest smoothness indicator
        beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)  -q(ii-1))**2
        beta2 = (q(ii)  -q(ii-1))**2 + (q(ii+1)-q(ii)  )**2
        beta3 = (q(ii+1)-q(ii)  )**2 + (q(ii+2)-q(ii+1))**2
        !     
        !# compute ratio of largest to smallest smoothness indicator.
        !     
        lambda(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
        !
        !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if(lambda(i).gt.maxlambda) then
           limit(i) = 1

           ! monotonize the cell edge values
           qmx = max(q(ilft),q(ibar))
           qmn = min(q(ilft),q(ibar))
           ql = MAX(qmn,MIN(qmx,ql))

           qmx = max(q(irgt),q(ibar))
           qmn = min(q(irgt),q(ibar))
           qr = MAX(qmn,MIN(qmx,qr))

           ! if qbar outside of [ql,qr], use piecewise constant,
           !   otherwise use parabola
           if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
              a0 = qbar ! piecewise constant reconstruction
              a1 = 0.d0
              a2 = 0.d0
              a3 = 0.d0 
           else
              ! define coefficients of parabola
              a0 = ql
              a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
              a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
              a3 = 0.d0 
              if(a2.ne.0.d0) then
                 ! check for extrema of the parabola within the interval
                 x1 = -a1/(2.d0*a2)
                 if(ABS(x1-0.5d0).lt.0.5d0) then
                    ! the parabola has an extrema in the interval, so that
                    !   we're going to modify the parabola to remove the extrema
                    if(ABS(qbar-ql).le.ABS(qbar-qr)) then
                       ! extrema is on the ql side of the interval from ql to qr, 
                       !   so that we bound parabola using ql to remove under/overshoot
                       a0 = ql
                       a1 = 0.d0
                       a2 = -3.d0*ql+3.d0*qbar
                       a3 = 0.d0
                    else
                       ! extrema is on the qr side of the interval from ql to qr, 
                       !   so that we bound parabola using qr to remove under/overshoot
                       a0 = 3.d0*qbar-2.d0*qr
                       a1 = 6.d0*qr-6.d0*qbar
                       a2 = -3.d0*qr+3.d0*qbar
                       a3 = 0.d0
                    end if
                 end if
              end if
           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qup(i) = qbar
        qcorr(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0)) - qbar

        fxup(i) = fxup(i) + rhoufrac*qbar
        fxcorr(i) = rhoufrac*qcorr(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

     imin = 0
     imax = N

     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(1) to account for possible mean gradient of scalar
        fxup(0) = dt*rhou(0)*MAX(0.d0,qup(1) - fbc(1))
        fxcorr(0) = dt*rhou(0)*qcorr(1)
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fxcorr(0) = dt*fxleft - fxup(0) 
        imin = 1 ! no limiting of left boundary flux
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(2) to account for possible mean gradient of scalar
        fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
        fxcorr(N) = dt*rhou(N)*qcorr(N-1)
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fxcorr(N) = dt*fxright - fxup(N)
        imax = N-1 ! no limiting of right boundary flux
     end if

     ! limit correction fluxes at faces to prevent the occurrence of 
     !   negative scalar values
     tmp_fxcorr = fxcorr(-1)
     do i = 0,N
        if(fxcorr(i).ge.0.d0) then
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i) + (fxup(i-1) - fxup(i))
           fluxout = max(0.d0,fxcorr(i)) - min(0.d0,tmp_fxcorr)
           tmp_fxcorr = fxcorr(i)
        else
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i+1) + (fxup(i) - fxup(i+1))
           fluxout = max(0.d0,fxcorr(i+1)) - min(0.d0,fxcorr(i))
           tmp_fxcorr = fxcorr(i)
        end if
        if(fluxout.gt.rhoqtd) then
           fxcorr(i) = fxcorr(i)*MAX(0.d0,MIN(1.d0,rhoqtd/(fluxout+EPSILON(rhoqtd))))
        end if
     end do

     !# update solution, eliminating any remaining negative values
     !============================================================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + (fx(i-1) - fx(i)) !MAX(0.d0,
        q(i) = rhoq(i)/rhoprime(i)
     end do

     lambdavec(0:N) = lambda(0:N)
     limitpts(0:N) = limit(0:N)

  else
     
     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -1-nmaxcfl,nmaxcfl
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = 0,N
        ! initialize qface and fdbeta to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+nmaxcfl+1)) &
                                          + fac2*(q(i+nmaxcfl-1) + q(i+nmaxcfl+2))

        if(rhou(i).ge.0.d0) then
           fx(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i
           do k = 1,nmaxcfl
              if(rhoufrac.lt.rhoprime(ii)) EXIT
              fx(i) = fx(i) + rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
              ii = ii - 1 ! change index of cell being added to flux
           end do
           !
           !# compute upwind and PCM correction flux
           !     
           cr = rhoufrac/rhoprime(ii)
           dq = fac3*(q(ii-1)-q(ii+1)) + fac4*(q(ii-2)-q(ii+2)) ! reverse sign
           ql = qface(ii) ! note: ql and qr flipped here.
           qr = qface(ii-1)
           qbar = q(ii)
           ilft = ii+1
           ibar = ii
           irgt = ii-1
        else
           fx(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i+1
           do k = 1,nmaxcfl
              if(-rhoufrac.lt.rhoprime(ii)) EXIT
              fx(i) = fx(i) - rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
              ii = ii + 1 ! change index of cell being added to flux
           end do
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhoufrac/rhoprime(ii))
           dq = fac3*(q(ii+1)-q(ii-1)) + fac4*(q(ii+2)-q(ii-2))
           ql = qface(ii-1)
           qr = qface(ii)
           qbar = q(ii)
           ilft = ii-1
           ibar = ii
           irgt = ii+1
        end if
        !     
        !# compute ratio of largest to smallest smoothness indicator
        beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)-q(ii-1))**2
        beta2 = (q(ii)-q(ii-1))**2 + (q(ii)-q(ii+1))**2
        beta3 = (q(ii+1)-q(ii))**2 + (q(ii+2)-q(ii+1))**2
        !     
        !# compute ratio of largest to smallest smoothness indicator.
        !     
        lambdavec(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
        !
        !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        ! if lambdavec (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if((lambdavec(i).gt.maxlambda)) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           qmx = max(q(ilft),q(ibar))
           qmn = min(q(ilft),q(ibar))
           ql = MAX(qmn,MIN(qmx,ql))

           qmx = max(q(irgt),q(ibar))
           qmn = min(q(irgt),q(ibar))
           qr = MAX(qmn,MIN(qmx,qr))

           ! if qbar outside of [ql,qr], use piecewise constant,
           !   otherwise use parabola
           if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
              a0 = qbar ! piecewise constant reconstruction
              a1 = 0.d0
              a2 = 0.d0
              a3 = 0.d0 
           else
              ! define coefficients of parabola
              a0 = ql
              a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
              a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
              a3 = 0.d0 
              if(a2.ne.0.d0) then
                 ! check for extrema of the parabola within the interval
                 x1 = -a1/(2.d0*a2)
                 if(ABS(x1-0.5d0).lt.0.5d0) then
                    ! the parabola has an extrema in the interval, so that
                    !   we're going to modify the parabola to remove the extrema
                    if(ABS(qbar-ql).le.ABS(qbar-qr)) then
                       ! extrema is on the ql side of the interval from ql to qr, 
                       !   so that we bound parabola using ql to remove under/overshoot
                       a0 = ql
                       a1 = 0.d0
                       a2 = -3.d0*ql+3.d0*qbar
                       a3 = 0.d0
                    else
                       ! extrema is on the qr side of the interval from ql to qr, 
                       !   so that we bound parabola using qr to remove under/overshoot
                       a0 = 3.d0*qbar-2.d0*qr
                       a1 = 6.d0*qr-6.d0*qbar
                       a2 = -3.d0*qr+3.d0*qbar
                       a3 = 0.d0
                    end if
                 end if
              end if
           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = fx(i) + rhoufrac*qavg(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
     
     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        !   fbc accounts for mean gradient
        fx(0) = dt*(rhou(0)*qavg(1) - fbc(1))
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = dt*fxleft 
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient
        !   fbc accounts for mean gradient
        fx(N) = dt*rhou(N)*qavg(N-1) !+ fbc(2))
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = dt*fxright 
     end if

     !# update solution
     !=================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime(i)
     end do

  end if ! doposlimit

  fx(:) = fx(:)/dt
end subroutine pcmsweep_pmod_select_semilagrangian
! ============================================================
subroutine ppmsweep_fct_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !# inputs
  logical, intent(in) :: doposlimit
  integer, intent(in) :: N, npad, nmaxcfl, bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: dt, scale
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !# in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhoprime(1-npad-nmaxcfl:N+npad+nmaxcfl)

  !# outputs
  real(kind=8) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, imin, imax, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) ::  crface, crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl), qtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-2:N+2), rhoprime_new(-1:N+2)
  real(kind=8) :: beta1, beta2, beta3, tmpeps, lambda(-2:N+2), rhoqtd(-1:N+2)
  real(kind=8) :: secnd(-2-nmaxcfl:N+3+nmaxcfl), tmp_fxcorr(-2:N+2)
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: qup(-2:N+2), qcorr(-2:N+2)
  integer :: monlimit(-3:N+3), itmp, ib, ic
  real(kind=8) :: sdir
  integer :: idir, iup, idn, ii, iu, iu2, ishift(-2:N+2)
  real(kind=8) :: cr, dq, ql, qr, qbar, a0, a1, a2, a3
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16, &
       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  monlimit(:) = 0
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad-nmaxcfl:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad+nmaxcfl) = q(N)
  !     
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  !# compute ppm estimate for flux
  !===============================
  !     
  do i = -1-nmaxcfl,-2+nmaxcfl
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do
  
  do i = -2,N+2
     !# qface holds interpolant of q at faces.
     qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+1+nmaxcfl)) &
          + fac2*(q(i-1+nmaxcfl) + q(i+2+nmaxcfl))

     if(rhou(i).ge.0.d0) then
        fxup(i) = 0.d0
        rhoufrac = rhou(i)*dt
        iup = i
        ii = iup
        do k = 1,nmaxcfl
!!$           if(rhoufrac.lt.rho(ii)) EXIT
!!$           fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
!!$           rhoufrac = rhoufrac - rho(ii) ! integer cfl mass flux
           if(rhoufrac.lt.rhoprime(ii)) EXIT
           fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
           rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
           ii = ii - 1 ! change index of cell being added to flux
        end do
        !
        !# compute upwind and PCM correction flux
        !     
        ishift(i) = iup - ii                ! size of integer shift at face i
        fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
!!$        cr = rhoufrac/rho(ii)
        cr = rhoufrac/rhoprime(ii)
        ql = qface(ii) ! note: ql and qr flipped here.
        qr = qface(ii-1)
        qbar = q(ii)
     else
        fxup(i) = 0.d0
        rhoufrac = rhou(i)*dt
        iup = i+1
        ii = iup
        do k = 1,nmaxcfl
!!$           if(-rhoufrac.lt.rho(ii)) EXIT
!!$           fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
!!$           rhoufrac = rhoufrac + rho(ii) ! integer cfl mass flux
           if(-rhoufrac.lt.rhoprime(ii)) EXIT
           fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
           rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
           ii = ii + 1 ! change index of cell being added to flux
        end do
        !     
        !# compute upwind and PCM correction flux
        !     
        ishift(i) = iup - ii                ! size of integer shift at face i
        fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
!!$        cr = ABS(rhoufrac/rho(ii))
        cr = ABS(rhoufrac/rhoprime(ii))
        ql = qface(ii-1)
        qr = qface(ii)
        qbar = q(ii)
     end if
     !     
     !# compute ratio of largest to smallest smoothness indicator
     beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)-q(ii-1))**2
     beta2 = (q(ii)-q(ii-1))**2 + (q(ii)-q(ii+1))**2
     beta3 = (q(ii+1)-q(ii))**2 + (q(ii+2)-q(ii+1))**2
     !     
     !# compute ratio of largest to smallest smoothness indicator.
     !     
     lambda(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
     !
     !# compute PPM flux by integrating across reconstruction.
     qup(i) = qbar
     qcorr(i) = (1.d0-cr)*(ql - qbar - cr*(ql - 2.d0*qbar + qr))
     fxcorr(i) = rhoufrac*qcorr(i)
     tmp_fxcorr(i) = fxcorr(i)
  end do
  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  imin = 0
  imax = N

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).or.(bcleft.eq.3)) then
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(1) to account for possible mean gradient of scalar
        fxup(0) = dt*rhou(0)*(qup(1) - fbc(1))
        fxcorr(0) = dt*rhou(0)*qcorr(1)

        lambda(0) = lambda(1)
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fxcorr(0) = dt*fxleft - fxup(0) 
        imin = 1 ! no limiting of left boundary flux
     end if

     ! if open or fixed flux bcs at the left end, specify q and fxup
     !   so that cells -1 and 0 match cell 1 (essentially constraining 
     !   any flux correction to be applied only using values from the
     !   domain interior
     fxup(-2:-1) = fxup(0)
     fxcorr(-2:-1) = fxcorr(0)
     tmp_fxcorr(-2:0) = fxcorr(-2:0)
     q(-1:0) = q(1)
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).or.(bcright.eq.3)) then
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(2) to account for possible mean gradient of scalar
        fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
        fxcorr(N) = dt*rhou(N)*qcorr(N-1)

        lambda(N) = lambda(N-1)
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fxcorr(N) = dt*fxright - fxup(N)
        imax = N-1 ! no limiting of right boundary flux
     end if

     ! if open or fixed flux bcs at the left end, specify q and fxup
     !   so that cells -1 and 0 match cell 1 (essentially constraining 
     !   any flux correction to be applied only using values from the
     !   domain interior)
     fxup(N+1:N+2) = fxup(N)
     fxcorr(N+1:N+2) = fxcorr(N)
     tmp_fxcorr(N:N+2) = fxcorr(N:N+2)
     q(N+1:N+2) = q(N)
  end if


  !
  !# update rhoprime -- kept to ensure mass consistency
  rhoprime_new(-1:N+2) = rhoprime(-1:N+2) + dt*(rhou(-2:N+1) - rhou(-1:N+2))
  rhoqtd(-1:N+2) = rhoprime(-1:N+2)*q(-1:N+2) + fxup(-2:N+1) - fxup(-1:N+2)
  qtd(-1:N+2) = rhoqtd(-1:N+2)/rhoprime_new(-1:N+2)

  tmpeps = eps2*scale ! a small parameter, used to prevent division by zero

  !============= selective, positive flux correction =============
  do i = imin,imax

     ! perform monotonic flux correction if lambda>lambda_max at this face.
     if(lambda(i).gt.maxlambda) then
        monlimit(i) = 1

        if(tmp_fxcorr(i).ge.0.d0) then

           ii = i - ishift(i)

           !check whether flux out of cell i will induce new minimum in q
           qmn = min(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)

           !check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
           fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)

        else

           ii = i - ishift(i)

           !check whether flux into cell i will induce new maximum in q 
           qmx = max(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
           fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
           ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)

           !check whether flux out of cell i+1 will induce new min in q
           qmn = min(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
           fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
           ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)

        end if

        !# limit flux correction 
        fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)

     end if

  end do

  if(doposlimit) then

     do i = imin,imax
        ! perform flux correction for non-negativity in all locations
        if(tmp_fxcorr(i).ge.0.d0) then
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoqtd(i)
        else
           fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
           ratiom = rhoqtd(i+1)
        end if
        if(fluxout.gt.ratiom) then
           fxcorr(i) = tmp_fxcorr(i)*max(0.d0,min(1.d0, &
                ratiom/(tmpeps+fluxout), &      ! ratiom for positivity
                fxcorr(i)/(tmp_fxcorr(i)+tmpeps))) ! ratiom for selective monotonicity
        end if
     end do

     !# update solution, eliminating any remaining negative values
     !============================================================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = rhoqtd(i) + fxcorr(i-1) - fxcorr(i) !max(0.d0,rhoq(i) + fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do

  else

     !# update solution
     !=================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do
  end if ! doposlimit

  fx(:) = fx(:)/dt
  limitpts(:) = monlimit(0:N)
  lambdavec(:) = lambda(0:N)

end subroutine ppmsweep_fct_select_semilagrangian
!
! ============================================================
! ============================================================
subroutine ppmsweep_pmod_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !# inputs
  logical, intent(in) :: doposlimit
  integer, intent(in) :: N, npad, nmaxcfl, bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: dt, scale
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !# in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhoprime(1-npad-nmaxcfl:N+npad+nmaxcfl)

  !# outputs
  real(kind=8) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, imin, imax, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: limit(-1:N+1), lambda(-1:N+1)
  real(kind=8) :: tmpeps, qmn, qmx, qavg(0:N)
  real(kind=8) :: ql, qr, qbar, a0, a1, a2, cr, qbarl, qbarr
  real(kind=8) :: discr, x1, x2, x1dis, x2dis
  real(kind=8) ::  crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: beta1, beta2, beta3
  real(kind=8) :: sdir, rhoqtd, fluxout, tmp_fxcorr
  real(kind=8) :: fxup(-1:N+1), fxcorr(-1:N+1)
  real(kind=8) :: qup(-1:N+1), qcorr(-1:N+1)
  integer :: idir, iup, idn, ii, iu, iu2, ishift(-2:N+2)
  logical :: poscorr
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16, &
       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
!!$  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, &
!!$!!$       eps2=1.d-16, epshybrid=1.d-4, maxlambda=20.d0
!!$       eps2=1.d-16, epshybrid=1.d-8, maxlambda=100.d0
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     q(1-npad-nmaxcfl:0) = q(1)
     rhoq(0) = rhoq(1)
  end if
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     q(N+1:N+npad+nmaxcfl) = q(N)
     rhoq(N+1) = rhoq(N)
  end if
  !     
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  if(doposlimit) then
     
     ! INITIALIZE MONLIMIT DIAGNOSTIC
     limit(:) = 0

     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -2-nmaxcfl,nmaxcfl
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = -1,N+1
        ! initialize qface and fdbeta to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+nmaxcfl+1)) &
                            + fac2*(q(i+nmaxcfl-1) + q(i+nmaxcfl+2))

        if(rhou(i).ge.0.d0) then
           fxup(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i
           do k = 1,nmaxcfl
!!$              if(rhoufrac.lt.rho(ii)) EXIT
!!$              fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
!!$              rhoufrac = rhoufrac - rho(ii) ! integer cfl mass flux
              if(rhoufrac.lt.rhoprime(ii)) EXIT
              fxup(i) = fxup(i) + rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
              ii = ii - 1 ! change index of cell being added to flux
           end do
           !
           !# compute upwind and PCM correction flux
           !     
!!$           cr = rhoufrac/rho(ii)
           cr = rhoufrac/rhoprime(ii)
           fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
           ql = qface(ii) ! note: ql and qr flipped here.
           qr = qface(ii-1)
           qbar = q(ii)
        else
           fxup(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i+1
           do k = 1,nmaxcfl
!!$              if(-rhoufrac.lt.rho(ii)) EXIT
!!$              fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
!!$              rhoufrac = rhoufrac + rho(ii) ! integer cfl mass flux
              if(-rhoufrac.lt.rhoprime(ii)) EXIT
              fxup(i) = fxup(i) - rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
              ii = ii + 1 ! change index of cell being added to flux
           end do
           !     
           !# compute upwind and PCM correction flux
           !     
           fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
!!$           cr = ABS(rhoufrac/rho(ii))
           cr = ABS(rhoufrac/rhoprime(ii))
           ql = qface(ii-1)
           qr = qface(ii)
           qbar = q(ii)
        end if
        !     
        !# compute ratio of largest to smallest smoothness indicator
        beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)  -q(ii-1))**2
        beta2 = (q(ii)  -q(ii-1))**2 + (q(ii+1)-q(ii)  )**2
        beta3 = (q(ii+1)-q(ii)  )**2 + (q(ii+2)-q(ii+1))**2
        !     
        !# compute ratio of largest to smallest smoothness indicator.
        !     
        lambda(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
        !
        ! define coefficients of parabola
        a0 =       ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 =  3.d0*ql + 3.d0*qr - 6.d0*qbar

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if(lambda(i).gt.maxlambda) then
           limit(i) = 1

           if(rhou(i).ge.0.d0) then
              qbarl = q(ii+1)
              qbarr = q(ii-1)
           else
              qbarl = q(ii-1)
              qbarr = q(ii+1)
           end if

           ! monotonize the cell edge values
           qmx = max(qbarl,qbar)
           qmn = min(qbarl,qbar)
           ql = MAX(qmn,MIN(qmx,ql))

           qmx = max(qbarr,qbar)
           qmn = min(qbarr,qbar)
           qr = MAX(qmn,MIN(qmx,qr))

           ! if qbar outside of [ql,qr], use piecewise constant,
           !   otherwise use parabola
           if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
              a0 = qbar ! piecewise constant reconstruction
              a1 = 0.d0
              a2 = 0.d0
           else
              ! define coefficients of parabola
              a0 = ql
              a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
              a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
              ! the following check is equivalent to 0 < -a1/(2*a2) < 1
              !   which is the condition for having an extrema in the interval
              if(abs(-a1-a2).lt.abs(a2)) then
                 ! the parabola has an extrema in the interval, so that
                 !   we're going to modify the parabola to remove the extrema
                 if(ABS(qbar-ql).le.ABS(qbar-qr)) then
                    ! extrema is on the ql side of the interval from ql to qr, 
                    !   so that we bound parabola using ql to remove under/overshoot
                    a0 = ql
                    a1 = 0.d0
                    a2 = -3.d0*ql+3.d0*qbar
                 else
                    ! extrema is on the qr side of the interval from ql to qr, 
                    !   so that we bound parabola using qr to remove under/overshoot
                    a0 = 3.d0*qbar-2.d0*qr
                    a1 = 6.d0*qr-6.d0*qbar
                    a2 = -3.d0*qr+3.d0*qbar
                 end if
              end if
           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qup(i) = qbar
        qcorr(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0) - qbar
        fxcorr(i) = rhoufrac*qcorr(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

     imin = 0
     imax = N

     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(1) to account for possible mean gradient of scalar
        fxup(0) = dt*rhou(0)*MAX(0.d0,qup(1) - fbc(1))
        fxcorr(0) = dt*rhou(0)*qcorr(1)
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fxcorr(0) = dt*fxleft - fxup(0) 
        imin = 1 ! no limiting of left boundary flux
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(2) to account for possible mean gradient of scalar
        fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
        fxcorr(N) = dt*rhou(N)*qcorr(N-1)
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fxcorr(N) = dt*fxright - fxup(N)
        imax = N-1 ! no limiting of right boundary flux
     end if

     ! limit correction fluxes at faces to prevent the occurrence of 
     !   negative scalar values
     tmp_fxcorr = fxcorr(-1)
     do i = 0,N
        if(fxcorr(i).ge.0.d0) then
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i) + (fxup(i-1) - fxup(i))
           fluxout = max(0.d0,fxcorr(i)) - min(0.d0,tmp_fxcorr)
           tmp_fxcorr = fxcorr(i)
        else
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i+1) + (fxup(i) - fxup(i+1))
           fluxout = max(0.d0,fxcorr(i+1)) - min(0.d0,fxcorr(i))
           tmp_fxcorr = fxcorr(i)
        end if
        if(fluxout.gt.rhoqtd) then
           fxcorr(i) = fxcorr(i)*MAX(0.d0,MIN(1.d0,rhoqtd/(fluxout+EPSILON(rhoqtd))))
        end if
     end do

     !# update solution, eliminating any remaining negative values
     !============================================================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + (fx(i-1) - fx(i)) !MAX(0.d0,
        q(i) = rhoq(i)/rhoprime(i)
     end do

     lambdavec(0:N) = lambda(0:N)
     limitpts(0:N) = limit(0:N)

  else
     
     ! INITIALIZE MONLIMIT DIAGNOSTIC
     limitpts(:) = 0

     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -1-nmaxcfl,nmaxcfl
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = 0,N
        ! initialize qface and fdbeta to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+nmaxcfl+1)) &
                                          + fac2*(q(i+nmaxcfl-1) + q(i+nmaxcfl+2))

        if(rhou(i).ge.0.d0) then
           fx(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i
           do k = 1,nmaxcfl
              if(rhoufrac.lt.rhoprime(ii)) EXIT
              fx(i) = fx(i) + rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac - rhoprime(ii) ! integer cfl mass flux
              ii = ii - 1 ! change index of cell being added to flux
           end do
           !
           !# compute upwind and PCM correction flux
           !     
           cr = rhoufrac/rhoprime(ii)
           ql = qface(ii) ! note: ql and qr flipped here.
           qr = qface(ii-1)
           qbar = q(ii)
           qbarl = q(ii+1)
           qbarr = q(ii-1)
        else
           fx(i) = 0.d0
           rhoufrac = rhou(i)*dt
           ii = i+1
           do k = 1,nmaxcfl
              if(-rhoufrac.lt.rhoprime(ii)) EXIT
              fx(i) = fx(i) - rhoprime(ii)*q(ii) ! integer cfl part
              rhoufrac = rhoufrac + rhoprime(ii) ! integer cfl mass flux
              ii = ii + 1 ! change index of cell being added to flux
           end do
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhoufrac/rhoprime(ii))
           ql = qface(ii-1)
           qr = qface(ii)
           qbar = q(ii)
           qbarl = q(ii-1)
           qbarr = q(ii+1)
        end if
        !     
        !# compute ratio of largest to smallest smoothness indicator
        beta1 = (q(ii-1)-q(ii-2))**2 + (q(ii)-q(ii-1))**2
        beta2 = (q(ii)-q(ii-1))**2 + (q(ii)-q(ii+1))**2
        beta3 = (q(ii+1)-q(ii))**2 + (q(ii+2)-q(ii+1))**2
        !     
        !# compute ratio of largest to smallest smoothness indicator.
        !     
        lambdavec(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
        !
        ! define coefficients of parabola
        a0 = ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

        ! if lambdavec (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if((lambdavec(i).gt.maxlambda)) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           qmx = max(qbarl,qbar)
           qmn = min(qbarl,qbar)
           ql = MAX(qmn,MIN(qmx,ql))

           qmx = max(qbarr,qbar)
           qmn = min(qbarr,qbar)
           qr = MAX(qmn,MIN(qmx,qr))

           ! if qbar outside of [ql,qr], use piecewise constant,
           !   otherwise use parabola
           if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
              a0 = qbar ! piecewise constant reconstruction
              a1 = 0.d0
              a2 = 0.d0
              !
              ! the following check is equivalent to 0 < -a1/(2*a2) < 1
              !   which is the condition for having an extrema in the interval
           elseif(abs(-a1-a2).lt.abs(a2)) then
              ! the parabola has an extrema in the interval, so that
              !   we're going to modify the parabola to remove the extrema
              if(ABS(qbar-ql).le.ABS(qbar-qr)) then
                 ! extrema is on the ql side of the interval from ql to qr, 
                 !   so that we bound parabola using ql to remove under/overshoot
                 a0 = ql
                 a1 = 0.d0
                 a2 = -3.d0*ql+3.d0*qbar
              else
                 ! extrema is on the qr side of the interval from ql to qr, 
                 !   so that we bound parabola using qr to remove under/overshoot
                 a0 = 3.d0*qbar-2.d0*qr
                 a1 = 6.d0*qr-6.d0*qbar
                 a2 = -3.d0*qr+3.d0*qbar
              end if
           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = fx(i) + rhoufrac*qavg(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
     
     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        !   fbc accounts for mean gradient
        fx(0) = dt*(rhou(0)*qavg(1) - fbc(1))
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = dt*fxleft 
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient
        !   fbc accounts for mean gradient
        fx(N) = dt*rhou(N)*qavg(N-1) !+ fbc(2))
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = dt*fxright 
     end if

     !# update solution
     !=================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime(i)
     end do

  end if ! doposlimit

  fx(:) = fx(:)/dt
end subroutine ppmsweep_pmod_select_semilagrangian
