! ============================================================
subroutine ppmsweep_select(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 


  !# local variables
  integer i, imin, imax, bcleft, bcright
  real(kind=8) :: fxleft, fxright
  logical doselimit, domonlimit

  !# local variables defined at cell faces
  real(kind=8) ::  crface, qface(-2:N+2), qtd(-1:N+2), rhoprime_new(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), rhoqtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1), tmp_fxcorr(-1:N+1)
  real(kind=8) :: beta1, beta2, beta3, tmpeps, lambda(-1:N+1), secnd(-2:N+3)
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: fd(-3:N+3), fdsq(-3:N+3), fdbeta(-3:N+3)
  integer :: monlimit(-2:N+2), itmp, ib, ic
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  monlimit(:) = 0
  !     
  !# compute ppm estimate for flux
  !===============================
  !     
  do i = -2,-1
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do

  if(rhou(-2).ge.0.d0) then
     qup(-2) = q(-2)
  else
     qup(-2) = q(-1)
  end if
  fxup(-2) = dt*rhou(-2)*qup(-2)
  
  !# apply monotonic limiting selectively -- only at faces at or
  !     adjacent to those which are found to be less than smooth using the
  !     ratio of the largest to smallest WENO smoothness ratios as an
  !     indicator. 
  
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  fdbeta(-2:0) = (q(-3:-1)-q(-2:0))**2 + (q(-2:0)-q(-1:1))**2

  do i = -1,N+1
     ! initialize qface, secnd and monlimit to the right of the current face.
     !  - qface(i+1) holds an interpolated scalar value at face i+1
     !  - secnd(i+2) is an estimate of the second derivative of the scalar 
     !       using scalar values from cells i+1, i+2 and i+3
     !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
     qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
     fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

     if(rhou(i).ge.0.d0) then
        !     
        !# compute upwind and PPM correction flux
        !     
        crface = rhou(i)*dt/rho(i) 
        qup(i)   = q(i)
        qcorr(i) = (1.d0-crface) &
             *( qface(i) - q(i) - crface*(qface(i-1) -2.d0*q(i) +qface(i)) )
        !     
        !# compute ratio of largest to smallest smoothness indicator
        lambda(i) = MAX(fdbeta(i-1),fdbeta(i),fdbeta(i+1)) &
                     /(tmpeps+MIN(fdbeta(i-1),fdbeta(i),fdbeta(i+1)))

     else
        !     
        !# compute upwind and PPM correction flux
        !     
        crface = rhou(i)*dt/rho(i+1) 
        qup(i)   = q(i+1)
        qcorr(i) = (1.d0+crface) &
             *(qface(i) - q(i+1) + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
        !     
        !# compute ratio of largest to smallest smoothness indicator
        lambda(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
                     /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
     end if
     !
     !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
     fxup(i) = dt*rhou(i)*qup(i)
     fxcorr(i) = dt*rhou(i)*qcorr(i)
     tmp_fxcorr(i) = fxcorr(i)
     !
     !# update rhoprime -- kept to ensure mass consistency
     rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
     !
     !# compute transported, diffused upwind solution
     rhoqtd(i) = rhoprime(i)*q(i) + fxup(i-1) - fxup(i)
     qtd(i) = rhoqtd(i)/rhoprime_new(i)

  end do

  if(rhou(N+2).ge.0.d0) then
     qup(N+2) = q(N+2)
  else
     qup(N+2) = q(N+3)
  end if
  fxup(N+2) = dt*rhou(N+2)*qup(N+2)
  rhoprime_new(N+2) = rhoprime(N+2) + dt*(rhou(N+1) - rhou(N+2))
  rhoqtd(N+2) = rhoprime(N+2)*q(N+2) + fxup(N+2-1) - fxup(N+2)
  qtd(N+2) = rhoqtd(N+2)/rhoprime_new(N+2)

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  imin = 0
  imax = N

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
  
  !============== LEFT BOUNDARY ==============
  if(((bcleft.eq.1).and.(rhou(0).gt.0.d0)).or.(bcleft.eq.3)) then
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        !   fbc accounts for mean gradient
        if(doposlimit) then
           fxup(0) = dt*rhou(0)*MAX(0.d0,(qup(1) - fbc(1)))
        else
           fxup(0) = dt*rhou(0)*(qup(1) - fbc(1))
        end if
        fxcorr(0) = dt*rhou(0)*qcorr(1)
        tmp_fxcorr(0) = fxcorr(0)
     else
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = dt*fxleft 
        fxcorr(0) = fx(0) - fxup(0)
        tmp_fxcorr(0) = fxcorr(0)

        imin = 1
     end if

     ! update transported, diffused solution in cell 1
     rhoqtd(1) = rhoprime(1)*q(1) + fxup(0) - fxup(1)
     qtd(1) = rhoqtd(1)/rhoprime_new(1)

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(-1:0) = q(1)
     qtd(-1:0) = qtd(1)

  elseif((bcleft.eq.1).and.(rhou(0).lt.0.d0)) then

     ! re-set q(0) and qtd(1) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(-1:0) = q(1)
     qtd(-1:0) = qtd(1)
 
  end if

  !============== RIGHT BOUNDARY ==============
  if(((bcright.eq.1).and.(rhou(N).lt.0.d0)).or.(bcright.eq.3)) then
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient
        !   fbc accounts for mean gradient
        if(doposlimit) then
           fxup(N) = dt*rhou(N)*MAX(0.d0,qup(N-1) + fbc(2))
        else
           fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
        end if
        fxcorr(N) = dt*rhou(N)*qcorr(N-1)
        tmp_fxcorr(N) = fxcorr(N)
     else
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = dt*fxright 
        fxcorr(N) = fx(N) - fxup(N)
        tmp_fxcorr(N) = fxcorr(N)

        imax = N-1
     end if

     ! update transported, diffused solution in cell 1
     rhoqtd(N) = rhoprime(N)*q(N) + fxup(N-1) - fxup(N)
     qtd(N) = rhoqtd(N)/rhoprime_new(N)

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(N+1) = q(N)
     qtd(N+1) = qtd(N)

  elseif((bcright.eq.1).and.(rhou(N).gt.0.d0)) then

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(N+1) = q(N)
     qtd(N+1) = qtd(N)

  end if

  tmpeps = eps2*scale ! a small parameter, used to prevent division by zero

  !============= selective, positive flux correction =============
  do i = imin,imax

     ! perform monotonic flux correction if lambda>lambda_max at this face.
     if(lambda(i).gt.maxlambda) then
        monlimit(i) = 1

        if(tmp_fxcorr(i).ge.0.d0) then

           !check whether flux out of cell i will induce new minimum in q
           qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)

           !check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)

        else

           !check whether flux into cell i will induce new maximum in q 
           qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
           ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)

           !check whether flux out of cell i+1 will induce new min in q
           qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
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
!!$        rhoq(i) = rhoq(i) + fx(i-1) - fx(i) !max(0.d0,rhoq(i) + fx(i-1) - fx(i))
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


end subroutine ppmsweep_select
!
! ============================================================
! ============================================================
! ============================================================
subroutine pcmsweep_select(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambda,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambda(0:N) 


  !# local variables
  integer i, imin, imax, bcleft, bcright
  real(kind=8) :: fxleft, fxright
  logical doselimit, domonlimit, poscorr

  !# local variables defined at cell faces
  real(kind=8) ::  crface, qface(-1:N+1), fdbeta(-1:N+2), qavg(0:N)
  real(kind=8) :: tmpeps, qmn, qmx
  real(kind=8) :: ql, qr, qbar, a0, a1, a2, a3, cr, dq
  real(kind=8) :: discr, x1, x2, x1dis, x2dis
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
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
!!$  if((bcleft.eq.1).or.(bcleft.eq.3)) then
!!$     imin = 3
!!$  else
!!$     imin = -1
!!$  end if
!!$  if((bcright.eq.1).or.(bcright.eq.3)) then
!!$     imax = N-3
!!$  else
!!$     imax = N+1
!!$  end if
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  limitpts(:) = 0
  
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  if(doposlimit) then
     
     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -1,0
        qface(i) = MAX(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do
     fdbeta(-1:1) = (q(-2:0)-q(-1:1))**2 + (q(-1:1)-q(0:2))**2

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        !  - secnd(i+2) is an estimate of the second derivative of the scalar 
        !       using scalar values from cells i+1, i+2 and i+3
        !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
        qface(i+1) = MAX(0.d0,fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3)))
        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = rhou(i)*dt/rho(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambda(i) = MAX(fdbeta(i-1),fdbeta(i),fdbeta(i+1)) &
                /(tmpeps+MIN(fdbeta(i-1),fdbeta(i),fdbeta(i+1)))

        else
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhou(i)*dt/rho(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambda(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
                /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
        end if
        !
        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if(lambda(i).gt.maxlambda) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           if(rhou(i).ge.0.d0) then
              qmx = max(q(i-1),q(i),q(i+1))
              qmn = min(q(i-1),q(i),q(i+1))
           else
              qmx = max(q(i),q(i+1),q(i+2))
              qmn = min(q(i),q(i+1),q(i+2))
           end if

           ql = MAX(0.d0,qmn,MIN(qmx,ql))
           qr = MAX(0.d0,qmn,MIN(qmx,qr))

           ! reduce order, kill cubic term
           a3 = 0.d0

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
              if (qr.gt.ql) then
                 if (3.d0*qbar .lt. 2.d0*ql + qr) then
                    !bound parabola from below to remove undershoot
                    a0 = ql
                    a1 = 0.d0
                    a2 = -3.d0*ql+3.d0*qbar
                 elseif (3.d0*qbar .gt. ql + 2.d0*qr) then
                    !bound parabola from above to remove overshoot
                    a0 = 3.d0*qbar-2.d0*qr
                    a1 = 6.d0*qr-6.d0*qbar
                    a2 = -3.d0*qr+3.d0*qbar
                 end if
              else
                 if (3.d0*qbar .gt. 2.d0*ql + qr) then
                    !bound parabola from above to remove overshoot
                    a0 = ql
                    a1 = 0.d0
                    a2 = -3.d0*ql+3.d0*qbar
                 elseif (3.d0*qbar .lt. ql + 2.d0*qr) then
                    !bound parabola from below to remove undershoot
                    a0 = 3.d0*qbar-2.d0*qr
                    a1 = 6.d0*qr-6.d0*qbar
                    a2 = -3.d0*qr+3.d0*qbar
                 end if
              end if

           end if
        end if

        ! find location of extrema and test whether value is negative there
        discr = a2**2 - 3.d0*a1*a3
        if(discr.ge.0.d0) then
           ! two extrema in general -- find distance from center of [0,1] interval
           x1 = (-a2+sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           x2 = (-a2-sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           x1dis = ABS(x1-0.5d0)
           x2dis = ABS(x2-0.5d0)

           if (((x1dis.lt.0.5d0).AND.(a0+x1*(a1+x1*(a2+x1*a3)).lt.0.d0)) .OR. &
                ((x2dis.lt.0.5d0).AND.(a0+x2*(a1+x2*(a2+x2*a3)).lt.0.d0))) then
              ! one of the extrema is negative, use piecewise
              ! constant reconstruction
              a0 = qbar
              a1 = 0.d0
              a2 = 0.d0
              a3 = 0.d0
           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = rhou(i)*qavg(i)
     end do

  else
     
     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     qface(-1:0) = fac1*(q(-1:0) + q(0:1)) + fac2*(q(-2:-1) + q(1:2))
     fdbeta(-1:1) = (q(-2:0)-q(-1:1))**2 + (q(-1:1)-q(0:2))**2

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        !  - secnd(i+2) is an estimate of the second derivative of the scalar 
        !       using scalar values from cells i+1, i+2 and i+3
        !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = rhou(i)*dt/rho(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambda(i) = MAX(fdbeta(i-1),fdbeta(i),fdbeta(i+1)) &
                /(tmpeps+MIN(fdbeta(i-1),fdbeta(i),fdbeta(i+1)))

        else
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhou(i)*dt/rho(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambda(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
                /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
        end if
        !
        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if((lambda(i).gt.maxlambda)) then !.OR.(i.lt.imin).OR.(i.gt.imax)) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           if(rhou(i).ge.0.d0) then
              qmx = max(q(i-1),q(i),q(i+1))
              qmn = min(q(i-1),q(i),q(i+1))
           else
              qmx = max(q(i),q(i+1),q(i+2))
              qmn = min(q(i),q(i+1),q(i+2))
           end if

           ql = MAX(qmn,MIN(qmx,ql))
           qr = MAX(qmn,MIN(qmx,qr))

           ! reduce order, kill cubic term
           a3 = 0.d0

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
              if (qr.gt.ql) then
                 if (3.d0*qbar .lt. 2.d0*ql + qr) then
                    !bound parabola from below to remove undershoot
                    a0 = ql
                    a1 = 0.d0
                    a2 = -3.d0*ql+3.d0*qbar
                 elseif (3.d0*qbar .gt. ql + 2.d0*qr) then
                    !bound parabola from above to remove overshoot
                    a0 = 3.d0*qbar-2.d0*qr
                    a1 = 6.d0*qr-6.d0*qbar
                    a2 = -3.d0*qr+3.d0*qbar
                 end if
              else
                 if (3.d0*qbar .gt. 2.d0*ql + qr) then
                    !bound parabola from above to remove overshoot
                    a0 = ql
                    a1 = 0.d0
                    a2 = -3.d0*ql+3.d0*qbar
                 elseif (3.d0*qbar .lt. ql + 2.d0*qr) then
                    !bound parabola from below to remove undershoot
                    a0 = 3.d0*qbar-2.d0*qr
                    a1 = 6.d0*qr-6.d0*qbar
                    a2 = -3.d0*qr+3.d0*qbar
                 end if
              end if

           end if
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = rhou(i)*qavg(i)
     end do

  end if

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     !   fbc accounts for mean gradient
     fx(0) = rhou(0)*qavg(1)! - fbc(1))
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary by modifying flux correction.
     fx(0) = fxleft 
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     ! inflow at right, open boundary -- impose zero scalar gradient
     !   fbc accounts for mean gradient
     fx(N) = rhou(N)*qavg(N-1) !+ fbc(2))
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fx(N) = fxright 
  end if

  if(doposlimit) then
     !# update solution, eliminating any remaining negative values
     !============================================================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i)) !MAX(0.d0,
        q(i) = rhoq(i)/rhoprime(i)
     end do

  else

     !# update solution
     !=================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do
  end if ! doposlimit

!!$     write(*,*) bcleft, bcright, fxleft, fxright
!!$     
!!$     i=0
!!$     write(*,990) i, rhoprime(i), q(i), rhou(i), qavg(i), fx(i), rho(i), rhou(i)*dt/rho(i), rhoq(i)
!!$     do i = 1,N
!!$        write(*,990) i, rhoprime(i), q(i), rhou(i), qavg(i), fx(i), rho(i), rhou(i)*dt/rho(i), rhoq(i),  fx(i) - fx(i-1)
!!$     end do
!!$     990 format(i4, 12e12.4)
!!$     STOP 'in pcmsweep_select'
!!$
end subroutine pcmsweep_select
!!$!
!!$! ============================================================
!!$! ============================================================
!!$! ============================================================
!!$subroutine pcmsweep_select(rhoq,q,rhou,rho,rhoprime,fx,dt, &
!!$     N,npad,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
!!$     lambdavec,limitpts)
!!$  !     
!!$  implicit none
!!$
!!$  !     # inputs
!!$  logical doposlimit
!!$  integer, intent(in) :: N
!!$  integer, intent(in) :: npad
!!$  integer, intent(in) :: bctype(2)
!!$  real(kind=8), intent(in) :: rhou(-2:N+2)
!!$  real(kind=8), intent(in) :: rho(1-npad:N+npad)
!!$  real(kind=8), intent(in) :: dt
!!$  real(kind=8), intent(in) :: fbc(2)
!!$  real(kind=8), intent(in) :: scale
!!$  integer, intent(in) :: nselpad
!!$  real(kind=8), intent(in) :: maxlambda
!!$  real(kind=8), intent(in) :: epshybrid
!!$
!!$
!!$  !     # in/outputs
!!$  real(kind=8), intent(inout) :: rhoq(0:N+1)
!!$  real(kind=8), intent(inout) :: q(1-npad:N+npad)
!!$  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)
!!$
!!$  !     # outputs
!!$  real(kind=8), intent(out) ::  fx(0:N) 
!!$  real(kind=8), intent(out) ::  limitpts(0:N) 
!!$  real(kind=8), intent(out) ::  lambdavec(0:N) 
!!$
!!$
!!$  !# local variables
!!$  integer i, imin, imax, bcleft, bcright
!!$  real(kind=8) :: fxleft, fxright
!!$  logical doselimit, domonlimit
!!$
!!$  !# local variables defined at cell faces
!!$  real(kind=8) ::  crface, qface(-2:N+2), qtd(-1:N+2), rhoprime_new(-1:N+2)
!!$  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), rhoqtd(-1:N+2)
!!$  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1), tmp_fxcorr(-1:N+1)
!!$  real(kind=8) :: beta1, beta2, beta3, tmpeps, lambda(-1:N+1), secnd(-2:N+3)
!!$  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
!!$  real(kind=8) :: fd(-3:N+3), fdsq(-3:N+3), fdbeta(-3:N+3)
!!$  real(kind=8) :: ql, qr, qbar, a0, a1, a2, a3, cr, dq
!!$  integer :: monlimit(-2:N+2), itmp, ib, ic
!!$  !# local parameters
!!$  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16, &
!!$       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
!!$  !
!!$  ! set up boundary condition variables
!!$  bcleft = bctype(1)
!!$  bcright = bctype(2)
!!$  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
!!$  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
!!$  !
!!$  ! INITIALIZE MONLIMIT DIAGNOSTIC
!!$  monlimit(:) = 0
!!$  !     
!!$  !# compute ppm estimate for flux
!!$  !===============================
!!$  !     
!!$  do i = -2,-1
!!$     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
!!$  end do
!!$
!!$  if(rhou(-2).ge.0.d0) then
!!$     qup(-2) = q(-2)
!!$  else
!!$     qup(-2) = q(-1)
!!$  end if
!!$  fxup(-2) = dt*rhou(-2)*qup(-2)
!!$  
!!$  !# apply monotonic limiting selectively -- only at faces at or
!!$  !     adjacent to those which are found to be less than smooth using the
!!$  !     ratio of the largest to smallest WENO smoothness ratios as an
!!$  !     indicator. 
!!$  
!!$  !# re-scale epshybrid by the square of scale.  This will prevent 
!!$  !#   epshybrid from overwhelming the smoothness indicators when 
!!$  !#   q is small.
!!$  tmpeps = epshybrid*scale*scale
!!$
!!$  fdbeta(-2:0) = (q(-3:-1)-q(-2:0))**2 + (q(-2:0)-q(-1:1))**2
!!$
!!$  do i = -1,N+1
!!$     ! initialize qface, secnd and monlimit to the right of the current face.
!!$     !  - qface(i+1) holds an interpolated scalar value at face i+1
!!$     !  - secnd(i+2) is an estimate of the second derivative of the scalar 
!!$     !       using scalar values from cells i+1, i+2 and i+3
!!$     !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
!!$     qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
!!$     fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2
!!$
!!$     if(rhou(i).ge.0.d0) then
!!$        !     
!!$        !# compute upwind and PCM correction flux
!!$        !     
!!$        cr = rhou(i)*dt/rho(i)
!!$        dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
!!$        ql = qface(i) ! note: ql and qr flipped here.
!!$        qr = qface(i-1)
!!$        qbar = q(i)
!!$!!$        crface = rhou(i)*dt/rho(i) 
!!$!!$        qup(i)   = q(i)
!!$!!$        qcorr(i) = (1.d0-crface) &
!!$!!$             *( qface(i) - q(i) - crface*(qface(i-1) -2.d0*q(i) +qface(i)) )
!!$        !     
!!$        !# compute ratio of largest to smallest smoothness indicator
!!$        lambda(i) = MAX(fdbeta(i-1),fdbeta(i),fdbeta(i+1)) &
!!$             /(tmpeps+MIN(fdbeta(i-1),fdbeta(i),fdbeta(i+1)))
!!$
!!$     else
!!$        !     
!!$        !# compute upwind and PCM correction flux
!!$        !     
!!$        cr = ABS(rhou(i)*dt/rho(i+1))
!!$        dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
!!$        ql = qface(i)
!!$        qr = qface(i+1)
!!$        qbar = q(i+1)
!!$!!$        crface = rhou(i)*dt/rho(i+1) 
!!$!!$        qup(i)   = q(i+1)
!!$!!$        qcorr(i) = (1.d0+crface) &
!!$!!$             *(qface(i) - q(i+1) + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
!!$        !     
!!$        !# compute ratio of largest to smallest smoothness indicator
!!$        lambda(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
!!$             /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
!!$     end if
!!$     !
!!$     !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
!!$     ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
!!$     a0 =       ql                                  ! constant term
!!$     a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
!!$     a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
!!$     a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term
!!$
!!$     !# compute PCM flux by integrating across reconstruction.
!!$     qup(i) = qbar
!!$     qcorr(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0)) - qbar
!!$     fxup(i) = dt*rhou(i)*qup(i)
!!$     fxcorr(i) = dt*rhou(i)*qcorr(i)
!!$     tmp_fxcorr(i) = fxcorr(i)
!!$     !
!!$     !# update rhoprime -- kept to ensure mass consistency
!!$     rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
!!$     !
!!$     !# compute transported, diffused upwind solution
!!$     rhoqtd(i) = rhoprime(i)*q(i) + fxup(i-1) - fxup(i)
!!$     qtd(i) = rhoqtd(i)/rhoprime_new(i)
!!$
!!$  end do
!!$
!!$  if(rhou(N+2).ge.0.d0) then
!!$     qup(N+2) = q(N+2)
!!$  else
!!$     qup(N+2) = q(N+3)
!!$  end if
!!$  fxup(N+2) = dt*rhou(N+2)*qup(N+2)
!!$  rhoprime_new(N+2) = rhoprime(N+2) + dt*(rhou(N+1) - rhou(N+2))
!!$  rhoqtd(N+2) = rhoprime(N+2)*q(N+2) + fxup(N+2-1) - fxup(N+2)
!!$  qtd(N+2) = rhoqtd(N+2)/rhoprime_new(N+2)
!!$
!!$  !===============================================================
!!$  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
!!$
!!$  imin = 0
!!$  imax = N
!!$
!!$  !===============================================================
!!$  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
!!$  
!!$  !============== LEFT BOUNDARY ==============
!!$  if(((bcleft.eq.1).and.(rhou(0).gt.0.d0)).or.(bcleft.eq.3)) then
!!$     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
!!$        ! inflow at left, open boundary -- impose zero scalar gradient there
!!$        !   fbc accounts for mean gradient
!!$        fxup(0) = dt*rhou(0)*(qup(1) - fbc(1))
!!$        fxcorr(0) = dt*rhou(0)*qcorr(1)
!!$        tmp_fxcorr(0) = fxcorr(0)
!!$     else
!!$        ! specify fixed flux at left boundary by modifying flux correction.
!!$        fx(0) = dt*fxleft 
!!$        fxcorr(0) = fx(0) - fxup(0)
!!$        tmp_fxcorr(0) = fxcorr(0)
!!$
!!$        imin = 1
!!$     end if
!!$
!!$     ! update transported, diffused solution in cell 1
!!$     rhoqtd(1) = rhoprime(1)*q(1) + fxup(0) - fxup(1)
!!$     qtd(1) = rhoqtd(1)/rhoprime_new(1)
!!$
!!$     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
!!$     ! be used in min/max computations while avoiding including
!!$     ! extrapolated scalar values in those min/max computations.
!!$     q(0) = q(1)
!!$     qtd(0) = qtd(1)
!!$
!!$  elseif((bcleft.eq.1).and.(rhou(0).lt.0.d0)) then
!!$
!!$     ! re-set q(0) and qtd(1) to values at cell 1, so that they can
!!$     ! be used in min/max computations while avoiding including
!!$     ! extrapolated scalar values in those min/max computations.
!!$     q(0) = q(1)
!!$     qtd(0) = qtd(1)
!!$ 
!!$  end if
!!$
!!$  !============== RIGHT BOUNDARY ==============
!!$  if(((bcright.eq.1).and.(rhou(N).lt.0.d0)).or.(bcright.eq.3)) then
!!$     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
!!$        ! inflow at right, open boundary -- impose zero scalar gradient
!!$        !   fbc accounts for mean gradient
!!$        fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2))
!!$        fxcorr(N) = dt*rhou(N)*qcorr(N-1)
!!$        tmp_fxcorr(N) = fxcorr(N)
!!$     else
!!$        ! specify fixed flux at right boundary by modifying flux correction.
!!$        fx(N) = dt*fxright 
!!$        fxcorr(N) = fx(N) - fxup(N)
!!$        tmp_fxcorr(N) = fxcorr(N)
!!$
!!$        imax = N-1
!!$     end if
!!$
!!$     ! update transported, diffused solution in cell 1
!!$     rhoqtd(N) = rhoprime(N)*q(N) + fxup(N-1) - fxup(N)
!!$     qtd(N) = rhoqtd(N)/rhoprime_new(N)
!!$
!!$     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
!!$     ! be used in min/max computations while avoiding including
!!$     ! extrapolated scalar values in those min/max computations.
!!$     q(N+1) = q(N)
!!$     qtd(N+1) = qtd(N)
!!$
!!$  elseif((bcright.eq.1).and.(rhou(N).gt.0.d0)) then
!!$
!!$     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
!!$     ! be used in min/max computations while avoiding including
!!$     ! extrapolated scalar values in those min/max computations.
!!$     q(N+1) = q(N)
!!$     qtd(N+1) = qtd(N)
!!$
!!$  end if
!!$
!!$  tmpeps = eps2*scale ! a small parameter, used to prevent division by zero
!!$
!!$  !============= selective, positive flux correction =============
!!$  do i = imin,imax
!!$
!!$     ! perform monotonic flux correction if lambda>lambda_max at this face.
!!$     if(lambda(i).gt.maxlambda) then
!!$        monlimit(i) = 1
!!$
!!$        if(tmp_fxcorr(i).ge.0.d0) then
!!$
!!$           !check whether flux out of cell i will induce new minimum in q
!!$           qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
!!$           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
!!$           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)
!!$
!!$           !check whether flux into cell i+1 will induce new maximum in q
!!$           qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
!!$           fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
!!$           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)
!!$
!!$        else
!!$
!!$           !check whether flux into cell i will induce new maximum in q 
!!$           qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
!!$           fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
!!$           ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)
!!$
!!$           !check whether flux out of cell i+1 will induce new min in q
!!$           qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
!!$           fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
!!$           ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)
!!$
!!$        end if
!!$
!!$        !# limit flux correction 
!!$        fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)
!!$
!!$     end if
!!$
!!$  end do
!!$
!!$  if(doposlimit) then
!!$
!!$     do i = imin,imax
!!$        ! perform flux correction for non-negativity in all locations
!!$        if(tmp_fxcorr(i).ge.0.d0) then
!!$           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
!!$           ratiom = rhoqtd(i)
!!$        else
!!$           fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
!!$           ratiom = rhoqtd(i+1)
!!$        end if
!!$        if(fluxout.gt.ratiom) then
!!$           fxcorr(i) = tmp_fxcorr(i)*max(0.d0,min(1.d0, &
!!$                ratiom/(tmpeps+fluxout), &      ! ratiom for positivity
!!$                fxcorr(i)/(tmp_fxcorr(i)+tmpeps))) ! ratiom for selective monotonicity
!!$        end if
!!$     end do
!!$
!!$     !# update solution, eliminating any remaining negative values
!!$     !============================================================
!!$     fx(0) = fxup(0) + fxcorr(0)
!!$     do i = 1,N
!!$        fx(i) = fxup(i) + fxcorr(i)
!!$        rhoq(i) = rhoqtd(i) + fxcorr(i-1) - fxcorr(i) !max(0.d0,rhoq(i) + fx(i-1) - fx(i))
!!$!!$        rhoq(i) = rhoq(i) + fx(i-1) - fx(i) !max(0.d0,rhoq(i) + fx(i-1) - fx(i))
!!$        q(i) = rhoq(i)/rhoprime_new(i)
!!$        rhoprime(i) = rhoprime_new(i)
!!$     end do
!!$
!!$  else
!!$
!!$     !# update solution
!!$     !=================
!!$     fx(0) = fxup(0) + fxcorr(0)
!!$     do i = 1,N
!!$        fx(i) = fxup(i) + fxcorr(i)
!!$        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
!!$        q(i) = rhoq(i)/rhoprime_new(i)
!!$        rhoprime(i) = rhoprime_new(i)
!!$     end do
!!$  end if ! doposlimit
!!$
!!$  fx(:) = fx(:)/dt
!!$  limitpts(:) = monlimit(0:N)
!!$  lambdavec(:) = lambda(0:N)
!!$
!!$
!!$end subroutine pcmsweep_select
!
! ============================================================
! ============================================================
subroutine ppmsweep_select_wenobeta(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambdavec,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale
  integer, intent(in) :: nselpad
  real(kind=8), intent(in) :: maxlambda
  real(kind=8), intent(in) :: epshybrid


  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 


  !# local variables
  integer i, imin, imax, bcleft, bcright
  real(kind=8) :: fxleft, fxright
  logical doselimit, domonlimit

  !# local variables defined at cell faces
  real(kind=8) ::  crface, qface(-2:N+2), qtd(-1:N+2), rhoprime_new(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), rhoqtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1), tmp_fxcorr(-1:N+1)
  real(kind=8) :: beta1, beta2, beta3, tmpeps, lambda(-1:N+1), secnd(-2:N+3)
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: fd(-3:N+3), fdsq(-3:N+3), fdbeta(-3:N+3)
  integer :: monlimit(-2:N+2), itmp, ib, ic
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
!!$  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, &
!!$!!$       eps2=1.d-16, epshybrid=1.d-4, maxlambda=20.d0
!!$       eps2=1.d-16, epshybrid=1.d-8, maxlambda=100.d0
  !

  doselimit = .true.
  domonlimit = .false.

  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !     
  !# compute ppm estimate for flux
  !===============================
  !     
  do i = -2,-1
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do

  if(rhou(-2).ge.0.d0) then
     qup(-2) = q(-2)
  else
     qup(-2) = q(-1)
  end if
  fxup(-2) = dt*rhou(-2)*qup(-2)
  
  if(doselimit) then
     !# apply monotonic limiting selectively -- only at faces at or
     !     adjacent to those which are found to be less than smooth using the
     !     ratio of the largest to smallest WENO smoothness ratios as an
     !     indicator. 

     !# re-scale epshybrid by the square of scale.  This will prevent 
     !#   epshybrid from overwhelming the smoothness indicators when 
     !#   q is small.
     tmpeps = epshybrid*scale*scale

     ! itmp is the number of faces where the WENO smoothness indicators
     itmp = 0 

     do i = -2,0
        secnd(i) = (13.d0/12.d0)*(q(i-1) - 2.d0*q(i) + q(i+1))**2
     end do
!!$     fdbeta(-2:0) = (q(-3:-1)-q(-2:0))**2 + (q(-2:0)-q(-1:1))**2

     monlimit(-2:1) = 0

     do i = -1,N+1
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        !  - secnd(i+2) is an estimate of the second derivative of the scalar 
        !       using scalar values from cells i+1, i+2 and i+3
        !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        monlimit(i+1) = 0

        secnd(i+2) = (13.d0/12.d0)*(q(i+1) - 2.d0*q(i+2) + q(i+3))**2
!!$        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

        ! estimate cfl number at face i
        !  - note that an average rho is used to remove any bias due to 
        !      density gradients.
        crface = 2.d0*rhou(i)*dt/(rho(i)+rho(i+1)) 


        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PPM correction flux
           !     
           qup(i)   = q(i)
           qcorr(i) = (1.d0-crface) &
                *( qface(i) - q(i) - crface*(qface(i-1) -2.d0*q(i) +qface(i)) )
           !     
           !# compute WENO5 smoothness indicators
           !     
!!$           beta1 = fdbeta(i-1)
!!$           beta2 = fdbeta(i)
!!$           beta3 = fdbeta(i+1)
!!$
           beta1 = secnd(i-1) + 0.25d0*(q(i-2) - 4.d0*q(i-1) + 3.d0*q(i))**2
           beta2 = secnd(i)   + 0.25d0*(q(i-1) - q(i+1))**2
           beta3 = secnd(i+1) + 0.25d0*(3.d0*q(i) - 4.d0*q(i+1) + q(i+2))**2

        else
           !     
           !# compute upwind and PPM correction flux
           !     
           qup(i)   = q(i+1)
           qcorr(i) = (1.d0+crface) &
                *(qface(i) - q(i+1) + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
           !     
           !# compute WENO5 smoothness indicators
           !     
!!$           beta1 = fdbeta(i+2)
!!$           beta2 = fdbeta(i+1)
!!$           beta3 = fdbeta(i)
!!$
           beta1 = secnd(i+2) + 0.25d0*(3.d0*q(i+1) - 4.d0*q(i+2) + q(i+3))**2
           beta2 = secnd(i+1) + 0.25d0*(q(i) - q(i+2))**2
           beta3 = secnd(i)   + 0.25d0*(q(i-1) - 4.d0*q(i) + 3.d0*q(i+1))**2

        end if
        !
        !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
        fxup(i) = dt*rhou(i)*qup(i)
        fxcorr(i) = dt*rhou(i)*qcorr(i)
        tmp_fxcorr(i) = fxcorr(i)
        !
        !# update rhoprime -- kept to ensure mass consistency
        rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        !
        !# compute transported, diffused upwind solution
        rhoqtd(i) = rhoprime(i)*q(i) + fxup(i-1) - fxup(i)
        !     
        !# apply monotonic limiting at faces at/adjacent to ones where
        !#   ratio of largest to smallest WENO smoothness indicator exceeds 
        !#   maxlamba -- inspired by the hybrid WENO of Hill & Pullin (JCP 2004) 
        !     

        lambda(i)=MAX(beta1,beta2,beta3)/(tmpeps+MIN(beta1,beta2,beta3))
        if(lambda(i).gt.maxlambda) then
!!$        if(MAX(beta1,beta2,beta3) &
!!$             .gt.(tmpeps+maxlambda*MIN(beta1,beta2,beta3))) then
           itmp = itmp + 1
           monlimit(i-1:i+1) = 1
        end if

     end do

  else

     !# non-selective monotonic limiting (either everywhere or nowhere)
     if(domonlimit) then
        lambda(:) = 0.
        monlimit(:) = 1
        itmp = N+1
     else
        lambda(:) = 0.
        monlimit(:) = 0
        itmp = 0
     end if

     do i = -1,N+1
        ! initialize qface
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))

        ! estimate cfl number at face i
        !  - note that an average rho is used to remove any bias due to 
        !      density gradients.
        crface = 2.d0*rhou(i)*dt/(rho(i)+rho(i+1)) 

        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PPM correction flux
           !     
           qup(i)   = q(i)
           qcorr(i) = (1.d0-crface) &
                *( qface(i) - q(i) - crface*(qface(i-1) -2.d0*q(i) +qface(i)) )
        else
           !     
           !# compute upwind and PPM correction flux
           !     
           qup(i)   = q(i+1)
           qcorr(i) = &
                (1.d0+crface)*(qface(i) - q(i+1) &
                               + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
        end if
        !
        !# set up upwind and correction flux !!!! NOTE: SCALED BY DT !!!!!!
        fxup(i) = dt*rhou(i)*qup(i)
        fxcorr(i) = dt*rhou(i)*qcorr(i)
        tmp_fxcorr(i) = fxcorr(i)
        !
        !# update rhoprime -- kept to ensure mass consistency
        rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        !
        !# compute transported, diffused upwind solution
        rhoqtd(i) = rhoprime(i)*q(i) + fxup(i-1) - fxup(i)
     end do

  end if

  if(rhou(N+2).ge.0.d0) then
     qup(N+2) = q(N+2)
  else
     qup(N+2) = q(N+3)
  end if
  fxup(N+2) = dt*rhou(N+2)*qup(N+2)
  rhoprime_new(N+2) = rhoprime(N+2) + dt*(rhou(N+1) - rhou(N+2))
  rhoqtd(N+2) = rhoprime(N+2)*q(N+2) + fxup(N+2-1) - fxup(N+2)

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  imin = 0
  imax = N

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     fxup(0) = dt*rhou(0)*(qup(1) - fbc(1)) !fbc accounts for mean gradient
     fxcorr(0) = dt*rhou(0)*qcorr(1)
     tmp_fxcorr(0) = fxcorr(0)
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary
     fxcorr(0) = dt*fxleft - fxup(0)
     imin = 1 ! no limiting of left boundary flux
     tmp_fxcorr(0) = fxcorr(0)
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     ! inflow at right, open boundary -- impose zero scalar gradient there
     fxup(N) = dt*rhou(N)*(qup(N-1) + fbc(2)) !fbc accounts for mean gradient
     fxcorr(N) = dt*rhou(N)*qcorr(N-1)
     tmp_fxcorr(N) = fxcorr(N)
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fxcorr(N) = dt*fxleft - fxup(N)
     tmp_fxcorr(N) = fxcorr(N)
     imax = N-1 ! no limiting of right boundary flux
  end if

  tmpeps = eps2*scale ! a small parameter, used to prevent division by zero

  if(doposlimit) then

     if(itmp.gt.0) then 

        !====LIMIT FOR MONOTONICITY AND POSITIVITY IN SOME PLACES ========
!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of correction flux
        qtd(-1:2) = rhoqtd(-1:2)/rhoprime_new(-1:2) ! scalar after upwind step

        do i = imin,imax

!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of correction flux
           qtd(i+2) = rhoqtd(i+2)/rhoprime_new(i+2) ! scalar after upwind step

           if(monlimit(i).gt.0) then

              if(tmp_fxcorr(i).ge.0.d0) then

                 !check whether flux out of cell i will induce new minimum in q
                 qmn = max(0.d0, &
                      min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1)))
                 fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
                 ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)

                 !check whether flux into cell i+1 will induce new maximum in q
                 qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
                 fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
                 ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)

              else

                 !check whether flux into cell i will induce new maximum in q 
                 qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
                 fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
                 ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)

                 !check whether flux out of cell i+1 will induce new min in q
                 qmn = max(0.d0, &
                      min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2)))
                 fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
                 ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)

              end if

              !# limit flux correction 
              fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)

           end if

        end do

     end if

     ! check whether limiting for positivity is necessary
!!$     if(2.d0*MAXVAL(ABS(fxcorr)).gt. &
!!$          MINVAL(rhoqtd(-1:N+2))) then

        !==========================================
        !===LIMIT FOR NON-NEGATIVITY EVERYWHERE ===

!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of flux correction
        do i = imin,imax
!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of flux correction
           if(tmp_fxcorr(i).ge.0.d0) then
              fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
              ratiom = rhoqtd(i)
           else
              fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
              ratiom = rhoqtd(i+1)
           end if

           if(fluxout.gt.ratiom) then
              !# limit flux correction for positivity
              ratiom = ratiom/(fluxout + EPSILON(1.d0))
              fxcorr(i) = max(0.d0,min(1.d0,ratiom))*tmp_fxcorr(i)
           end if

        end do

!!$     end if

     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        !# update solution, eliminating any remaining negative values
        !============================================================
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = max(0.d0,rhoq(i) + fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do

  else

     ! NO LIMITING FOR POSITIVITY

     if(itmp.gt.0) then 
        !======================================================================
        !============ LIMIT ONLY FOR MONOTONICITY, ONLY IN SOME PLACES ========

!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of correction flux
        qtd(-1:2) = rhoqtd(-1:2)/rhoprime_new(-1:2) ! scalar after upwind step

        do i = imin,imax

!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of correction flux
           qtd(i+2) = rhoqtd(i+2)/rhoprime_new(i+2) ! scalar after upwind step

           if(monlimit(i).gt.0) then

              if(tmp_fxcorr(i).ge.0.d0) then

                 !check whether flux out of cell i will induce new minimum in q
                 qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
                 fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
                 ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)

                 !check whether flux into cell i+1 will induce new maximum in q
                 qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
                 fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
                 ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)

              else

                 !check whether flux into cell i will induce new maximum in q 
                 qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
                 fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
                 ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)

                 !check whether flux out of cell i+1 will induce new min in q
                 qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
                 fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
                 ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)

              end if

              !# limit flux correction 
              fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)

           end if ! monlimit(i).gt.0
        end do ! i = imin,imax
     end if ! itmp.gt.0

     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        !# update solution
        !=================
        fx(i) = fxup(i) + fxcorr(i)
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do
  end if

  fx(:) = fx(:)/dt
  limitpts(:) = monlimit(0:N)
  lambdavec(:) = lambda(0:N)


end subroutine ppmsweep_select_wenobeta
!
! ============================================================
subroutine ppmsweep(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,domonlimit,doposlimit,lambdavec,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit, domonlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)

  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: cr, ql, qr, qbar, qavg(-1:N+1)
  real(kind=8) ::  qface(-2:N+2), rhoprime_new(-1:N+2), rhoqtd(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1), tmp_fxcorr(-1:N+1)
  real(kind=8) ::  fluxout, ratiom
  
  integer :: monlimit(-2:N+2), itmp
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
  !
  if(domonlimit) then
     call ppmsweep_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
          N,npad,bctype,fbc,domonlimit,doposlimit,lambdavec,limitpts) 
     return
  end if
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)

  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)
  !
  if(.not.doposlimit) then
     ! =========== no flux correction ================
     do i = -1,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = 0,N
        !# qface holds interpolant of q at faces.
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))

        if(rhou(i).ge.0.d0) then
           !# compute upwind and PPM correction flux
           cr = rhou(i)*dt/rho(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rho(i+1))
           ql = qface(i)
           qbar = q(i+1)
           qr = qface(i+1)
        end if
        qavg(i) = qbar + (1.d0-cr)*(ql - qbar - cr*(qr -2.d0*qbar +ql))
        fx(i) = rhou(i)*qavg(i)
     end do
     
     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
     
     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        fx(0) = rhou(0)*(qavg(1) - fbc(1)) !fbc accounts for mean gradient
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = fxleft 
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient there
        fx(N) = rhou(N)*(qavg(N-1) + fbc(2)) !fbc accounts for mean gradient
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = fxright
     end if

     !# update solution using total flux
     !==================================
     !     
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do

     ! =========== END OF ROUTINE FOR NO FLUX CORRECTION ===============

  else

     ! ============== flux correction for non-negativity ===============
     do i = -2,1
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     ! compute upwind flux at i=-2
     if(rhou(-2).ge.0.d0) then
        qup(-2) = q(-2)
     else
        qup(-2) = q(-1)
     end if
     fxup(-2) = rhou(-2)*qup(-2)

     ! compute upwind and correction fluxes at two faces up to left
     ! face of this domain 
     do i = -1,0
        if(rhou(i).ge.0.d0) then
           !# compute upwind and PPM correction flux
           cr = rhou(i)*dt/rho(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rho(i+1))
           ql = qface(i)
           qbar = q(i+1)
           qr = qface(i+1)
        end if
        ! upwind and correction scalar values
        qup(i) = qbar
        qcorr(i) = (1.d0-cr)*(ql - qbar - cr*(qr -2.d0*qbar +ql))
        ! upwind and correction scalar fluxes
        fxup(i) = rhou(i)*qbar
        fxcorr(i) = rhou(i)*qcorr(i)
        tmp_fxcorr(i) = fxcorr(i) ! copy of correction flux for use in fct
        ! updated mass (after applying mass fluxes only in this direction)
        rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        ! transported, diffused scalar mass and mixing ratio
        rhoqtd(i) = rhoprime(i)*q(i) + dt*(fxup(i-1) - fxup(i))
     end do
     
     ! compute upwind and correction fluxes across domain
     ! -- correct fluxes for non-negativity if necessary
     do i = 1,N+1
        !# qface holds interpolant of q at faces.
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))

        if(rhou(i).ge.0.d0) then
           !# compute upwind and PPM correction flux
           cr = rhou(i)*dt/rho(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rho(i+1))
           ql = qface(i)
           qbar = q(i+1)
           qr = qface(i+1)
        end if
        ! upwind and correction scalar values
        qup(i) = qbar
        qcorr(i) = (1.d0-cr)*(ql - qbar - cr*(qr -2.d0*qbar +ql))
        ! upwind and correction scalar fluxes
        fxup(i) = rhou(i)*qbar
        fxcorr(i) = rhou(i)*qcorr(i)
        tmp_fxcorr(i) = fxcorr(i) ! copy of correction flux for use in fct
        ! updated mass (after applying mass fluxes only in this direction)
        rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        ! transported, diffused scalar mass and mixing ratio
        rhoqtd(i) = rhoprime(i)*q(i) + dt*(fxup(i-1) - fxup(i))

        ! apply flux correction at face i-1 to preserve non-negativity
        if(fxcorr(i-1).ge.0.d0) then
           !# check whether flux out of cell i-1 will evacuate that cell
           fluxout = max(0.d0,tmp_fxcorr(i-1)) - min(0.d0,tmp_fxcorr(i-2))
           ratiom = rhoqtd(i-1)/dt/(fluxout+eps2)
        else
           !# check whether flux out of cell i will evacuate that cell
           fluxout = - min(0.d0,tmp_fxcorr(i-1)) + max(0.d0,tmp_fxcorr(i))
           ratiom = rhoqtd(i)/dt/(fluxout+eps2)
        end if
        fx(i-1) = fxup(i-1) + max(0.d0,min(1.d0,ratiom))*fxcorr(i-1)
     end do
     
     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
     
     !============== LEFT BOUNDARY ==============
     if(((bcleft.eq.1).and.(rhou(0).gt.0.d0)).or.(bcleft.eq.3)) then
        if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
           ! inflow at left, open boundary -- impose zero scalar gradient there
           !   fbc accounts for mean gradient
           fx(0) = max(0.d0,rhou(0)*(qup(1) + qcorr(1) - fbc(1)))
           fxup(0) = rhou(0)*(qup(1) - fbc(1))
           tmp_fxcorr(0) = fx(0) - fxup(0)
        else
           ! specify fixed flux at left boundary by modifying flux correction.
           fx(0) = fxleft 
           fxcorr(0) = fx(0) - fxup(0)
           tmp_fxcorr(0) = fxcorr(0)
        end if

        ! if flux at face 1 is out of cell 1, recompute flux
        ! correction, since flux at face 0 has changed.
        if(tmp_fxcorr(1).ge.0.d0) then
           rhoqtd(1) = rhoprime(1)*q(1) + dt*(fxup(0) - fxup(1))
           !# check whether flux out of cell i-1 will evacuate that cell
           fluxout = max(0.d0,tmp_fxcorr(1)) - min(0.d0,tmp_fxcorr(0))
           ratiom = rhoqtd(1)/dt/(fluxout+eps2)
           fx(1) = fxup(1) + max(0.d0,min(1.d0,ratiom))*tmp_fxcorr(1)
        end if
     end if

     !============== RIGHT BOUNDARY ==============
     if(((bcright.eq.1).and.(rhou(N).lt.0.d0)).or.(bcright.eq.3)) then
        if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
           ! inflow at right, open boundary -- impose zero scalar gradient
           !   fbc accounts for mean gradient
           fx(N) = rhou(N)*max(0.d0,qup(N-1) + qcorr(N-1) + fbc(2)) 
           fxup(N) = rhou(N)*(qup(N-1) + fbc(2))
           tmp_fxcorr(N) = fx(N) - fxup(N)
        else
           ! specify fixed flux at right boundary by modifying flux correction.
           fx(N) = fxright 
           fxcorr(N) = fx(N) - fxup(N)
           tmp_fxcorr(N) = fxcorr(N)
        end if

        ! if flux at face N-1 is out of cell N, recompute flux
        ! correction, since flux at face N has changed.
        if(tmp_fxcorr(N-1).lt.0.d0) then
           rhoqtd(N) = rhoprime(N)*q(N) + dt*(fxup(N-1) - fxup(N))
           !# check whether flux out of cell i-1 will evacuate that cell
           fluxout = max(0.d0,tmp_fxcorr(N-1)) - min(0.d0,tmp_fxcorr(N))
           ratiom = rhoqtd(N)/dt/(fluxout+eps2)
           fx(N-1) = fxup(N-1) + max(0.d0,min(1.d0,ratiom))*tmp_fxcorr(N-1)
        end if
     end if

     !# update solution using total flux
     !==================================
     !     
     do i = 1,N
        rhoq(i) = max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
        rhoprime(i) = rhoprime_new(i)
        q(i) = rhoq(i)/rhoprime_new(i)
     end do

     ! ============ end of routine for positive flux correction ==============

  end if

  limitpts(:) = 0.
  lambdavec(:) = 0.

end subroutine ppmsweep


!
! ============================================================
subroutine ppmsweep_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,domonlimit,doposlimit,lambdavec,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit, domonlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)

  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, ii, bcleft, bcright, imin, imax
  real(kind=8) :: fxleft, fxright, tmpeps

  !# local variables defined at cell faces
  real(kind=8) :: cr, ql, qr, qbar, qavg(-1:N+1)
  real(kind=8) ::  qface(-2:N+2), rhoprime_new(-1:N+2), rhoqtd(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), qtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1), tmp_fxcorr(-1:N+1)
  real(kind=8) ::  fluxin, fluxout, ratiom, ratiop, qmn, qmx
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !
  !  FILL GHOST CELLS IF INFLOW AT BOUNDARIES
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad) = q(N)

  ! ============== flux correction for monotonicity ===============
  do i = -2,-1
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do

  ! compute upwind flux at i=-2
  if(rhou(-2).ge.0.d0) then
     qup(-2) = q(-2)
  else
     qup(-2) = q(-1)
  end if
  fxup(-2) = rhou(-2)*qup(-2)

  ! compute upwind and correction fluxes at two faces up to left
  ! face of this domain 
  do i = -1,N+1
     qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
     if(rhou(i).ge.0.d0) then
        !# compute upwind and PPM correction flux
        cr = rhou(i)*dt/rho(i)
        ql = qface(i)
        qbar = q(i)
        qr = qface(i-1)
     else
        !# compute upwind and PPM correction flux
        cr = ABS(rhou(i)*dt/rho(i+1))
        ql = qface(i)
        qbar = q(i+1)
        qr = qface(i+1)
     end if
     ! upwind and correction scalar values
     qup(i) = qbar
     qcorr(i) = (1.d0-cr)*(ql - qbar - cr*(qr -2.d0*qbar +ql))
     ! upwind and correction scalar fluxes
     fxup(i) = rhou(i)*qbar
     fxcorr(i) = rhou(i)*qcorr(i)
     tmp_fxcorr(i) = fxcorr(i) ! copy of correction flux for use in fct
     ! updated mass (after applying mass fluxes only in this direction)
     rhoprime_new(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
     ! transported, diffused scalar mass and mixing ratio
     rhoqtd(i) = rhoprime(i)*q(i) + dt*(fxup(i-1) - fxup(i))
     qtd(i) = rhoqtd(i)/rhoprime_new(i)
  end do

  ! compute upwind flux at i=N+2
  if(rhou(N+2).ge.0.d0) then
     qup(N+2) = q(N+2)
  else
     qup(N+2) = q(N+3)
  end if
  fxup(N+2) = rhou(N+2)*qup(N+2)
  ! updated mass (after applying mass fluxes only in this direction)
  rhoprime_new(N+2) = rhoprime(N+2) + dt*(rhou(N+1) - rhou(N+2))
  ! transported, diffused scalar mass and mixing ratio
  rhoqtd(N+2) = rhoprime(N+2)*q(N+2) + dt*(fxup(N+1) - fxup(N+2))
  qtd(N+2) = rhoqtd(N+2)/rhoprime_new(N+2)

  imin = 0
  imax = N

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
  
  !============== LEFT BOUNDARY ==============
  if(((bcleft.eq.1).and.(rhou(0).gt.0.d0)).or.(bcleft.eq.3)) then
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        !   fbc accounts for mean gradient
        fx(0) = rhou(0)*(qup(1) + qcorr(1) - fbc(1)) 
        fxup(0) = rhou(0)*(qup(1) - fbc(1))
        tmp_fxcorr(0) = fx(0) - fxup(0)
     else
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = fxleft 
        fxcorr(0) = fx(0) - fxup(0)
        tmp_fxcorr(0) = fxcorr(0)

        imin = 1
     end if

     ! update transported, diffused solution in cell 1
     rhoqtd(1) = rhoprime(1)*q(1) + dt*(fxup(0) - fxup(1))
     qtd(1) = rhoqtd(1)/rhoprime_new(1)

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(-1:0) = q(1)
     qtd(-1:0) = qtd(1)

  elseif((bcleft.eq.1).and.(rhou(0).lt.0.d0)) then

     ! re-set q(0) and qtd(1) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(-1:0) = q(1)
     qtd(-1:0) = qtd(1)

  end if

  !============== RIGHT BOUNDARY ==============
  if(((bcright.eq.1).and.(rhou(N).lt.0.d0)).or.(bcright.eq.3)) then
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient
        !   fbc accounts for mean gradient
        fx(N) = rhou(N)*max(0.d0,qup(N-1) + qcorr(N-1) + fbc(2)) 
        fxup(N) = rhou(N)*(qup(N-1) + fbc(2))
        tmp_fxcorr(N) = fx(N) - fxup(N)
     else
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = fxright 
        fxcorr(N) = fx(N) - fxup(N)
        tmp_fxcorr(N) = fxcorr(N)

        imax = N-1
     end if

     ! update transported, diffused solution in cell 1
     rhoqtd(N) = rhoprime(N)*q(N) + dt*(fxup(N-1) - fxup(N))
     qtd(N) = rhoqtd(N)/rhoprime_new(N)

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(N+1:N+2) = q(N)
     qtd(N+1:N+2) = qtd(N)

  elseif((bcright.eq.1).and.(rhou(N).gt.0.d0)) then

     ! re-set q(0) and qtd(0) to values at cell 1, so that they can
     ! be used in min/max computations while avoiding including
     ! extrapolated scalar values in those min/max computations.
     q(N+1:N+2) = q(N)
     qtd(N+1:N+2) = qtd(N)

  end if

  ! -- correct fluxes for monotonicity if necessary
  do i = imin,imax
     ! apply flux correction at face i-1 to preserve monotonicity
     if(fxcorr(i).ge.0.d0) then

        ! check whether flux out of cell i (i-1) will induce new minimum in q
        qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
        fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
        ratiom = (rhoqtd(i)-rhoprime_new(i)*qmn)/dt/(fluxout+eps2)

        ! check whether flux into cell i+1 (i) will induce new maximum in q
        qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
        fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
        ratiop = (rhoprime_new(i+1)*qmx-rhoqtd(i+1))/dt/(fluxin+eps2)

     else

        ! check whether flux into cell i (i-1) will induce new maximum in q 
        qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
        fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
        ratiop = (rhoprime_new(i)*qmx-rhoqtd(i))/dt/(fluxin+eps2)

        !check whether flux out of cell i+1 (i) will induce new minimum in q
        qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
        fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
        ratiom = (rhoqtd(i+1)-rhoprime_new(i+1)*qmn)/dt/(fluxout+eps2)

     end if

     !# add limited flux correction 
     fx(i) = fxup(i) + max(0.d0,min(1.d0,ratiom,ratiop))*fxcorr(i)
  end do

  if(doposlimit) then
     !# update solution, eliminating any remaining negative values
     !============================================================
     do i = 1,N
        rhoq(i) = max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
        rhoprime(i) = rhoprime_new(i)
        q(i) = rhoq(i)/rhoprime_new(i)
     end do
  else
     !# update solution using total flux
     !==================================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime_new(i)
        q(i) = rhoq(i)/rhoprime_new(i)
     end do

  end if

  rhoprime(1:N) = rhoprime_new(1:N)
  limitpts(:) = 1.
  lambdavec(:) = 0.

end subroutine ppmsweep_monotonic


!
! ============================================================
subroutine lwsweep(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,domonlimit,doposlimit)
  !     
  implicit none

  !     # inputs
  logical domonlimit, doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)

  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8) ::  fx(0:N) 

  !     # local variables
  integer i, bcleft, bcright
  real(kind=8) :: fxleft, fxright, tmpeps

  !     # local variables defined at cell faces
  real(kind=8) ::  crface, qface(-2:N+2), qtd(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), fxcorr(-1:N+1) 
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: rhoprime_new(-1:N+2)
  integer :: monlimit(-2:N+2), itmp
  !     # local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
  !
  monlimit(:) = 0
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)
  !     
  !     # compute ppm estimate for flux
  !     ===============================
  !     
  do i = -2,-1
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do

  if(rhou(-2).ge.0.d0) then
     qup(-2) = q(-2)
  else
     qup(-2) = q(-1)
  end if

  do i = -1,N+1
     !     # qface holds interpolant of q at faces.
     qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
     crface = 2.d0*rhou(i)*dt/(rho(i)+rho(i+1))

     if(rhou(i).ge.0.d0) then
        !     
        !     # compute upwind and PPM correction flux
        !     
        qup(i)   = q(i)
        qcorr(i) = (1.d0-crface) &
             *( q(i+1) - q(i))/2.d0
!!$             *( qface(i) - q(i) - crface*(qface(i-1) -2.d0*q(i) +qface(i)) )
     else
        !     
        !     # compute upwind and PPM correction flux
        !     
        qup(i)   = q(i+1)
        qcorr(i) = (1.d0+crface) &
             *(q(i) - q(i+1))/2.d0
!!$             *(qface(i) - q(i+1) + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
     end if
  end do

  if(rhou(N+2).ge.0.d0) then
     qup(N+2) = q(N+2)
  else
     qup(N+2) = q(N+3)
  end if

  !     # fix boundary conditions (for open bcs) if necessary
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     qup(0) = qup(1) - fbc(1)
     qcorr(0) = qcorr(1)
  end if
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     qup(N) = qup(N-1) + fbc(2)
     qcorr(N) = qcorr(N-1)
  end if

  !     # update rhoprime -- kept to ensure mass consistency
  rhoprime_new(-1:N+2) = rhoprime(-1:N+2) + dt*(rhou(-2:N+1) - rhou(-1:N+2))

  if(.not.domonlimit) then

     ! define full PPM flux
     fx(0:N) = rhou(0:N)*(qup(0:N) + qcorr(0:N))

  else

     !!!!!!!!!! MONOTONIC LIMITING EVERYWHERE !!!!!!!!!!!!!!!
     ! DEFINE TRANSPORTED AND DIFFUSED SOLUTION AND ADD ONLY THAT
     ! FRACTION OF CORRECTION FLUX TO THE UPWIND FLUX THAT WILL ALLOW LOCAL
     ! MONOTONICITY TO BE PRESERVED.
     !     
     !     # find new scalar mass using only upwind flux
     qtd(-1:N+2) = &
          (rhoprime(-1:N+2)*q(-1:N+2) &
          + dt*(rhou(-2:N+1)*qup(-2:N+1) - rhou(-1:N+2)*qup(-1:N+2))) &
          /rhoprime_new(-1:N+2)
     !     
     !     # set up total flux (just upwind now) and correction flux
     fx(0:N) = rhou(0:N)*qup(0:N)
     fxcorr(-1:N+1) = rhou(-1:N+1)*qcorr(-1:N+1)

     !     # if fixed flux bcs, fix flux correction
     !     ========================================
     if(bcleft.eq.3) fxcorr(0) = fxleft - fx(0)
     if(bcright.eq.3) fxcorr(N) = fxright - fx(N)

     tmpeps = eps2

     do i = 0,N
        if(fxcorr(i).ge.0.d0) then

           ! check whether flux out of cell i will induce new minimum in q
           qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxout = max(0.d0,fxcorr(i)) - min(0.d0,fxcorr(i-1))
           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/dt/(fluxout+tmpeps)

           ! check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxin = max(0.d0,fxcorr(i)) - min(0.d0,fxcorr(i+1))
           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/dt/(fluxin+tmpeps)

        else

           ! check whether flux into cell i will induce new maximum in q 
           qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxin = - min(0.d0,fxcorr(i)) + max(0.d0,fxcorr(i-1))
           ratiop = rhoprime_new(i)*(qmx-qtd(i))/dt/(fluxin+tmpeps)

           !check whether flux out of cell i+1 will induce new minimum in q
           qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxout = - min(0.d0,fxcorr(i)) + max(0.d0,fxcorr(i+1))
           ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/dt/(fluxout+tmpeps)

        end if

        !     # add limited flux correction 
        fx(i) = fx(i) + max(0.d0,min(1.d0,ratiom,ratiop))*fxcorr(i)

     end do

  end if

  !     # if fixed flux bcs, fix flux
  !     =============================
  if(bcleft.eq.3) fx(0) = fxleft
  if(bcright.eq.3) fx(N) = fxright

  if(.not.domonlimit.and.doposlimit) then
     
     ! eliminate negative values in absence of limiting for monotonicity
     tmpeps = eps2

     do i = 0,N
        if(fx(i).ge.0.d0) then

           !     # check whether flux out of cell i will evacuate that cell
           fluxout = max(0.d0,fx(i)) - min(0.d0,fx(i-1))
           ratiom = rhoq(i)/dt/(fluxout+tmpeps)

        else

           !     # check whether flux out of cell i+1 will evacuate that cell
           fluxout = - min(0.d0,fx(i)) + max(0.d0,fx(i+1))
           ratiom = rhoq(i+1)/dt/(fluxout+tmpeps)

        end if

        !     # limit flux 
        fx(i) = max(0.d0,min(1.d0,ratiom))*fx(i)

     end do

  end if

  if(doposlimit) then
     !     # update solution, eliminating any (small) remaining negative values
     !     ====================================================================
     !     
     do i = 1,N
        rhoq(i) = max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
        q(i) = rhoq(i)/rhoprime_new(i)
     end do

  else

     !     # update solution using total flux
     !     ==================================
     !     
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime_new(i)
     end do

  end if

  rhoprime(1:N) = rhoprime_new(1:N)

end subroutine lwsweep

! ============================================================
subroutine wenosweep(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,doposlimit,scale,lambdavec,limitpts)
  !     
  implicit none

  !     # inputs
  logical doposlimit
  integer, intent(in) :: N
  integer, intent(in) :: npad
  integer, intent(in) :: bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad:N+npad)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: scale

  !     # in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad:N+npad)
  real(kind=8), intent(inout) :: rhoprime(1-npad:N+npad)

  !     # outputs
  real(kind=8), intent(out) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 


  !# local variables
  integer i, imin, imax, bcleft, bcright, rkstep
  real(kind=8) :: fxleft, fxright, tmp
  real(kind=8) :: rhoq0(1:N), rhoprime0(1:N), xflx(0:N), secnd(-1:N+2)
  real(kind=8) :: qface(0:N), wgt1,wgt2,wgt3, beta1,beta2,beta3
  real(kind=8) :: rka(3), rkb(3), rkc(3), q1,q2,q3

  !# local parameters
  real, parameter :: epsweno = 1.d-6
  integer, parameter :: pweno = 2
  !
  ! set up boundary condition variables
  bcleft = bctype(1)
  bcright = bctype(2)
  fxleft = fbc(1) ! fixed flux at left boundary (used if bcleft==3)
  fxright = fbc(2)  ! fixed flux at right boundary (used if bcright==3)

  ! set up parameters for third order SSP RK scheme
  rka(1) = 1.d0
  rka(2) = 0.75d0
  rka(3) = 1.d0/3.d0

  rkb(1:3) = 1.d0 - rka(1:3)

  rkc(1) = 1.d0
  rkc(2) = 0.25d0
  rkc(3) = 2.d0/3.d0
  !
  ! save a copy of the initial rhoq.
  rhoq0(1:N) = rhoq(1:N)
  rhoprime0(1:N) = rhoprime(1:N)
  fx(:) = 0.d0
  !
  !# loop over three SSP Runge-Kutta substeps
  !==========================================
  do rkstep = 1,3

     ! assume periodic bcs for ghost cells
     q(-2:0) = q(N-2:N)
     q(N+1:N+3) = q(1:3)

     secnd(-1:1) = (13./12.)*(q(-2:0) - 2.*q(-1:1) + q(0:2))**2        
     !
     !# compute weno fluxes
     !===============================
     do i = 0,N
        
        secnd(i+2) = (13./12.)*(q(i+1) - 2.*q(i+2) + q(i+3))**2

        if(rhou(i).ge.0.d0) then

!     # compute smoothness indicators for each stencil.
            beta1 = secnd(i-1) &
                 + 0.25d0*(q(i-2) - 4.*q(i-1) + 3.*q(i))**2
            beta2 = secnd(i)   + 0.25d0*(q(i-1) - q(i+1))**2
            beta3 = secnd(i+1) &
                 + 0.25d0*(3.*q(i) - 4.*q(i+1) + q(i+2))**2

!     # compute interpolated value using each stencil
            q1 = ( 2.*q(i-2) - 7.*q(i-1) + 11.*q(i)  )/6.
            q2 = (     -q(i-1) + 5.*q(i)   +  2.*q(i+1))/6.
            q3 = ( 2.*q(i)   + 5.*q(i+1) -       q(i+2))/6.

         else
!     # compute smoothness indicators for each stencil.
            beta1 = secnd(i+2) &
                 + 0.25d0*(3.*q(i+1) - 4.*q(i+2) + q(i+3))**2
            beta2 = secnd(i+1) + 0.25d0*(q(i) - q(i+2))**2
            beta3 = secnd(i)   &
                 + 0.25d0*(q(i-1) - 4.*q(i) + 3.*q(i+1))**2

!     # compute interpolated value using each stencil
            q1 = (11.*q(i+1) - 7.*q(i+2) +  2.*q(i+3))/6.
            q2 = ( 2.*q(i)   + 5.*q(i+1) -       q(i+2))/6.
            q3 = (     -q(i-1) + 5.*q(i)   +  2.*q(i+1))/6.

         end if

!     # modify ideal weights according to smoothness of corresponding stencil
         wgt1 = 0.1d0*(epsweno + beta1)**(-pweno)
         wgt2 = 0.6d0*(epsweno + beta2)**(-pweno)
         wgt3 = 0.3d0*(epsweno + beta3)**(-pweno)
         tmp = (wgt1+wgt2+wgt3)**(-1) 

!     # compute interpolated value using weighted combination of stencils
         qface(i) = tmp*(wgt1*q1 + wgt2*q2 + wgt3*q3)
         xflx(i) = rhou(i)*qface(i)

      end do

      !============== LEFT BOUNDARY ==============
      if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
         ! inflow at left, open boundary -- impose zero scalar gradient there
         xflx(0) = rhou(0)*(qface(1)-fbc(1)) !fbc accounts for mean gradient
      elseif(bcleft.eq.3) then
         ! specify fixed flux at left boundary
         xflx(0) = fxleft
         imin = 1 ! no limiting of left boundary flux
      end if

      !============== RIGHT BOUNDARY ==============
      if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
         ! inflow at right, open boundary -- impose zero scalar gradient there
         xflx(N) = rhou(N)*(qface(N-1)+fbc(2)) !fbc accounts for mean gradient
      elseif(bcright.eq.3) then
         ! specify fixed flux at right boundary by modifying flux correction.
         xflx(N) = fxleft
         imax = N-1 ! no limiting of right boundary flux
      end if

      ! update scalar mass rhoq
      rhoq(1:N) = &
             rka(rkstep)*rhoq0(1:N) &
           + rkb(rkstep)*rhoq(1:N) &
           + rkc(rkstep)*dt*(xflx(0:N-1) - xflx(1:N))

      ! update synthetic density rhoprime
      rhoprime(1:N) = &
             rka(rkstep)*rhoprime0(1:N) &
           + rkb(rkstep)*rhoprime(1:N) &
           + rkc(rkstep)*dt*(rhou(0:N-1) - rhou(1:N))

      ! update scalar mixing ratio q and increment flux fx.
      q(1:N) = rhoq(1:N)/rhoprime(1:N)
      fx(0:N) = fx(0:N) + xflx(0:N)

   end do

  limitpts(:) = 0.d0
  lambdavec(:) = 0.d0


end subroutine wenosweep
!
