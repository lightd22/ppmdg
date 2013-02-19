! ============================================================
subroutine ppmsweep_pmod_select(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
  integer i, imin, imax, bcleft, bcright, iup, ilft, irgt, ibar
  real(kind=8) :: fxleft, fxright, tmpeps, fluxout, rhoqtd, tmp_fxcorr
  real(kind=8) :: limit(-1:N+1), lambda(-1:N+1)
  real(kind=8) :: fxup(-1:N+1), fxcorr(-1:N+1)
  real(kind=8) :: qup(-1:N+1), qcorr(-1:N+1)
  logical doselimit, domonlimit

  !# local variables defined at cell faces
  real(kind=8) :: cr, ql, qr, qbar, qbarl, qbarr, qavg(0:N), a0, a1, a2
  real(kind=8) :: qmn, qmx, qface(-2:N+2), x1, fdbeta(-2:N+3)

  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0
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
  !# compute ppm estimate for flux
  !===============================
  !     
  !# apply monotonic limiting selectively -- only at faces at or
  !     adjacent to those which are found to be less than smooth using the
  !     ratio of the largest to smallest WENO smoothness ratios as an
  !     indicator. 
  
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  if(doposlimit) then
     
     !
     ! INITIALIZE MONLIMIT DIAGNOSTIC
     limit(:) = 0

     do i = -2,-1
!!$     do i = -1,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do
     do i = -2,0
        fdbeta(i) = (q(i-1)-q(i))**2 + (q(i)-q(i+1))**2
     end do

     do i = -1,N+1
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = rhou(i)*dt/rhoprime(i)
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
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambda(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
                /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
        end if
        !
        ! define coefficients of parabola
        a0 = ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !  exceeds maxlambda (a threshold value), modify the parabola to
        !  ensure that the construction does not introduce a new extrema
        !  within this grid cell.
        if(lambda(i).gt.maxlambda) then
           limit(i) = 1

           if(rhou(i).ge.0.d0) then
              qbarl = q(i+1)
              qbarr = q(i-1)
           else
              qbarl = q(i)
              qbarr = q(i+2)
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
              !
              ! re-define coefficients of parabola (in case ql/qr have changed)
              a0 = ql
              a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
              a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
              !
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

        fxup(i) = rhou(i)*qbar
        fxcorr(i) = rhou(i)*qcorr(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

     imin = 0
     imax = N

     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(1) to account for possible mean gradient of scalar
        fxup(0) = rhou(0)*MAX(0.d0,qup(1) - fbc(1))
        fxcorr(0) = rhou(0)*qcorr(1)
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fxcorr(0) = fxleft - fxup(0) 
        imin = 1 ! no limiting of left boundary flux
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right boundary, impose no scalar gradient across cell next to 
        !  boundary.  Use fbc(2) to account for possible mean gradient of scalar
        fxup(N) = rhou(N)*(qup(N-1) + fbc(2))
        fxcorr(N) = rhou(N)*qcorr(N-1)
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fxcorr(N) = fxright - fxup(N)
        imax = N-1 ! no limiting of right boundary flux
     end if

     ! limit correction fluxes at faces to prevent the occurrence of 
     !   negative scalar values
     tmp_fxcorr = fxcorr(-1)
     do i = 0,N
        if(fxcorr(i).ge.0.d0) then
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i) + dt*(fxup(i-1) - fxup(i))
           fluxout = max(0.d0,fxcorr(i)) - min(0.d0,tmp_fxcorr)
           tmp_fxcorr = fxcorr(i)
        else
           ! ensure correction flux out of cell i+1 will not empty cell
           rhoqtd = rhoq(i+1) + dt*(fxup(i) - fxup(i+1))
           fluxout = max(0.d0,fxcorr(i+1)) - min(0.d0,fxcorr(i))
           tmp_fxcorr = fxcorr(i)
        end if
        if(dt*fluxout.gt.rhoqtd) then
           fxcorr(i) = fxcorr(i) &
                *MAX(0.d0,MIN(1.d0,rhoqtd/(dt*fluxout+EPSILON(rhoqtd))))
        end if
     end do

     !# update solution, eliminating any remaining negative values
     !============================================================
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        fx(i) = fxup(i) + fxcorr(i)
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i)) !MAX(0.d0,
        q(i) = rhoq(i)/rhoprime(i)
     end do

     lambdavec(0:N) = lambda(0:N)
     limitpts(0:N) = limit(0:N)

  else

     !
     ! INITIALIZE MONLIMIT DIAGNOSTIC
     limitpts(:) = 0

     do i = -1,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do
     fdbeta(-1:1) = (q(-2:0)-q(-1:1))**2 + (q(-1:1)-q(0:2))**2

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

        if(rhou(i).ge.0.d0) then
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
           qbarl = q(i+1)
           qbarr = q(i-1)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambdavec(i) = MAX(fdbeta(i-1),fdbeta(i),fdbeta(i+1)) &
                /(tmpeps+MIN(fdbeta(i-1),fdbeta(i),fdbeta(i+1)))

        else
           !     
           !# compute upwind and PCM correction flux
           !     
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           qbarl = q(i)
           qbarr = q(i+2)
           !     
           !# compute ratio of largest to smallest smoothness indicator
           lambdavec(i) = MAX(fdbeta(i+2),fdbeta(i+1),fdbeta(i)) &
                /(tmpeps+MIN(fdbeta(i+2),fdbeta(i+1),fdbeta(i)))
        end if
        !
        ! define coefficients of parabola
        a0 = ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !  exceeds maxlambda (a threshold value), modify the parabola to
        !  ensure that the construction does not introduce a new extrema
        !  within this grid cell.
        if(lambdavec(i).gt.maxlambda) then
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
           else
              !
              ! re-define coefficients of parabola (in case ql/qr have changed)
              a0 = ql
              a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
              a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
              !
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
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = rhou(i)*qavg(i)
     end do

     !===============================================================
     !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
     
     !============== LEFT BOUNDARY ==============
     if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        !   fbc accounts for mean gradient
        fx(0) = rhou(0)*(qavg(1) - fbc(1))
     elseif(bcleft.eq.3) then
        ! specify fixed flux at left boundary by modifying flux correction.
        fx(0) = fxleft 
     end if

     !============== RIGHT BOUNDARY ==============
     if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
        ! inflow at right, open boundary -- impose zero scalar gradient
        !   fbc accounts for mean gradient
        fx(N) = rhou(N)*(qavg(N-1) + fbc(2))
     elseif(bcright.eq.3) then
        ! specify fixed flux at right boundary by modifying flux correction.
        fx(N) = fxright 
     end if

     !# update solution
     !=================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do

  end if ! doposlimit

end subroutine ppmsweep_pmod_select
!
! ============================================================
! ============================================================
! ============================================================
subroutine ppmsweep_pmod(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
  real(kind=8) :: cr, ql, qr, qbar, qavg(0:N), a0, a1, a2
  real(kind=8) :: qface(-1:N+1), x1
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0
  !
  if(domonlimit) then
     call ppmsweep_pmod_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  limitpts(:) = 0
  lambdavec(:) = 0.d0
  !     
  !# compute ppm estimate for flux
  !===============================
  if(doposlimit) then

     do i = -1,0
        qface(i) = MAX(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = MAX(0.d0,fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3)))

        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
        end if
        !
        ! define coefficients of parabola
        a0 = ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

        ! find location of extrema and test whether value is negative there
        ! one extrema, located at
        x1 = -a1/(2.d0*a2 + EPSILON(1.d0))
        if((ABS(x1-0.5d0).lt.0.5d0).AND.(a0+x1*(a1+x1*a2).lt.0.d0)) then
           ! the extrema is within the cell and negative, 
           !   --> use piecewise constant reconstruction
           a0 = qbar
           a1 = 0.d0
           a2 = 0.d0
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = rhou(i)*qavg(i)
     end do

  else

     do i = -1,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))

        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
        end if
        !
        ! define coefficients of parabola
        a0 = ql
        a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
        a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = rhou(i)*qavg(i)
     end do

  end if

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
  
  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     !   fbc accounts for mean gradient
     fx(0) = rhou(0)*(qavg(1) - fbc(1))
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary by modifying flux correction.
     fx(0) = fxleft 
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     ! inflow at right, open boundary -- impose zero scalar gradient
     !   fbc accounts for mean gradient
     fx(N) = rhou(N)*(qavg(N-1) + fbc(2))
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fx(N) = fxright 
  end if

  if(doposlimit) then
     !# update solution, eliminating any remaining negative values
     !============================================================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i)) !MAX(0.d0,
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do

  else

     !# update solution
     !=================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do
  end if ! doposlimit
end subroutine ppmsweep_pmod
!
! ============================================================
subroutine ppmsweep_pmod_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
  integer i, ii, bcleft, bcright, iup, imin, imax
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: cr, ql, qr, qbar, qavg(0:N), a0, a1, a2
  real(kind=8) :: qmn, qmx, qface(-1:N+1), x1
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0
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
  limitpts(:) = 1
  lambdavec(:) = 0.d0
  !     
  !# compute ppm estimate for flux
  !===============================
  !     
  !# apply polynomial modification everywhere so that the polynomial 
  !    reconstruction is monotonic.  Additional modification to
  !    ensure a positive reconstuction is added when doposlimit=.true.
  if(doposlimit) then

     do i = -1,0
        qface(i) = MAX(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = MAX(0.d0,fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3)))

        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
           iup = i
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           iup = i+1
        end if
        ! monotonize the cell edge values
        qmx = max(q(iup-1),q(iup),q(iup+1))
        qmn = min(q(iup-1),q(iup),q(iup+1))

        ql = MAX(0.d0,qmn,MIN(qmx,ql))
        qr = MAX(0.d0,qmn,MIN(qmx,qr))

        ! modify the parabola to ensure that the construction does not 
        !  introduce a new extrema  within this grid cell.

        ! if qbar outside of [ql,qr], use piecewise constant,
        !   otherwise use parabola
        if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
           a0 = qbar ! piecewise constant reconstruction
           a1 = 0.d0
           a2 = 0.d0
        else
           !
           ! re-define coefficients of parabola
           a0 = ql
           a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
           a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

           ! check to see if the parabola has an extrema outside ql or qr.
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

        ! find location of extrema and test whether value is negative there
        ! one extrema, located at
        x1 = -a1/(2.d0*a2 + EPSILON(1.d0))
        if((ABS(x1-0.5d0).lt.0.5d0).AND.(a0+x1*(a1+x1*a2).lt.0.d0)) then
           ! the extrema is within the cell and negative, 
           !   --> use piecewise constant reconstruction
           a0 = qbar
           a1 = 0.d0
           a2 = 0.d0
        end if

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = rhou(i)*qavg(i)
     end do

  else

     do i = -1,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do

     do i = 0,N
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))

        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
           iup = i
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
           iup = i+1
        end if
        ! Modify the parabola to ensure that the construction does
        !   not introduce a new extrema within this grid cell.

        ! monotonize the cell edge values
        qmx = max(q(iup-1),q(iup),q(iup+1))
        qmn = min(q(iup-1),q(iup),q(iup+1))

        ql = MAX(qmn,MIN(qmx,ql))
        qr = MAX(qmn,MIN(qmx,qr))

        ! if qbar outside of [ql,qr], use piecewise constant,
        !   otherwise use parabola
        if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
           a0 = qbar ! piecewise constant reconstruction
           a1 = 0.d0
           a2 = 0.d0
        else
           !
           ! re-define coefficients of parabola
           a0 = ql
           a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
           a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar

           ! check to see if the parabola has an extrema outside ql or qr.
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

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*a2/3.d0)
        fx(i) = rhou(i)*qavg(i)
     end do

  end if

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========
  
  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     !   fbc accounts for mean gradient
     fx(0) = rhou(0)*(qavg(1) - fbc(1))
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary by modifying flux correction.
     fx(0) = fxleft 
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     ! inflow at right, open boundary -- impose zero scalar gradient
     !   fbc accounts for mean gradient
     fx(N) = rhou(N)*(qavg(N-1) + fbc(2))
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fx(N) = fxright 
  end if

  if(doposlimit) then
     !# update solution, eliminating any remaining negative values
     !============================================================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i)) !MAX(0.d0,
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do

  else

     !# update solution
     !=================
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do
  end if ! doposlimit

end subroutine ppmsweep_pmod_monotonic
