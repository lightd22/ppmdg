! ============================================================
subroutine ppmsweep_fct_select(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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


end subroutine ppmsweep_fct_select
!
! ============================================================
! ============================================================
! ============================================================
subroutine ppmsweep_fct(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
     call ppmsweep_fct_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
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
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
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
           cr = rhou(i)*dt/rhoprime(i)
           ql = qface(i)
           qbar = q(i)
           qr = qface(i-1)
        else
           !# compute upwind and PPM correction flux
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
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

end subroutine ppmsweep_fct


!
! ============================================================
subroutine ppmsweep_fct_monotonic(rhoq,q,rhou,rho,rhoprime,fx,dt, &
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
        cr = rhou(i)*dt/rhoprime(i)
        ql = qface(i)
        qbar = q(i)
        qr = qface(i-1)
     else
        !# compute upwind and PPM correction flux
        cr = ABS(rhou(i)*dt/rhoprime(i+1))
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

end subroutine ppmsweep_fct_monotonic
