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

!!$     do i = -3,0
!!$ !!$         secnd(i) = (13.d0/12.d0)*(q(i-1) - 2.d0*q(i) + q(i+1))**2
!!$        fd(i) = q(i)-q(i+1)
!!$        fdsq(i) = fd(i)*fd(i)
!!$     end do
     fdbeta(-2:0) = (q(-3:-1)-q(-2:0))**2 + (q(-2:0)-q(-1:1))**2

     monlimit(-2:1) = 0

     do i = -1,N+1
        ! initialize qface, secnd and monlimit to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        !  - secnd(i+2) is an estimate of the second derivative of the scalar 
        !       using scalar values from cells i+1, i+2 and i+3
        !  - monlimit(i+1) indicates whether to limit monotonically at face i+1
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        monlimit(i+1) = 0

!!$        secnd(i+2) = (13.d0/12.d0)*(q(i+1) - 2.d0*q(i+2) + q(i+3))**2
!!$        fd(i+2) = q(i+2)-q(i+3)
!!$        fdsq(i+2) = fd(i+2)*fd(i+2)
        fdbeta(i+2) = (q(i+1)-q(i+2))**2 + (q(i+2)-q(i+3))**2

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
           beta1 = fdbeta(i-1)
           beta2 = fdbeta(i)
           beta3 = fdbeta(i+1)

!!$           beta1 = fdsq(i-2) + 5.*fdsq(i-1) - 4.*fd(i-2)*fd(i-1)
!!$           beta2 = fdsq(i-1) + fdsq(i)
!!$           beta3 = fdsq(i+1) + 5.*fdsq(i) - 4.*fd(i)*fd(i+1)
!!$
!!$           beta1 = secnd(i-1) + 0.25d0*(q(i-2) - 4.d0*q(i-1) + 3.d0*q(i))**2
!!$           beta2 = secnd(i)   + 0.25d0*(q(i-1) - q(i+1))**2
!!$           beta3 = secnd(i+1) + 0.25d0*(3.d0*q(i) - 4.d0*q(i+1) + q(i+2))**2

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
           beta1 = fdbeta(i+2)
           beta2 = fdbeta(i+1)
           beta3 = fdbeta(i)

!!$           beta1 = fdsq(i+2) + 5.*fdsq(i+1) - 4.*fd(i+1)*fd(i+2)
!!$           beta2 = fdsq(i+1) + fdsq(i)
!!$           beta3 = fdsq(i-1) + 5.*fdsq(i) - 4.*fd(i)*fd(i-1)
!!$
!!$           beta1 = secnd(i+2) + 0.25d0*(3.d0*q(i+1) - 4.d0*q(i+2) + q(i+3))**2
!!$           beta2 = secnd(i+1) + 0.25d0*(q(i) - q(i+2))**2
!!$           beta3 = secnd(i)   + 0.25d0*(q(i-1) - 4.d0*q(i) + 3.d0*q(i+1))**2
!!$

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
           monlimit(i-nselpad:i+nselpad) = 1
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


end subroutine ppmsweep_select
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
  real(kind=8) :: fxleft, fxright, tmpeps

  !# local variables defined at cell faces
  real(kind=8) ::  crface, qface(-2:N+2), qtd(-1:N+2), rhoprime_new(-1:N+2)
  real(kind=8) ::  qup(-2:N+2), qcorr(-1:N+1), rhoqtd(-1:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-1:N+1) 
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  integer :: monlimit(-2:N+2), itmp
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
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
  fxup(-2) = rhou(-2)*qup(-2)

  do i = -1,N+1
     !# qface holds interpolant of q at faces.
     qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
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
        qcorr(i) = (1.d0+crface) &
             *(qface(i) - q(i+1) + crface*(qface(i) -2.d0*q(i+1) +qface(i+1)))
     end if
     fxup(i) = rhou(i)*qup(i)
     fxcorr(i) = rhou(i)*qcorr(i)
  end do

  if(rhou(N+2).ge.0.d0) then
     qup(N+2) = q(N+2)
  else
     qup(N+2) = q(N+3)
  end if
  fxup(N+2) = rhou(N+2)*qup(N+2)

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !============== LEFT BOUNDARY ==============
  if(bcleft.eq.1) then
     if(rhou(0).gt.0.d0) then
        ! inflow at left, open boundary -- impose zero scalar gradient there
        fxup(-2:0) = rhou(0)*(qup(1) - fbc(1)) !fbc accounts for mean gradient
        fxcorr(-1:0) = rhou(0)*qcorr(1)
     else
        ! outflow at left, open boundary -- no gradient flux
        fxup(-2:-1) = fxup(0)
        fxcorr(-1) = fxcorr(0)
     end if
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary by modifying flux correction.
     fxcorr(0) = fxleft - fxup(0) 
  end if

  !============== RIGHT BOUNDARY ==============
  if(bcright.eq.1) then
     if(rhou(N).lt.0.d0) then
        ! inflow at right, open boundary -- impose zero scalar gradient there
        !fbc accounts for mean gradient
        fxup(N:N+2) = rhou(N)*(qup(N-1) + fbc(2)) 
        fxcorr(N:N+1) = rhou(N)*qcorr(N-1)
     else
        fxup(N+1:N+2) = fxup(N)
        fxcorr(N+1) = fxcorr(N)
     end if
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fxcorr(N) = fxright - fxup(N)
  end if

  !# update rhoprime -- kept to ensure mass consistency
  rhoprime_new(-1:N+2) = rhoprime(-1:N+2) + dt*(rhou(-2:N+1) - rhou(-1:N+2))

  if(.not.domonlimit) then

     ! define full PPM flux
     fx(0:N) = fxup(0:N) + fxcorr(0:N)

  else

     !!!!!!!!!! MONOTONIC LIMITING EVERYWHERE !!!!!!!!!!!!!!!
     ! DEFINE TRANSPORTED AND DIFFUSED SOLUTION AND ADD ONLY THAT
     ! FRACTION OF CORRECTION FLUX TO THE UPWIND FLUX THAT WILL ALLOW LOCAL
     ! MONOTONICITY TO BE PRESERVED.
     !     
     !# find new scalar mass using only upwind flux
     rhoqtd(-1:N+2) = rhoprime(-1:N+2)*q(-1:N+2) &
          + dt*(fxup(-2:N+1) - fxup(-1:N+2))
     qtd(-1:N+2) = rhoqtd(-1:N+2)/rhoprime_new(-1:N+2)

     tmpeps = eps2 ! a small parameter, used to prevent division by zero

     do i = 0,N
        if(fxcorr(i).ge.0.d0) then

           ! check whether flux out of cell i will induce new minimum in q
           qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxout = max(0.d0,fxcorr(i)) - min(0.d0,fxcorr(i-1))
           ratiom = (rhoqtd(i)-rhoprime_new(i)*qmn)/dt/(fluxout+tmpeps)

           ! check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxin = max(0.d0,fxcorr(i)) - min(0.d0,fxcorr(i+1))
           ratiop = (rhoprime_new(i+1)*qmx-rhoqtd(i+1))/dt/(fluxin+tmpeps)

        else

           ! check whether flux into cell i will induce new maximum in q 
           qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxin = - min(0.d0,fxcorr(i)) + max(0.d0,fxcorr(i-1))
           ratiop = (rhoprime_new(i)*qmx-rhoqtd(i))/dt/(fluxin+tmpeps)

           !check whether flux out of cell i+1 will induce new minimum in q
           qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxout = - min(0.d0,fxcorr(i)) + max(0.d0,fxcorr(i+1))
           ratiom = (rhoqtd(i+1)-rhoprime_new(i+1)*qmn)/dt/(fluxout+tmpeps)

        end if

        !# add limited flux correction 
        fx(i) = fxup(i) + max(0.d0,min(1.d0,ratiom,ratiop))*fxcorr(i)

     end do

  end if

  !# if fixed flux bcs, fix flux
  !=============================
  if(bcleft.eq.3) fx(0) = fxleft
  if(bcright.eq.3) fx(N) = fxright

  if(domonlimit) then
     ! FLUXES HAVE BEEN CORRECTED FOR MONOTONICITY
     if(doposlimit) then

        !# update solution, eliminating any (small) remaining negative values
        !====================================================================
        !     
        do i = 1,N
           rhoq(i) = max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
           q(i) = rhoq(i)/rhoprime_new(i)
        end do

     else

        !# update solution using total flux
        !==================================
        !     
        do i = 1,N
           rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
           q(i) = rhoq(i)/rhoprime_new(i)
        end do

     end if

  else
     ! NO LIMITING FOR MONOTONICITY.

     if(doposlimit) then

        ! NEED TO BE CAREFUL WITH LIMITING FOR NON-NEGATIVITY
        
        ! eliminate negative values in absence of limiting for monotonicity
        tmpeps = eps2

        rhoqtd(-1:N+2) = rhoprime(-1:N+2)*q(-1:N+2) &
             + dt*(rhou(-2:N+1)*qup(-2:N+1) - rhou(-1:N+2)*qup(-1:N+2))

        ! limit flux at left boundary to preserve positivity.
        if(fxcorr(0).ge.0.d0) then
           !# check whether flux out of cell 0 will evacuate that cell
           fluxout = max(0.d0,fxcorr(0)) - min(0.d0,fxcorr(0-1))
           ratiom = rhoqtd(0)/dt/(fluxout+tmpeps)
        else
           !# check whether flux out of cell 1 will evacuate that cell
           fluxout = - min(0.d0,fxcorr(0)) + max(0.d0,fxcorr(0+1))
           ratiom = rhoqtd(1)/dt/(fluxout+tmpeps)
        end if

        !# limit flux 
        fx(0) = fxup(0) + max(0.d0,min(1.d0,ratiom))*fxcorr(0)


        ! limit the rest of the fluxes and update scalar values.
        do i = 1,N
!!$           fxcorr(i+1) = fx(i+1) ! use fxcorr to hold the uncorrected flux
           if(fxcorr(i).ge.0.d0) then

              !# check whether flux out of cell i will evacuate that cell
              fluxout = max(0.d0,fxcorr(i)) - min(0.d0,fxcorr(i-1))
              ratiom = rhoqtd(i)/dt/(fluxout+tmpeps)

           else

              !# check whether flux out of cell i+1 will evacuate that cell
              fluxout = - min(0.d0,fxcorr(i)) + max(0.d0,fxcorr(i+1))
              ratiom = rhoqtd(i+1)/dt/(fluxout+tmpeps)

           end if

           !# limit flux 
           fx(i) = fxup(i) + max(0.d0,min(1.d0,ratiom))*fxcorr(i)

           !# update solution, eliminating any remaining negative values
           !============================================================
           rhoq(i) = max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
           q(i) = rhoq(i)/rhoprime_new(i)
        end do

     else

        !# update solution using total flux
        !==================================
        !     
        do i = 1,N
           rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
           q(i) = rhoq(i)/rhoprime_new(i)
        end do

     end if

  end if

  rhoprime(1:N) = rhoprime_new(1:N)
  if(domonlimit) then
     limitpts(:) = 1.
  else
     limitpts(:) = 0.
  end if
  lambdavec(:) = 0.

end subroutine ppmsweep


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
