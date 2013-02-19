! ============================================================
subroutine ppmsweep_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
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
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2=1.d-16
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
     !     
     !# compute upwind and PPM correction flux using semi-Lagrangian
     !     technique described in Skamarock (2006)
     sdir = sign(1.d0,rhou(i))
     idir = nint(sdir)

     iup = i - min(0,idir) ! first point upwind from face

     rhoufrac = dt*rhou(i)
     fxup(i) = 0.d0
     ii = iup

     do k = 1,nmaxcfl
        if(sdir*rhoufrac.lt.rhoprime(ii)) EXIT
        fxup(i) = fxup(i) + sdir*rhoprime(ii)*q(ii) ! integer cfl part
        rhoufrac = rhoufrac - sdir*rhoprime(ii) ! integer cfl mass flux
        ii = ii - idir ! change index of cell being added to flux
     end do

     ishift(i) = iup - ii                ! size of integer shift at face i
     fxup(i) = fxup(i) + rhoufrac*q(ii)  ! total (integer + frac) upwind flux
     crfrac = rhoufrac/rhoprime(ii)      ! fractional cfl number

     if (rhou(i).eq.0.d0) then
        qup(i) = q(ii)
     else
        qup(i) = fxup(i)/(dt*rhou(i))
     end if

     iu  = ii + min(0,idir)
     iu2 = ii - max(0,idir)
     ! compute PPM correction to upwind flux -- this is only
     !   associated with the fractional cfl.  The upwind flux due to
     !   the integer part of the cfl is exact.
     qcorr(i) = (1.d0-sdir*crfrac)*( qface(iu) - q(ii) &
                                    - sdir*crfrac*(qface(iu) -2.d0*q(ii) +qface(iu2)))
     fxcorr(i) = rhoufrac*qcorr(i)
     tmp_fxcorr(i) = fxcorr(i)
     !     
     !# compute WENO5 smoothness indicators
     !     
     beta1 = (q(ii-idir)-q(ii-idir-1))**2 + (q(ii-idir)-q(ii-idir+1))**2
     beta2 = (q(ii)-q(ii-1))**2 + (q(ii)-q(ii+1))**2
     beta3 = (q(ii+idir)-q(ii+idir-1))**2 + (q(ii+idir)-q(ii+idir+1))**2
     !     
     !# compute ratio of largest to smallest smoothness indicator.
     !     
     lambda(i) = MAX(beta1,beta2,beta3)/(tmpeps + MIN(beta1,beta2,beta3))
     !     
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

!!$  if(MINVAL(qtd(1:N)).lt.0.d0) then
!!$     write(*,*) dt, MINVAL(rhoqtd(1:N)), MINVAL(qtd(1:N))
!!$     write(*,*)
!!$
!!$     do i = 1,N
!!$        write(*,990) i, rhoprime(i), rhoprime_new(i), q(i), qtd(i), &
!!$             rhou(i), fxup(i), rho(i), rhou(i)*dt/rho(i), rhoq(i), fxup(i) - fxup(i-1), &
!!$             rhoq(i), rhoprime(i)*q(i)
!!$     end do
!!$     990 format(i4, 12e12.4)
!!$     STOP 'in pmmsweep_select_semilagrangian'
!!$  end if

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

!!$  if(doposlimit) then
!!$
!!$     if(itmp.gt.0) then 
!!$
!!$        !====LIMIT FOR MONOTONICITY AND POSITIVITY IN SOME PLACES ========
!!$!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of correction flux
!!$        qtd(-1:2) = rhoqtd(-1:2)/rhoprime_new(-1:2) ! scalar after upwind step
!!$
!!$        do i = imin,imax
!!$
!!$!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of correction flux
!!$           qtd(i+2) = rhoqtd(i+2)/rhoprime_new(i+2) ! scalar after upwind step
!!$
!!$           if(monlimit(i).gt.0) then
!!$
!!$              if(tmp_fxcorr(i).ge.0.d0) then
!!$
!!$                 ii = i-ishift(i)
!!$
!!$                 !check whether flux out of cell i will induce new minimum in q
!!$                 qmn = max(0.d0, &
!!$                      min(q(ii-1),qtd(i-1), &
!!$                          q(ii),  qtd(i), &
!!$                          q(ii+1),qtd(i+1)))
!!$                 fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
!!$                 ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)
!!$
!!$                 !check whether flux into cell i+1 will induce new maximum in q
!!$                 qmx = max(q(ii),  qtd(i), &
!!$                           q(ii+1),qtd(i+1), &
!!$                           q(ii+2),qtd(i+2))
!!$                 fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
!!$                 ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)
!!$
!!$              else
!!$
!!$                 ii = i-ishift(i)
!!$
!!$                 !check whether flux into cell i will induce new maximum in q 
!!$                 qmx = max(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
!!$                 fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
!!$                 ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)
!!$
!!$                 !check whether flux out of cell i+1 will induce new min in q
!!$                 qmn = max(0.d0, &
!!$                      min(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2)))
!!$                 fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
!!$                 ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)
!!$
!!$              end if
!!$
!!$              !# limit flux correction 
!!$              fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)
!!$
!!$           end if
!!$
!!$        end do
!!$
!!$     end if
!!$
!!$        !==========================================
!!$        !===LIMIT FOR NON-NEGATIVITY EVERYWHERE ===
!!$
!!$!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of flux correction
!!$        do i = imin,imax
!!$!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of flux correction
!!$           if(tmp_fxcorr(i).ge.0.d0) then
!!$              fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
!!$              ratiom = rhoqtd(i)
!!$           else
!!$              fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
!!$              ratiom = rhoqtd(i+1)
!!$           end if
!!$
!!$           if(fluxout.gt.ratiom) then
!!$              !# limit flux correction for positivity
!!$              ratiom = ratiom/(fluxout + EPSILON(1.d0))
!!$              fxcorr(i) = max(0.d0,min(1.d0,ratiom))*tmp_fxcorr(i)
!!$           end if
!!$
!!$        end do
!!$
!!$     fx(0) = fxup(0) + fxcorr(0)
!!$     do i = 1,N
!!$        !# update solution, eliminating any remaining negative values
!!$        !============================================================
!!$        fx(i) = fxup(i) + fxcorr(i)
!!$        rhoq(i) = max(0.d0,rhoq(i) + fx(i-1) - fx(i))
!!$        q(i) = rhoq(i)/rhoprime_new(i)
!!$        rhoprime(i) = rhoprime_new(i)
!!$     end do
!!$
!!$  else
!!$
!!$     ! NO LIMITING FOR POSITIVITY
!!$
!!$     if(itmp.gt.0) then 
!!$        !======================================================================
!!$        !============ LIMIT ONLY FOR MONOTONICITY, ONLY IN SOME PLACES ========
!!$
!!$        tmp_fxcorr(-1:1) = fxcorr(-1:1) ! unlimited copy of correction flux
!!$        qtd(-1:2) = rhoqtd(-1:2)/rhoprime_new(-1:2) ! scalar after upwind step
!!$
!!$        do i = imin,imax
!!$
!!$           tmp_fxcorr(i+1) = fxcorr(i+1) ! unlimited copy of correction flux
!!$           qtd(i+2) = rhoqtd(i+2)/rhoprime_new(i+2) ! scalar after upwind step
!!$
!!$           if(monlimit(i).gt.0) then
!!$
!!$              if(tmp_fxcorr(i).ge.0.d0) then
!!$
!!$                 ii = i-ishift(i)
!!$
!!$                 !check whether flux out of cell i will induce new minimum in q
!!$                 qmn = min(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
!!$                 fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
!!$                 ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(tmpeps + fluxout)
!!$
!!$                 !check whether flux into cell i+1 will induce new maximum in q
!!$                 qmx = max(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
!!$                 fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
!!$                 ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(tmpeps + fluxin)
!!$
!!$              else
!!$
!!$                 ii = i-ishift(i)
!!$
!!$                 !check whether flux into cell i will induce new maximum in q 
!!$                 qmx = max(q(ii-1),qtd(i-1),q(ii),qtd(i),q(ii+1),qtd(i+1))
!!$                 fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
!!$                 ratiop = rhoprime_new(i)*(qmx-qtd(i))/(tmpeps + fluxin)
!!$
!!$                 !check whether flux out of cell i+1 will induce new min in q
!!$                 qmn = min(q(ii),qtd(i),q(ii+1),qtd(i+1),q(ii+2),qtd(i+2))
!!$                 fluxout = -min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
!!$                 ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(tmpeps + fluxout)
!!$
!!$              end if
!!$
!!$              !# limit flux correction 
!!$              fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)
!!$
!!$           end if ! monlimit(i).gt.0
!!$        end do ! i = imin,imax
!!$     end if ! itmp.gt.0
!!$
!!$     fx(0) = fxup(0) + fxcorr(0)
!!$     do i = 1,N
!!$        !# update solution
!!$        !=================
!!$        fx(i) = fxup(i) + fxcorr(i)
!!$        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
!!$        q(i) = rhoq(i)/rhoprime_new(i)
!!$        rhoprime(i) = rhoprime_new(i)
!!$     end do
!!$  end if

  fx(:) = fx(:)/dt
  limitpts(:) = monlimit(0:N)
  lambdavec(:) = lambda(0:N)

end subroutine ppmsweep_select_semilagrangian
!
! ============================================================
! ============================================================
subroutine pcmsweep_select_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,doposlimit,scale,nselpad,maxlambda,epshybrid, &
     lambda,limitpts)
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
  real(kind=8), intent(out) ::  lambda(0:N) 

  !# local variables
  integer i, imin, imax, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: tmpeps, qmn, qmx, qavg(0:N)
  real(kind=8) :: ql, qr, qbar, a0, a1, a2, a3, cr, dq
  real(kind=8) :: discr, x1, x2, x1dis, x2dis
  real(kind=8) ::  crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: beta1, beta2, beta3
  real(kind=8) :: sdir
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
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) q(1-npad-nmaxcfl:0) = q(1)
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) q(N+1:N+npad+nmaxcfl) = q(N)
  !
  ! INITIALIZE MONLIMIT DIAGNOSTIC
  limitpts(:) = 0
  !     
  !# re-scale epshybrid by the square of scale.  This will prevent 
  !#   epshybrid from overwhelming the smoothness indicators when 
  !#   q is small.
  tmpeps = epshybrid*scale*scale

  if(doposlimit) then
     
     !# apply order reduction selectively -- only at faces where the ratio of
     !     the largest to smallest smoothness indicators exceeds a threshold
     do i = -1-nmaxcfl,nmaxcfl
        qface(i) = MAX(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do

     do i = 0,N
        ! initialize qface and fdbeta to the right of the current face.
        !  - qface(i+1) holds an interpolated scalar value at face i+1
        qface(i+nmaxcfl) = MAX(0.d0,fac1*(q(i+nmaxcfl) + q(i+nmaxcfl+1)) &
                                          + fac2*(q(i+nmaxcfl-1) + q(i+nmaxcfl+2)))

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

        ! find location of extrema and test whether value is negative there
        poscorr = .false.
        discr = a2**2 - 3.d0*a1*a3
        if(discr.ge.0.d0) then
           ! two extrema in general -- find position
           x1 = (-a2+sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           x2 = (-a2-sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           ! ...and distance from center of [0,1] interval
           x1dis = ABS(x1-0.5d0)
           x2dis = ABS(x2-0.5d0)

           if (((x1dis.lt.0.5d0).AND.(a0+x1*(a1+x1*(a2+x1*a3)).lt.0.d0)) .OR. &
                ((x2dis.lt.0.5d0).AND.(a0+x2*(a1+x2*(a2+x2*a3)).lt.0.d0))) then
              ! if one of the extrema is within the interval and is negative, 
              !   then use piecewise constant reconstruction
              a0 = qbar
              a1 = 0.d0
              a2 = 0.d0
              a3 = 0.d0
              poscorr = .true.
           end if
        end if

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if((lambda(i).gt.maxlambda).and.(.NOT.poscorr)) then !.OR.(i.lt.imin).OR.(i.gt.imax)) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           qmx = max(q(ii-1),q(ii),q(ii+1))
           qmn = min(q(ii-1),q(ii),q(ii+1))

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

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = fx(i) + rhoufrac*qavg(i)
     end do

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

        ! if lambda (ratio of largest to smallest smoothness indicators)
        !   exceeds maxlambda (a threshold value), apply order reduction
        !   to the cubic polynomial in this cell.
        if((lambda(i).gt.maxlambda)) then !.OR.(i.lt.imin).OR.(i.gt.imax)) then
           limitpts(i) = 1

           ! monotonize the cell edge values
           qmx = max(q(ii-1),q(ii),q(ii+1))
           qmn = min(q(ii-1),q(ii),q(ii+1))

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
        fx(i) = fx(i) + rhoufrac*qavg(i)
     end do

  end if

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     !   fbc accounts for mean gradient
     fx(0) = dt*rhou(0)*qavg(1)! - fbc(1))
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

  if(doposlimit) then
     !# update solution, eliminating any remaining negative values
     !============================================================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i) !MAX(0.d0,
!!$        if(rhoq(i).lt.0.d0) then
!!$           if(rhoq(i)/rhoprime(i).gt.-1.d-10) then
!!$              rhoq(i) = 0.d0
!!$           else
!!$              write(*,*) rhoq(i), rhoprime(i), rhoq(i)/rhoprime(i)
!!$              write(*,*)
!!$              do ii = 1,i
!!$                 write(*,990) ii, rhoq(ii), rhoprime(ii), rhoq(ii)/rhoprime(ii), qavg(ii), fx(ii), rhou(ii)*dt/rho(ii), fx(ii-1)-fx(ii)
!!$              end do
!!$990           format(i4, 12e12.4)
!!$              STOP 'in pcmsweep_select'
!!$           end if
!!$        end if
        q(i) = rhoq(i)/rhoprime(i)
     end do

  else

     !# update solution
     !=================
     do i = 1,N
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1) - rhou(i))
        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        q(i) = rhoq(i)/rhoprime(i)
     end do
  end if ! doposlimit

  fx(:) = fx(:)/dt
end subroutine pcmsweep_select_semilagrangian
!
! ============================================================
subroutine ppmsweep_semilagrangian(rhoq,q,rho,rhoprime,rhou,fx,dt, &
     N,npad,nmaxcfl,bctype,fbc,domonlimit,doposlimit,lambdavec,limitpts)
  !     
  implicit none

  !# inputs
  logical, intent(in) :: domonlimit, doposlimit
  integer, intent(in) :: N, npad, nmaxcfl, bctype(2)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: fbc(2)
  real(kind=8), intent(in) :: dt
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)


  !# in/outputs
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhoprime(1-npad-nmaxcfl:N+npad+nmaxcfl)

  !# outputs
  real(kind=8) ::  fx(0:N) 
  real(kind=8), intent(out) ::  limitpts(0:N) 
  real(kind=8), intent(out) ::  lambdavec(0:N) 

  !# local variables
  integer i, k, nshift, bcleft, bcright
  real(kind=8) :: fxleft, fxright, tmpeps

  !# local variables defined at cell faces
  real(kind=8) ::  crface, crfrac, rhoufrac
  real(kind=8) ::  qface(-npad-nmaxcfl:N+npad+nmaxcfl), qtd(-1:N+2)
  real(kind=8) :: rhoprime_new(-1:N+2), tmp_fxcorr(-2:N+2)
  real(kind=8) ::  fxup(-2:N+2), fxcorr(-2:N+2), rhoqtd(-1:N+2)
  real(kind=8) :: fluxin, fluxout, ratiop, ratiom, qmn, qmx
  real(kind=8) :: sdir
  integer :: idir, iup, idn, ii, iu, iu2
  !# local parameters
  real(kind=8), parameter :: fac1=7.d0/12.d0, fac2=-1.d0/12.d0, eps2&
       &=1.d-16
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
  do i = -2-nmaxcfl,-3+nmaxcfl
     qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
  end do

  do i = -2,N+2
     !# qface holds interpolant of q at faces.
     qface(i+nmaxcfl) = fac1*(q(i+nmaxcfl) + q(i+1+nmaxcfl)) &
                          + fac2*(q(i-1+nmaxcfl) + q(i+2+nmaxcfl))

     ! find cfl number using this new velocity
     crface = 2.d0*dt*rhou(i)/(rho(i)+rho(i+1))

     sdir = sign(1.d0,rhou(i))
     idir = nint(sdir)

     iup = i - min(0,idir) ! first point upwind from face
     idn = i + max(0,idir) ! first point downwind from face

     nshift = floor(abs(crface)) ! integer part of cfl
     
     rhoufrac = dt*rhou(i)
     fxup(i) = 0.d0
     ii = iup
     do k = 1,nshift
        fxup(i) = fxup(i) + sdir*rho(ii)*q(ii) ! integer cfl part
        rhoufrac = rhoufrac - sdir*rho(ii) ! integer cfl mass flux
        ii = ii - idir ! change index of cell being added to flux
     end do
     fxup(i) = fxup(i) + rhoufrac*q(ii)
     crfrac = rhoufrac/rho(ii)
     
     iu  = ii + min(0,idir)
     iu2 = ii - max(0,idir)
     ! compute PPM correction to upwind flux -- this is only
     !   associated with the fractional cfl.  The upwind flux due to
     !   the integer part of the cfl is exact.
     fxcorr(i) = rhoufrac*(1.d0-sdir*crfrac) & 
          *( qface(iu) - q(ii) &
          - sdir*crfrac*(qface(iu) -2.d0*q(ii) +qface(iu2)))

  end do

  !===============================================================
  !======== DEAL WITH BOUNDARY CONDITIONS, IF NECESSARY ==========

  !============== LEFT BOUNDARY ==============
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     ! inflow at left, open boundary -- impose zero scalar gradient there
     crface = 2.d0*dt*rhou(0)/(rho(0)+rho(1))
     ! compute integer part of cfl and mass flux associated with
     !   fractional cfl
     nshift = floor(crface) 

     rhoufrac = dt*rhou(0)
     fxup(0) = 0.d0
     do k = 1,nshift
        fxup(0) = fxup(0) + rho(2-k)*q(2-k) ! integer cfl part
        rhoufrac = rhoufrac - rho(2-k) ! integer cfl part of mass flux
     end do
     fxup(0) = fxup(0) + rhoufrac*q(1-nshift)
     crfrac = rhoufrac/rho(1-nshift)
     
     ! compute PPM correction to upwind flux -- this is only
     !   associated with the fractional cfl.  The upwind flux due to
     !   the integer part of the cfl is exact.
     fxcorr(0) = rhoufrac*(1.d0-crfrac) & 
          *( qface(1-nshift) - q(1-nshift) &
          - crfrac*(qface(-nshift) -2.d0*q(1-nshift) +qface(1-nshift)))
     
     ! account for affect of mean gradient of scalar on inflow flux
     fxup(0) = fxup(0) - rhoufrac*fbc(1)
  elseif(bcleft.eq.3) then
     ! specify fixed flux at left boundary by modifying flux correction.
     fxcorr(0) = dt*fxleft - fxup(0) 
  end if

  !============== RIGHT BOUNDARY ==============
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     ! inflow at right, open boundary -- impose zero scalar gradient there
     crface = 2.d0*dt*rhou(N)/(rho(N)+rho(N+1))
     ! compute integer part of cfl and mass flux associated with
     !   fractional cfl
     nshift = floor(abs(crface))

     rhoufrac = dt*rhou(N)
     fxup(N) = 0.d0
     do k = 1,nshift
        fxup(N) = fxup(N) - rho(N-1+k)*q(N-1+k) ! integer cfl part
        rhoufrac = rhoufrac + rho(N-1+k) ! integer cfl mass flux
     end do
     fxup(N) = fxup(N) + rhoufrac*q(N+nshift) ! fractional cfl contribution
     crfrac = rhoufrac/rho(N+nshift)

     ! compute PPM correction to upwind flux -- this is only
     !   associated with the fractional cfl.  The upwind flux due to
     !   the integer part of the cfl is exact.
     fxcorr(N) = rhoufrac*(1.d0+crfrac)*(qface(N-1+nshift) - q(N+nshift) &
          +crfrac*(qface(N-1+nshift)-2.d0*q(N+nshift) +qface(N+nshift)))

     ! account for affect of mean gradient of scalar on inflow flux
     fxup(N) = fxup(N) + rhoufrac*fbc(2)
  elseif(bcright.eq.3) then
     ! specify fixed flux at right boundary by modifying flux correction.
     fxcorr(N) = dt*fxright - fxup(N)
  end if

  !# update rhoprime -- kept to ensure mass consistency
  rhoprime_new(-1:N+2) = rhoprime(-1:N+2) + dt*(rhou(-2:N+1) - rhou(-1:N+2))
  rhoqtd(-1:N+2) = rhoprime(-1:N+2)*q(-1:N+2) + fxup(-2:N+1) - fxup(-1:N+2)

  tmpeps = eps2*dt ! a small parameter, used to prevent division by zero

  if(domonlimit) then

     !!!!!!!!!! MONOTONIC LIMITING EVERYWHERE !!!!!!!!!!!!!!!
     !# find new scalar mass using only upwind flux 
     !   NOTE: fxup and fxcorr have already been multiplied by dt.
     qtd(-1:N+2) = (rhoprime(-1:N+2)*q(-1:N+2) &
                     + fxup(-2:N+1)-fxup(-1:N+2))/rhoprime_new(-1:N+2)

     tmp_fxcorr(-1:1) = fxcorr(-1:1)

     do i = 0,N
        tmp_fxcorr(i+1) = fxcorr(i+1)

        if(tmp_fxcorr(i).ge.0.d0) then

           ! check whether flux out of cell i will induce new minimum in q
           qmn = min(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
           ratiom = rhoprime_new(i)*(qtd(i)-qmn)/(fluxout+tmpeps)

           ! check whether flux into cell i+1 will induce new maximum in q
           qmx = max(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxin = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i+1))
           ratiop = rhoprime_new(i+1)*(qmx-qtd(i+1))/(fluxin+tmpeps)

        else

           ! check whether flux into cell i will induce new maximum in q 
           qmx = max(q(i-1),qtd(i-1),q(i),qtd(i),q(i+1),qtd(i+1))
           fluxin = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i-1))
           ratiop = rhoprime_new(i)*(qmx-qtd(i))/(fluxin+tmpeps)

           !check whether flux out of cell i+1 will induce new minimum in q
           qmn = min(q(i),qtd(i),q(i+1),qtd(i+1),q(i+2),qtd(i+2))
           fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
           ratiom = rhoprime_new(i+1)*(qtd(i+1)-qmn)/(fluxout+tmpeps)

        end if

        !# add only that portion of the flux correction that will
        ! guarantee monotonicity to be preserved locally.
        fxcorr(i) = max(0.d0,min(1.d0,ratiom,ratiop))*tmp_fxcorr(i)

     end do

  end if

  !# if fixed flux bcs, fix flux
  !=============================
  if(bcleft.eq.3) fx(0) = fxleft*dt
  if(bcright.eq.3) fx(N) = fxright*dt
  
  rhoprime(-1:N+2) = rhoprime_new(-1:N+2)

  if((.not.domonlimit).and.(doposlimit)) then

     ! compute transported/diffused value
     qtd(-1:N+2) = &
          (rhoprime(-1:N+2)*q(-1:N+2) + fxup(-2:N+1) - fxup(-1:N+2)) &
          /rhoprime_new(-1:N+2)

     ! check whether limiting for positivity is necessary
     if(2.d0*MAXVAL(ABS(fxcorr)).gt. &
          MINVAL(rhoprime_new(-1:N+2)*qtd(-1:N+2))) then

        !==========================================
        !===LIMIT FOR NON-NEGATIVITY EVERYWHERE ===

        tmp_fxcorr(-1:1) = fxcorr(-1:1)
        do i = 0,N
           tmp_fxcorr(i+1) = fxcorr(i+1) ! copy fxcorr into new array
                                         ! so that unlimited copy is around

           if(tmp_fxcorr(i).ge.0.d0) then
              fluxout = max(0.d0,tmp_fxcorr(i)) - min(0.d0,tmp_fxcorr(i-1))
              ratiom = rhoprime_new(i)*qtd(i)/(tmpeps + fluxout)
           else
              fluxout = - min(0.d0,tmp_fxcorr(i)) + max(0.d0,tmp_fxcorr(i+1))
              ratiom = rhoprime_new(i+1)*qtd(i+1)/(tmpeps + fluxout)
           end if

           !# limit flux correction for positivity if necessary (if ratiom<1)
           fxcorr(i) = max(0.d0,min(1.d0,ratiom))*tmp_fxcorr(i)
        end do

     end if

  end if

  if(doposlimit) then
     fx(0) = fxup(0) + fxcorr(0)
     do i = 1,N
        !# update solution, eliminating any remaining negative values
        !============================================================
        fx(i) = fxup(i) + fxcorr(i)
!!$        rhoq(i) = rhoq(i) + fx(i-1) - fx(i)
        rhoq(i) = max(0.d0,rhoq(i) + fx(i-1) - fx(i))
        q(i) = rhoq(i)/rhoprime_new(i)
        rhoprime(i) = rhoprime_new(i)
     end do
  else
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
  if(domonlimit) then
     limitpts(:) = 1.
  else
     limitpts(:) = 0.
  end if
  lambdavec(:) = 0.

end subroutine ppmsweep_semilagrangian
