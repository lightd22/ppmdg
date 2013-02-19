
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

     if(bcleft.eq.1.and.rhou(0).ge.0.d0) then
        q(-2:0) = q(1) !inflow at open bc
     elseif(bcleft.eq.2) then
        q(-2:0) = q(N-2:N) ! periodic
     elseif((bcleft.eq.3).or.(bcleft.eq.1.and.rhou(0).lt.0.d0)) then
        ! outflow at open bc or fixed flux
        q(0)  = q(1) + (q(1) - q(2))
        q(-1) = q(1) + 2.d0*(q(1) - q(2))
        q(-2) = q(1) + 3.d0*(q(1) - q(2))
     end if

     if(bcright.eq.1.and.rhou(N).lt.0.d0) then
        q(N+1:N+3) = q(N) ! inflow at open bc
     elseif(bcright.eq.2) then
        q(N+1:N+3) = q(1:3) ! periodic
     elseif((bcright.eq.3).or.(bcright.eq.1.and.rhou(0).ge.0.d0)) then
        ! outflow at open bc or fixed flux
        q(N+1) = q(N) + (q(N) - q(N-1))
        q(N+2) = q(N) + 2.d0*(q(N) - q(N-1))
        q(N+3) = q(N) + 3.d0*(q(N) - q(N-1))
     end if

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
