! ============================================================
subroutine pcmsweep(rhoq,q,rhou,rho,rhoprime,fx,dt, &
     N,npad,bctype,fbc,domonlimit,doposlimit,lambdavec,limitpts)
  !     
  use pcm
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
  integer i, ii, bcleft, bcright
  real(kind=8) :: fxleft, fxright

  !# local variables defined at cell faces
  real(kind=8) :: qface(-3:N+3), cr, dq, p(4), qavg(0:N), qq(6)
  real(kind=8) :: ql, qr, qbar, a0, a1, a2, a3
  real(kind=8) :: discr, x1, x2, x1dis, x2dis
  integer :: monlimit(-3:N+3)

  !# local parameters
  real(kind=8), parameter :: &
       fac1=7.d0/12.d0,  fac2=-1.d0/12.d0, &
       fac3=34.d0/48.d0, fac4=-5.d0/48.d0
  !
  monlimit(:) = 0
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
  ! CHOOSE VARIETY OF PCM (MONONTONIC OR NON, POSITIVE OR NON)
  if((.not.domonlimit).and.(.not.doposlimit)) then
     !=========== plain PCM, no positivity of monotonicity filters =========
     do i = -2,0
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
     end do
     !     
     !# compute pcm flux across domain
     !================================
     do i = 0,N
        !# qface holds interpolant of q at faces.
        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
        !     
        !# compute PCM reconstruction of scalar field in upwind cell
        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
!!$           p = pcm_coef(qface(i),qface(i-1),q(i),dq,doposlimit)
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
!!$           p = pcm_coef(qface(i),qface(i+1),q(i+1),dq,doposlimit) 
        end if

        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = rhou(i)*qavg(i)
     end do

  elseif(doposlimit.and.(.not.domonlimit)) then
     !=========== PCM with positivity filter (not monotonicity) =========
     do i = -2,0
        qface(i) = max(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
     end do
     !     
     !# compute pcm flux across domain
     !================================
     do i = 0,N
        !# qface holds interpolant of q at faces.
        qface(i+1) = max(0.d0,fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3)))
        !     
        !# compute PCM reconstruction of scalar field in upwind cell
        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
           ql = qface(i) ! note: ql and qr flipped here.
           qr = qface(i-1)
           qbar = q(i)
!!$           p = pcm_coef(qface(i),qface(i-1),q(i),dq,doposlimit)
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
!!$           p = pcm_coef(qface(i),qface(i+1),q(i+1),dq,doposlimit) 
        end if

        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        discr = a2**2 - 3.d0*a1*a3

        if(discr.ge.0.d0) then
           ! two extrema in general
           x1 = (-a2+sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           x2 = (-a2-sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           if (((abs(x1-0.5d0).lt.0.5d0).AND. &
                (a0+x1*(a1+x1*(a2+x1*a3)).lt.0.d0)) .OR. &
                ((abs(x2-0.5d0).lt.0.5d0).AND. &
                (a0+x2*(a1+x2*(a2+x2*a3)).lt.0.d0))) then
              ! an extrema is both in the interval and negative
              ! --> USE PIECEWISE CONSTANT RECONSTRUCTION
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

  elseif(domonlimit.and.(.not.doposlimit)) then

     do i = -3,N+3
!!$        qface(i) = pcm_face_value_monotonic(q(i-2:i+3),doposlimit)
        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
        if ((qface(i)-q(i))*(q(i+1)-qface(i)).lt.0.d0) then 
           ! possible grid scale monotonicity violation
           !   -- make further check to confirm monotonicity violation
           if (((q(i) - q(i-1))*(q(i+2) - q(i+1)) .ge. 0.d0) .OR. &
                ((q(i) - q(i-1))*(q(i-1) - q(i-2)) .le. 0.d0) .OR. &
                ((q(i+2) - q(i+1))*(q(i+3) - q(i+2)) .le. 0.d0) .OR. &
                ((qface(i)   - q(i))*(q(i) - q(i-1)) .le. 0.d0)) then
              monlimit(i) = 1
              ! use either left or right cell avg value 
              !   -- whichever is closer to the current estimate
              if (abs(qface(i) - q(i+1)) .ge. abs(qface(i) - q(i))) then
                 qface(i) = q(i)
              else
                 qface(i) = q(i+1)
              end if
           end if
        end if
     end do
     !     
     !# compute pcm flux across domain
     !================================
     do i = 0,N
!!$        !# qface holds interpolant of q at faces.
!!$        qface(i+3) = pcm_face_value_monotonic(q(i+1:i+7),doposlimit)
        !     
        !# compute PCM reconstruction of scalar field in upwind cell
        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
           qq = qface(i+2:i-3:-1)
           ql = qface(i)
           qr = qface(i-1)
           qbar = q(i)
!!$           p = pcm_coef_monotonic(qq,q(i),dq,doposlimit) ! ql and qr flipped
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
           qq = qface(i-2:i+3)
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
!!$           p = pcm_coef_monotonic(qq,q(i+1),dq,doposlimit) 
        end if

        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        discr = a2**2 - 3.d0*a1*a3

        if(discr.ge.0.d0) then
           ! two extrema in general
           x1 = (-a2+sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
           x2 = (-a2-sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))

           ! check if both extrema are within interval
           if ((abs(x1-0.5d0).lt.0.5d0).AND.(abs(x2-0.5d0).lt.0.5d0)) then
              monlimit(i) = monlimit(i) + 2 
              ! bad extrema, zero out cubic term
              a3 = 0.d0 
              ! reduce order, if qbar outside of [ql,qr], use
              ! piecewise constant, otherwise use parabola
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
           elseif ((abs(x1-0.5d0).lt.0.5d0).OR.(abs(x2-0.5d0).lt.0.5d0)) then
              ! exactly one extrema within the interval
              ! make a check for subgrid monotonicity violations
              if((qq(3)-qq(2))*(qq(5)-qq(4)).ge.0.d0 & !bad extrema
                   .OR.(qq(2)-qq(1))*(qq(3)-qq(2)).le.0.d0 &
                   .OR.(qq(5)-qq(4))*(qq(6)-qq(5)).le.0.d0 &
                   .OR.(qq(3)-qq(2))*a1.le.0.d0) then 
                 monlimit(i) = monlimit(i) + 2
                 ! bad extrema, zero out cubic term
                 a3 = 0.d0 
                 ! reduce order, if qbar outside of [ql,qr], use
                 ! piecewise constant, otherwise use parabola
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

           end if ! checks for extrema within interval

        end if ! check for discr of first derivative polynomial > 0
        ! ========== END MONOTONICITY FILTER =====================

         !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = rhou(i)*qavg(i)

     end do

  else ! Both monotonicity and positivity filters in use.

     do i = -3,N+3
!!$        qface(i) = pcm_face_value_monotonic(q(i-2:i+3),doposlimit)
        qface(i) = max(0.d0,fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2)))
        if ((qface(i)-q(i))*(q(i+1)-qface(i)).lt.0.d0) then 
           ! possible grid scale monotonicity violation
           !   -- make further check to confirm monotonicity violation
           if (((q(i) - q(i-1))*(q(i+2) - q(i+1)) .ge. 0.d0) .OR. &
                ((q(i) - q(i-1))*(q(i-1) - q(i-2)) .le. 0.d0) .OR. &
                ((q(i+2) - q(i+1))*(q(i+3) - q(i+2)) .le. 0.d0) .OR. &
                ((qface(i)  - q(i))*(q(i) - q(i-1)) .le. 0.d0)) then
              monlimit(i) = 1
              ! use either left or right cell avg value 
              !   -- whichever is closer to the current estimate
              if (abs(qface(i)- q(i+1)) .ge. abs(qface(i)- q(i))) then
                 qface(i)= max(0.d0,q(i))
              else
                 qface(i)= max(0.d0,q(i+1))
              end if
           end if
        end if
     end do
     !     
     !# compute pcm flux across domain
     !================================
     do i = 0,N
!!$        !# qface holds interpolant of q at faces.
!!$        qface(i+3) = pcm_face_value_monotonic(q(i+1:i+7),doposlimit)
        !     
        !# compute PCM reconstruction of scalar field in upwind cell
        if(rhou(i).ge.0.d0) then
           cr = rhou(i)*dt/rhoprime(i)
           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
!!$           qq = qface(i+2:i-3:-1)
           ql = qface(i)
           qr = qface(i-1)
           qbar = q(i)
!!$           p = pcm_coef_monotonic(qq,q(i),dq,doposlimit) ! ql and qr flipped
        else
           cr = ABS(rhou(i)*dt/rhoprime(i+1))
           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
!!$            qq = qface(i-2:i+3)
           ql = qface(i)
           qr = qface(i+1)
           qbar = q(i+1)
!!$           p = pcm_coef_monotonic(qq,q(i+1),dq,doposlimit) 
        end if

        ! this is equation 13 in Zerroukat et al QJRMS 128:2801 (2002)
        a0 =       ql                                  ! constant term
        a1 = -6.d0*ql + 6.d0*qbar           - 2.d0*dq  ! linear term
        a2 =  9.d0*ql - 6.d0*qbar - 3.d0*qr + 6.d0*dq  ! quadratic term
        a3 = -4.d0*ql             + 4.d0*qr - 4.d0*dq  ! cubic term

        if(a1*(a1+2.d0*a2+3.d0*a3).lt.0.d0) then

           if(rhou(i).ge.0.d0) then
              do ii = 1,6
                 qq(ii) = qface(i+3-ii)
              end do
           else
              qq = qface(i-2:i+3)
           end if

           ! one extrema in interval.
           !   check for subgrid-scale violation
           if((qq(3)-qq(2))*(qq(5)-qq(4)).ge.0.d0 & !bad extrema
                .OR.(qq(2)-qq(1))*(qq(3)-qq(2)).le.0.d0 &
                .OR.(qq(5)-qq(4))*(qq(6)-qq(5)).le.0.d0 &
                .OR.(qq(3)-qq(2))*a1.le.0.d0) then 
              monlimit(i) = monlimit(i) + 2
              ! bad extrema, reduce order
              ! start by zeroing out cubic term
              a3 = 0.d0
              ! reduce order, if qbar outside of [ql,qr], use
              ! piecewise constant, otherwise use parabola
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
           
        else
           ! there are either zero or two extrema in the interval
           ! check whether there are two extrema
           if((abs(-a2/(3.d0*a3)-0.5d0).le.0.5d0) &
                .AND. (a1*(a1 - a2**2/(3.d0*a3)).lt.0.d0)) then
              ! there are two extrema in the interval
              !   this is a subgrid-scale violation of monotonicity
              !   reduce order
              monlimit(i) = monlimit(i) + 2
              ! bad extrema, reduce order
              ! start by zeroing out cubic term
              a3 = 0.d0
              ! reduce order, if qbar outside of [ql,qr], use
              ! piecewise constant, otherwise use parabola
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
           
        end if

        ! check for negative values in the interval and eliminate
        if(a3.ne.0.d0) then

           discr = a2**2 - 3.d0*a1*a3

           if(discr.ge.0.d0) then
              ! two extrema in general -- find distance from center of [0,1] interval
              x1 = (-a2+sqrt(discr))/(3.d0*a3)
              x2 = (-a2-sqrt(discr))/(3.d0*a3)
              x1dis = ABS(x1-0.5d0)
              x2dis = ABS(x2-0.5d0)

              if (((x1dis.lt.0.5d0).AND.(a0+x1*(a1+x1*(a2+x1*a3)).lt.0.d0)) .OR. &
                   ((x2dis.lt.0.5d0).AND.(a0+x2*(a1+x2*(a2+x2*a3)).lt.0.d0))) then
                 ! one of the extrema ia negative, use piecewise
                 ! constant reconstruction
                 a0 = qbar
                 a1 = 0.d0
                 a2 = 0.d0
                 a3 = 0.d0
              end if
           end if

        elseif(a2.ne.0.d0) then

           x1 = -a1/(2.d0*a2)
           if((ABS(x1-0.5d0).lt.0.5d0).AND.(a0+x1*(a1+x1*a2).lt.0.d0)) then
              a0 = qbar
              a1 = 0.d0
              a2 = 0.d0
              a3 = 0.d0
           end if
        end if

!!$
!!$        discr = a2**2 - 3.d0*a1*a3
!!$
!!$        if(discr.ge.0.d0) then
!!$           ! two extrema in general -- find distance from center of [0,1] interval
!!$           x1 = (-a2+sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
!!$           x2 = (-a2-sqrt(discr))/(3.d0*a3 + EPSILON(1.d0))
!!$           x1dis = ABS(x1-0.5d0)
!!$           x2dis = ABS(x2-0.5d0)
!!$
!!$           if (((x1dis.lt.0.5d0).AND.(a0+x1*(a1+x1*(a2+x1*a3)).lt.0.d0)) .OR. &
!!$                ((x2dis.lt.0.5d0).AND.(a0+x2*(a1+x2*(a2+x2*a3)).lt.0.d0))) then
!!$              ! one of the extrema ia negative, use piecewise
!!$              ! constant reconstruction
!!$              a0 = qbar
!!$              a1 = 0.d0
!!$              a2 = 0.d0
!!$              a3 = 0.d0
!!$
!!$           elseif (MAX(x1dis,x2dis).lt.0.5d0) then
!!$              monlimit(i) = monlimit(i) + 2
!!$              ! both extrema are within interval -->
!!$              ! reduce order, start by zeroing out cubic term
!!$              a3 = 0.d0
!!$              ! reduce order, if qbar outside of [ql,qr], use
!!$              ! piecewise constant, otherwise use parabola
!!$              if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
!!$                 a0 = qbar ! piecewise constant reconstruction
!!$                 a1 = 0.d0
!!$                 a2 = 0.d0
!!$              else
!!$                 ! define coefficients of parabola
!!$                 a0 = ql
!!$                 a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
!!$                 a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
!!$                 if (qr.gt.ql) then
!!$                    if (3.d0*qbar .lt. 2.d0*ql + qr) then
!!$                       !bound parabola from below to remove undershoot
!!$                       a0 = ql
!!$                       a1 = 0.d0
!!$                       a2 = -3.d0*ql+3.d0*qbar
!!$                    elseif (3.d0*qbar .gt. ql + 2.d0*qr) then
!!$                       !bound parabola from above to remove overshoot
!!$                       a0 = 3.d0*qbar-2.d0*qr
!!$                       a1 = 6.d0*qr-6.d0*qbar
!!$                       a2 = -3.d0*qr+3.d0*qbar
!!$                    end if
!!$                 else
!!$                    if (3.d0*qbar .gt. 2.d0*ql + qr) then
!!$                       !bound parabola from above to remove overshoot
!!$                       a0 = ql
!!$                       a1 = 0.d0
!!$                       a2 = -3.d0*ql+3.d0*qbar
!!$                    elseif (3.d0*qbar .lt. ql + 2.d0*qr) then
!!$                       !bound parabola from below to remove undershoot
!!$                       a0 = 3.d0*qbar-2.d0*qr
!!$                       a1 = 6.d0*qr-6.d0*qbar
!!$                       a2 = -3.d0*qr+3.d0*qbar
!!$                    end if
!!$                 end if
!!$              end if
!!$           elseif (MIN(x1dis,x2dis).lt.0.5d0) then
!!$              ! exactly one extrema within the interval
!!$              ! make a check for subgrid monotonicity violations
!!$              if((qq(3)-qq(2))*(qq(5)-qq(4)).ge.0.d0 & !bad extrema
!!$                   .OR.(qq(2)-qq(1))*(qq(3)-qq(2)).le.0.d0 &
!!$                   .OR.(qq(5)-qq(4))*(qq(6)-qq(5)).le.0.d0 &
!!$                   .OR.(qq(3)-qq(2))*a1.le.0.d0) then 
!!$                 monlimit(i) = monlimit(i) + 2
!!$                 ! bad extrema, reduce order
!!$                 ! start by zeroing out cubic term
!!$                 a3 = 0.d0
!!$                 ! reduce order, if qbar outside of [ql,qr], use
!!$                 ! piecewise constant, otherwise use parabola
!!$                 if ((qbar - ql)*(qbar - qr) .gt. 0.d0) then
!!$                    a0 = qbar ! piecewise constant reconstruction
!!$                    a1 = 0.d0
!!$                    a2 = 0.d0
!!$                 else
!!$                    ! define coefficients of parabola
!!$                    a0 = ql
!!$                    a1 = -4.d0*ql - 2.d0*qr + 6.d0*qbar
!!$                    a2 = 3.d0*ql + 3.d0*qr - 6.d0*qbar
!!$                    if (qr.gt.ql) then
!!$                       if (3.d0*qbar .lt. 2.d0*ql + qr) then
!!$                          !bound parabola from below to remove undershoot
!!$                          a0 = ql
!!$                          a1 = 0.d0
!!$                          a2 = -3.d0*ql+3.d0*qbar
!!$                       elseif (3.d0*qbar .gt. ql + 2.d0*qr) then
!!$                          !bound parabola from above to remove overshoot
!!$                          a0 = 3.d0*qbar-2.d0*qr
!!$                          a1 = 6.d0*qr-6.d0*qbar
!!$                          a2 = -3.d0*qr+3.d0*qbar
!!$                       end if
!!$                    else
!!$                       if (3.d0*qbar .gt. 2.d0*ql + qr) then
!!$                          !bound parabola from above to remove overshoot
!!$                          a0 = ql
!!$                          a1 = 0.d0
!!$                          a2 = -3.d0*ql+3.d0*qbar
!!$                       elseif (3.d0*qbar .lt. ql + 2.d0*qr) then
!!$                          !bound parabola from below to remove undershoot
!!$                          a0 = 3.d0*qbar-2.d0*qr
!!$                          a1 = 6.d0*qr-6.d0*qbar
!!$                          a2 = -3.d0*qr+3.d0*qbar
!!$                       end if
!!$                    end if
!!$                 end if
!!$              end if
!!$
!!$           end if ! checks for extrema within interval
!!$
!!$        end if ! check for discr of first derivative polynomial > 0
        ! ========== END MONOTONICITY FILTER =====================

         !# compute PCM flux by integrating across reconstruction.
        qavg(i) = a0 + cr*(a1/2.d0 + cr*(a2/3.d0 + cr*a3/4.d0))
        fx(i) = rhou(i)*qavg(i)

     end do

!!$
!!$     do i = -2,0
!!$        qface(i) = pcm_face_value(q(i-1:i+2),doposlimit)
!!$!!$        qface(i) = fac1*(q(i) + q(i+1)) + fac2*(q(i-1) + q(i+2))
!!$     end do
!!$     !     
!!$     !# compute pcm flux across domain
!!$     !================================
!!$     do i = 0,N
!!$        !# qface holds interpolant of q at faces.
!!$        qface(i+1) = pcm_face_value(q(i:i+3),doposlimit)
!!$!!$        qface(i+1) = fac1*(q(i+1) + q(i+2)) + fac2*(q(i) + q(i+3))
!!$        !     
!!$        !# compute PCM reconstruction of scalar field in upwind cell
!!$        if(rhou(i).ge.0.d0) then
!!$           cr = rhou(i)*dt/rhoprime(i)
!!$           dq = fac3*(q(i-1)-q(i+1)) + fac4*(q(i-2)-q(i+2)) ! reverse sign
!!$           p = pcm_coef(qface(i),qface(i-1),q(i),dq,doposlimit) ! ql and qr flipped
!!$        else
!!$           cr = ABS(rhou(i)*dt/rhoprime(i+1))
!!$           dq = fac3*(q(i+2)-q(i)) + fac4*(q(i+3)-q(i-1))
!!$           p = pcm_coef(qface(i),qface(i+1),q(i+1),dq,doposlimit) 
!!$        end if
!!$        !# compute PCM flux by integrating across reconstruction.
!!$        qavg(i) = p(1) + cr*(p(2)/2.d0 + cr*(p(3)/3.d0 + cr*p(4)/4.d0))
!!$        fx(i) = rhou(i)*qavg(i)
!!$
!!$     end do
!!$
  end if

  ! fix bcs at left end if necessary
  if((bcleft.eq.1).and.(rhou(0).gt.0.d0)) then
     fx(0) = rhou(0)*qavg(1)
  elseif(bcleft.eq.3) then
     fx(0) = fxleft
  end if

  ! fix bcs at right end if necessary
  if((bcright.eq.1).and.(rhou(N).lt.0.d0)) then
     fx(N) = rhou(N)*qavg(N-1)
  elseif(bcright.eq.3) then
     fx(N) = fxright
  end if


  !# UPDATE SCALAR MIXING RATIO
  !==================================
  if(doposlimit) then
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i)) !max(0.d0,rhoq(i) + dt*(fx(i-1) - fx(i)))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1)-rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do
  else
     do i = 1,N
        rhoq(i) = rhoq(i) + dt*(fx(i-1) - fx(i))
        rhoprime(i) = rhoprime(i) + dt*(rhou(i-1)-rhou(i))
        q(i) = rhoq(i)/rhoprime(i)
     end do
  end if

  lambdavec(:) = 0.
  limitpts(:) = monlimit(0:N)
end subroutine pcmsweep

