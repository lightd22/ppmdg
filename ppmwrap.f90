subroutine ppmwrap(rhoq,q,rhou,rho,rhop,flx,dt, &
     N,npad,nmaxcfl,bctype,fbc, &
     dosemilagr,doselect,domonotonic,dopositive,dopcm,doweno, &
     scale,nmethod,lambdamax,epslambda,lambda,monlimit, &
     DGu, num_elem,elemdx,nnodes,DGCmat,DGCmatINV,DGnodes,DGwghts,DGDmat)

  implicit none

  integer, intent(in) :: N, npad, nmaxcfl, bctype(2), nmethod
  real(kind=8), intent(inout) :: rhoq(0:N+1)
  real(kind=8), intent(inout) :: q(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(in) :: rhou(-2:N+2)
  real(kind=8), intent(in) :: rho(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(inout) :: rhop(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8), intent(out) :: flx(0:N)
  real(kind=8), intent(in) :: dt, fbc(2)

  logical, intent(in) :: dosemilagr,doselect,domonotonic,dopositive, &
       dopcm,doweno

  real(kind=8), intent(in) :: scale, lambdamax, epslambda
  real(kind=8), intent(out) :: lambda(0:N), monlimit(0:N)

  real(kind=8) :: q_new(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: rhop_new(1-npad-nmaxcfl:N+npad+nmaxcfl)
  integer :: npad2, nselpad

  integer, parameter :: ns = 1
  real(kind=8) :: scalar_old(ns,1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: scalar_new(ns,1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: u(1-npad-nmaxcfl:N+npad+nmaxcfl+1)
  real(kind=8) :: r_old(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: r_new(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: dx(1-npad-nmaxcfl:N+npad+nmaxcfl)
  real(kind=8) :: dt_real
  integer :: limiter, its, ite, idata_s, idata_e, ims, ime

  ! DG Parameters
  INTEGER, INTENT(IN) :: num_elem,nnodes
  REAL(KIND=8), INTENT(IN) :: elemdx
  REAL(KIND=8), DIMENSION(0:nnodes) :: DGnodes, DGwghts
  REAL(KIND=8), DIMENSION(0:nnodes,0:nnodes), INTENT(IN) :: DGCmat, DGCmatINV,DGDmat
  REAL(KIND=8), DIMENSION(1:N), INTENT(IN) :: DGu


  nselpad = 0

  if(dosemilagr) then

     select case(nmethod)
!!$     case(1:4) ! PPM FCT (nolimit, poslimit, mnlimit, posmonl)
!!$        call ppmsweep_fct_semilagrangian(rhoq,q,rho,rhop,rhou, &
!!$             flx, dt,N,npad,nmaxcfl,bctype,fbc, &
!!$             domonotonic,dopositive,lambda,monlimit)
     case(5:6) ! PPM FCT (selimit, poselim)
        call ppmsweep_fct_select_semilagrangian(rhoq,q,rho,rhop,rhou, &
             flx, dt,N,npad,nmaxcfl,bctype,fbc, &
             dopositive,scale,nselpad,lambdamax,epslambda,lambda,monlimit)

!!$     case(21:24) ! PPM Polynomial Modification (nolimit, poslimit, mnlimit, posmonl)
!!$        call ppmsweep_pmod_semilagrangian(rhoq,q,rho,rhop,rhou, &
!!$             flx, dt,N,npad,nmaxcfl,bctype,fbc, &
!!$             domonotonic,dopositive,lambda,monlimit)
     case(25:26) ! PPM Polynomial Modification (selimit, poselim)
        call ppmsweep_pmod_select_semilagrangian(rhoq,q,rho,rhop,rhou, &
             flx, dt,N,npad,nmaxcfl,bctype,fbc, &
             dopositive,scale,nselpad,lambdamax,epslambda,lambda,monlimit)

     CASE(81:82) !PPM PMOD CLEAN (HOPEFULLY SPEEDY IMPLEMENTATION)

       npad2 = npad + nmaxcfl
             call ppmsweep_pmod_select_clean_sl(q,rhou(0:N),rhop, &
                  flx, dt,N,npad2,nmaxcfl,bctype,fbc,dopositive,scale)
             rhoq(1:N) = rhop(1:N)*q(1:N)
             lambda(:) = 0.
             monlimit(:) = 0.
     case default
        write(*,*) 'nmethod = ', nmethod
        STOP 'in ppmwrap.f90: no semi-Lagrangian method with this nmethod'
     end select
        
  else

     npad2 = npad + nmaxcfl

     select case(nmethod)
     case(1:4) ! PPM FCT (nolimit, poslimit, mnlimit, posmonl)
        call ppmsweep_fct(rhoq,q,rhou,rho,rhop, &
             flx, dt,N,npad2,bctype,fbc, &
             domonotonic,dopositive,lambda,monlimit)
     case(5:6) ! PPM FCT (selimit, poselim)
        call ppmsweep_fct_select(rhoq,q,rhou,rho,rhop, &
             flx, dt,N,npad2,bctype,fbc, &
             dopositive,scale,nselpad,lambdamax,epslambda,lambda,monlimit)

     case(21:24) ! PPM Polynomial Modification (nolimit, poslimit, mnlimit, posmonl)
        call ppmsweep_pmod(rhoq,q,rhou,rho,rhop, &
             flx, dt,N,npad2,bctype,fbc, &
             domonotonic,dopositive,lambda,monlimit)
     case(25:26) ! PPM Polynomial Modification (selimit, poselim)
        call ppmsweep_pmod_select(rhoq,q,rhou,rho,rhop, &
             flx, dt,N,npad2,bctype,fbc, &
             dopositive,scale,nselpad,lambdamax,epslambda,lambda,monlimit)

     CASE(81:82) !PPM PMOD CLEAN (HOPEFULLY SPEEDY IMPLEMENTATION)
             call ppmsweep_pmod_select_clean(q,rhou(0:N),rhop(0:N+1), &
                  flx, dt,N,npad2,bctype,fbc,dopositive,scale)
             rhoq(1:N) = rhop(1:N)*q(1:N)
             lambda(:) = 0.
             monlimit(:) = 0.
     case(99) ! PPM/DG Hybrid, no limiting, ADDED BY DEVIN
        CALL dgsweep(rhoq,rhop,num_elem,elemdx,nnodes,DGnodes,DGwghts,DGu,N,DGCmat,DGCmatINV,DGDmat,npad2,dt)
     case default
        write(*,*) 'nmethod = ', nmethod
        STOP 'in ppmwrap.f90: no CFL<1 method with this nmethod'
     end select
        
  end if

end subroutine ppmwrap

