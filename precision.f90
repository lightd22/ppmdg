! module which defines kinds of single/double precision variables.
module precision
  implicit none
  integer, parameter       :: single            & ! single precision kind.
                               = selected_real_kind(6,35) 
  integer, parameter       :: double            & ! double precision kind.
                               = selected_real_kind(12,150) 
end module precision
