module precision
   use, intrinsic :: iso_fortran_env
   implicit none
   
   private

   integer, parameter, public :: si = int16
   integer, parameter, public :: ri = int32
   integer, parameter, public :: li = int64
   integer, parameter, public :: wi = ri

   integer, parameter, public :: sp = real32
   integer, parameter, public :: dp = real64
   integer, parameter, public :: wp = dp

   integer, parameter, public :: wl = 4

end module precision
