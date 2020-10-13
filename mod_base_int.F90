module base_int
   use precision
   implicit none
  
   type, abstract :: int_base
      private

      real(kind=wp), public, allocatable :: dt(:)
      
   contains
      procedure(explicit_int), deferred, pass(self) :: integrate
   end type int_base

   abstract interface
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine explicit_int(self,dv,f,q)
         use precision
         import int_base
         implicit none
         class(int_base), intent(inout) :: self
         real(kind=wp), intent(in) :: dv(:),f(:,:)
         real(kind=wp), intent(inout) :: q(:,:)
      end subroutine explicit_int
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   end interface

   public :: int_base
   
end module base_int
  
