module base_int
   use precision
   implicit none
  
   type, abstract :: int_base
      private

      real(kind=wp), public, allocatable :: dt(:)
      
      real(kind=wp), public, allocatable :: q0(:,:)
      
   contains
      procedure(explicit_int), deferred, pass(self) :: integrate
     
      procedure, public, pass(self) :: begin => begin_integration
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

contains
   !**********************************************************************
   subroutine begin_integration(self,q)
      implicit none

      class(int_base), intent(inout) :: self

      real(kind=wp), intent(in) :: q(:,:)

      integer(kind=wi) :: n1,n2,i,j

      n1 = size(q,1)
      n2 = size(q,2)
      
      do j=1,n2
         do i=1,n1
            self%q0(i,j) = q(i,j)
         enddo
      enddo

      return
   end subroutine begin_integration
   !**********************************************************************
end module base_int
  
