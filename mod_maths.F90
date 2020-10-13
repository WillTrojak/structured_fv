module maths
   use precision
   implicit none

   private

   public :: l2norm

   interface l2norm
      module procedure l2norm_dp, l2norm_sp
   end interface l2norm
   
contains
   !**********************************************************************
   function l2norm_dp(x) result(l)
      implicit none

      real(kind=dp), intent(in) :: x(:)

      real(kind=dp) :: l
      
      integer(kind=wi) :: i,n

      n = size(x)

      l = 0e0_dp
      do i=1,n
         l = l + x(i)*x(i)
      enddo
      l = sqrt(l)

      return
   end function l2norm_dp
   !**********************************************************************
   function l2norm_sp(x) result(l)
      implicit none

      real(kind=sp), intent(in) :: x(:)

      real(kind=sp) :: l
      
      integer(kind=wi) :: i,n

      n = size(x)

      l = 0e0_sp
      do i=1,n
         l = l + x(i)*x(i)
      enddo
      l = sqrt(l)

      return
   end function l2norm_sp
   !**********************************************************************
end module maths
