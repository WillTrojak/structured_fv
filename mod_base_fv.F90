module base_fv
   use precision
   use grid, only : domain
   implicit none

   private

   type, abstract :: fv_base
      private
      
      type(domain), public :: omega
      
      real(kind=wp), public, allocatable :: q(:,:)
      real(kind=wp), public, allocatable :: qb(:,:)
      real(kind=wp), public, allocatable :: df(:,:)
      
   contains
      procedure(prim_to_cons), public, deferred, pass(self) :: cons
      procedure(cons_to_prim), public, deferred, pass(self) :: prim
      procedure(flux), public, deferred, pass(self) :: iflux
      procedure(flux), public, deferred, pass(self) :: bflux
      procedure(flux), public, deferred, pass(self) :: residual
      procedure(flux), public, deferred, pass(self) :: bcs
      procedure(riemann_flux), deferred, pass(self) :: rsolve
      procedure(point_flux_x), deferred, pass(self) :: pflux_x
      procedure(initial_cond), deferred, pass(self) :: init_sol
      
      procedure, public, nopass :: transform_to
      procedure, public, nopass :: transform_from
   end type fv_base

   abstract interface
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function prim_to_cons(self,qp) result(qc)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(in) :: self
         real(kind=wp), intent(in) :: qp(:)
         real(kind=wp) :: qc(size(qp))
       end function prim_to_cons
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
      function cons_to_prim(self,qc) result(qp)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(in) :: self
         real(kind=wp), intent(in) :: qc(:)
         real(kind=wp) :: qp(size(qc))
       end function cons_to_prim
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       
      subroutine flux(self)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(inout) :: self
      end subroutine flux
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function riemann_flux(self,ql,qr) result(fc)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(in) :: self
         real(kind=wp), intent(in) :: ql(:),qr(:)
         real(kind=wp) :: fc(size(ql))
      end function riemann_flux
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function point_flux_x(self,q) result(fx)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(in) :: self
         real(kind=wp), intent(in) :: q(:)
         real(kind=wp) :: fx(size(q))
      end function point_flux_x
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine initial_cond(self)
         use precision
         import fv_base
         implicit none
         class(fv_base), intent(inout) :: self
      end subroutine initial_cond
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   end interface

   public :: fv_base

contains
   !**********************************************************************   
   function transform_to(n,u) result(t)
      implicit none
      
      real(kind=wp), intent(in) :: n(:),u(:) 

      real(kind=wp) :: t(size(u))

      t(1) = u(1)
      t(2) =  n(1)*u(2) + n(2)*u(3)
      t(3) = -n(2)*u(2) + n(1)*u(3)
      t(4) = u(4)
      
      return
   end function transform_to
   !**********************************************************************
   function transform_from(n,t) result(u)
      implicit none

      real(kind=wp), intent(in) :: n(:),t(:) 

      real(kind=wp) :: u(size(t))
      
      u(1) = t(1)
      u(2) = n(1)*t(2) - n(2)*t(3)
      u(3) = n(2)*t(2) + n(1)*t(3)
      u(4) = t(4)

      return
   end function transform_from
   !**********************************************************************
end module base_fv
