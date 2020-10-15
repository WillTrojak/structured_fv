module erk_int
   use precision
   use base_int
   implicit none

   type, extends(int_base) :: int_erk
      private

      integer(kind=wi) :: nrk,ne,nv,ir
      
      real(kind=wp), allocatable :: k(:,:,:)

      real(kind=wp), allocatable :: ab(:,:)

      logical(kind=wl), public :: stage_complete
      logical(kind=wl) :: initd = .false.
   contains
      generic, public :: init => init_int_erk
      procedure :: init_int_erk
      procedure, public :: integrate => explicit_rk
   end type int_erk

contains
   !**********************************************************************
   subroutine init_int_erk(self,int_type,ne,nv)
      implicit none

      class(int_erk), intent(inout) :: self

      character(*), intent(in) :: int_type

      integer(kind=wi), intent(in) :: ne,nv

      if(self%initd)then
         print *,'already initiated'
         return
      endif
      
      select case(int_type)
      case('rk3')
         self%nrk = 3
         allocate(self%ab(self%nrk+1,self%nrk))
         self%ab(1,:) = [0e0_wp,  0e0_wp, 0e0_wp]
         self%ab(2,:) = [0.5e0_wp,  0e0_wp, 0e0_wp]
         self%ab(3,:) = [-1e0_wp,  2e0_wp, 0e0_wp]
         self%ab(4,:) = [1e0_wp, 4e0_wp, 1e0_wp]/6e0_wp
      case('ssp-rk3')
         self%nrk = 3
         allocate(self%ab(self%nrk+1,self%nrk))
         
         self%ab(1,:) = [0e0_wp,  0e0_wp, 0e0_wp]
         self%ab(2,:) = [1e0_wp,  0e0_wp, 0e0_wp]
         self%ab(3,:) = [0.25e0_wp, 0.25e0_wp, 0e0_wp]
         self%ab(4,:) = [1e0_wp, 1e0_wp, 4e0_wp]/6e0_wp
      case('rk4')
         self%nrk = 4
         allocate(self%ab(self%nrk+1,self%nrk))
         
         self%ab(1,:) = [  0e0_wp,   0e0_wp, 0e0_wp, 0e0_wp]
         self%ab(2,:) = [0.5e0_wp,   0e0_wp, 0e0_wp, 0e0_wp]
         self%ab(3,:) = [  0e0_wp, 0.5e0_wp, 0e0_wp, 0e0_wp]
         self%ab(4,:) = [  0e0_wp,   0e0_wp, 1e0_wp, 0e0_wp]
         self%ab(5,:) = [1e0_wp, 2e0_wp, 2e0_wp, 1e0_wp]/6e0_wp
      case default
         print *,'ERK type not recognised:',int_type
         stop
      end select

      allocate(self%q0(ne,nv))
      allocate(self%k(self%nrk,ne,nv))
      allocate(self%dt(nv))

      self%ir = 1
      self%stage_complete = .true.
      
      self%initd = .true.
      
      return
   end subroutine init_int_erk
   !**********************************************************************
   subroutine explicit_rk(self,dv,f,q)
      implicit none

      class(int_erk), intent(inout) :: self

      real(kind=wp), intent(in) :: dv(:),f(:,:)

      real(kind=wp), intent(inout) :: q(:,:)

      integer(kind=wi) :: iv,ie,is
      
      real(kind=wp) :: eps
      
      do iv=1,self%nv
         do ie=1,self%ne
            self%k(self%ir,ie,iv) = q(ie,iv)/dv(iv)
            
            eps = 0e0_wp
            do is=1,self%ir
               eps = eps + self%ab(self%ir+1,is)*self%k(is,ie,iv)
            enddo
            
            q(ie,iv) = self%q0(ie,iv) - self%dt(iv)*eps
         enddo
      enddo

      self%ir = self%ir + 1
      self%stage_complete = .false.
      if(self%ir .eq. self%nrk + 1) then
         self%ir = 1
         self%stage_complete = .true.
      endif
      
      return
   end subroutine explicit_rk
   !**********************************************************************
end module erk_int
