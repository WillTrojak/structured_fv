module euler_fv
   use precision
   use base_fv
   implicit none

   type, extends(fv_base) :: fv_euler
      private
      integer(kind=wp), public :: ne

      real(kind=wp), private :: gamma, gm1, rgm1
      
      logical(kind=wl) :: initd = .false.
   contains
      generic, public :: init => init_fv_euler
      procedure, pass(self) :: init_fv_euler
      
      procedure, public, pass(self) :: cons => euler_prim_to_cons
      procedure, public, pass(self) :: prim => euler_cons_to_prim
      procedure, public, pass(self) :: iflux => euler_interior_flux
      procedure, public, pass(self) :: bflux => euler_boundary_flux
      procedure, public, pass(self) :: rsolve => rusanov_flux
      procedure, public, pass(self) :: pflux_x => euler_point_flux_x
      
      procedure, public, pass(self) :: consts => compute_consts
   end type fv_euler
   
contains
   !**********************************************************************
   subroutine init_fv_euler(self,omega,ne,gamma)
      use grid 
      implicit none

      class(fv_euler), intent(inout) :: self

      type(domain) :: omega

      integer(kind=wi), intent(in) :: ne

      real(kind=wp), intent(in) :: gamma
      
      if(self%initd)then
         print *,'already initiated'
         return
      endif

      self%omega = omega
      
      self%ne = ne

      call self%consts(gamma)
      
      allocate(self%q(ne,self%omega%nv))
      allocate(self%df(ne,self%omega%nv))
         
      self%initd = .true.
      
      return
   end subroutine init_fv_euler
   !**********************************************************************
   subroutine compute_consts(self,gamma)
      implicit none

      class(fv_euler), intent(inout) :: self
      
      real(kind=wp), intent(in) :: gamma

      self%gamma = gamma
      self%gm1 = self%gamma - 1e0_wp
      self%rgm1 = 1e0_wp/self%gm1
      
      return
   end subroutine compute_consts
   !**********************************************************************
   function euler_prim_to_cons(self,qp) result(qc)
      implicit none

      class(fv_euler), intent(in) :: self
      
      real(kind=wp), intent(in) :: qp(:) 
      real(kind=wp) :: qc(size(qp))
      
      qc(1) = qp(1)
      qc(2) = qp(1)*qp(2)
      qc(3) = qp(1)*qp(3)
      qc(4) = 0.5e0_wp*qp(1)*(qp(2)*qp(2) + qc(3)*qc(3)) + qp(4)*self%rgm1
      
      return
   end function euler_prim_to_cons
   !**********************************************************************
   function euler_cons_to_prim(self,qc) result(qp)
      implicit none

      class(fv_euler), intent(in) :: self
      
      real(kind=wp), intent(in) :: qc(:) 
      real(kind=wp) :: qp(size(qc))

      real(kind=wp) :: rrho
      
      rrho = 1e0_wp/qc(1)
            
      qp(1) = qc(1)
      qp(2) = qc(2)*rrho
      qp(3) = qc(3)*rrho
      qp(4) = self%gm1*(qc(4) - 0.5e0_wp*(qc(2)*qc(2) + qc(3)*qc(3))*rrho)
      
      return
   end function euler_cons_to_prim
   !**********************************************************************
   function euler_point_flux_x(self,q) result(fx)
      implicit none

      class(fv_euler), intent(in) :: self
      
      real(kind=wp), intent(in) :: q(:)

      real(kind=wp) :: fx(size(q))

      real(kind=wp) :: rrho,p
      
      rrho = 1e0_wp/q(1)
      p = self%gm1*(q(4) - 0.5e0_wp*(q(2)*q(2) + q(3)*q(3))*rrho)
      
      fx(1) = q(2)
      fx(2) = q(2)*q(2)*rrho + p
      fx(3) = q(2)*q(3)*rrho
      fx(4) = q(2)*(q(4) + p)*rrho
      
      return
   end function euler_point_flux_x
   !**********************************************************************
   function rusanov_flux(self,ql,qr) result(fc)
      implicit none

      class(fv_euler), intent(in) :: self
      
      real(kind=wp), intent(in) :: ql(:),qr(:)

      real(kind=wp) :: fc(size(ql))

      real(kind=wp) :: hl,fl(size(ql)),hr,fr(size(qr))
      real(kind=wp) :: rrhosqrt,ub,hb,ab,s
            
      hl = self%gm1*(ql(4) - 0.5e0_wp*(ql(2)*ql(2) + ql(3)*ql(3))/ql(1))
      hr = self%gm1*(qr(4) - 0.5e0_wp*(qr(2)*qr(2) + qr(3)*qr(3))/qr(1))
      ! Roe average velocity and wave speed
      rrhosqrt = 1e0_wp/(sqrt(ql(1)) + sqrt(qr(1)))
      ub = (sqrt(qr(1))*ql(2) + sqrt(ql(1))*qr(2))*rrhosqrt*rrhosqrt
      hb = (sqrt(ql(1))*hl + sqrt(qr(1))*hr)*rrhosqrt
      ab = sqrt(self%gm1*(hb - 0.5e0_wp*ub*ub))

      s = abs(ub) + ab

      fl = self%pflux_x(ql)
      fr = self%pflux_x(qr)
      
      fc(1) = 0.5e0_wp*(fl(1) + fr(1)) - 0.5e0_wp*s*(qr(1) - ql(1)) 
      fc(2) = 0.5e0_wp*(fl(2) + fr(2)) - 0.5e0_wp*s*(qr(2) - ql(2))
      fc(3) = 0.5e0_wp*(fl(3) + fr(3)) - 0.5e0_wp*s*(qr(3) - ql(3))
      fc(4) = 0.5e0_wp*(fl(4) + fr(4)) - 0.5e0_wp*s*(qr(4) - ql(4))
      
      return
   end function rusanov_flux
   !**********************************************************************
   subroutine euler_interior_flux(self)
      use maths, only : l2norm
      implicit none

      class(fv_euler), intent(inout) :: self

      integer(kind=wi) :: ie,el,er
      
      real(kind=wp) :: qlt(self%ne),qrt(self%ne)
      real(kind=wp) :: fct(self%ne),fc(self%ne)
      real(kind=wp) :: n(self%omega%nd),nu(self%omega%nd),dl
      
      ! Interior fluxes
      do ie=1,self%omega%nfi
         el = self%omega%face_i(3,ie) 
         er = self%omega%face_i(4,ie)

         n(:) = self%omega%n_i(:,ie)
         dl = l2norm(n)
         nu(:) = n(:)/dl

         qlt = self%transform_to(nu,self%q(:,el))
         qrt = self%transform_to(nu,self%q(:,er))

         fct = self%rsolve(qlt,qrt)
         fc = self%transform_from(nu,fct)
         self%df(:,el) =  fc(:)*dl
         self%df(:,er) = -fc(:)*dl
      enddo
      
      return
   end subroutine euler_interior_flux
   !**********************************************************************
   subroutine euler_boundary_flux(self)
      use maths, only : l2norm
      implicit none

      class(fv_euler), intent(inout) :: self

      integer(kind=wi) :: ie,el,er
      
      real(kind=wp) :: qlt(self%ne),qrt(self%ne)
      real(kind=wp) :: fct(self%ne),fc(self%ne)
      real(kind=wp) :: n(self%omega%nd),nu(self%omega%nd),dl
      
      do ie=1,self%omega%nfb
         el = self%omega%face_b(3,ie) 
         er = self%omega%face_b(4,ie)
         
         n(:) = self%omega%n_b(:,ie)
         dl = l2norm(n)
         nu(:) = n(:)/dl

         qlt = self%transform_to(nu,self%q(:,el))
         qrt = self%transform_to(nu,self%qb(:,er))

         fct = self%rsolve(qlt,qrt)
         fc = self%transform_from(nu,fct)
         self%df(:,el) =  fc(:)*dl
      enddo
      
      return
   end subroutine euler_boundary_flux
   !**********************************************************************
   subroutine boundary_value(self)
      use maths, only : l2norm
      implicit none
      
      class(fv_euler), intent(inout) :: self

      integer(kind=wi) :: i,j,nf,nl,nb
      real(kind=wp) :: nu(self%omega%nd),rvdn
      real(kind=wp) :: r,ma,p
      
      nf = 0
      do j=1,self%omega%nj
         do i=1,self%omega%ni
            ! South
            if((j .eq. 1) .and. (i .ne. self%omega%ni))then
               nf = nf + 1
               nl = self%omega%face_b(3,nf)
               nb = self%omega%face_b(4,nf) 

               nu(:) = self%omega%n_b(:,nf)
               nu(:) = nu(:)/l2norm(nu)
               
               self%qb(1,nb) = self%q(1,nl)
               self%qb(4,nb) = self%q(4,nl)

               rvdn = self%q(2,nl)*nu(1) + self%q(3,nl)*nu(2)
               self%qb(2,nb) = self%q(2,nl) - rvdn*nu(1)
               self%qb(3,nb) = self%q(3,nl) - rvdn*nu(1)
               
            endif

            ! North
            if((j .eq. self%omega%nj) .and. (i .ne. self%omega%ni))then
               nf = nf + 1
               nl = self%omega%face_b(3,nf)
               nb = self%omega%face_b(4,nf) 

               nu(:) = self%omega%n_b(:,nf)
               nu(:) = nu(:)/l2norm(nu)
               
               self%qb(1,nb) = self%q(1,nl)
               self%qb(4,nb) = self%q(4,nl)

               rvdn = self%q(2,nl)*nu(1) + self%q(3,nl)*nu(2)
               self%qb(2,nb) = self%q(2,nl) - rvdn*nu(1)
               self%qb(3,nb) = self%q(3,nl) - rvdn*nu(1)               
            endif

            ! West
            if((i .ne. 1) .and. (j .ne. self%omega%nj))then
               nf = nf + 1
               nl = self%omega%face_b(3,nf)
               nb = self%omega%face_b(4,nf)

               self%qb(1,nb) = self%q(1,nl)
               self%qb(2,nb) = self%q(2,nl)
               self%qb(3,nb) = self%q(3,nl)
               self%qb(4,nb) = self%q(4,nl)              
            endif

            ! East
            if((i .ne. self%omega%ni) .and. (j .ne. self%omega%nj))then
               nf = nf + 1
               nl = self%omega%face_b(3,nf)
               nb = self%omega%face_b(4,nf)

               r = 1e0_wp
               ma = 2e0_wp
               p = 1e0_wp
               self%qb(1,nb) = r
               self%qb(2,nb) = ma*sqrt(r*self%gamma*p)
               self%qb(3,nb) = 0e0_wp
               self%qb(4,nb) = 0.5e0_wp*self%qb(2,nb)*self%qb(2,nb)/r + p*self%rgm1
            endif
         enddo
      enddo
      
      return
   end subroutine boundary_value
   !**********************************************************************
end module euler_fv
 
