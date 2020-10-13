module grid
   use precision
   implicit none
  
   private
   
   type :: domain
      integer(kind=wi) :: ierr
      integer(kind=wi) :: nd,nv,nfi,nfb
      integer(kind=wi) :: ni,nj

      real(kind=wp), allocatable :: x(:,:)
      real(kind=wp), allocatable :: dv(:)
      
      integer(kind=wi), allocatable :: face_i(:,:),face_b(:,:)
      real(kind=wp), allocatable :: n_i(:,:),n_b(:,:)
   contains
      generic :: init => init_domain
      procedure, private :: init_domain
   end type domain
   
   public :: domain
   
contains
   !**********************************************************************
   subroutine init_domain(self,nd,ni,nj)
      implicit none

      class(domain) :: self
      
      integer(kind=wi), intent(in) :: nd,ni,nj
      
      self%ierr = 0
      
      self%nd = nd
      self%nv = (ni-1)*(nj-1)
      self%nfi = (ni-1)*(nj-2) + (ni-2)*(nj-1)
      self%nfb = 2*(nj-1) + 2*(ni-1)
      
      self%ni = ni; self%nj = nj

      allocate(self%x(nd,ni*nj))
      allocate(self%dv(ni*nj))
      
      allocate(self%face_i(4,self%nfi))
      allocate(self%n_i(nd,self%nfi))
      
      allocate(self%face_b(4,self%nfb))
      allocate(self%n_b(nd,self%nfb))
      
      self%x = grid_coords(nd,ni,nj)
      self%dv = volumes(ni,nj,self%x)
      
      self%face_i = interior_faces(ni,nj,self%nfi)
      self%n_i = normals(self%face_i,self%x)
      
      self%face_b = boundary_faces(ni,nj,self%nfb)
      self%n_b = normals(self%face_b,self%x)
      
      return
   end subroutine init_domain
   !**********************************************************************
   function grid_coords(nd,ni,nj) result(x)
      implicit none

      integer(kind=wi), intent(in) :: nd,ni,nj

      real(kind=wp) :: x(nd,ni*nj)

      integer(kind=wi) :: i,j,m

      do j=1,nj
         do i=1,ni
            m = (j-1)*ni + i
            x(1,m) = real(i-1,kind=wp)/real(ni-1,kind=wp)
            x(2,m) = real(j-1,kind=wp)/real(nj-1,kind=wp)
         enddo
      enddo

      return
   end function grid_coords
   !**********************************************************************
   function volumes(ni,nj,x) result(dv)
      implicit none

      integer(kind=wi), intent(in) :: ni,nj

      real(kind=wp), intent(in) :: x(:,:)

      real(kind=wp) :: dv((ni-1)*(nj-1))
      
      integer(kind=wi) :: i,j,m

      real(kind=wp) :: a(size(x,1)),b(size(x,1))
      
      do j=1,nj-1
         do i=1,ni-1
            m = (j-1)*(ni-1) + i
            
            a(:) = x(:,m + (ni-1) + 1) - x(:,m + 0)
            b(:) = x(:,m + (ni-1) + 0) - x(:,m + 1)

            dv(m) = abs(a(1)*b(2) - b(1)*a(2))
         enddo
      enddo
      
      return
   end function volumes
   !**********************************************************************
   function interior_faces(ni,nj,nfi) result(face_i)
      implicit none

      integer(kind=wi), intent(in) :: ni,nj,nfi

      integer(kind=wi) :: face_i(4,nfi)

      integer(kind=wi) :: i,j,nf

      nf = 0
      do j=1,nj
         do i=1,ni
            if((j .ne. 1) .and. (j .ne. nj) .and. (i .ne. ni))then
               nf = nf + 1
               face_i(1,nf) = (j - 1)*ni + i; 
               face_i(2,nf) = (j - 1)*ni + i + 1
               face_i(3,nf) = (j + 0)*(ni-1) + i ! Left volume 
               face_i(4,nf) = (j - 1)*(ni-1) + i ! Right volume
            endif 
            if((i .ne. 1) .and. (i .ne. ni) .and. (j .ne. nj))then
               nf = nf + 1
               face_i(1,nf) = (j - 1)*ni + i
               face_i(2,nf) = (j + 0)*ni + i               
               face_i(3,nf) = (j - 1)*(ni-1) + i - 1 ! Left volume 
               face_i(4,nf) = (j - 1)*(ni-1) + i     ! Right volume
            endif
         enddo
      enddo
      
   end function interior_faces
   !**********************************************************************
   function boundary_faces(ni,nj,nfb) result(face_b)
      implicit none
     
      integer(kind=wi), intent(in) :: ni,nj,nfb

      integer(kind=wi) :: face_b(4,nfb)

      integer(kind=wi) :: i,j,nf

      nf = 0
      do j=1,nj
         do i=1,ni
            if((j .eq. 1) .and. (i .ne. ni))then
               nf = nf + 1
               face_b(1,nf) = (j - 1)*ni + i; 
               face_b(2,nf) = (j - 1)*ni + i + 1
               face_b(3,nf) = (j + 0)*(ni-1) + i ! Left volume 
               face_b(4,nf) = nf ! Ghost boundary point
            endif
            if((j .eq. nj) .and. (i .ne. ni))then
               nf = nf + 1
               face_b(1,nf) = (j - 1)*ni + i; 
               face_b(2,nf) = (j - 1)*ni + i + 1
               face_b(3,nf) = (j - 1)*(ni-1) + i ! Left volume 
               face_b(4,nf) = nf ! Ghost boundary point
            endif
            if((i .eq. 1) .and. (j .ne. nj))then
               nf = nf + 1
               face_b(1,nf) = (j - 1)*ni + i
               face_b(2,nf) = (j + 0)*ni + i               
               face_b(3,nf) = (j - 1)*(ni-1) + i ! Left volume 
               face_b(4,nf) = nf ! Ghost boundary point
            endif
            if((i .eq. ni) .and. (j .ne. nj))then
               nf = nf + 1
               face_b(1,nf) = (j - 1)*ni + i
               face_b(2,nf) = (j + 0)*ni + i               
               face_b(3,nf) = (j - 1)*(ni-1) + i - 1 ! Left volume 
               face_b(4,nf) = nf ! Ghost boundary point
            endif
         enddo
      enddo
      
      return
   end function boundary_faces
   !**********************************************************************
   function normals(face,x) result(n)
      implicit none

      integer(kind=wi), intent(in) :: face(:,:)

      real(kind=wp), intent(in) :: x(:,:)

      real(kind=wp) :: n(size(x,1),size(face,2))
      
      integer(kind=wi) :: nf,ie

      real(kind=wp) :: xe(size(x,1)),nt(size(x,1)),l,a
      
      nf = size(face,2)

      do ie=1,nf
         xe = x(:,face(2,ie)) - x(:,face(1,ie))

         a = xe(2)/xe(1)
         nt(1) = sqrt(1e0_wp/(1e0_wp + a*a))
         nt(2) = -nt(1)/a
         if((xe(1)*nt(2) - xe(2)*nt(1)) .gt. 0e0)then
            nt(1) = -nt(1)
            nt(2) = -nt(2)
         endif

         l = sqrt(xe(1)*xe(1) + xe(2)*xe(2))
         
         n(1,ie) = nt(1)*l
         n(1,ie) = nt(2)*l
      enddo
            
      return
   end function normals
   !**********************************************************************
end module grid
