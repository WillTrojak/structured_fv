program main
   use grid
   use euler_fv
   implicit none

   type(domain) :: test

   type(fv_euler) :: test_spatial
   
   call test%init(2,3,3)

   call test_spatial%init(test,4,1.4e0_wp)
   
end program main
