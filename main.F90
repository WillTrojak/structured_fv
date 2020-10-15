program main
   use grid
   use euler_fv
   use gmsh_mesh
   use erk_int
   implicit none

   type(domain) :: mesh

   type(fv_euler) :: space

   type(mesh_gmsh) :: test_mesh

   type(int_erk) :: time

   integer(kind=wi) :: it
   
   call mesh%init(2,3,3)
   
   call space%init(mesh,4,1.4e0_wp)
   
   call time%init('rk3',4,mesh%nv)

   call space%init_sol()
   
   time%dt = 1e-3_wp

   do it=1,10      
      call space%residual()
      call time%integrate(space%omega%dv,space%df,space%q)
      
      do while (.not. time%stage_complete)
         call space%residual()
         call time%integrate(space%omega%dv,space%df,space%q)
      enddo
   enddo
   
end program main
