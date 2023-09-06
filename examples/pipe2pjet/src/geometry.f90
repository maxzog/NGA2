!> Various definitions and tools for initializing NGA2 config
module pjet_geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: pjet_geometry_init
   
contains
    
   !> Initialization of problem geometry
   subroutine pjet_geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use parallel,    only: group
      implicit none
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,nx,ny,nz
         logical :: stretch_x,stretch_y,stretch_z
         real(WP) :: Lx,Ly,Lz
         real(WP) :: D
         real(WP) :: xtilde,ytilde,ztilde
         real(WP) :: alpha,minAlpha,maxAlpha,r
         real(WP), dimension(:), allocatable :: x
         real(WP), dimension(:), allocatable :: y
         real(WP), dimension(:), allocatable :: z
         integer, dimension(3) :: partition
         ! Read in grid definition
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         call param_read('Lz',Lz)
         call param_read('nx',nx)
         call param_read('ny',ny)
         call param_read('nz',nz)
         call param_read('Pipe diameter', D)

         call param_read('Stretch x', stretch_x, default=.false.)
         call param_read('Stretch y', stretch_y, default=.false.)
         call param_read('Stretch z', stretch_z, default=.false.)

         call param_read('Stretch factor', r, default=2.0_WP)
         allocate(x(nx+1))
         allocate(y(ny+1))
         allocate(z(nz+1))
         ! Create simple rectilinear grid
         do i=1,nx+1
            if (stretch_x) then
              xtilde = 2.0_WP * real(i-1,WP) / real(nx-1,WP) - 1.0_WP
               alpha = 0.5_WP * (1.0_WP+sinh(r*xtilde))
               minAlpha = 0.5_WP * (1.0_WP + sinh(-r))
               maxAlpha = 0.5_WP * (1.0_WP + sinh(r))
               alpha = (alpha-minAlpha)/(maxAlpha-minAlpha)
               x(i) = Lx * alpha + D - !!!! TODO
            else
               x(i)=real(i-1,WP)/real(nx,WP)*Lx
            end if
         end do
         do i=1,ny+1
            if (stretch_y) then
               ytilde = 2.0_WP * real(i-1,WP) / real(ny-1,WP) - 1.0_WP
               alpha = 0.5_WP * (1.0_WP+sinh(r*ytilde))
               minAlpha = 0.5_WP * (1.0_WP + sinh(-r))
               maxAlpha = 0.5_WP * (1.0_WP + sinh(r))
               alpha = (alpha-minAlpha)/(maxAlpha-minAlpha)
               y(i) = Ly * alpha - 0.5_WP * Ly
            else
               y(i)=real(i-1,WP)/real(ny,WP)*Ly - 0.5_WP * Ly
            end if
         end do
         do i=1,nz+1
            if (stretch_z) then
               ztilde = 2.0_WP * real(i-1,WP) / real(nz-1,WP) - 1.0_WP
               alpha = 0.5_WP * (1.0_WP+sinh(r*ztilde))
               minAlpha = 0.5_WP * (1.0_WP + sinh(-r))
               maxAlpha = 0.5_WP * (1.0_WP + sinh(r))
               alpha = (alpha-minAlpha)/(maxAlpha-minAlpha)
               z(i) = Lz * alpha - 0.5_WP * Lz
            else
               z(i)=real(i-1,WP)/real(nz,WP)*Lz - 0.5_WP * Lz
            end if
         end do
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='jet')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         cfg%VF=1.0_WP
      end block create_grid
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      ! Create masks for this config
      create_walls: block
         cfg%VF=1.0_WP
      end block create_walls
      
   end subroutine pjet_geometry_init
   
   
end module pjet_geometry
