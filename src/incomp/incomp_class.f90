!> Incompressible flow solver class:
!> Provides support for various BC, RHS calculation,
!> implicit solver, and pressure solution
!> Assumes constant viscosity and density.
module incomp_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: incomp
   
   !> Incompressible solver object definition
   type :: incomp
      
      ! This is our config
      class(config), pointer :: cfg                       !< This is the config the solver is build for
      
      ! This is the name of the solver
      character(len=str_medium) :: name='UNNAMED_INCOMP'  !< Solver name (default=UNNAMED_INCOMP)
      
      ! Constant property fluid
      real(WP) :: rho                                     !< This is our constant fluid density
      real(WP) :: visc                                    !< These is our constant fluid dynamic viscosity
      
      ! Boundary condition list
      !type(bcond), pointer :: bc=NULL()                   !<
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: itp_u_xm,itp_u_y ,itp_u_z   !< Interpolation for U
      real(WP), dimension(:,:,:,:), allocatable :: itp_v_x ,itp_v_ym,itp_v_z   !< Interpolation for V
      real(WP), dimension(:,:,:,:), allocatable :: itp_w_x ,itp_w_y ,itp_w_zm  !< Interpolation for W
      real(WP), dimension(:,:,:,:), allocatable :: divu_x  ,divu_y  ,divu_z    !< Divergence for U-cell
      real(WP), dimension(:,:,:,:), allocatable :: divv_x  ,divv_y  ,divv_z    !< Divergence for V-cell
      real(WP), dimension(:,:,:,:), allocatable :: divp_x  ,divp_y  ,divp_z    !< Divergence for P-cell
      real(WP), dimension(:,:,:,:), allocatable :: grad_u_x,grad_u_y,grad_u_z  !< Velocity gradient for U
      real(WP), dimension(:,:,:,:), allocatable :: grad_v_x,grad_v_y,grad_v_z  !< Velocity gradient for V
      real(WP), dimension(:,:,:,:), allocatable :: grad_w_x,grad_w_y,grad_w_z  !< Velocity gradient for W
      
      ! Flow variables
      real(WP), dimension(:,:,:), allocatable :: U        !< U velocity array
      real(WP), dimension(:,:,:), allocatable :: V        !< V velocity array
      real(WP), dimension(:,:,:), allocatable :: W        !< W velocity array
      real(WP), dimension(:,:,:), allocatable :: P        !< Pressure array
      
      ! MAYBE SOME WORK ARRAYS TOO?
      
   contains
      procedure :: print=>incomp_print                    !< Output solver to the screen
   end type incomp
   
   
   !> Declare incompressible solver constructor
   interface incomp
      procedure constructor
   end interface incomp
   
contains
   
   
   !> Default constructor for incompressible flow solver
   function constructor(cfg,name) result(self)
      implicit none
      type(incomp) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: name
      
      ! Set the name for the iterator
      self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg
      
      ! Prepare metrics
      
      
      ! Allocate flow variables
      allocate(self%U(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%U=0.0_WP
      allocate(self%V(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%V=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%W=0.0_WP
      allocate(self%P(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%P=0.0_WP
      
   end function constructor
   
   
   
   !> Print out info for incompressible flow solver
   subroutine incomp_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(incomp), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Incompressible solver [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
         write(output_unit,'(" >   density = ",es12.5)') this%rho
         write(output_unit,'(" > viscosity = ",es12.5)') this%visc
      end if
      
   end subroutine incomp_print
   
   
end module incomp_class