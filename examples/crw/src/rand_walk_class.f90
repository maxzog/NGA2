!> Basic Lagrangian random walk solver class:
!> Provides support for Lagrangian-transported objects
module rand_walk_class
  use precision,      only: WP
  use string,         only: str_medium
  use config_class,   only: config
  use ddadi_class,    only: ddadi
  use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
  implicit none
  private

  ! Expose type/constructor/methods
  public :: lpt

  !> Memory adaptation parameter
  real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
  real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor

  real(WP), parameter :: C_poz=1.0_WP         !< Pozorski & Apte (2009) model constant

  !> I/O chunk size to read at a time
  integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing

  !> Basic particle object definition
  type :: part
     !> MPI_INTEGER8 data
     integer(kind=8) :: id                !< Particle ID
     !> MPI_DOUBLE_PRECISION data
     real(WP) :: d                        !< Particle diameter
     real(WP), dimension(3) :: pos        !< Particle center coordinates
     real(WP), dimension(3) :: vel        !< Velocity of particle
     real(WP), dimension(3) :: us         !< Fluid veloicty seen
     real(WP), dimension(3) :: dW         !< Wiener increment
     real(WP) :: tau_eddy                 !< Time step size for the particle
     real(WP) :: taup                     !< Time step size for the particle
     real(WP) :: dt                       !< Time step size for the particle
     !> MPI_INTEGER data
     integer , dimension(3) :: ind        !< Index of cell containing particle center
     integer  :: flag                     !< Control parameter (0=normal, 1=done->will be removed)
  end type part
  !> Number of blocks, block length, and block types in a particle
  integer, parameter                         :: part_nblock=3
  integer           , dimension(part_nblock) :: part_lblock=[1,16,4]
  type(MPI_Datatype), dimension(part_nblock) :: part_tblock=[MPI_INTEGER8,MPI_DOUBLE_PRECISION,MPI_INTEGER]
  !> MPI_PART derived datatype and size
  type(MPI_Datatype) :: MPI_PART
  integer :: MPI_PART_SIZE

  !> Lagrangian particle tracking solver object definition
  type :: lpt

     ! This is our underlying config
     class(config), pointer :: cfg                       !< This is the config the solver is build for

     type(ddadi) :: implicit                             !< Implicit solver for filtering

     ! This is the name of the solver
     character(len=str_medium) :: name='UNNAMED_LPT'     !< Solver name (default=UNNAMED_LPT)

     ! Particle data
     integer :: np                                       !< Global number of particles
     integer :: np_                                      !< Local number of particles
     integer, dimension(:), allocatable :: np_proc       !< Number of particles on each processor
     type(part), dimension(:), allocatable :: p          !< Array of particles of type part

     ! Overlap particle (i.e., ghost) data
     integer :: ng_                                      !< Local number of ghosts
     type(part), dimension(:), allocatable :: g          !< Array of ghosts of type part

     ! CFL numbers
     real(WP) :: CFLp_x,CFLp_y,CFLp_z                    !< CFL numbers

     ! Particle density
     real(WP) :: rho                                     !< Density of particle

     ! Gravitational acceleration
     real(WP), dimension(3) :: gravity=0.0_WP            !< Acceleration of gravity

     ! Solver parameters
     real(WP) :: nstep=1                                 !< Number of substeps (default=1)

     ! Injection parameters
     real(WP) :: mfr                                     !< Mass flow rate for particle injection
     real(WP), dimension(3) :: inj_pos                   !< Center location to inject particles
     real(WP), dimension(3) :: inj_vel                   !< Celocity assigned during injection
     real(WP) :: inj_dmean                               !< Mean diameter assigned during injection
     real(WP) :: inj_dsd                                 !< STD diameter assigned during injection
     real(WP) :: inj_dmin                                !< Min diameter assigned during injection
     real(WP) :: inj_dmax                                !< Max diameter assigned during injection
     real(WP) :: inj_dshift                              !< Diameter shift assigned during injection
     real(WP) :: inj_D=0.0_WP                            !< Diameter to inject particles within

     ! Monitoring info
     real(WP) :: dmin,dmax,dmean,dvar                    !< Diameter info
     real(WP) :: Umin,Umax,Umean,Uvar                    !< U velocity info
     real(WP) :: Vmin,Vmax,Vmean,Vvar                    !< V velocity info
     real(WP) :: Wmin,Wmax,Wmean,Wvar                    !< W velocity info
     real(WP) :: Usmin,Usmax,Usmean,Usvar                !< U velocity seen info
     real(WP) :: Vsmin,Vsmax,Vsmean,Vsvar                !< V velocity seen info
     real(WP) :: Wsmin,Wsmax,Wsmean,Wsvar                !< W velocity seen info
     integer  :: np_new,np_out                           !< Number of new and removed particles

     ! Particle volume fraction
     real(WP), dimension(:,:,:), allocatable :: VF       !< Particle volume fraction, cell-centered

     ! Filtering operation
     logical :: implicit_filter                          !< Solve implicitly
     real(WP) :: filter_width                            !< Characteristic filter width
     real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z !< Divergence operator
     real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z !< Gradient operator

   contains
     procedure :: update_partmesh                        !< Update a partmesh object using current particles
     procedure :: advance                                !< Step forward the particle ODEs - discrete random walk
     procedure :: advance_crw                            !< Step forward the particle ODEs - continuous random walk
     procedure :: advance_crw_fede                       !< Step forward the particle ODEs - continuous random walk with corrections
     procedure :: get_rhs                                !< Compute rhs of particle odes
     procedure :: get_diffusion                          !< Compute diffusion coefficients
     procedure :: get_drift                              !< Compute drift coefficients
     procedure :: resize                                 !< Resize particle array to given size
     procedure :: resize_ghost                           !< Resize ghost array to given size
     procedure :: recycle                                !< Recycle particle array by removing flagged particles
     procedure :: sync                                   !< Synchronize particles across interprocessor boundaries
     procedure :: share                                  !< Share particles across interprocessor boundaries
     procedure :: read                                   !< Parallel read particles from file
     procedure :: write                                  !< Parallel write particles to file
     procedure :: get_max                                !< Extract various monitoring data
     procedure :: get_cfl                                !< Calculate maximum CFL
     procedure :: filter                                 !< Apply volume filtering to field
     procedure :: inject                                 !< Inject particles at a prescribed boundary
  end type lpt


  !> Declare lpt solver constructor
  interface lpt
     procedure constructor
  end interface lpt

contains


  !> Default constructor for lpt solver
  function constructor(cfg,name) result(self)
    implicit none
    type(lpt) :: self
    class(config), target, intent(in) :: cfg
    character(len=*), optional :: name
    integer :: i,j,k,l

    ! Set the name for the solver
    if (present(name)) self%name=trim(adjustl(name))

    ! Point to pgrid object
    self%cfg=>cfg

    ! Allocate variables
    allocate(self%np_proc(1:self%cfg%nproc)); self%np_proc=0
    self%np_=0; self%np=0
    call self%resize(0)
    self%np_new=0; self%np_out=0

    ! Initialize MPI derived datatype for a particle
    call prepare_mpi_part()

    ! Allocate finite volume divergence operators
    allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
    ! Create divergence operator to cell center [xm,ym,zm]
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
             self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
             self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
          end do
       end do
    end do

    ! Allocate finite difference velocity gradient operators
    allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
    allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
    allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
    ! Create gradient coefficients to cell faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
             self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
             self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
          end do
       end do
    end do

!!    ! Allocate finite difference density interpolation coefficients
!!    allocate(self%itpr_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
!!    allocate(self%itpr_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
!!    allocate(self%itpr_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
!!    ! Create density interpolation coefficients to cell face
!!    do k=self%cfg%kmin_,self%cfg%kmax_+1
!!       do j=self%cfg%jmin_,self%cfg%jmax_+1
!!          do i=self%cfg%imin_,self%cfg%imax_+1
!!             self%itpr_x(:,i,j,k)=self%cfg%dxmi(i)*[self%cfg%xm(i)-self%cfg%x(i),self%cfg%x(i)-self%cfg%xm(i-1)] !< Linear interpolation in x from [xm,ym,zm] to [x,ym,zm]
!!             self%itpr_y(:,i,j,k)=self%cfg%dymi(j)*[self%cfg%ym(j)-self%cfg%y(j),self%cfg%y(j)-self%cfg%ym(j-1)] !< Linear interpolation in y from [xm,ym,zm] to [xm,y,zm]
!!             self%itpr_z(:,i,j,k)=self%cfg%dzmi(k)*[self%cfg%zm(k)-self%cfg%z(k),self%cfg%z(k)-self%cfg%zm(k-1)] !< Linear interpolation in z from [xm,ym,zm] to [xm,ym,z]
!!          end do
!!       end do
!!    end do
    
    ! Loop over the domain and zero divergence in walls
    do k=self%cfg%kmin_,self%cfg%kmax_
       do j=self%cfg%jmin_,self%cfg%jmax_
          do i=self%cfg%imin_,self%cfg%imax_
             if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                self%div_x(:,i,j,k)=0.0_WP
                self%div_y(:,i,j,k)=0.0_WP
                self%div_z(:,i,j,k)=0.0_WP
             end if
          end do
       end do
    end do

    ! Zero out gradient to wall faces
    do k=self%cfg%kmin_,self%cfg%kmax_+1
       do j=self%cfg%jmin_,self%cfg%jmax_+1
          do i=self%cfg%imin_,self%cfg%imax_+1
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
             if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
          end do
       end do
    end do

    ! Adjust metrics to account for lower dimensionality
    if (self%cfg%nx.eq.1) then
       self%div_x=0.0_WP
       self%grd_x=0.0_WP
    end if
    if (self%cfg%ny.eq.1) then
       self%div_y=0.0_WP
       self%grd_y=0.0_WP
    end if
    if (self%cfg%nz.eq.1) then
       self%div_z=0.0_WP
       self%grd_z=0.0_WP
    end if

    ! Create implicit solver object for filtering
    self%implicit=ddadi(cfg=self%cfg,name='Filter',nst=7)
    self%implicit%stc(1,:)=[ 0, 0, 0]
    self%implicit%stc(2,:)=[+1, 0, 0]
    self%implicit%stc(3,:)=[-1, 0, 0]
    self%implicit%stc(4,:)=[ 0,+1, 0]
    self%implicit%stc(5,:)=[ 0,-1, 0]
    self%implicit%stc(6,:)=[ 0, 0,+1]
    self%implicit%stc(7,:)=[ 0, 0,-1]
    call self%implicit%init()

    ! Set filter width to zero by default
    self%filter_width=0.0_WP

    ! Solve implicitly by default
    self%implicit_filter=.true.

!!    ! Generate a wall distance/norm function
!!    allocate(self%Wdist(  self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
!!    allocate(self%Wnorm(3,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
!!    ! First pass to set correct sign
!!    do k=self%cfg%kmino_,self%cfg%kmaxo_
!!       do j=self%cfg%jmino_,self%cfg%jmaxo_
!!          do i=self%cfg%imino_,self%cfg%imaxo_
!!             if (self%cfg%VF(i,j,k).eq.0.0_WP) then
!!                self%Wdist(i,j,k)=-sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
!!             else
!!                self%Wdist(i,j,k)=+sqrt(self%cfg%xL**2+self%cfg%yL**2+self%cfg%zL**2)
!!             end if
!!             self%Wnorm(:,i,j,k)=0.0_WP
!!          end do
!!       end do
!!    end do
!!    ! Second pass to compute local distance (get 2 closest cells)
!!    do l=1,2
!!       do k=self%cfg%kmino_,self%cfg%kmaxo_
!!          do j=self%cfg%jmino_,self%cfg%jmaxo_
!!             do i=self%cfg%imino_+l,self%cfg%imaxo_
!!                if (self%Wdist(i,j,k)*self%Wdist(i-l,j,k).lt.0.0_WP) then
!!                   ! There is a wall at x(i)
!!                   if (abs(self%cfg%xm(i  )-self%cfg%x(i-l+1)).lt.abs(self%Wdist(i  ,j,k))) then
!!                      self%Wdist(i  ,j,k)=sign(self%cfg%xm(i  )-self%cfg%x(i-l+1),self%Wdist(i  ,j,k))
!!                      self%Wnorm(:,i  ,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
!!                   end if
!!                   if (abs(self%cfg%xm(i-l)-self%cfg%x(i)).lt.abs(self%Wdist(i-l,j,k))) then
!!                      self%Wdist(i-l,j,k)=sign(self%cfg%xm(i-l)-self%cfg%x(i),self%Wdist(i-l,j,k))
!!                      self%Wnorm(:,i-l,j,k)=[self%cfg%VF(i,j,k)-self%cfg%VF(i-l,j,k),0.0_WP,0.0_WP]
!!                   end if
!!                end if
!!             end do
!!          end do
!!       end do
!!       do k=self%cfg%kmino_,self%cfg%kmaxo_
!!          do j=self%cfg%jmino_+l,self%cfg%jmaxo_
!!             do i=self%cfg%imino_,self%cfg%imaxo_
!!                if (self%Wdist(i,j,k)*self%Wdist(i,j-l,k).lt.0.0_WP) then
!!                   ! There is a wall at y(j)
!!                   if (abs(self%cfg%ym(j  )-self%cfg%y(j-l+1)).lt.abs(self%Wdist(i,j  ,k))) then
!!                      self%Wdist(i,j  ,k)=sign(self%cfg%ym(j  )-self%cfg%y(j-l+1),self%Wdist(i,j  ,k))
!!                      self%Wnorm(:,i,j  ,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
!!                   end if
!!                   if (abs(self%cfg%ym(j-l)-self%cfg%y(j)).lt.abs(self%Wdist(i,j-l,k))) then
!!                      self%Wdist(i,j-l,k)=sign(self%cfg%ym(j-l)-self%cfg%y(j),self%Wdist(i,j-l,k))
!!                      self%Wnorm(:,i,j-l,k)=[0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j-l,k),0.0_WP]
!!                   end if
!!                end if
!!             end do
!!          end do
!!       end do
!!       do k=self%cfg%kmino_+l,self%cfg%kmaxo_
!!          do j=self%cfg%jmino_,self%cfg%jmaxo_
!!             do i=self%cfg%imino_,self%cfg%imaxo_
!!                if (self%Wdist(i,j,k)*self%Wdist(i,j,k-l).lt.0.0_WP) then
!!                   ! There is a wall at z(k)
!!                   if (abs(self%cfg%zm(k  )-self%cfg%z(k-l+1)).lt.abs(self%Wdist(i,j,k  ))) then
!!                      self%Wdist(i,j,k  )=sign(self%cfg%zm(k  )-self%cfg%z(k-l+1),self%Wdist(i,j,k  ))
!!                      self%Wnorm(:,i,j,k  )=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
!!                   end if
!!                   if (abs(self%cfg%zm(k-l)-self%cfg%z(k)).lt.abs(self%Wdist(i,j,k-l))) then
!!                      self%Wdist(i,j,k-l)=sign(self%cfg%zm(k-l)-self%cfg%z(k),self%Wdist(i,j,k-l))
!!                      self%Wnorm(:,i,j,k-l)=[0.0_WP,0.0_WP,self%cfg%VF(i,j,k)-self%cfg%VF(i,j,k-l)]
!!                   end if
!!                end if
!!             end do
!!          end do
!!       end do
!!    end do
!!    call self%cfg%sync(self%Wdist)
!!    call self%cfg%sync(self%Wnorm)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (self%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end function constructor

  subroutine get_drift(this, p, rho, sgs_visc, U, V, W, a)
   use random, only: random_normal
   implicit none
   class(lpt), intent(inout) :: this
   type(part), intent(inout) :: p
   real(WP), intent(inout) :: a
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho           !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: sgs_visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(3) :: fvel
   real(WP) :: fvisc,frho,delta,tke_sgs,sig_sgs,tau_crwi

   !> Interpolate values
   frho =this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=rho     ,    bc='n')
   fvisc=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=sgs_visc,    bc='n')
   delta=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=this%cfg%vol,bc='n')
   fvel=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=U,V=V,W=W)
   delta=delta**(0.333333_WP)

   !> Compute a
   tke_sgs = (fvisc/frho/delta/0.067_WP)**2    ! SGS tke
   sig_sgs = sqrt(0.6666_WP*tke_sgs)           ! SGS velocity
   tau_crwi= sig_sgs/delta/C_poz               ! Inverse SGS timescale 
   a = tau_crwi

  end subroutine get_drift

  subroutine get_diffusion(this, p, rho, sgs_visc, U, V, W, b)
   use random, only: random_normal
   implicit none
   class(lpt), intent(inout) :: this
   type(part), intent(inout) :: p
   real(WP), intent(inout) :: b
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho           !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: sgs_visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W             !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:
   real(WP), dimension(3) :: fvel
   real(WP) :: fvisc,frho,delta,tke_sgs,sig_sgs,tau_crwi
   real(WP) :: Tl,Tcr,Le,cross

   !> Interpolate values
   frho =this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=rho     ,    bc='n')
   fvisc=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=sgs_visc,    bc='n')
   delta=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=this%cfg%vol,bc='n')
   fvel=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=U,V=V,W=W)
   delta=delta**(0.333333_WP)

   !> Compute b
   tke_sgs = (fvisc/frho/delta/0.067_WP)**2    ! SGS tke
   sig_sgs = sqrt(0.6666_WP*tke_sgs)           ! SGS velocity
   b = sig_sgs

   !> If tau_eddy has passed update dW
   p%tau_eddy = p%tau_eddy - p%dt
   if (p%tau_eddy.gt.0.0_WP) return
  
   tau_crwi= sig_sgs/delta/C_poz               ! Inverse SGS timescale
   if (tau_crwi.gt.0.0_WP) then
      Tl = 1.0_WP / tau_crwi                        ! Subgrid lagrangian timescale
      Le = 2.0_WP * Tl * sqrt(0.66667_WP * tke_sgs) ! eddy lengthscale
      cross = Le / (p%taup * norm2(fvel - p%vel))
      if (cross.lt.1.0_WP) then                     ! Check if crossing time can be computed
         Tcr = -p%taup * log(1.0_WP - cross)        ! eddy crossing time
      else
         Tcr = 2.0_WP * Tl                          ! Set to eddy lifetime otherwise
      end if
      p%tau_eddy = min(2.0_WP * Tl, Tcr)
      p%dW = [random_normal(m=0.0_WP,sd=1.0_WP), &
              random_normal(m=0.0_WP,sd=1.0_WP), &
              random_normal(m=0.0_WP,sd=1.0_WP)]
   else
      p%tau_eddy = 0.0_WP
   end if

  end subroutine get_diffusion
  

  !> Advance the particle equations by a specified time step dt
  !  using a discrete random walk model
  !> p%id=0 => no coll, no solve
  !> p%id=-1=> no coll, no move
  subroutine advance(this,dt,U,V,W,rho,visc,sgs_visc)
    use mpi_f08, only : MPI_SUM,MPI_INTEGER
    use mathtools, only: Pi
    use random, only: random_normal
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: sgs_visc  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(3) :: acc
    integer :: i,ierr
    real(WP) :: mydt,dt_done
    real(WP) :: b 
    type(part) :: myp,pold

    ! Zero out number of particles removed
    this%np_out=0

    ! Advance the equations
    do i=1,this%np_
       ! Avoid particles with id=0
       if (this%p(i)%id.eq.0) cycle
       ! Create local copy of particle
       myp=this%p(i)
       ! Time-integrate until dt_done=dt
       dt_done=0.0_WP
       do while (dt_done.lt.dt)
          ! Decide the timestep size
          mydt=min(myp%dt,dt-dt_done)
          ! Remember the particle
          pold=myp

          ! Advance with Euler prediction
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+0.5_WP*mydt*myp%vel
          myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity)

          ! Get diffusion coefficient
          call this%get_diffusion(p=myp,rho=rho,sgs_visc=sgs_visc,U=U,V=V,W=W,b=b)
         
          ! Correct with midpoint rule
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+mydt*myp%vel
          myp%vel=pold%vel+mydt*(acc+this%gravity)

          myp%Us(1)=b*myp%dW(1)
          myp%Us(2)=b*myp%dW(2)
          myp%Us(3)=b*myp%dW(3)

          ! Relocalize
          myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
          ! Increment
          dt_done=dt_done+mydt
       end do
       ! Correct the position to take into account periodicity
       if (this%cfg%xper) myp%pos(1)=this%cfg%x(this%cfg%imin)+modulo(myp%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
       if (this%cfg%yper) myp%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(myp%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
       if (this%cfg%zper) myp%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(myp%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
       ! Handle particles that have left the domain
       if (myp%pos(1).lt.this%cfg%x(this%cfg%imin).or.myp%pos(1).gt.this%cfg%x(this%cfg%imax+1)) myp%flag=1
       if (myp%pos(2).lt.this%cfg%y(this%cfg%jmin).or.myp%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) myp%flag=1
       if (myp%pos(3).lt.this%cfg%z(this%cfg%kmin).or.myp%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) myp%flag=1
       ! Relocalize the particle
       myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
       ! Count number of particles removed
       if (myp%flag.eq.1) this%np_out=this%np_out+1
       ! Copy back to particle
       if (myp%id.ne.-1) this%p(i)=myp
    end do

    ! Communicate particles
    call this%sync()

    ! Sum up particles removed
    call MPI_ALLREDUCE(this%np_out,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=i

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine advance

  !> Advance the particle equations by a specified time step dt
  !  using a continuous random walk model
  !> p%id=0 => no coll, no solve
  !> p%id=-1=> no coll, no move
  subroutine advance_crw(this,dt,U,V,W,rho,visc,sgs_visc)
    use mpi_f08, only : MPI_SUM,MPI_INTEGER
    use mathtools, only: Pi
    use random, only: random_normal
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: sgs_visc  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(3) :: acc
    integer :: i,ierr
    real(WP) :: rmydt,mydt,dt_done
    real(WP) :: a,b    ! Coefficients of OU process 
    type(part) :: myp,pold

    ! Zero out number of particles removed
    this%np_out=0

    ! Advance the equations
    do i=1,this%np_
       ! Avoid particles with id=0
       if (this%p(i)%id.eq.0) cycle
       ! Create local copy of particle
       myp=this%p(i)
       ! Time-integrate until dt_done=dt
       dt_done=0.0_WP
       do while (dt_done.lt.dt)
          ! Decide the timestep size
          mydt=min(myp%dt,dt-dt_done)
          rmydt=sqrt(mydt)
          ! Remember the particle
          pold=myp

          ! Advance with Euler prediction
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+0.5_WP*mydt*myp%vel
          myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity)

          ! Get diffusion coefficient
          call this%get_diffusion(p=myp,rho=rho,sgs_visc=sgs_visc,U=U,V=V,W=W,b=b)
          call this%get_drift(p=myp,rho=rho,sgs_visc=sgs_visc,U=U,V=V,W=W,a=a)
         
          ! Compute Wiener increment
          myp%dW = [random_normal(m=0.0_WP,sd=1.0_WP), &
                    random_normal(m=0.0_WP,sd=1.0_WP), &
                    random_normal(m=0.0_WP,sd=1.0_WP)]
         
          ! Correct with midpoint rule
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+mydt*myp%vel
          myp%vel=pold%vel+mydt*(acc+this%gravity)

          myp%Us(1)=(1.0_WP - a*mydt)*myp%Us(1) + b*myp%dW(1)*rmydt
          myp%Us(2)=(1.0_WP - a*mydt)*myp%Us(1) + b*myp%dW(2)*rmydt
          myp%Us(3)=(1.0_WP - a*mydt)*myp%Us(1) + b*myp%dW(3)*rmydt

          ! Relocalize
          myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
          ! Increment
          dt_done=dt_done+mydt
       end do
       ! Correct the position to take into account periodicity
       if (this%cfg%xper) myp%pos(1)=this%cfg%x(this%cfg%imin)+modulo(myp%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
       if (this%cfg%yper) myp%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(myp%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
       if (this%cfg%zper) myp%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(myp%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
       ! Handle particles that have left the domain
       if (myp%pos(1).lt.this%cfg%x(this%cfg%imin).or.myp%pos(1).gt.this%cfg%x(this%cfg%imax+1)) myp%flag=1
       if (myp%pos(2).lt.this%cfg%y(this%cfg%jmin).or.myp%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) myp%flag=1
       if (myp%pos(3).lt.this%cfg%z(this%cfg%kmin).or.myp%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) myp%flag=1
       ! Relocalize the particle
       myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
       ! Count number of particles removed
       if (myp%flag.eq.1) this%np_out=this%np_out+1
       ! Copy back to particle
       if (myp%id.ne.-1) this%p(i)=myp
    end do

    ! Communicate particles
    call this%sync()

    ! Sum up particles removed
    call MPI_ALLREDUCE(this%np_out,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=i

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine advance_crw

  !> Advance the particle equations by a specified time step dt
  !  using a continuous random walk model with well-mixed preserving drift
  !> p%id=0 => no coll, no solve
  !> p%id=-1=> no coll, no move
  subroutine advance_crw_fede(this,dt,U,V,W,rho,visc,sgs_visc,gradu,dtaurdx,dtaurdy,dtaurdz)
    use mpi_f08, only : MPI_SUM,MPI_INTEGER
    use mathtools, only: Pi
    use random, only: random_normal
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(inout) :: dt  !< Timestep size over which to advance
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: sgs_visc  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dtaurdx   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dtaurdy   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: dtaurdz   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(1:3,1:3,this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: gradu !< Velocity gradient
    real(WP), dimension(3) :: acc
    integer :: i,ierr
    real(WP) :: rmydt,mydt,dt_done
    real(WP) :: a,b    ! Coefficients of OU process 
    type(part) :: myp,pold

    ! Zero out number of particles removed
    this%np_out=0

    ! Advance the equations
    do i=1,this%np_
       ! Avoid particles with id=0
       if (this%p(i)%id.eq.0) cycle
       ! Create local copy of particle
       myp=this%p(i)
       ! Time-integrate until dt_done=dt
       dt_done=0.0_WP
       do while (dt_done.lt.dt)
          ! Decide the timestep size
          mydt=min(myp%dt,dt-dt_done)
          rmydt=sqrt(mydt)
          ! Remember the particle
          pold=myp

          ! Advance with Euler prediction
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+0.5_WP*mydt*myp%vel
          myp%vel=pold%vel+0.5_WP*mydt*(acc+this%gravity)

          ! Get diffusion coefficient
          call this%get_diffusion(p=myp,rho=rho,sgs_visc=sgs_visc,U=U,V=V,W=W,b=b)
          call this%get_drift(p=myp,rho=rho,sgs_visc=sgs_visc,U=U,V=V,W=W,a=a)
         
          ! Compute Wiener increment
          myp%dW = [random_normal(m=0.0_WP,sd=1.0_WP), &
                    random_normal(m=0.0_WP,sd=1.0_WP), &
                    random_normal(m=0.0_WP,sd=1.0_WP)]
         
          ! Correct with midpoint rule
          call this%get_rhs(U=U,V=V,W=W,rho=rho,visc=visc,p=myp,acc=acc,opt_dt=myp%dt)
          myp%pos=pold%pos+mydt*myp%vel
          myp%vel=pold%vel+mydt*(acc+this%gravity)

          !> Interpolate divergence of SGS stress tensor to particle location
          tau = this%cfg%get_velocity(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),U=dtaurdx,V=dtaurdy,W=dtaurdz)    

          !> Interpolate the velocity gradient tensor to the particle location
          gu11 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(1,1,:,:,:),bc='n')
          gu12 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(1,2,:,:,:),bc='n')
          gu13 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(1,3,:,:,:),bc='n')
          gu21 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(2,1,:,:,:),bc='n')
          gu22 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(2,2,:,:,:),bc='n')
          gu23 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(2,3,:,:,:),bc='n')
          gu31 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(3,1,:,:,:),bc='n')
          gu32 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(3,2,:,:,:),bc='n')
          gu33 = this%cfg%get_scalar(pos=myp%pos,i0=myp%ind(1),j0=myp%ind(2),k0=myp%ind(3),S=gradu(3,3,:,:,:),bc='n')
          
          gux = [gu11,gu21,gu31]
          guy = [gu12,gu22,gu32]
          guz = [gu13,gu23,gu33]

          myp%Us(1) = (-dot_product(myp%Us(:),gux) + tau(1))*mydt + (1.0_WP - a*mydt)*myp%Us(1) + b*myp%dW(1)*rmydt
          myp%Us(2) = (-dot_product(myp%Us(:),guy) + tau(2))*mydt + (1.0_WP - a*mydt)*myp%Us(2) + b*myp%dW(2)*rmydt
          myp%Us(3) = (-dot_product(myp%Us(:),guz) + tau(3))*mydt + (1.0_WP - a*mydt)*myp%Us(3) + b*myp%dW(3)*rmydt

          ! Relocalize
          myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
          ! Increment
          dt_done=dt_done+mydt
       end do
       ! Correct the position to take into account periodicity
       if (this%cfg%xper) myp%pos(1)=this%cfg%x(this%cfg%imin)+modulo(myp%pos(1)-this%cfg%x(this%cfg%imin),this%cfg%xL)
       if (this%cfg%yper) myp%pos(2)=this%cfg%y(this%cfg%jmin)+modulo(myp%pos(2)-this%cfg%y(this%cfg%jmin),this%cfg%yL)
       if (this%cfg%zper) myp%pos(3)=this%cfg%z(this%cfg%kmin)+modulo(myp%pos(3)-this%cfg%z(this%cfg%kmin),this%cfg%zL)
       ! Handle particles that have left the domain
       if (myp%pos(1).lt.this%cfg%x(this%cfg%imin).or.myp%pos(1).gt.this%cfg%x(this%cfg%imax+1)) myp%flag=1
       if (myp%pos(2).lt.this%cfg%y(this%cfg%jmin).or.myp%pos(2).gt.this%cfg%y(this%cfg%jmax+1)) myp%flag=1
       if (myp%pos(3).lt.this%cfg%z(this%cfg%kmin).or.myp%pos(3).gt.this%cfg%z(this%cfg%kmax+1)) myp%flag=1
       ! Relocalize the particle
       myp%ind=this%cfg%get_ijk_global(myp%pos,myp%ind)
       ! Count number of particles removed
       if (myp%flag.eq.1) this%np_out=this%np_out+1
       ! Copy back to particle
       if (myp%id.ne.-1) this%p(i)=myp
    end do

    ! Communicate particles
    call this%sync()

    ! Sum up particles removed
    call MPI_ALLREDUCE(this%np_out,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_out=i

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
         if (verbose.gt.1) write(output_unit,'(a)') trim(message)
         if (verbose.gt.0) call log(message)
      end if
    end block logging

  end subroutine advance_crw_fede

  !> Calculate RHS of the particle ODEs
  subroutine get_rhs(this,U,V,W,rho,visc,p,acc,opt_dt)
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: rho       !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: visc      !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    type(part), intent(inout) :: p
    real(WP), dimension(3), intent(out) :: acc
    real(WP), intent(out) :: opt_dt
    real(WP) :: fvisc,frho
    real(WP), dimension(3) :: fvel

    ! Interpolate fluid quantities to particle location
    interpolate: block
      ! Interpolate the fluid phase velocity to the particle location
      fvel=this%cfg%get_velocity(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),U=U,V=V,W=W)
      fvel=fvel+p%Us
      ! Interpolate the fluid phase viscosity to the particle location
      fvisc=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=visc,bc='n')
      fvisc=fvisc+epsilon(1.0_WP)
      ! Interpolate the fluid phase density to the particle location
      frho=this%cfg%get_scalar(pos=p%pos,i0=p%ind(1),j0=p%ind(2),k0=p%ind(3),S=rho,bc='n')
    end block interpolate

    ! Compute acceleration due to drag
    compute_drag: block
      real(WP) :: Re,tau,corr
      ! Particle Reynolds number
      Re=frho*norm2(p%vel-fvel)*p%d/fvisc+epsilon(1.0_WP)
      ! Schiller-Naumann correction
      corr=1.0_WP+0.15_WP*Re**(0.687_WP)
      ! Particle response time
      tau=this%rho*p%d**2/(18.0_WP*fvisc*corr)
      p%taup=tau
      ! Return acceleration and optimal timestep size
      acc=(fvel-p%vel)/tau
      opt_dt=tau/real(this%nstep,WP)
    end block compute_drag

  end subroutine get_rhs
  

  !> Laplacian filtering operation
  subroutine filter(this,A)
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP) :: filter_coeff
    integer :: i,j,k,n,nstep
    real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

    ! Return without filtering if filter width is zero
    if (this%filter_width.le.0.0_WP) return

    ! Recompute filter coefficient
    filter_coeff=max(this%filter_width**2-this%cfg%min_meshsize**2,0.0_WP)/(16.0_WP*log(2.0_WP))
    if (filter_coeff.le.0.0_WP) return

    if (this%implicit_filter) then  !< Apply filter implicitly
       if (.not.this%implicit%setup_done) then
          ! Prepare diffusive operator (only need to do this once)
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   this%implicit%opr(1,i,j,k)=1.0_WP-(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x(-1,i+1,j,k)+&
                   &                                  this%div_x( 0,i,j,k)*filter_coeff*this%grd_x( 0,i  ,j,k)+&
                   &                                  this%div_y(+1,i,j,k)*filter_coeff*this%grd_y(-1,i,j+1,k)+&
                   &                                  this%div_y( 0,i,j,k)*filter_coeff*this%grd_y( 0,i,j  ,k)+&
                   &                                  this%div_z(+1,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k+1)+&
                   &                                  this%div_z( 0,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k  ))
                   this%implicit%opr(2,i,j,k)=      -(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x( 0,i+1,j,k))
                   this%implicit%opr(3,i,j,k)=      -(this%div_x( 0,i,j,k)*filter_coeff*this%grd_x(-1,i  ,j,k))
                   this%implicit%opr(4,i,j,k)=      -(this%div_y(+1,i,j,k)*filter_coeff*this%grd_y( 0,i,j+1,k))
                   this%implicit%opr(5,i,j,k)=      -(this%div_y( 0,i,j,k)*filter_coeff*this%grd_y(-1,i,j  ,k))
                   this%implicit%opr(6,i,j,k)=      -(this%div_z(+1,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k+1))
                   this%implicit%opr(7,i,j,k)=      -(this%div_z( 0,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k  ))
                end do
             end do
          end do
       end if
       ! Solve the linear system
       call this%implicit%setup()
       this%implicit%rhs=A
       this%implicit%sol=0.0_WP
       call this%implicit%solve()
       A=this%implicit%sol
       call this%cfg%sync(A)
       
    else  !< Apply filter explicitly
       ! Allocate flux arrays
       allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
       nstep=ceiling(6.0_WP*filter_coeff/this%cfg%min_meshsize**2)
       filter_coeff=filter_coeff/real(nstep,WP)
       do n=1,nstep
          ! Diffusive flux of A
          do k=this%cfg%kmin_,this%cfg%kmax_+1
             do j=this%cfg%jmin_,this%cfg%jmax_+1
                do i=this%cfg%imin_,this%cfg%imax_+1
                   FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
                   FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
                   FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
                end do
             end do
          end do
          ! Divergence of fluxes
          do k=this%cfg%kmin_,this%cfg%kmax_
             do j=this%cfg%jmin_,this%cfg%jmax_
                do i=this%cfg%imin_,this%cfg%imax_
                   A(i,j,k)=A(i,j,k)+sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
                end do
             end do
          end do
          ! Sync A
          call this%cfg%sync(A)
       end do
       ! Deallocate flux arrays
       deallocate(FX,FY,FZ)
    end if

  end subroutine filter


  !> Inject particles from a prescribed location with given mass flowrate
  !> Requires injection parameters to be set beforehand
  subroutine inject(this,dt,U,V,W)
    use mpi_f08
    use parallel, only: MPI_REAL_WP
    use mathtools, only: Pi
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: U         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: V         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: W         !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
    real(WP), intent(inout) :: dt                  !< Timestep size over which to advance
    real(WP) :: inj_min(3),inj_max(3)              !< Min/max extents of injection
    real(WP) :: Mgoal,Madded,Mtmp,buf              !< Mass flow rate parameters
    real(WP), save :: previous_error=0.0_WP        !< Store mass left over from previous timestep
    integer(kind=8) :: maxid_,maxid                !< Keep track of maximum particle id
    integer :: i,j,np0_,np2,np_tmp,count,ierr
    
    ! Initial number of particles
    np0_=this%np_
    this%np_new=0

    ! Get the particle mass that should be added to the system
    Mgoal  = this%mfr*dt+previous_error
    Madded = 0.0_WP

    ! Determine id to assign to particle
    maxid_=0
    do i=1,this%np_
       maxid_=max(maxid_,this%p(i)%id)
    end do
    call MPI_ALLREDUCE(maxid_,maxid,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr)

    ! Add new particles until desired mass is achieved
    do while (Madded.lt.Mgoal)

       if (this%cfg%amRoot) then
          ! Initialize parameters
          Mtmp = 0.0_WP
          np_tmp = 0
          ! Loop while the added volume is not sufficient
          do while (Mtmp.lt.Mgoal-Madded)
             ! Increment counter
             np_tmp=np_tmp+1
             count = np0_+np_tmp
             ! Create space for new particle
             call this%resize(count)
             ! Generate a diameter
             this%p(count)%d=get_diameter()
             ! Set various parameters for the particle
             this%p(count)%id    =maxid+int(np_tmp,8)
             this%p(count)%dt    =0.0_WP
             this%p(count)%Us    =0.0_WP
             this%p(count)%dW    =0.0_WP
             this%p(count)%tau_eddy=0.0_WP
             this%p(count)%taup  =0.0_WP
             ! Give a position at the injector to the particle
             this%p(count)%pos=get_position()
             ! Localize the particle
             this%p(count)%ind(1)=this%cfg%imin; this%p(count)%ind(2)=this%cfg%jmin; this%p(count)%ind(3)=this%cfg%kmin
             this%p(count)%ind=this%cfg%get_ijk_global(this%p(count)%pos,this%p(count)%ind)
             ! Give it a velocity
             this%p(count)%vel=this%cfg%get_velocity(pos=this%p(count)%pos,i0=this%p(count)%ind(1),j0=this%p(count)%ind(2), & 
                                                     k0=this%p(count)%ind(3),U=U,V=V,W=W)
             ! Make it an "official" particle
             this%p(count)%flag=0
             ! Update the added mass for the timestep
             Mtmp = Mtmp + this%rho*Pi/6.0_WP*this%p(count)%d**3
          end do
       end if
       ! Communicate particles
       call this%sync()
       ! Loop through newly created particles
       buf=0.0_WP
       do i=np0_+1,this%np_
          ! Remove if out of bounds
          if (this%cfg%VF(this%p(i)%ind(1),this%p(i)%ind(2),this%p(i)%ind(3)).le.0.0_WP) this%p(i)%flag=1
          if (this%p(i)%flag.eq.0) then
             ! Update the added mass for the timestep
             buf = buf + this%rho*Pi/6.0_WP*this%p(i)%d**3
             ! Update the max particle id
             maxid = max(maxid,this%p(i)%id)
             ! Increment counter
             this%np_new=this%np_new+1
          end if
       end do
       ! Total mass added
       call MPI_ALLREDUCE(buf,Mtmp,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); Madded=Madded+Mtmp
       ! Clean up particles
       call this%recycle()
       ! Update initial npart
       np0_=this%np_
       ! Maximum particle id
       call MPI_ALLREDUCE(maxid,maxid_,1,MPI_INTEGER8,MPI_MAX,this%cfg%comm,ierr); maxid=maxid_
    end do

    ! Remember the error
    previous_error = Mgoal-Madded

    ! Sum up injected particles
    call MPI_ALLREDUCE(this%np_new,i,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr); this%np_new=i


  contains

    ! Compute particle diameter
    function get_diameter() result(dp)
      use random, only: random_lognormal
      implicit none
      real(WP) :: dp
      dp=random_lognormal(m=this%inj_dmean-this%inj_dshift,sd=this%inj_dsd)+this%inj_dshift
      do while (dp.gt.this%inj_dmax+epsilon(1.0_WP).or.dp.lt.this%inj_dmin-epsilon(1.0_WP))
         dp=random_lognormal(m=this%inj_dmean-this%inj_dshift,sd=this%inj_dsd)+this%inj_dshift
      end do
    end function get_diameter

    ! Position for bulk injection of particles
    function get_position() result(pos)
      use random, only: random_uniform
      use mathtools, only: twoPi
      implicit none
      real(WP), dimension(3) :: pos
      real(WP) :: rand,r,theta
      ! Set x position
      pos(1) = this%inj_pos(1)
      ! Set in y & z
      if (this%inj_D.gt.0.0_WP) then
         ! Random y & z position within a circular region
         if (this%cfg%nz.eq.1) then
            pos(2)=random_uniform(lo=this%inj_pos(2)-0.5_WP*this%inj_D,hi=this%inj_pos(3)+0.5_WP*this%inj_D)
            pos(3) = this%cfg%zm(this%cfg%kmin_)
         else
            rand=random_uniform(lo=0.0_WP,hi=1.0_WP)
            r=0.5_WP*this%inj_D*sqrt(rand) !< sqrt(rand) avoids accumulation near the center
            call random_number(rand)
            theta=random_uniform(lo=0.0_WP,hi=twoPi)
            pos(2) = this%inj_pos(2)+r*sin(theta)
            pos(3) = this%inj_pos(3)+r*cos(theta)
         end if
      else
         ! Random y & z position across domain width
         pos(2)=random_uniform(lo=this%cfg%y(this%cfg%jmin),hi=this%cfg%y(this%cfg%jmax+1))
         if (this%cfg%nz.eq.1) then
            pos(3) = this%cfg%zm(this%cfg%kmin_)
         else
            pos(3)=random_uniform(lo=this%cfg%z(this%cfg%kmin),hi=this%cfg%z(this%cfg%kmax+1))
         end if
      end if
      
    end function get_position

  end subroutine inject
  
  
  !> Calculate the CFL
  subroutine get_cfl(this,dt,cflc,cfl)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
    use parallel, only: MPI_REAL_WP
    implicit none
    class(lpt), intent(inout) :: this
    real(WP), intent(in)  :: dt
    real(WP), intent(out) :: cflc
    real(WP), optional :: cfl
    integer :: i,ierr
    real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,my_CFL_col

    ! Set the CFLs to zero
    my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP; my_CFL_col=0.0_WP
    do i=1,this%np_
       my_CFLp_x=max(my_CFLp_x,abs(this%p(i)%vel(1))*this%cfg%dxi(this%p(i)%ind(1)))
       my_CFLp_y=max(my_CFLp_y,abs(this%p(i)%vel(2))*this%cfg%dyi(this%p(i)%ind(2)))
       my_CFLp_z=max(my_CFLp_z,abs(this%p(i)%vel(3))*this%cfg%dzi(this%p(i)%ind(3)))
       my_CFL_col=max(my_CFL_col,sqrt(sum(this%p(i)%vel**2))/this%p(i)%d)
    end do
    my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt

    ! Get the parallel max
    call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
    call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

    ! Return the maximum convective CFL
    cflc=max(this%CFLp_x,this%CFLp_y,this%CFLp_z)
  end subroutine get_cfl

  !> Extract various monitoring data from particle field
  subroutine get_max(this)
    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM
    use parallel, only: MPI_REAL_WP
    implicit none
    class(lpt), intent(inout) :: this
    real(WP) :: buf,safe_np
    integer :: i,j,k,ierr

    ! Create safe np
    safe_np=real(max(this%np,1),WP)

    ! Diameter and velocity min/max/mean
    this%dmin=huge(1.0_WP); this%dmax=-huge(1.0_WP); this%dmean=0.0_WP
    this%Umin=huge(1.0_WP); this%Umax=-huge(1.0_WP); this%Umean=0.0_WP
    this%Vmin=huge(1.0_WP); this%Vmax=-huge(1.0_WP); this%Vmean=0.0_WP
    this%Wmin=huge(1.0_WP); this%Wmax=-huge(1.0_WP); this%Wmean=0.0_WP
    this%Usmin=huge(1.0_WP); this%Usmax=-huge(1.0_WP); this%Usmean=0.0_WP
    this%Vsmin=huge(1.0_WP); this%Vsmax=-huge(1.0_WP); this%Vsmean=0.0_WP
    this%Wsmin=huge(1.0_WP); this%Wsmax=-huge(1.0_WP); this%Wsmean=0.0_WP
    do i=1,this%np_
       this%dmin=min(this%dmin,this%p(i)%d     ); this%dmax=max(this%dmax,this%p(i)%d     ); this%dmean=this%dmean+this%p(i)%d
       this%Umin=min(this%Umin,this%p(i)%vel(1)); this%Umax=max(this%Umax,this%p(i)%vel(1)); this%Umean=this%Umean+this%p(i)%vel(1)
       this%Vmin=min(this%Vmin,this%p(i)%vel(2)); this%Vmax=max(this%Vmax,this%p(i)%vel(2)); this%Vmean=this%Vmean+this%p(i)%vel(2)
       this%Wmin=min(this%Wmin,this%p(i)%vel(3)); this%Wmax=max(this%Wmax,this%p(i)%vel(3)); this%Wmean=this%Wmean+this%p(i)%vel(3)
       this%Usmin=min(this%Usmin,this%p(i)%us(1)); this%Usmax=max(this%Usmax,this%p(i)%us(1)); this%Usmean=this%Usmean+this%p(i)%us(1)
       this%Vsmin=min(this%Vsmin,this%p(i)%us(2)); this%Vsmax=max(this%Vsmax,this%p(i)%us(2)); this%Vsmean=this%Vsmean+this%p(i)%us(2)
       this%Wsmin=min(this%Wsmin,this%p(i)%us(3)); this%Wsmax=max(this%Wsmax,this%p(i)%us(3)); this%Wsmean=this%Wsmean+this%p(i)%us(3)
    end do
    call MPI_ALLREDUCE(this%dmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%dmin =buf
    call MPI_ALLREDUCE(this%dmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%dmax =buf
    call MPI_ALLREDUCE(this%dmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dmean=buf/safe_np
    call MPI_ALLREDUCE(this%Umin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Umin =buf
    call MPI_ALLREDUCE(this%Umax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Umax =buf
    call MPI_ALLREDUCE(this%Umean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Umean=buf/safe_np
    call MPI_ALLREDUCE(this%Vmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vmin =buf
    call MPI_ALLREDUCE(this%Vmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vmax =buf
    call MPI_ALLREDUCE(this%Vmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vmean=buf/safe_np
    call MPI_ALLREDUCE(this%Wmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wmin =buf
    call MPI_ALLREDUCE(this%Wmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wmax =buf
    call MPI_ALLREDUCE(this%Wmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wmean=buf/safe_np
    call MPI_ALLREDUCE(this%Usmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Usmin =buf
    call MPI_ALLREDUCE(this%Usmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Usmax =buf
    call MPI_ALLREDUCE(this%Usmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Usmean=buf/safe_np
    call MPI_ALLREDUCE(this%Vsmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Vsmin =buf
    call MPI_ALLREDUCE(this%Vsmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Vsmax =buf
    call MPI_ALLREDUCE(this%Vsmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vsmean=buf/safe_np
    call MPI_ALLREDUCE(this%Wsmin ,buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%Wsmin =buf
    call MPI_ALLREDUCE(this%Wsmax ,buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%Wsmax =buf
    call MPI_ALLREDUCE(this%Wsmean,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wsmean=buf/safe_np

    ! Diameter and velocity variance
    this%dvar=0.0_WP
    this%Uvar=0.0_WP
    this%Vvar=0.0_WP
    this%Wvar=0.0_WP
    this%Usvar=0.0_WP
    this%Vsvar=0.0_WP
    this%Wsvar=0.0_WP
    do i=1,this%np_
       this%dvar=this%dvar+(this%p(i)%d     -this%dmean)**2.0_WP
       this%Uvar=this%Uvar+(this%p(i)%vel(1)-this%Umean)**2.0_WP
       this%Vvar=this%Vvar+(this%p(i)%vel(2)-this%Vmean)**2.0_WP
       this%Wvar=this%Wvar+(this%p(i)%vel(3)-this%Wmean)**2.0_WP
       this%Usvar=this%Usvar+(this%p(i)%Us(1)-this%Usmean)**2.0_WP
       this%Vsvar=this%Vsvar+(this%p(i)%Us(2)-this%Vsmean)**2.0_WP
       this%Wsvar=this%Wsvar+(this%p(i)%Us(3)-this%Wsmean)**2.0_WP
    end do
    call MPI_ALLREDUCE(this%dvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%dvar=buf/safe_np
    call MPI_ALLREDUCE(this%Uvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Uvar=buf/safe_np
    call MPI_ALLREDUCE(this%Vvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vvar=buf/safe_np
    call MPI_ALLREDUCE(this%Wvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wvar=buf/safe_np
    call MPI_ALLREDUCE(this%Usvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Usvar=buf/safe_np
    call MPI_ALLREDUCE(this%Vsvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Vsvar=buf/safe_np
    call MPI_ALLREDUCE(this%Wsvar,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); this%Wsvar=buf/safe_np

  end subroutine get_max


  !> Update particle mesh using our current particles
  subroutine update_partmesh(this,pmesh)
    use partmesh_class, only: partmesh
    implicit none
    class(lpt), intent(inout) :: this
    class(partmesh), intent(inout) :: pmesh
    integer :: i
    ! Reset particle mesh storage
    call pmesh%reset()
    if (this%np_.gt.0) then
      ! Copy particle info
      call pmesh%set_size(this%np_)
      do i=1,this%np_
         pmesh%pos(:,i)=this%p(i)%pos
      end do
    end if
    ! Root adds a particle if there are none
    if (this%np.eq.0.and.this%cfg%amRoot) then
      call pmesh%set_size(1)
      pmesh%pos(1,1)=this%cfg%x(this%cfg%imin)
      pmesh%pos(2,1)=this%cfg%y(this%cfg%jmin)
      pmesh%pos(3,1)=this%cfg%z(this%cfg%kmin)
    end if
  end subroutine update_partmesh


  !> Creation of the MPI datatype for particle
  subroutine prepare_mpi_part()
    use mpi_f08
    use messager, only: die
    implicit none
    integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
    integer(MPI_ADDRESS_KIND) :: lb,extent
    type(MPI_Datatype) :: MPI_PART_TMP
    integer :: i,mysize,ierr
    ! Prepare the displacement array
    disp(1)=0
    do i=2,part_nblock
       call MPI_Type_size(part_tblock(i-1),mysize,ierr)
       disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(part_lblock(i-1),MPI_ADDRESS_KIND)
    end do
    ! Create and commit the new type
    call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
    call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
    call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
    call MPI_Type_commit(MPI_PART,ierr)
    ! If a problem was encountered, say it
    if (ierr.ne.0) call die('[lpt prepare_mpi_part] MPI Particle type creation failed')
    ! Get the size of this type
    call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
  end subroutine prepare_mpi_part


  !> Synchronize particle arrays across processors
  subroutine sync(this)
    use mpi_f08
    implicit none
    class(lpt), intent(inout) :: this
    integer, dimension(0:this%cfg%nproc-1) :: nsend_proc,nrecv_proc
    integer, dimension(0:this%cfg%nproc-1) :: nsend_disp,nrecv_disp
    integer :: n,prank,ierr
    type(part), dimension(:), allocatable :: buf_send
    ! Recycle first to minimize communication load
    call this%recycle()
    ! Prepare information about what to send
    nsend_proc=0
    do n=1,this%np_
       prank=this%cfg%get_rank(this%p(n)%ind)
       nsend_proc(prank)=nsend_proc(prank)+1
    end do
    nsend_proc(this%cfg%rank)=0
    ! Inform processors of what they will receive
    call MPI_ALLtoALL(nsend_proc,1,MPI_INTEGER,nrecv_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    ! Prepare displacements for all-to-all
    nsend_disp(0)=0
    nrecv_disp(0)=this%np_   !< Directly add particles at the end of main array
    do n=1,this%cfg%nproc-1
       nsend_disp(n)=nsend_disp(n-1)+nsend_proc(n-1)
       nrecv_disp(n)=nrecv_disp(n-1)+nrecv_proc(n-1)
    end do
    ! Allocate buffer to send particles
    allocate(buf_send(sum(nsend_proc)))
    ! Pack the particles in the send buffer
    nsend_proc=0
    do n=1,this%np_
       ! Get the rank
       prank=this%cfg%get_rank(this%p(n)%ind)
       ! Skip particles still inside
       if (prank.eq.this%cfg%rank) cycle
       ! Pack up for sending
       nsend_proc(prank)=nsend_proc(prank)+1
       buf_send(nsend_disp(prank)+nsend_proc(prank))=this%p(n)
       ! Flag particle for removal
       this%p(n)%flag=1
    end do
    ! Allocate buffer for receiving particles
    call this%resize(this%np_+sum(nrecv_proc))
    ! Perform communication
    call MPI_ALLtoALLv(buf_send,nsend_proc,nsend_disp,MPI_PART,this%p,nrecv_proc,nrecv_disp,MPI_PART,this%cfg%comm,ierr)
    ! Deallocate buffer
    deallocate(buf_send)
    ! Recycle to remove duplicate particles
    call this%recycle()
  end subroutine sync


   !> Share particles across processor boundaries
   subroutine share(this,nover)
      use mpi_f08
      use messager, only: warn,die
      implicit none
      class(lpt), intent(inout) :: this
      integer, optional :: nover
      type(part), dimension(:), allocatable :: tosend
      type(part), dimension(:), allocatable :: torecv
      integer :: n,no,nsend,nrecv
      type(MPI_Status) :: status
      integer :: isrc,idst,ierr
      
      ! Check overlap size
      if (present(nover)) then
         no=nover
         if (no.gt.this%cfg%no) then
            call warn('[rand_walk_class share] Specified overlap is larger than that of cfg - reducing no')
            no=this%cfg%no
         else if (no.le.0) then
            call die('[rand_walk_class share] Specified overlap cannot be less or equal to zero')
         end if
      else
         no=1
      end if
      
      ! Clean up ghost array
      call this%resize_ghost(n=0); this%ng_=0
      
      ! Share ghost particles in -x (no ghosts are sent here)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).lt.this%cfg%imin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).lt.this%cfg%imin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).lt.this%cfg%imin+no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)+this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)+this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles in +x (no ghosts are sent here)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).gt.this%cfg%imax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(1).gt.this%cfg%imax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%xper.and.tosend(nsend)%ind(1).gt.this%cfg%imax-no) then
               tosend(nsend)%pos(1)=tosend(nsend)%pos(1)-this%cfg%xL
               tosend(nsend)%ind(1)=tosend(nsend)%ind(1)-this%cfg%nx
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,0,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles in -y (ghosts need to be sent now)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
            end if
         end if
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(2).lt.this%cfg%jmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).lt.this%cfg%jmin+no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)+this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)+this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles in +y (ghosts need to be sent now - but not newly received ghosts!)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
            end if
         end if
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(2).gt.this%cfg%jmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%yper.and.tosend(nsend)%ind(2).gt.this%cfg%jmax-no) then
               tosend(nsend)%pos(2)=tosend(nsend)%pos(2)-this%cfg%yL
               tosend(nsend)%ind(2)=tosend(nsend)%ind(2)-this%cfg%ny
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,1,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles in -z (ghosts need to be sent now)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).lt.this%cfg%kmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
            end if
         end if
      end do
      do n=1,this%ng_
         if (this%g(n)%ind(3).lt.this%cfg%kmin_+no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).lt.this%cfg%kmin+no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)+this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)+this%cfg%nz
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,2,-1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
      ! Share ghost particles in +z (ghosts need to be sent now - but not newly received ghosts!)
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) nsend=nsend+1
      end do
      allocate(tosend(nsend))
      nsend=0
      do n=1,this%np_
         if (this%p(n)%ind(3).gt.this%cfg%kmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%p(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
            end if
         end if
      end do
      do n=1,this%ng_-nrecv
         if (this%g(n)%ind(3).gt.this%cfg%kmax_-no) then
            nsend=nsend+1
            tosend(nsend)=this%g(n)
            if (this%cfg%zper.and.tosend(nsend)%ind(3).gt.this%cfg%kmax-no) then
               tosend(nsend)%pos(3)=tosend(nsend)%pos(3)-this%cfg%zL
               tosend(nsend)%ind(3)=tosend(nsend)%ind(3)-this%cfg%nz
            end if
         end if
      end do
      nrecv=0
      call MPI_CART_SHIFT(this%cfg%comm,2,+1,isrc,idst,ierr)
      call MPI_SENDRECV(nsend,1,MPI_INTEGER,idst,0,nrecv,1,MPI_INTEGER,isrc,0,this%cfg%comm,status,ierr)
      allocate(torecv(nrecv))
      call MPI_SENDRECV(tosend,nsend,MPI_PART,idst,0,torecv,nrecv,MPI_PART,isrc,0,this%cfg%comm,status,ierr)
      call this%resize_ghost(this%ng_+nrecv)
      this%g(this%ng_+1:this%ng_+nrecv)=torecv
      this%ng_=this%ng_+nrecv
      if (allocated(tosend)) deallocate(tosend)
      if (allocated(torecv)) deallocate(torecv)
      
   end subroutine share
   
   
   !> Adaptation of particle array size
   subroutine resize(this,n)
      implicit none
      class(lpt), intent(inout) :: this
      integer, intent(in) :: n
      type(part), dimension(:), allocatable :: tmp
      integer :: size_now,size_new
      ! Resize particle array to size n
      if (.not.allocated(this%p)) then
         ! Allocate directly to size n
         allocate(this%p(n))
         this%p(1:n)%flag=1
      else
         ! Update from a non-zero size to another non-zero size
         size_now=size(this%p,dim=1)
         if (n.gt.size_now) then
            size_new=max(n,int(real(size_now,WP)*coeff_up))
            allocate(tmp(size_new))
            tmp(1:size_now)=this%p
            tmp(size_now+1:)%flag=1
            call move_alloc(tmp,this%p)
         else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
            allocate(tmp(n))
            tmp(1:n)=this%p(1:n)
            call move_alloc(tmp,this%p)
         end if
      end if
   end subroutine resize


  !> Adaptation of ghost array size
  subroutine resize_ghost(this,n)
    implicit none
    class(lpt), intent(inout) :: this
    integer, intent(in) :: n
    type(part), dimension(:), allocatable :: tmp
    integer :: size_now,size_new
    ! Resize ghost array to size n
    if (.not.allocated(this%g)) then
       ! Allocate directly to size n
       allocate(this%g(n))
       this%g(1:n)%flag=1
    else
       ! Update from a non-zero size to another non-zero size
       size_now=size(this%g,dim=1)
       if (n.gt.size_now) then
          size_new=max(n,int(real(size_now,WP)*coeff_up))
          allocate(tmp(size_new))
          tmp(1:size_now)=this%g
          tmp(size_now+1:)%flag=1
          call move_alloc(tmp,this%g)
       else if (n.lt.int(real(size_now,WP)*coeff_dn)) then
          allocate(tmp(n))
          tmp(1:n)=this%g(1:n)
          call move_alloc(tmp,this%g)
       end if
    end if
  end subroutine resize_ghost


  !> Clean-up of particle array by removing flag=1 particles
  subroutine recycle(this)
    implicit none
    class(lpt), intent(inout) :: this
    integer :: new_size,i,ierr
    ! Compact all active particles at the beginning of the array
    new_size=0
    if (allocated(this%p)) then
       do i=1,size(this%p,dim=1)
          if (this%p(i)%flag.ne.1) then
             new_size=new_size+1
             if (i.ne.new_size) then
                this%p(new_size)=this%p(i)
                this%p(i)%flag=1
             end if
          end if
       end do
    end if
    ! Resize to new size
    call this%resize(new_size)
    ! Update number of particles
    this%np_=new_size
    call MPI_ALLGATHER(this%np_,1,MPI_INTEGER,this%np_proc,1,MPI_INTEGER,this%cfg%comm,ierr)
    this%np=sum(this%np_proc)
  end subroutine recycle


  !> Parallel write particles to file
  subroutine write(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(lpt), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: i,ierr,iunit

    ! Root serial-writes the file header
    if (this%cfg%amRoot) then
       ! Open the file
       open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
       if (ierr.ne.0) call die('[lpt write] Problem encountered while serial-opening data file: '//trim(filename))
       ! Number of particles and particle object size
       write(iunit) this%np,MPI_PART_SIZE
       ! Done with the header
       close(iunit)
    end if

    ! The rest is done in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[lpt write] Problem encountered while parallel-opening data file: '//trim(filename))

    ! Get current position
    call MPI_FILE_GET_POSITION(ifile,offset,ierr)

    ! Compute the offset and write
    do i=1,this%cfg%rank
       offset=offset+int(this%np_proc(i),MPI_OFFSET_KIND)*int(MPI_PART_SIZE,MPI_OFFSET_KIND)
    end do
    if (this%np_.gt.0) call MPI_FILE_WRITE_AT(ifile,offset,this%p,this%np_,MPI_PART,status,ierr)

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Wrote ",i0," particles to file [",a,"] on partitioned grid [",a,"]")') this%np,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine write


  !> Parallel read particles to file
  subroutine read(this,filename)
    use mpi_f08
    use messager, only: die
    use parallel, only: info_mpiio
    implicit none
    class(lpt), intent(inout) :: this
    character(len=*), intent(in) :: filename
    type(MPI_File) :: ifile
    type(MPI_Status):: status
    integer(kind=MPI_OFFSET_KIND) :: offset,header_offset
    integer :: i,j,ierr,npadd,psize,nchunk,cnt
    integer, dimension(:,:), allocatable :: ppp

    ! First open the file in parallel
    call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
    if (ierr.ne.0) call die('[lpt read] Problem encountered while reading data file: '//trim(filename))

    ! Read file header first
    call MPI_FILE_READ_ALL(ifile,npadd,1,MPI_INTEGER,status,ierr)
    call MPI_FILE_READ_ALL(ifile,psize,1,MPI_INTEGER,status,ierr)

    ! Remember current position
    call MPI_FILE_GET_POSITION(ifile,header_offset,ierr)

    ! Check compatibility of particle type
    if (psize.ne.MPI_PART_SIZE) call die('[lpt read] Particle type unreadable')

    ! Naively share reading task among all processors
    nchunk=int(npadd/(this%cfg%nproc*part_chunk_size))+1
    allocate(ppp(this%cfg%nproc,nchunk))
    ppp=int(npadd/(this%cfg%nproc*nchunk))
    cnt=0
    out:do j=1,nchunk
       do i=1,this%cfg%nproc
          cnt=cnt+1
          if (cnt.gt.mod(npadd,this%cfg%nproc*nchunk)) exit out
          ppp(i,j)=ppp(i,j)+1
       end do
    end do out

    ! Read by chunk
    do j=1,nchunk
       ! Find offset
       offset=header_offset+int(MPI_PART_SIZE,MPI_OFFSET_KIND)*int(sum(ppp(1:this%cfg%rank,:))+sum(ppp(this%cfg%rank+1,1:j-1)),MPI_OFFSET_KIND)
       ! Resize particle array
       call this%resize(this%np_+ppp(this%cfg%rank+1,j))
       ! Read this file
       call MPI_FILE_READ_AT(ifile,offset,this%p(this%np_+1:this%np_+ppp(this%cfg%rank+1,j)),ppp(this%cfg%rank+1,j),MPI_PART,status,ierr)
       ! Most general case: relocate every droplet
       do i=this%np_+1,this%np_+ppp(this%cfg%rank+1,j)
          this%p(i)%ind=this%cfg%get_ijk_global(this%p(i)%pos,this%p(i)%ind)
       end do
       ! Exchange all that
       call this%sync()
    end do

    ! Close the file
    call MPI_FILE_CLOSE(ifile,ierr)

    ! Log/screen output
    logging: block
      use, intrinsic :: iso_fortran_env, only: output_unit
      use param,    only: verbose
      use messager, only: log
      use string,   only: str_long
      character(len=str_long) :: message
      if (this%cfg%amRoot) then
         write(message,'("Read ",i0," particles from file [",a,"] on partitioned grid [",a,"]")') npadd,trim(filename),trim(this%cfg%name)
         if (verbose.gt.2) write(output_unit,'(a)') trim(message)
         if (verbose.gt.1) call log(message)
      end if
    end block logging

  end subroutine read


end module rand_walk_class