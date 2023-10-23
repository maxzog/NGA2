!> Various definitions and tools for running an NGA2 simulation
module pjet_class
   use precision,         only: WP
   use parallel,          only: rank, amRoot, group, comm
   use hypre_str_class,   only: hypre_str
   use fft3d_class,       only: fft3d
   use fft2d_class,       only: fft2d
   use config_class,      only: config
   use ibconfig_class,    only: ibconfig
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use datafile_class,    only: datafile
   use string,            only: str_medium
   use coupler_class,     only: coupler
   use mpi_f08,           only: MPI_GROUP, MPI_UNDEFINED
   implicit none
   private

   public :: pjetdomain
   
   type :: pjetdomain
      character(len=str_medium) :: desc
      type(config),      public :: cfg
      type(incomp),      public :: fs
      type(hypre_str),   public :: ps
      type(ddadi),       public :: vs
      
      type(monitor)     :: cflfile,mfile,simfile 
      type(timetracker) :: time
      type(datafile)    :: df
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out
      type(event)    :: ens_evt
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
      real(WP), dimension(:,:,:,:), allocatable :: vort
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
      real(WP) :: visc,mfr_target,mfr,bforce,cfl
      
      integer :: pipeRank
      logical :: inpipeGrp
   contains
      procedure :: init
      procedure :: geometry_init
      procedure :: setup_monitors
      procedure :: write_monitors
      procedure :: step
      procedure :: final
   end type pjetdomain
   
   !> Jet parameters
   real(WP) :: Djet        ! Jet diameter
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: stk,tau_eta,visc
   logical  :: use_sgs
   integer  :: sgs_type

 contains

   !> Function that localizes the jet at -x
   function jet_loc(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     real(WP) :: radius
     isIn = .false.
     ! Jet in yz plane
     radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
     if (radius.le.0.5_WP*Djet.and.i.eq.pg%imin) isIn=.true.
     !if (i.eq.pg%imin) isIn=.true.
   end function jet_loc

  !> Function that localizes the coflow at -x
   function coflow_loc(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     real(WP) :: radius
     isIn=.false.
     ! Coflow in yz plane
     radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
     if (radius.gt.0.5_WP*Djet.and.i.eq.pg%imin) isIn=.true.
   end function coflow_loc

   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain
  
   !> Function that localizes the bottom (y-) of the domain
   function bottom_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn 
     isIn=.false.
     if (j.eq.pg%jmin) isIn=.true.
   end function bottom_of_domain

   !> Function that localizes the top (y+) of the domain
   function top_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (j.eq.pg%jmax+1) isIn=.true.
   end function top_of_domain

   !> Function that localizes the front (z+) of the domain
   function front_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmax+1) isIn=.true.
   end function front_of_domain
   
   !> Function that localizes the back (z-) of the domain
   function back_of_domain(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn
     isIn=.false.
     if (k.eq.pg%kmin) isIn=.true.
   end function back_of_domain
   
   !> Monitor init
   subroutine  setup_monitors(this)
      implicit none
      class(pjetdomain), intent(inout) :: this

      ! Create simulation monitor
      this%mfile = monitor(this%fs%cfg%amRoot, trim(this%desc) // '_simulation')
      call this%mfile%add_column(this%time%n,'Timestep number')
      call this%mfile%add_column(this%time%t,'Time')
      call this%mfile%add_column(this%time%dt,'Timestep size')
      call this%mfile%add_column(this%time%cfl,'Maximum CFL')
      call this%mfile%add_column(this%fs%Umax,'Umax')
      call this%mfile%add_column(this%fs%Vmax,'Vmax')
      call this%mfile%add_column(this%fs%Wmax,'Wmax')
      call this%mfile%add_column(this%fs%Pmax,'Pmax')
      call this%mfile%add_column(this%fs%divmax,'Maximum divergence')
      call this%mfile%add_column(this%fs%psolv%it,'Pressure iteration')
      call this%mfile%add_column(this%fs%psolv%rerr,'Pressure error')
      call this%mfile%write()

      ! Create CFL monitor
      this%cflfile = monitor(this%fs%cfg%amRoot, trim(this%desc) // '_cfl')
      call this%cflfile%add_column(this%time%n,'Timestep number')
      call this%cflfile%add_column(this%time%t,'Time')
      call this%cflfile%add_column(this%fs%CFLc_x,'Convective xCFL')
      call this%cflfile%add_column(this%fs%CFLc_y,'Convective yCFL')
      call this%cflfile%add_column(this%fs%CFLc_z,'Convective zCFL')
      call this%cflfile%add_column(this%fs%CFLv_x,'Viscous xCFL')
      call this%cflfile%add_column(this%fs%CFLv_y,'Viscous yCFL')
      call this%cflfile%add_column(this%fs%CFLv_z,'Viscous zCFL')
      call this%cflfile%write()
   end subroutine setup_monitors

   subroutine write_monitors(this)
      implicit none
      class(pjetdomain), intent(inout) :: this
  
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
   end subroutine write_monitors

   !> Initializes pjet domain stretched grid
   subroutine geometry_init(this,pjetpart)
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use parallel,    only: group
      implicit none
      class(pjetdomain), intent(inout) :: this
      type(sgrid) :: grid
      integer, dimension(3), intent(in) :: pjetpart
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,nx,ny,nz
         logical :: stretch_x,stretch_y,stretch_z
         real(WP) :: Lx,Ly,Lz
         real(WP) :: xtilde,ytilde,ztilde
         real(WP) :: alpha,minAlpha,maxAlpha,r
         real(WP), dimension(:), allocatable :: x
         real(WP), dimension(:), allocatable :: y
         real(WP), dimension(:), allocatable :: z
         ! Read in grid definition
         call param_read('Jet diameter',Djet,default=0.1_WP)
         call param_read('Jet Lx',Lx)
         call param_read('Jet Ly',Ly)
         call param_read('Jet Lz',Lz)
         call param_read('Jet nx',nx)
         call param_read('Jet ny',ny)
         call param_read('Jet nz',nz)

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
               x(i) = Lx * alpha
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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='jet')
         deallocate(x, y, z)
         ! Create partitioned grid
         this%cfg=config(grp=group,decomp=pjetpart,grid=grid)
         this%cfg%VF=1.0_WP
      end block create_grid
   end subroutine geometry_init
   
   !> Initialization of pjet solver
   subroutine init(this)
      use param, only: param_read
      use messager, only: log
      implicit none
      class(pjetdomain), intent(inout) :: this

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker

      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: smg,pcg_pfmg,pfmg,gmres_smg
         use incomp_class,    only: clipped_neumann, dirichlet, slip
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='NS solver')
         ! Add BCs
         call this%fs%add_bcond(name='jet',   type=dirichlet,locator=jet_loc,face='x',dir=-1,canCorrect=.false. )
         call this%fs%add_bcond(name='coflow',type=dirichlet,locator=coflow_loc,face='x',dir=-1,canCorrect=.false. )
         call this%fs%add_bcond(name='right', type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true.)
         call this%fs%add_bcond(name='bottom',type=slip, locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call this%fs%add_bcond(name='top',   type=slip, locator=top_of_domain,face='y',dir=+1,canCorrect=.false.)
         call this%fs%add_bcond(name='front', type=slip, locator=front_of_domain,face='z',dir=+1,canCorrect=.false.)
         call this%fs%add_bcond(name='back',  type=slip, locator=back_of_domain,face='z',dir=-1,canCorrect=.false.)
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); this%fs%visc=visc
         ! Assign constant density
         call param_read('Density',this%fs%rho)
         ! Prepare and configure solvers
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         this%ps=hypre_str(cfg=this%cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',this%ps%maxit)
         call param_read('Pressure tolerance',this%ps%rcvg)
         call param_read('Implicit iteration', this%vs%maxit)
         call param_read('Implicit tolerance', this%vs%rcvg)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps, implicit_solver=this%vs)
      end block create_and_initialize_flow_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%vort(1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays

      ! Prepare initial velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         use random,   only: random_normal
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         use mathtools,only: Pi
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: r

         this%fs%U = 0.0_WP
         this%fs%V = 0.0_WP
         this%fs%W = 0.0_WP

         call this%fs%get_mfr()
         call this%fs%correct_mfr()
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         call this%fs%get_vorticity(this%vort)

         this%desc = 'pjet'
      end block initialize_velocity
      
      create_ensight: block
         this%ens_out=ensight(cfg=this%cfg,name='pjet')
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period', this%ens_evt%tper)
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('divergence',this%fs%div)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight     

   end subroutine init

   !> Time integrate our problem
   subroutine step(this)
      use parallel,       only: parallel_time
      implicit none
      class(pjetdomain), intent(inout) :: this
      integer :: it
      
      call this%fs%get_cfl(this%time%dt,this%time%cfl)
      call this%time%adjust_dt()
      call this%time%increment()

      ! Remember old velocity
      this%fs%Uold=this%fs%U
      this%fs%Vold=this%fs%V
      this%fs%Wold=this%fs%W
      
      ! Perform sub-iterations
      do it=1,this%time%itmax
         ! Build mid-time velocity
         this%fs%U=0.5_WP*(this%fs%U+this%fs%Uold)
         this%fs%V=0.5_WP*(this%fs%V+this%fs%Vold)
         this%fs%W=0.5_WP*(this%fs%W+this%fs%Wold)
         
         ! Explicit calculation of drho*u/dt from NS
         call this%fs%get_dmomdt(this%resU,this%resV,this%resW)

         ! Assemble explicit residual
         this%resU=-2.0_WP*(this%fs%rho*this%fs%U-this%fs%rho*this%fs%Uold)+this%time%dt*this%resU
         this%resV=-2.0_WP*(this%fs%rho*this%fs%V-this%fs%rho*this%fs%Vold)+this%time%dt*this%resV
         this%resW=-2.0_WP*(this%fs%rho*this%fs%W-this%fs%rho*this%fs%Wold)+this%time%dt*this%resW

         ! Implicit velocity solve
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)

         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW

         call this%fs%apply_bcond(this%time%t,this%time%dt)

         ! Solve Poisson equation
         call this%fs%get_mfr()
         call this%fs%correct_mfr()
         call this%fs%get_div()
         this%fs%psolv%rhs=-this%fs%cfg%vol*this%fs%div*this%fs%rho/this%time%dt
         this%fs%psolv%sol=0.0_WP
         call this%fs%psolv%solve()
         call this%fs%shift_p(this%fs%psolv%sol)
         
         ! Correct velocity
         call this%fs%get_pgrad(this%fs%psolv%sol,this%resU,this%resV,this%resW)
         this%fs%P=this%fs%P+this%fs%psolv%sol
         this%fs%U=this%fs%U-this%time%dt*this%resU/this%fs%rho
         this%fs%V=this%fs%V-this%time%dt*this%resV/this%fs%rho
         this%fs%W=this%fs%W-this%time%dt*this%resW/this%fs%rho
      end do

      call this%fs%get_div()
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      
      call this%write_monitors()
      call this%ens_out%write_data(this%time%t)
   end subroutine step
   
   !> Finalize the NGA2 simulation
   subroutine final(this)
      implicit none
      class(pjetdomain), intent(inout) :: this
      
   end subroutine final
   
end module pjet_class
