!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
!   use geometry,          only: cfg
   use parallel,          only: rank, amRoot
   use hypre_str_class,   only: hypre_str
   use fft3d_class,       only: fft3d
   use config_class,      only: config
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
   use flowrate,          only: qdat,tdat,dt_dat
   implicit none
   private
   
   !> TODO: add some more stuff here
   public :: simulation_init,simulation_run,simulation_final
   type(timetracker) :: time
   type(event)    :: ens_evt
   type(datafile) :: df

   ! General domain type definition
   type :: gendomain
     character(len=str_medium) :: desc
     type(config) :: cfg
     type(incomp) :: fs
     real(WP) :: cfl, meanU, meanV, meanW
     real(WP), dimension(:,:,:), allocatable :: resU, resV, resW, Ui, Vi, Wi
     real(WP), dimension(:,:,:,:), allocatable :: SR
     real(WP), dimension(:,:,:,:), allocatable :: vort
     real(WP), dimension(:,:,:,:,:), allocatable :: gradu
     type(monitor)  :: mfile, cflfile, simfile
     type(ensight)  :: ens_out     
   end type gendomain

   type, extends(gendomain) :: pipedomain
     type(fft3d) :: ps
     type(hypre_str) :: vs
     real(WP), dimension(:,:,:), allocatable :: G
     real(WP) :: visc,mfr_target,mfr,bforce
     real(WP) :: Qbulk
   end type pipedomain

   type, extends(gendomain) :: pjetdomain
     type(hypre_str)   :: ps
     type(ddadi)       :: vs
   end type pjetdomain

   type(pjetdomain) :: pjet
   type(pipedomain) :: pipe

   type(coupler) :: dom_cpl
    
   !> Fluid, forcing, and particle parameters
   real(WP) :: visc
   real(WP) :: stk,tau_eta
   logical  :: use_sgs
   integer  :: sgs_type

   !> Jet parameters
   real(WP) :: Djet        ! Jet diameter
   real(WP) :: Qt,Q0       ! Volumetric flow rates (time, max)
   real(WP) :: U0          ! Bulk velocity 
   real(WP) :: Uco         ! Coflow velocity 
   real(WP) :: omegaj,tauj ! Pulse frequency and decay rate
   real(WP) :: theta       ! Jet inlet momentum thickness

 contains
   
   !> Compute massflow rate
   function get_bodyforce_mfr(srcU) result(mfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(pipe%fs%cfg%imino_:,pipe%fs%cfg%jmino_:,pipe%fs%cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol,mfr
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=pipe%fs%cfg%kmin_,pipe%fs%cfg%kmax_
            do j=pipe%fs%cfg%jmin_,pipe%fs%cfg%jmax_
               do i=pipe%fs%cfg%imin_,pipe%fs%cfg%imax_
                  vol=pipe%fs%cfg%dxm(i)*pipe%fs%cfg%dy(j)*pipe%fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  !myRhoU=myRhoU+vol*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))*fs%rho+vol*srcU(i,j,k)
                  myRhoU=myRhoU+vol*(pipe%fs%rho*pipe%fs%U(i,j,k)+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=pipe%fs%cfg%kmin_,pipe%fs%cfg%kmax_
            do j=pipe%fs%cfg%jmin_,pipe%fs%cfg%jmax_
               do i=pipe%fs%cfg%imin_,pipe%fs%cfg%imax_
                  vol=pipe%fs%cfg%dxm(i)*pipe%fs%cfg%dy(j)*pipe%fs%cfg%dz(k)*get_VF(i,j,k,'U')
                  myUvol=myUvol+vol
                  !myRhoU=myRhoU+vol*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))*fs%rho
                  myRhoU=myRhoU+vol*pipe%fs%rho*pipe%fs%U(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,pipe%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,pipe%fs%cfg%comm,ierr); mfr=mfr/Uvol
   end function get_bodyforce_mfr

   !> Get volume fraction for direct forcing
   function get_VF(i,j,k,dir) result(VF)
      implicit none
      integer, intent(in)    :: i,j,k
      character(len=*)       :: dir
      real(WP)               :: VF
      real(WP)               :: r,eta,lam,delta,VFx,VFy,VFz
      real(WP), dimension(3) :: norm
      select case(trim(dir))
      case('U','u')
         delta=(pipe%cfg%dxm(i)*pipe%cfg%dy(j)*pipe%cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(pipe%cfg%ym(j)**2+pipe%cfg%zm(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=pipe%cfg%ym(j)/r; norm(3)=pipe%cfg%zm(k)/r
      case('V','v')
         delta=(pipe%cfg%dx(i)*pipe%cfg%dym(j)*pipe%cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(pipe%cfg%y(j)**2+pipe%cfg%zm(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=pipe%cfg%y(j)/r; norm(3)=pipe%cfg%zm(k)/r
      case('W','w')
         delta=(pipe%cfg%dx(i)*pipe%cfg%dy(j)*pipe%cfg%dzm(k))**(1.0_WP/3.0_WP)
         r=sqrt(pipe%cfg%ym(j)**2+pipe%cfg%z(k)**2)+epsilon(1.0_WP)
         norm(1)=0.0_WP; norm(2)=pipe%cfg%ym(j)/r; norm(3)=pipe%cfg%z(k)/r
      case default
         delta=(pipe%cfg%dx(i)*pipe%cfg%dy(j)*pipe%cfg%dz(k))**(1.0_WP/3.0_WP)
         r=sqrt(pipe%cfg%ym(j)**2+pipe%cfg%zm(k)**2)
         norm(1)=0.0_WP; norm(2)=pipe%cfg%ym(j)/r; norm(3)=pipe%cfg%zm(k)/r
      end select
      lam=sum(abs(norm)); eta=0.065_WP*(1.0_WP-lam**2)+0.39_WP
      VF=0.5_WP*(1.0_WP-tanh((r-0.5_WP*Djet)/(sqrt(2.0_WP)*lam*eta*delta+epsilon(1.0_WP))))
   end function get_VF
   
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
     !if (radius.le.0.5_WP*Djet.and.i.eq.pg%imin) isIn=.true.
     if (i.eq.pg%imin) isIn=.true.
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
  
   subroutine gendomain_allocate(d)
     implicit none
     class(gendomain), intent(inout) :: d
     integer :: imino_, imaxo_, jmino_, jmaxo_, kmino_, kmaxo_

     imino_ = d%cfg%imino_; imaxo_ = d%cfg%imaxo_;
     jmino_ = d%cfg%jmino_; jmaxo_ = d%cfg%jmaxo_;
     kmino_ = d%cfg%kmino_; kmaxo_ = d%cfg%kmaxo_;

     allocate(d%resU(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%resV(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%resW(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%Ui(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%Vi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%Wi(imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%SR(1:6,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%vort(1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
     allocate(d%gradu(1:3,1:3,imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_))
   end subroutine gendomain_allocate

   subroutine gendomain_deallocate(d)
      implicit none
      class(gendomain), intent(inout) :: d
      deallocate(d%resU, d%resV, d%resW)
      deallocate(d%Ui, d%Vi, d%Wi)
      deallocate(d%SR)
      deallocate(d%vort)
      deallocate(d%gradu)
    end subroutine gendomain_deallocate

     !> Monitor output shared between the two domains

   subroutine setup_gendom_monitors(d)
      implicit none
      class(gendomain), intent(inout) :: d

      ! Create simulation monitor
      d%mfile = monitor(d%fs%cfg%amRoot, trim(d%desc) // '_simulation')
      call d%mfile%add_column(time%n,'Timestep number')
      call d%mfile%add_column(time%t,'Time')
      call d%mfile%add_column(time%dt,'Timestep size')
      call d%mfile%add_column(time%cfl,'Maximum CFL')
      call d%mfile%add_column(d%fs%Umax,'Umax')
      call d%mfile%add_column(d%meanU,'Umean')
      call d%mfile%add_column(d%fs%Vmax,'Vmax')
      call d%mfile%add_column(d%meanV,'Vmean')
      call d%mfile%add_column(d%fs%Wmax,'Wmax')
      call d%mfile%add_column(d%meanW,'Wmean')
      call d%mfile%add_column(d%fs%Pmax,'Pmax')
      call d%mfile%add_column(d%fs%divmax,'Maximum divergence')
      call d%mfile%add_column(d%fs%psolv%it,'Pressure iteration')
      call d%mfile%add_column(d%fs%psolv%rerr,'Pressure error')
      call d%mfile%write()

      ! Create CFL monitor
      d%cflfile = monitor(d%fs%cfg%amRoot, trim(d%desc) // '_cfl')
      call d%cflfile%add_column(time%n,'Timestep number')
      call d%cflfile%add_column(time%t,'Time')
      call d%cflfile%add_column(d%fs%CFLc_x,'Convective xCFL')
      call d%cflfile%add_column(d%fs%CFLc_y,'Convective yCFL')
      call d%cflfile%add_column(d%fs%CFLc_z,'Convective zCFL')
      call d%cflfile%add_column(d%fs%CFLv_x,'Viscous xCFL')
      call d%cflfile%add_column(d%fs%CFLv_y,'Viscous yCFL')
      call d%cflfile%add_column(d%fs%CFLv_z,'Viscous zCFL')
      call d%cflfile%write()
   end subroutine setup_gendom_monitors

   subroutine write_gendom_monitors(d)
      implicit none
      class(gendomain), intent(inout) :: d
  
      call d%fs%get_max()
      call d%mfile%write()
      call d%cflfile%write()
  
    end subroutine write_gendom_monitors

   !> Initializes pjet domain stretched grid
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
         real(WP) :: xtilde,ytilde,ztilde
         real(WP) :: alpha,minAlpha,maxAlpha,r
         real(WP), dimension(:), allocatable :: x
         real(WP), dimension(:), allocatable :: y
         real(WP), dimension(:), allocatable :: z
         integer, dimension(3) :: partition
         ! Read in grid definition
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
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='jet')
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         pjet%cfg=config(grp=group,decomp=partition,grid=grid)
         pjet%cfg%VF=1.0_WP
      end block create_grid
   end subroutine pjet_geometry_init
   
   !> Initialization of pjet solver
   subroutine pjet_init
      use param, only: param_read
      use messager, only: log
      implicit none

      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: smg,pcg_pfmg,pfmg,gmres_smg
         use incomp_class,    only: clipped_neumann, dirichlet, slip
         ! Create flow solver
         pjet%fs=incomp(cfg=pjet%cfg,name='NS solver')
         ! Get jet diameter
         call param_read('Jet diameter',Djet,default=0.1_WP)
         theta = Djet / 20.0_WP
         call param_read('Flow rate',Q0,default=0.001_WP)
         call param_read('Pulse frequency',omegaj,default=1.0_WP)
         call param_read('Pulse decay time',tauj,default=0.2_WP)
         call param_read('Coflow velocity',Uco,default=0.46_WP)
         ! Add BCs
         call pjet%fs%add_bcond(name='jet',   type=dirichlet,locator=jet_loc,face='x',dir=-1,canCorrect=.false. )
         call pjet%fs%add_bcond(name='right', type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true.)
         call pjet%fs%add_bcond(name='bottom',type=slip, locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call pjet%fs%add_bcond(name='top',   type=slip, locator=top_of_domain,face='y',dir=+1,canCorrect=.false.)
         call pjet%fs%add_bcond(name='front', type=slip, locator=front_of_domain,face='z',dir=+1,canCorrect=.false.)
         call pjet%fs%add_bcond(name='back',  type=slip, locator=back_of_domain,face='z',dir=-1,canCorrect=.false.)
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); pjet%fs%visc=visc
         ! Assign constant density
         call param_read('Density',pjet%fs%rho)
         ! Prepare and configure pressure solver
         pjet%vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         pjet%ps=hypre_str(cfg=cfg,name='Pressure',method=smg,nst=7)
         call param_read('Pressure iteration',pjet%ps%maxit)
         call param_read('Pressure tolerance',pjet%ps%rcvg)
         ! Setup the solver
         call pjet%fs%setup(pressure_solver=pjet%ps, implicit_solver=pjet%vs)
      end block create_and_initialize_flow_solver
      
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
         
         Qt = qdat(1)
         U0 = Qt / (Pi * Djet**2 / 4.0_WP) * 1.667E-05_WP 
         Uco = 0.05_WP * U0
         
         pjet%fs%U = 0.0_WP
         pjet%fs%V = 0.0_WP
         pjet%fs%W = 0.0_WP

         ! pjet%fs%U=0.0_WP; pjet%fs%V=0.0_WP; pjet%fs%W=0.0_WP
         ! call pjet%fs%get_bcond('jet', mybc)
         ! do n=1,mybc%itr%no_
         !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
         !    r = norm2([pjet%fs%cfg%ym(j), pjet%fs%cfg%zm(k)])
         !    pjet%fs%U(i,j,k) = 0.5_WP * (U0 + Uco) + 0.5_WP * (U0 - Uco) * TANH(0.5_WP * (0.5_WP * Djet - r) / theta)
         ! end do

         call pjet%fs%get_mfr()
         call pjet%fs%correct_mfr()
         call pjet%fs%interp_vel(pjet%Ui,pjet%Vi,pjet%Wi)
         call pjet%fs%get_vorticity(pjet%vort)
      end block initialize_velocity

   end subroutine pjet_init

   !> Initialization of problem solver
   subroutine pipe_init
      use param, only: param_read
      implicit none

      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         pipe%fs=incomp(cfg=pipe%cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',pipe%fs%rho)
         call param_read('Dynamic viscosity',visc); pipe%fs%visc=visc
         ! Configure pressure solver
         pipe%ps=fft3d(cfg=pipe%cfg,name='Pressure',nst=7)
         !ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         !ps%maxlevel=14
         !call param_read('Pressure iteration',ps%maxit)
         !call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         pipe%vs=hypre_str(cfg=pipe%cfg,name='Velocity',method=pcg_pfmg,nst=7)
         call param_read('Implicit iteration',pipe%vs%maxit)
         call param_read('Implicit tolerance',pipe%vs%rcvg)
         ! Setup the solver
         call pipe%fs%setup(pressure_solver=pipe%ps,implicit_solver=pipe%vs)
      end block create_flow_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(pipe%G   (pipe%cfg%imino_:pipe%cfg%imaxo_,pipe%cfg%jmino_:pipe%cfg%jmaxo_,pipe%cfg%kmino_:pipe%cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k

         ! Handle restart/saves here
         pipeic: block
            character(len=str_medium) :: timestamp
            ! Create event for saving restart files
            ! Check if we are restarting
            call param_read(tag='Pipe IC File',val=timestamp,short='r',default='')
            ! Read the datafile
            df=datafile(pg=pipe%cfg,fdata='pipe_ic/data_'//trim(adjustl(timestamp)))
         end block pipeic

         call df%pullval(name='dt',val=time%dt)
         call df%pullvar(name='U',var=pipe%fs%U)
         call df%pullvar(name='V',var=pipe%fs%V)
         call df%pullvar(name='W',var=pipe%fs%W)
         call df%pullvar(name='P',var=pipe%fs%P)
         
         call pipe%fs%cfg%sync(pipe%fs%U)
         call pipe%fs%cfg%sync(pipe%fs%V)
         call pipe%fs%cfg%sync(pipe%fs%W)
         ! Compute cell-centered velocity
         call pipe%fs%interp_vel(pipe%Ui,pipe%Vi,pipe%Wi)
         ! Compute divergence
         call pipe%fs%get_div()
         call pipe%fs%get_vorticity(pipe%vort)
         ! Get target MFR and zero bodyforce
         pipe%mfr=get_bodyforce_mfr()
         pipe%mfr_target=pipe%mfr
         pipe%bforce=0.0_WP
      end block initialize_velocity
      
      ! Initialize IBM fields
      initialize_ibm: block
         integer :: i,j,k
         do k=pipe%fs%cfg%kmino_,pipe%fs%cfg%kmaxo_
            do j=pipe%fs%cfg%jmino_,pipe%fs%cfg%jmaxo_
               do i=pipe%fs%cfg%imino_,pipe%fs%cfg%imaxo_
                  pipe%G(i,j,k)=0.5_WP*Djet-sqrt(pipe%fs%cfg%ym(j)**2+pipe%fs%cfg%zm(k)**2)
               end do
            end do
         end do
      end block initialize_ibm
      
   end subroutine pipe_init
   
   !> Initialization of problem geometry
   subroutine pipe_geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx,overlap
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Pipe length',Lx)
         call param_read('Overlap', overlap)
         call param_read('Pipe ny',ny); allocate(y(ny+1))
         call param_read('Pipe nx',nx); allocate(x(nx+1))
         call param_read('Pipe nz',nz); allocate(z(nz+1))
         
         dx=Lx/real(nx,WP)
         no=6
         if (ny.gt.1) then
            Ly=Djet+real(2*no,WP)*Djet/real(ny-2*no,WP)
         else
            Ly=dx
         end if
         if (nz.gt.1) then
            Lz=Djet+real(2*no,WP)*Djet/real(ny-2*no,WP)
         else
            Lz=dx
         end if
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx + overlap
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='pipe')
      end block create_grid
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         pipe%cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      ! Create masks for this config
      create_walls: block
         integer :: i,j,k
         do k=pipe%cfg%kmin_,pipe%cfg%kmax_
            do j=pipe%cfg%jmin_,pipe%cfg%jmax_
               do i=pipe%cfg%imin_,pipe%cfg%imax_
                  pipe%cfg%VF(i,j,k)=max(get_VF(i,j,k,'SC'),epsilon(1.0_WP))
               end do
            end do
         end do
         call pipe%cfg%sync(pipe%cfg%VF)
         call pipe%cfg%calc_fluid_vol()
      end block create_walls
   end subroutine pipe_geometry_init 

   !> ENSIGHT ROUTINES 
   subroutine setup_ens(d)
      implicit none
      class(gendomain), intent(inout) :: d
      d%ens_out = ensight(cfg=d%cfg, name=trim(d%desc))
      call d%ens_out%add_vector('velocity', d%Ui, d%Vi, d%Wi)
      call d%ens_out%add_scalar('pressure', d%fs%P)
      call d%ens_out%add_scalar('vorticity', d%vort)
   end subroutine setup_ens

   subroutine write_ens(d)
      implicit none
      class(gendomain), intent(inout) :: d
      integer :: i
      call d%fs%interp_vel(d%Ui, d%Vi, d%Wi)
      call d%fs%get_vorticity(d%vort)
      call d%ens_out%write_data(time%t)
   end subroutine write_ens

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize domain geometries
      call pipe_geometry_init()
      call pjet_geometry_init()

      ! Allocate work arrays
      call gendomain_allocate(pipe)
      call gendomain_allocate(pjet)

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max iter',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Initialize solvers
      call pjet_init()
      call pipe_init()

      ! Set up domain couplers
      dom_cpl = coupler(src_grp=pipe%cfg%group, dst_grp=pjet%cfg%group,name='Domain Coupler')
      call dom_cpl%set_src(pipe%cfg)
      if (.true.) then    ! all tasks are in the destination group
         call dom_cpl%set_dst(pjet%cfg)
      end if
      call dom_cpl%initialize()

   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      use parallel,       only: parallel_time
      implicit none
      integer :: ii, it
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call pipe%fs%get_cfl(time%dt,pipe%cfl)
         call pjet%fs%get_cfl(time%dt,pjet%cfl)
         time%cfl = max(pipe%cfl, pjet%cfl)
         call time%adjust_dt()
         call time%increment()

         !! PIPE STEP
         ! Remember old velocity
         pipe%fs%Uold=pipe%fs%U
         pipe%fs%Vold=pipe%fs%V
         pipe%fs%Wold=pipe%fs%W
         
         ! Perform sub-iterations
         do it = 1, time%itmax
            
            ! Build mid-time velocity
            pipe%fs%U=0.5_WP*(pipe%fs%U+pipe%fs%Uold)
            pipe%fs%V=0.5_WP*(pipe%fs%V+pipe%fs%Vold)
            pipe%fs%W=0.5_WP*(pipe%fs%W+pipe%fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call pipe%fs%get_dmomdt(pipe%resU,pipe%resV,pipe%resW)
            
            ! Assemble explicit residual
            pipe%resU=-2.0_WP*(pipe%fs%rho*pipe%fs%U-pipe%fs%rho*pipe%fs%Uold)+time%dt*pipe%resU
            pipe%resV=-2.0_WP*(pipe%fs%rho*pipe%fs%V-pipe%fs%rho*pipe%fs%Vold)+time%dt*pipe%resV
            pipe%resW=-2.0_WP*(pipe%fs%rho*pipe%fs%W-pipe%fs%rho*pipe%fs%Wold)+time%dt*pipe%resW
            
            ! Add body forcing
            bodyforcing: block
               real(WP) :: mfr
               pipe%mfr=get_bodyforce_mfr(pipe%resU)
               pipe%bforce=(pipe%mfr_target-pipe%mfr)/time%dtmid
               pipe%resU=pipe%resU+time%dt*pipe%bforce
            end block bodyforcing

            ! Form implicit residuals
            call pipe%fs%solve_implicit(time%dt,pipe%resU,pipe%resV,pipe%resW)
            
            ! Apply these residuals
            pipe%fs%U=2.0_WP*pipe%fs%U-pipe%fs%Uold+pipe%resU
            pipe%fs%V=2.0_WP*pipe%fs%V-pipe%fs%Vold+pipe%resV
            pipe%fs%W=2.0_WP*pipe%fs%W-pipe%fs%Wold+pipe%resW
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
               integer :: i,j,k
               do k=pipe%fs%cfg%kmin_,pipe%fs%cfg%kmax_
                  do j=pipe%fs%cfg%jmin_,pipe%fs%cfg%jmax_
                     do i=pipe%fs%cfg%imin_,pipe%fs%cfg%imax_
                        pipe%fs%U(i,j,k)=get_VF(i,j,k,'U')*pipe%fs%U(i,j,k)
                        pipe%fs%V(i,j,k)=get_VF(i,j,k,'V')*pipe%fs%V(i,j,k)
                        pipe%fs%W(i,j,k)=get_VF(i,j,k,'W')*pipe%fs%W(i,j,k)
                     end do
                  end do
               end do
               call pipe%fs%cfg%sync(pipe%fs%U)
               call pipe%fs%cfg%sync(pipe%fs%V)
               call pipe%fs%cfg%sync(pipe%fs%W)
            end block ibforcing
           
            ! Apply other boundary conditions on the resulting fields
            call pipe%fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call pipe%fs%correct_mfr()
            call pipe%fs%get_div()
            pipe%fs%psolv%rhs=-pipe%fs%cfg%vol*pipe%fs%div*pipe%fs%rho/time%dt
            pipe%fs%psolv%sol=0.0_WP
            call pipe%fs%psolv%solve()
            call pipe%fs%shift_p(pipe%fs%psolv%sol)
            
            ! Correct velocity
            call pipe%fs%get_pgrad(pipe%fs%psolv%sol,pipe%resU,pipe%resV,pipe%resW)
            pipe%fs%P=pipe%fs%P+pipe%fs%psolv%sol
            pipe%fs%U=pipe%fs%U-time%dt*pipe%resU/pipe%fs%rho
            pipe%fs%V=pipe%fs%V-time%dt*pipe%resV/pipe%fs%rho
            pipe%fs%W=pipe%fs%W-time%dt*pipe%resW/pipe%fs%rho
        
         end do

         ! Inflow coupling
         call dom_cpl%push(pipe%fs%U)
         call dom_cpl%transfer()
         call dom_cpl%pull(pjet%resU)
         call dom_cpl%push(pipe%fs%V)
         call dom_cpl%transfer()
         call dom_cpl%pull(pjet%resV)
         call dom_cpl%push(pipe%fs%W)
         call dom_cpl%transfer()
         call dom_cpl%pull(pjet%resW)
         apply_inflow_vals: block
           if (pjet%cfg%imin .eq. pjet%cfg%imin_) then
             pjet%fs%U(pjet%cfg%imino_:pjet%cfg%imin_,:,:) = pjet%resU(pjet%cfg%imino_:pjet%cfg%imin_,:,:)
             pjet%fs%V(pjet%cfg%imino_:pjet%cfg%imin_-1,:,:) = pjet%resV(pjet%cfg%imino_:pjet%cfg%imin_-1,:,:)
             pjet%fs%W(pjet%cfg%imino_:pjet%cfg%imin_-1,:,:) = pjet%resW(pjet%cfg%imino_:pjet%cfg%imin_-1,:,:)
           end if
         end block apply_inflow_vals
         call pjet%fs%apply_bcond(time%t, time%dt)
         
         !! PJET STEP
         ! Remember old velocity
         pjet%fs%Uold=pjet%fs%U
         pjet%fs%Vold=pjet%fs%V
         pjet%fs%Wold=pjet%fs%W
         
         ! Perform sub-iterations
         do it = 1, time%itmax
            
            ! Build mid-time velocity
            pjet%fs%U=0.5_WP*(pjet%fs%U+pjet%fs%Uold)
            pjet%fs%V=0.5_WP*(pjet%fs%V+pjet%fs%Vold)
            pjet%fs%W=0.5_WP*(pjet%fs%W+pjet%fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call pjet%fs%get_dmomdt(pjet%resU,pjet%resV,pjet%resW)
            ! Assemble explicit residual
            pjet%resU=-2.0_WP*(pjet%fs%rho*pjet%fs%U-pjet%fs%rho*pjet%fs%Uold)+time%dt*pjet%resU
            pjet%resV=-2.0_WP*(pjet%fs%rho*pjet%fs%V-pjet%fs%rho*pjet%fs%Vold)+time%dt*pjet%resV
            pjet%resW=-2.0_WP*(pjet%fs%rho*pjet%fs%W-pjet%fs%rho*pjet%fs%Wold)+time%dt*pjet%resW
            ! Implicit velocity solve
            call pjet%fs%solve_implicit(time%dt,pjet%resU,pjet%resV,pjet%resW)

            ! Apply these residuals
            pjet%fs%U=2.0_WP*pjet%fs%U-pjet%fs%Uold+pjet%resU
            pjet%fs%V=2.0_WP*pjet%fs%V-pjet%fs%Vold+pjet%resV
            pjet%fs%W=2.0_WP*pjet%fs%W-pjet%fs%Wold+pjet%resW

            !> TODO: Time-varying boundary condition
            dirichlet_velocity: block
               use incomp_class, only: bcond
               use mathtools,    only: Pi
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,jn,jp
               real(WP) :: r

               jn = FLOOR(time%t/dt_dat) + 1
               jp = jn + 1;
               
               ! Linear interpolation of provided flow rate data
               if (jp.lt.1500) then
                  Qt = (qdat(jp) - qdat(jn)) / (tdat(jp) - tdat(jn)) * (time%t - tdat(jn)) + qdat(jn) 
                  U0 = Qt / (Pi * Djet**2 / 4.0_WP) * 1.667E-05_WP 
                  Uco = 0.05_WP * U0
               else
                  U0 = qdat(1500) / (Pi * Djet**2 / 4.0_WP) * 1.667E-05_WP 
                  Uco = 0.05_WP * U0
               end if
            end block dirichlet_velocity
             
            ! Apply other boundary conditions on the resulting fields
            call pjet%fs%apply_bcond(time%t,time%dt)

            ! Solve Poisson equation
            call pjet%fs%get_mfr()
            call pjet%fs%correct_mfr()
            call pjet%fs%get_div()
            pjet%fs%psolv%rhs=-pjet%fs%cfg%vol*pjet%fs%div*pjet%fs%rho/time%dt
            pjet%fs%psolv%sol=0.0_WP
            call pjet%fs%psolv%solve()
            call pjet%fs%shift_p(pjet%fs%psolv%sol)
            
            ! Correct velocity
            call pjet%fs%get_pgrad(pjet%fs%psolv%sol,pjet%resU,pjet%resV,pjet%resW)
            pjet%fs%P=pjet%fs%P+pjet%fs%psolv%sol
            pjet%fs%U=pjet%fs%U-time%dt*pjet%resU/pjet%fs%rho
            pjet%fs%V=pjet%fs%V-time%dt*pjet%resV/pjet%fs%rho
            pjet%fs%W=pjet%fs%W-time%dt*pjet%resW/pjet%fs%rho

         end do
         
         ! Recompute interpolated velocity and divergence
         call pjet%fs%interp_vel(pjet%Ui,pjet%Vi,pjet%Wi)
         call pjet%fs%get_div()
         
         !TODO: ENSIGHT AND MONITOR OUTPUTS

         call write_gendom_monitors(pipe)
         call write_gendom_monitors(pjet)

         if (ens_evt%occurs()) then
            call write_ens(pjet)
            call write_ens(pipe)
         end if
         
      end do
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      call gendomain_deallocate(pjet)
      call gendomain_deallocate(pipe)
   end subroutine simulation_final
   
end module simulation
