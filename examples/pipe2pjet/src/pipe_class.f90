!> Various definitions and tools for running an NGA2 simulation
module pipe_class
   use precision,         only: WP
!   use geometry,          only: cfg,D
   use ibconfig_class,    only: ibconfig
   use fft3d_class,       only: fft3d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use datafile_class,    only: datafile
   use string,            only: str_medium
   implicit none
   private

   public :: pipedomain
   
   type :: pipedomain
      character(len=str_medium) :: desc
      !> Get an an incompressible solver, pressure solver, and corresponding time tracker
      type(incomp),      public :: fs
      type(fft3d),       public :: ps
      type(ddadi),       public :: vs
      type(timetracker), public :: time
      type(ibconfig),    public :: cfg
      
      !> Ensight postprocessing
      type(ensight)  :: ens_out
      type(event)    :: ens_evt

      !> Provide a datafile and an event tracker for saving restarts
      type(event)    :: save_evt
      type(datafile) :: df
      ! type(partmesh) :: pmesh
      
      !> Simulation monitor file
      type(monitor) :: mfile,cflfile,simfile
      
      !> Work arrays
      real(WP), dimension(:,:,:,:,:), allocatable :: gradU
      real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
      real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
      real(WP) :: visc,mfr_target,mfr,bforce,cfl

      ! Stats
      real(WP), dimension (:), allocatable :: Uxr, Urr, Utr, Ux2, Ur2, Ut2,vol
      real(WP), dimension (:), allocatable :: Uxr_, Urr_, Utr_, Ux2_, Ur2_, Ut2_,vol_

      !> Event for calling post-processing script
      type(event) :: ppevt
   contains
      procedure :: get_bodyforce_mfr
      procedure :: setup_monitors
      procedure :: write_monitors
      procedure :: init
      procedure :: geometry_init
      procedure :: step
      procedure :: final
   end type pipedomain
    
   integer :: nr
   
   logical :: restarted
   
contains
   
   !> Compute massflow rate
   subroutine get_bodyforce_mfr(this, srcU, mfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      class(pipedomain), intent(inout) :: this
      real(WP), dimension(this%fs%cfg%imino_:,this%fs%cfg%jmino_:,this%fs%cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(inout) :: mfr
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  vol=this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*(this%fs%U(i,j,k)*this%fs%rho+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
            do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
               do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                  vol=this%fs%cfg%dxm(i)*this%fs%cfg%dy(j)*this%fs%cfg%dz(k)*sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*this%fs%U(i,j,k)*this%fs%rho
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,this%fs%cfg%comm,ierr); mfr=mfr/Uvol
   end subroutine get_bodyforce_mfr

   ! !> Subroutine for getting pipe velocity stats
   ! subroutine postproc_vel(this)
   !    use string,    only: str_medium
   !    use mpi_f08,   only: MPI_ALLREDUCE, MPI_SUM
   !    use parallel,  only: MPI_REAL_WP
   !    implicit none
   !    class(pipedomain), intent(inout) :: this
   !    integer :: iunit, ierr, i, j, k, rind
   !    real(WP) :: rnorm, rmax, dr
   !    real(WP), dimension(3) :: r, vel, rhat, rtheta
   !    character(len=str_medium) :: filename, timestamp

   !    rmax = 0.5_WP * Djet
   !    dr = rmax / real(nr,WP)

   !    ! GET MEAN PROFILES
   !    do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
   !       do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
   !          do i=this%fs%cfg%imin_,this%fs%cfg%imax_
   !             r = [0.0_WP, this%fs%cfg%ym(j), this%fs%cfg%zm(k)]    ! get radial position on axial slice
   !             rnorm = norm2(r)                          
   !             if (rnorm.gt.rmax) cycle    ! don't look at anything r > d/2
   !             rind = floor(rnorm/dr) + 1  
   !             rhat = r / rnorm

   !             ! assemble theta-wise unit vector
   !             rtheta = cross(rhat, [0.0_WP, 1.0_WP, 0.0_WP])  
   !             rtheta = cross(rtheta, rhat) / norm2(cross(rtheta, rhat))

   !             ! compute it all
   !             vel = [Ui(i,j,k), Vi(i,j,k), Wi(i,j,k)]
   !             vol_(rind)=vol_(rind)+this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Uxr_(rind) = Uxr_(rind) + Ui(i,j,k)*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Urr_(rind) = Urr_(rind) + dot_product(vel, rhat)*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Utr_(rind) = Utr_(rind) + dot_product(vel, rtheta)*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Ux2_(rind) = Ux2_(rind) + Ui(i,j,k)**2*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Ur2_(rind) = Ur2_(rind) + dot_product(vel, rhat)**2*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !             Ut2_(rind) = Ut2_(rind) + dot_product(vel, rtheta)**2*this%fs%cfg%vol(i,j,k)*this%fs%cfg%VF(i,j,k)*time%dt
   !          end do
   !       end do
   !    end do

   !    call MPI_ALLREDUCE(vol_, vol, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Uxr_, Uxr, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Urr_, Urr, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Utr_, Utr, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Ux2_, Ux2, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Ur2_, Ur2, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)
   !    call MPI_ALLREDUCE(Ut2_, Ut2, nr, MPI_REAL_WP, MPI_SUM, this%fs%cfg%comm, ierr)

   !    do i=1,nr
   !       if (vol(i).gt.0.0_WP) then
   !          Uxr(i)=Uxr(i)/vol(i)
   !          Urr(i)=Urr(i)/vol(i)
   !          Utr(i)=Utr(i)/vol(i)
   !          Ux2(i)=Ux2(i)/vol(i)
   !          Ur2(i)=Ur2(i)/vol(i)
   !          Ut2(i)=Ut2(i)/vol(i)
   !       end if
   !    end do

   !     ! If root, print it out
   !    if (this%fs%cfg%amRoot.and.ppevt%occurs()) then
   !       filename='rstat_'
   !       write(timestamp,'(es12.5)') time%t
   !       open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
   !       write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'r','Ux','Ux_rms','Ut_rms','Ur_rms'
   !       do i=1,nr
   !          write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') dr*real(i-1,WP),Uxr(i),Uxr(i)**2-Ux2(i),Utr(i)**2-Ut2(i),Urr(i)**2-Ur2(i)
   !       end do
   !       close(iunit)
   !    end if
      
   !    contains

   !       function cross(a, b) result(c)
   !          real(WP), dimension(3) :: c
   !          real(WP), dimension(3), intent(IN) :: a, b
   !          c(1) = a(2) * b(3) - a(3) * b(2)
   !          c(2) = a(3) * b(1) - a(1) * b(3)
   !          c(3) = a(1) * b(2) - a(2) * b(1)
   !       end function cross

   ! end subroutine postproc_vel   
   
   !> Monitor init
   subroutine  setup_monitors(this)
      implicit none
      class(pipedomain), intent(inout) :: this

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
      class(pipedomain), intent(inout) :: this
  
      call this%fs%get_max()
      call this%mfile%write()
      call this%cflfile%write()
   end subroutine write_monitors

   !> Initialization of problem geometry
   subroutine geometry_init(this, group, pipepart)
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use mpi_f08,     only: MPI_GROUP
      implicit none
      type(sgrid) :: grid
      type(MPI_GROUP), intent(in) :: group
      class(pipedomain), intent(inout) :: this
      integer, dimension(3), intent(in) :: pipepart
      real(WP) :: Djet
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,no
         real(WP) :: Lx,Ly,Lz,dx,overlap
         real(WP), dimension(:), allocatable :: x,y,z
       
         ! Read in grid definition
         call param_read('Jet diameter', Djet)
         call param_read('Pipe length',Lx)
         call param_read('Overlap', overlap)
         call param_read('Pipe nx',nx); allocate(x(nx+1))
         call param_read('Pipe ny',ny); allocate(y(ny+1))
         call param_read('Pipe nz',nz); allocate(z(nz+1))
          
         dx=Lx/real(nx,WP)
         no=4
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
            x(i)=real(i-1,WP)/real(nx,WP)*Lx - Lx + overlap
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.true.,zper=.true.,name='pipe')
         deallocate(x, y, z)
      end block create_grid
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         ! Create partitioned grid
         this%cfg=ibconfig(grp=group,decomp=pipepart,grid=grid)
      end block create_cfg
      
      ! Create masks for this config
      create_walls: block
         use ibconfig_class,  only: bigot, sharp
         integer :: i,j,k
         do k=this%cfg%kmino_,this%cfg%kmaxo_
            do j=this%cfg%jmino_,this%cfg%jmaxo_
               do i=this%cfg%imino_,this%cfg%imaxo_
                  this%cfg%Gib(i,j,k)=0.5_WP*Djet - sqrt(this%cfg%ym(j)**2 + this%cfg%zm(k)**2)
               end do
            end do
         end do
         call this%cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
         call this%cfg%sync(this%cfg%VF)
      end block create_walls
   end subroutine geometry_init 

   !> Initialization of problem solver
   subroutine init(this)
      use param, only: param_read
      implicit none
      class(pipedomain), intent(inout) :: this

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         use param, only: param_read
         this%time=timetracker(amRoot=this%cfg%amRoot)
         call param_read('Max timestep size',this%time%dtmax)
         call param_read('Max cfl number',this%time%cflmax)
         this%time%dt=this%time%dtmax
         this%time%itmax=2
      end block initialize_timetracker

      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
         ! Create flow solver
         this%fs=incomp(cfg=this%cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',this%fs%rho)
         call param_read('Dynamic viscosity',this%visc); this%fs%visc=this%visc
         ! Configure pressure solver
         this%ps=fft3d(cfg=this%cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         this%vs=ddadi(cfg=this%cfg,name='Velocity',nst=7)
         call param_read('Implicit iteration',this%vs%maxit)
         call param_read('Implicit tolerance',this%vs%rcvg)
         ! Setup the solver
         call this%fs%setup(pressure_solver=this%ps,implicit_solver=this%vs)
      end block create_flow_solver
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(this%gradU(1:3,1:3,this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))   
         allocate(this%resU(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resV(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%resW(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Ui  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Vi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
         allocate(this%Wi  (this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      end block allocate_work_arrays

      ! Handle restart/saves here
      restart_and_save: block
        character(len=str_medium) :: timestamp
        ! Check if we are restarting
        call param_read(tag='Pipe IC File',val=timestamp,short='r',default='')
        restarted=.false.; if (len_trim(timestamp).gt.0) restarted=.true.
        if (restarted) then
           ! If we are, read the name of the directory
           call param_read('Pipe IC File',timestamp,'r')
           ! Read the datafile
           this%df=datafile(pg=this%cfg,fdata='pipe_ic/data_'//trim(adjustl(timestamp)))
        else
           ! If we are not restarting, we will still need a datafile for saving restart files
           this%df=datafile(pg=this%cfg,filename=trim(this%cfg%name),nval=2,nvar=4)
           this%df%valname(1)='t'
           this%df%valname(2)='dt'
           this%df%varname(1)='U'
           this%df%varname(2)='V'
           this%df%varname(3)='W'
           this%df%varname(4)='P'
        end if
      end block restart_and_save
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Ubulk,amp
         ! Initial fields
         call param_read('Bulk velocity',Ubulk)
         this%fs%U=Ubulk; this%fs%V=0.0_WP; this%fs%W=0.0_WP; this%fs%P=0.0_WP
         ! For faster transition
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         if (restarted) then
            call this%df%pullvar(name='U',var=this%fs%U)
            call this%df%pullvar(name='V',var=this%fs%V)
            call this%df%pullvar(name='W',var=this%fs%W)
            call this%df%pullvar(name='P',var=this%fs%P)
         else
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%fs%U(i,j,k)=this%fs%U(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%zm(k)/this%fs%cfg%zL)*cos(8.0_WP*twoPi*this%fs%cfg%ym(j)/this%fs%cfg%yL)
                     this%fs%V(i,j,k)=this%fs%V(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
                     this%fs%W(i,j,k)=this%fs%W(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*this%fs%cfg%xm(i)/this%fs%cfg%xL)
                  end do
               end do
            end do
         end if
         call this%fs%cfg%sync(this%fs%U)
         call this%fs%cfg%sync(this%fs%V)
         call this%fs%cfg%sync(this%fs%W)
         ! Compute cell-centered velocity
         call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
         ! Compute divergence
         call this%fs%get_div()
         ! Get target MFR and zero bodyforce
         call get_bodyforce_mfr(this=this, mfr=this%mfr)
         this%mfr_target=this%mfr
         this%bforce=0.0_WP
         this%desc = 'pipe'
      end block initialize_velocity

      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         nr=this%fs%cfg%ny/2
         this%ppevt=event(time=this%time,name='Postproc output')
         call param_read('Postproc output period',this%ppevt%tper)
         ! Allocate
         allocate(this%vol(nr));  this%vol=0.0_WP
         allocate(this%Uxr(nr));  this%Uxr=0.0_WP
         allocate(this%Urr(nr));  this%Urr=0.0_WP
         allocate(this%Utr(nr));  this%Utr=0.0_WP
         allocate(this%Ux2(nr));  this%Ux2=0.0_WP
         allocate(this%Ur2(nr));  this%Ur2=0.0_WP
         allocate(this%Ut2(nr));  this%Ut2=0.0_WP
         allocate(this%vol_(nr)); this%vol_=0.0_WP
         allocate(this%Uxr_(nr)); this%Uxr_=0.0_WP
         allocate(this%Urr_(nr)); this%Urr_=0.0_WP
         allocate(this%Utr_(nr)); this%Utr_=0.0_WP
         allocate(this%Ux2_(nr)); this%Ux2_=0.0_WP
         allocate(this%Ur2_(nr)); this%Ur2_=0.0_WP
         allocate(this%Ut2_(nr)); this%Ut2_=0.0_WP
      end block create_postproc

      create_ensight: block
         this%ens_out=ensight(cfg=this%cfg,name='pipe')
         this%ens_evt=event(time=this%time,name='Ensight output')
         call param_read('Ensight output period', this%ens_evt%tper)
         call this%ens_out%add_vector('velocity',this%Ui,this%Vi,this%Wi)
         call this%ens_out%add_scalar('divergence',this%fs%div)
         call this%ens_out%add_scalar('pressure',this%fs%P)
         if (this%ens_evt%occurs()) call this%ens_out%write_data(this%time%t)
      end block create_ensight     
 
   end subroutine init
   
   
   !> Perform an NGA2 simulation
   subroutine step(this, dt)
      implicit none
      class(pipedomain), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: it
      this%time%dt=dt; call this%time%increment()
      
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
         
         ! Add body forcing
         bodyforcing: block
            call get_bodyforce_mfr(this=this,srcU=this%resU,mfr=this%mfr)
            this%bforce=(this%mfr_target-this%mfr)/this%time%dtmid
            this%resU=this%resU+this%time%dt*this%bforce
         end block bodyforcing

         ! Form implicit residuals
         call this%fs%solve_implicit(this%time%dt,this%resU,this%resV,this%resW)
         
         ! Apply these residuals
         this%fs%U=2.0_WP*this%fs%U-this%fs%Uold+this%resU
         this%fs%V=2.0_WP*this%fs%V-this%fs%Vold+this%resV
         this%fs%W=2.0_WP*this%fs%W-this%fs%Wold+this%resW
         
         ! Apply IB forcing to enforce BC at the pipe walls
         ibforcing: block
            integer :: i,j,k
            do k=this%fs%cfg%kmin_,this%fs%cfg%kmax_
               do j=this%fs%cfg%jmin_,this%fs%cfg%jmax_
                  do i=this%fs%cfg%imin_,this%fs%cfg%imax_
                     this%fs%U(i,j,k)=sum(this%fs%itpr_x(:,i,j,k)*this%cfg%VF(i-1:i,j,k))*this%fs%U(i,j,k)
                     this%fs%V(i,j,k)=sum(this%fs%itpr_y(:,i,j,k)*this%cfg%VF(i,j-1:j,k))*this%fs%V(i,j,k)
                     this%fs%W(i,j,k)=sum(this%fs%itpr_z(:,i,j,k)*this%cfg%VF(i,j,k-1:k))*this%fs%W(i,j,k)
                  end do
               end do
            end do
            call this%fs%cfg%sync(this%fs%U)
            call this%fs%cfg%sync(this%fs%V)
            call this%fs%cfg%sync(this%fs%W)
         end block ibforcing

         ! Apply other boundary conditions on the resulting fields
         call this%fs%apply_bcond(this%time%t,this%time%dt)
         
         ! Solve Poisson equation
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
      
      ! Recompute interpolated velocity and divergence
      call this%fs%interp_vel(this%Ui,this%Vi,this%Wi)
      call this%fs%get_div()

      call this%write_monitors()
      call this%ens_out%write_data(this%time%t)
   end subroutine step
   
   
   !> Finalize the NGA2 simulation
   subroutine final(this)
      implicit none
      class(pipedomain), intent(inout) :: this
      
      ! Deallocate work arrays
      ! deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      ! deallocate(Uxr, Urr, Utr, Ux2, Ur2, Ut2,vol,Uxr_, Urr_, Utr_, Ux2_, Ur2_, Ut2_,vol_)
   end subroutine final
   
end module pipe_class
