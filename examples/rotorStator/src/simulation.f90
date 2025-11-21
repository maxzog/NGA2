!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use timetracker_class, only: timetracker
   use sgsmodel_class,    only: sgsmodel
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use lpt_class,         only: lpt
   use partmesh_class,    only: partmesh
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(incomp),      public :: fs
   type(sgsmodel),    public :: sgs
   type(lpt),         public :: lp
   type(timetracker), public :: time
   type(partmesh),    public :: pmesh
   
   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,lpfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib,srcM
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:,:), allocatable :: SR

   !> Global
   real(WP) :: visc

   !> IB motion
   real(WP) :: rotor_velocity,y_shift
   real(WP), dimension(:,:,:), allocatable :: Gib0
   
   
contains
   
   
   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain
   
   
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
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Uib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcM(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Gib0(cfg%imin_:cfg%imax_,cfg%jmin:cfg%jmax,cfg%kmin_:cfg%kmax_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Create a flow solver with inflow-outflow
      create_flow_solver: block
         use incomp_class, only: dirichlet,clipped_neumann
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Evaluate IB velocity and mass source
      calc_ib_velocity: block
         use parallel, only: MPI_REAL_WP
         use mpi_f08,  only: MPI_IN_PLACE, MPI_ALLREDUCE, MPI_SUM
         integer :: i,j,k,ierr
         ! Get rotor velocity
         call param_read('Rotor velocity',rotor_velocity)
         Gib0=0.0_WP
         y_shift=0.0_WP
         ! Set IB velocities and store Gib0
         do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            do j=fs%cfg%jmino_,fs%cfg%jmaxo_
               do i=fs%cfg%imino_,fs%cfg%imaxo_
                  ! Velocity
                  Uib(i,j,k)=0.0_WP
                  if (cfg%xm(i).lt.0.0_WP) then
                     Vib(i,j,k)=rotor_velocity
                  else
                     Vib(i,j,k)=0.0_WP
                  end if
                  Wib(i,j,k)=0.0_WP
               end do
            end do
         end do
         ! Compute IB mass source
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  srcM(i,j,k)=fs%rho*(1.0_WP-cfg%VF(i,j,k))*(sum(fs%divp_x(:,i,j,k)*Uib(i:i+1,j,k))+&
                  &                                          sum(fs%divp_y(:,i,j,k)*Vib(i,j:j+1,k))+&
                  &                                          sum(fs%divp_z(:,i,j,k)*Wib(i,j,k:k+1)))

                  ! Signed distance
                  Gib0(i,j,k)=cfg%Gib(i,j,k) 
               end do
            end do
         end do
         call cfg%sync(srcM)
      end block calc_ib_velocity
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use random,       only: random_normal
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         real(WP) :: Uin
         ! Read inflow velocity
         call param_read('Inlet velocity',Uin)
         ! Make initial velocity field random to trigger transition
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=Uin   ,sd=0.01_WP*Uin)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=0.01_WP*Uin)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=0.01_WP*Uin)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Set inflow velocity
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=Uin
         end do
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr(src=srcM)
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         resU=srcM/fs%rho           !< Careful, we need to provide
         call fs%get_div(src=resU)  !< a volume source term to div
      end block initialize_velocity

      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.1_WP
      end block create_sgs

      ! Initialize our LPT
      initialize_lpt: block
        use random, only: random_uniform
        ! Create solver
        lp=lpt(cfg=cfg,name='LPT')
        ! Get particle density from the input
        call param_read('Particle density',lp%rho)
        ! Set gravity
        call param_read('Gravity',lp%gravity)
        ! Initialize with zero particles
        call lp%resize(0)
        ! Get initial particle volume fraction
        call lp%update_VF()
        ! Collision parameters
        lp%tau_col=25.0_WP*time%dt
        ! Set coefficient of restitution
        call param_read('Coefficient of restitution',lp%e_n)
        call param_read('Wall restitution',lp%e_w)
        call param_read('Friction coefficient',lp%mu_f)
        ! Injection parameters
        call param_read('Particle mass flow rate',lp%mfr)
        call param_read('Particle velocity',lp%inj_vel)
        call param_read('Particle mean diameter',lp%inj_dmean)
        call param_read('Particle standard deviation',lp%inj_dsd,default=0.0_WP)
        call param_read('Particle min diameter',lp%inj_dmin,default=tiny(1.0_WP))
        call param_read('Particle max diameter',lp%inj_dmax,default=huge(1.0_WP))
        call param_read('Particle diameter shift',lp%inj_dshift,default=0.0_WP)
        if (lp%inj_dsd.le.epsilon(1.0_WP)) then
           lp%inj_dmin=lp%inj_dmean
           lp%inj_dmax=lp%inj_dmean
        end if
        call param_read('Particle inject diameter',lp%inj_d)
        lp%inj_pos(1)=lp%cfg%x(lp%cfg%imin)+lp%inj_dmax
        lp%inj_pos(2:3)=0.0_WP
      end block initialize_lpt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
        integer :: i
        pmesh=partmesh(nvar=3,nvec=2,name='lpt')
        pmesh%varname(1)='id'
        pmesh%varname(2)='diameter'
        pmesh%varname(3)='delta'
        pmesh%vecname(1)='velocity'
        pmesh%vecname(2)='ang_vel'
        call lp%update_partmesh(pmesh)
        do i=1,lp%np_
           pmesh%var(1,i)=real(lp%p(i)%id,WP)
           pmesh%var(2,i)=lp%p(i)%d
           pmesh%var(3,i)=lp%p(i)%delta_n
           pmesh%vec(:,1,i)=lp%p(i)%vel
           pmesh%vec(:,2,i)=lp%p(i)%angVel
        end do
      end block create_pmesh
        
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='blade')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('Gib',cfg%Gib)
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Prepare some info about fields
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call lp%get_max()
         ! Create simulation monitor
         lpfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lpfile%add_column(time%n,'Timestep number')
         call lpfile%add_column(time%t,'Time')
         call lpfile%add_column(time%dt,'Timestep size')
         call lpfile%add_column(lp%np,'Particle number')
         call lpfile%add_column(lp%np_new,'Npart new')
         call lpfile%add_column(lp%np_out,'Npart removed')
         call lpfile%add_column(lp%ncol,'Particle collisions')
         call lpfile%add_column(lp%VFmax,'Max VF')
         call lpfile%add_column(lp%Umin,'Particle Umin')
         call lpfile%add_column(lp%Umax,'Particle Umax')
         call lpfile%add_column(lp%Vmin,'Particle Vmin')
         call lpfile%add_column(lp%Vmax,'Particle Vmax')
         call lpfile%add_column(lp%Wmin,'Particle Wmin')
         call lpfile%add_column(lp%Wmax,'Particle Wmax')
         call lpfile%add_column(lp%dmin,'Particle dmin')
         call lpfile%add_column(lp%dmax,'Particle dmax')
         call lpfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   subroutine move_IB()
      use ibconfig_class, only: sharp
      implicit none
      integer :: i,j,k
      integer :: j0
      real(WP) :: y0,yd,frac
      real(WP) :: y_period,yj,real_idx
      integer  :: j1

      y_shift = y_shift + rotor_velocity * time%dt
      y0 = cfg%ym(cfg%jmin)

      do k = cfg%kmin_, cfg%kmax_
        do j = cfg%jmin_, cfg%jmax_
          do i = cfg%imin_, cfg%imax_
         
            ! Only rot the rotor
            if (cfg%xm(i).gt.0.0_WP) cycle
            
            ! Physical coordinate of current cell
            yj = y0 + (real(j - cfg%jmin, WP)) * cfg%dy(j)
            
            ! Backtraced position (accounting for periodicity)
            yd = modulo(yj - y_shift - y0, cfg%yL) + y0

            ! Fractional index location in initial field
            real_idx = (yd - y0) / cfg%dy(j) + cfg%jmin
            j0 = floor(real_idx)
            frac = real_idx - real(j0, WP)

            ! Handle periodic boundaries
            j1 = j0 + 1
            if (j1 > cfg%jmax) j1 = cfg%jmin

            ! Linear interpolation
            cfg%Gib(i,j,k) = (1.0_WP - frac)*Gib0(i,j0,k) + frac*Gib0(i,j1,k)

          end do
        end do
      end do

      call cfg%sync(cfg%Gib)

      ! Recalculate normals and volume fraction
      call cfg%calculate_normal()
      call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
   end subroutine move_IB
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! rot the rotor
         call move_IB()

         ! Inject particles
         call lp%inject(dt=time%dt,avoid_overlap=.false.)

         call lp%collide(time%dt,cfg%Gib,cfg%Nib(1,:,:,:),cfg%Nib(2,:,:,:),cfg%Nib(3,:,:,:),Uib,Vib,Wib)
         resU=fs%rho; resV=visc
         call lp%advance(time%dt,fs%U,fs%V,fs%W,resU,resv)

         subgrid: block
            integer :: i
            ! Get terms
            call fs%interp_vel(Ui,Vi,Wi)
            call fs%get_strainrate(SR=SR)
            call fs%get_gradU(gradU=gradU)
            fs%visc=visc; resU=fs%rho

            ! DSM
            call sgs%get_visc(3,time%dt,resU,Ui,Vi,Wi,SR,gradU)
            where (sgs%visc.lt.-fs%visc)
               sgs%visc=-fs%visc
            end where
            fs%visc=fs%visc+sgs%visc
         end block subgrid

         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dt*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dt*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply direct IB forcing
            ibforcing: block
               integer :: i,j,k
               real(WP) :: vf,sd
               do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
                  ! U cell
                  vf=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                  fs%U(i,j,k)=vf*fs%U(i,j,k)+(1.0_WP-vf)*Uib(i,j,k)
                  ! V cell
                  vf=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
                  fs%V(i,j,k)=vf*fs%V(i,j,k)+(1.0_WP-vf)*Vib(i,j,k)
                  ! W cell
                  vf=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
                  fs%W(i,j,k)=vf*fs%W(i,j,k)+(1.0_WP-vf)*Wib(i,j,k)
                  ! IB mass source
                  srcM(i,j,k)=fs%rho*(1.0_WP-cfg%VF(i,j,k))*(sum(fs%divp_x(:,i,j,k)*Uib(i:i+1,j,k))+&
                  &                                          sum(fs%divp_y(:,i,j,k)*Vib(i,j:j+1,k))+&
                  &                                          sum(fs%divp_z(:,i,j,k)*Wib(i,j,k:k+1)))
               end do; end do; end do
               call fs%cfg%sync(fs%U)
               call fs%cfg%sync(fs%V)
               call fs%cfg%sync(fs%W)
               call cfg%sync(srcM)
            end block ibforcing

             ! ! Apply IB forcing to enforce BC at the walls
             ! ibforcing: block
             !   use ibconfig_class, only: VFhi,VFlo
             !   integer :: i,j,k
             !   real(WP) :: vf,vol,dudn,delta, sigma
             !   real(WP) :: Cslip = sqrt(2.0_WP) ! Wall model coefficient
             !   do k=fs%cfg%kmin_,fs%cfg%kmax_
             !      do j=fs%cfg%jmin_,fs%cfg%jmax_
             !         do i=fs%cfg%imin_,fs%cfg%imax_
             !            ! U cell
             !            if (fs%umask(i,j,k).eq.0) then
             !               ! Interpolate VF to face
             !               vf=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
             !               ! Skip if not IB cell
             !               !if (vf.lt.VFlo.or.vf.gt.VFhi) cycle
             !               ! Apply wall model
             !               vol=(fs%cfg%VF(i  ,j,k)*fs%cfg%vol(i  ,j,k)+&
             !                    &    fs%cfg%VF(i-1,j,k)*fs%cfg%vol(i-1,j,k))
             !               dudn=-(fs%cfg%VF(i  ,j,k)*fs%cfg%vol(i  ,j,k)*sum(gradU(:,1,i  ,j,k)*cfg%Nib(:,i  ,j,k))+&
             !                    &      fs%cfg%VF(i-1,j,k)*fs%cfg%vol(i-1,j,k)*sum(gradU(:,1,i-1,j,k)*cfg%Nib(:,i-1,j,k)))/vol
             !               delta=(vol)**(1.0_WP/3.0_WP)
             !               ! Apply IB forcing
             !               fs%U(i,j,k)=vf*fs%U(i,j,k)+(1.0_WP-vf)*Cslip*delta*dudn
             !               ! Store IB velocities
             !               Uib(i,j,k)=Cslip*dudn*delta
             !            end if
             !            ! V cell
             !            if (fs%vmask(i,j,k).eq.0) then
             !               ! Interpolate VF to face
             !               !if (vf.lt.VFlo.or.vf.gt.VFhi) cycle
             !               ! Skip if not IB cell
             !               vf=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))
             !               ! Apply wall model
             !               vol=(fs%cfg%VF(i,j  ,k)*fs%cfg%vol(i,j  ,k)+&
             !                    &    fs%cfg%VF(i,j-1,k)*fs%cfg%vol(i,j-1,k))
             !               dudn=-(fs%cfg%VF(i,j  ,k)*fs%cfg%vol(i,j  ,k)*sum(gradU(:,2,i,j  ,k)*cfg%Nib(:,i,j  ,k))+&
             !                    & fs%cfg%VF(i,j-1,k)*fs%cfg%vol(i,j-1,k)*sum(gradU(:,2,i,j-1,k)*cfg%Nib(:,i,j-1,k)))/vol
             !               delta=(vol)**(1.0_WP/3.0_WP)
             !               ! Apply IB forcing
             !               fs%V(i,j,k)=vf*fs%V(i,j,k)+(1.0_WP-vf)*Cslip*delta*dudn
             !               ! Store IB velocities
             !               Vib(i,j,k)=Cslip*dudn*delta
             !            
             !            end if
             !            ! W cell
             !            if (fs%wmask(i,j,k).eq.0) then
             !               ! Interpolate VF to face
             !               !if (vf.lt.VFlo.or.vf.gt.VFhi) cycle
             !               ! Skip if not IB cell
             !               vf=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))
             !               ! Apply wall model
             !               vol=(fs%cfg%VF(i,j,k  )*fs%cfg%vol(i,j,k  )+&
             !                    &    fs%cfg%VF(i,j,k-1)*fs%cfg%vol(i,j,k-1))
             !               dudn=-(fs%cfg%VF(i,j,k  )*fs%cfg%vol(i,j,k  )*sum(gradU(:,3,i,j,k  )*cfg%Nib(:,i,j,k  ))+&
             !                     &      fs%cfg%VF(i,j,k-1)*fs%cfg%vol(i,j,k-1)*sum(gradU(:,3,i,j,k-1)*cfg%Nib(:,i,j,k-1)))/vol
             !               delta=(vol)**(1.0_WP/3.0_WP)
             !               ! Apply IB forcing
             !               fs%W(i,j,k)=vf*fs%W(i,j,k)+(1.0_WP-vf)*Cslip*delta*dudn
             !               ! Store IB velocities
             !               Wib(i,j,k)=Cslip*dudn*delta
             !            end if
             !         end do
             !      end do
             !   end do
             !   call fs%cfg%sync(fs%U)
             !   call fs%cfg%sync(fs%V)
             !   call fs%cfg%sync(fs%W)
             ! end block ibforcing
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr(src=srcM)
            resU=srcM/fs%rho           !< Careful, we need to provide
            call fs%get_div(src=resU)  !< a volume source term to div
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho
            fs%V=fs%V-time%dt*resV/fs%rho
            fs%W=fs%W-time%dt*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         resU=srcM/fs%rho           !< Careful, we need to provide
         call fs%get_div(src=resU)  !< a volume source term to div
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              integer :: i
              call lp%update_partmesh(pmesh)
              do i=1,lp%np_
                 pmesh%var(1,i)=real(lp%p(i)%id,WP)
                 pmesh%var(2,i)=lp%p(i)%d
                 pmesh%var(3,i)=lp%p(i)%delta_n
                 pmesh%vec(:,1,i)=lp%p(i)%vel
                 pmesh%vec(:,2,i)=lp%p(i)%angVel
              end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
         call lp%get_max()
         call lpfile%write()
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,Uib,Vib,Wib,srcM,gradU)
      
   end subroutine simulation_final
   
   
end module simulation
