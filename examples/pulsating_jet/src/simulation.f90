!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   !use fourier3d_class,   only: fourier3d
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use sgsmodel_class,    only: sgsmodel
   use datafile_class,    only: datafile
   use string,            only: str_medium
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   type(hypre_str),   public :: ps
   !type(fourier3d),   public :: ps
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   logical :: restarted
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: visc,meanU,meanV,meanW
   real(WP) :: stk,tau_eta
   logical  :: use_sgs
   integer  :: sgs_type

   !> For monitoring
   ! real(WP) :: EPS,sgsEPSp
   ! real(WP) :: Re_L,Re_lambda
   ! real(WP) :: eta,ell
   ! real(WP) :: dx_eta,ell_Lx,Re_ratio,eps_ratio,tke_ratio,nondtime
   ! real(WP) :: ftvar,fvar,usvar,pvvar,varratio
   ! real(WP), dimension(3) :: fmean,pvmean,usmean

   ! !> Wallclock time for monitoring
   ! type :: timer
   !    real(WP) :: time_in
   !    real(WP) :: time
   !    real(WP) :: percent
   ! end type timer
   ! type(timer) :: wt_total,wt_vel,wt_pres,wt_lpt,wt_rest,wt_sgs,wt_stat,wt_force
   
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
   
   !> Function that localizes the jet inflow
   function jet_loc(pg,i,j,k) result(isIn)
     use pgrid_class, only: pgrid
     implicit none
     class(pgrid), intent(in) :: pg
     integer, intent(in) :: i,j,k
     logical :: isIn, inRadius
     real(WP) :: c_x, c_y  ! Center of jet
     real(WP) :: dist 
     isIn = .false.; inRadius = .false.
     c_x = 0.0_WP; c_y = 0.0_WP
     dist = sqrt((c_x - pg%x(i))**2 + (c_y - pg%y(j))**2) 
     !dist = abs(c_x - pg%x(i))
     if (dist.lt.0.1_WP) inRadius = .true.
     if (k.eq.pg%kmin.and.inRadius) isIn = .true. 
   end function jet_loc
   
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
         allocate(SR  (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradu(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max iter',time%nmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Handle restart/saves here
      restart_and_save: block
        character(len=str_medium) :: timestamp
        ! Create event for saving restart files
        save_evt=event(time,'Restart output')
        call param_read('Restart output period',save_evt%tper)
        ! Check if we are restarting
        call param_read(tag='Restart from',val=timestamp,short='r',default='')
        restarted=.false.; if (len_trim(timestamp).gt.0) restarted=.true.
        if (restarted) then
           ! If we are, read the name of the directory
           call param_read('Restart from',timestamp,'r')
           ! Read the datafile
           df=datafile(pg=cfg,fdata='restart/data_'//trim(adjustl(timestamp)))
        else
           ! If we are not restarting, we will still need a datafile for saving restart files
           df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=4)
           df%valname(1)='t'
           df%valname(2)='dt'
           df%varname(1)='U'
           df%varname(2)='V'
           df%varname(3)='W'
           df%varname(4)='P'
        end if
      end block restart_and_save

      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (restarted) then
            call df%pullval(name='t' ,val=time%t )
            call df%pullval(name='dt',val=time%dt)
            time%told=time%t-time%dt
         end if
      end block update_timetracker

      ! ! Initialize timers
      ! initialize_timers: block
      !    wt_total%time=0.0_WP; wt_total%percent=0.0_WP
      !    wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
      !    wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
      !    wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
      !    wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
      !    wt_sgs%time=0.0_WP;   wt_sgs%percent=0.0_WP
      !    wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
      !    wt_force%time=0.0_WP; wt_force%percent=0.0_WP
      ! end block initialize_timers
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg,pfmg,gmres_smg
         use incomp_class,    only: clipped_neumann, dirichlet, slip
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Add BCs
         call fs%add_bcond(name='left',  type=slip, locator=left_of_domain,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='right', type=slip, locator=right_of_domain,face='x',dir=+1,canCorrect=.false.)
         call fs%add_bcond(name='bottom',type=slip, locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top',   type=slip, locator=top_of_domain,face='y',dir=+1,canCorrect=.false. )
         call fs%add_bcond(name='front', type=clipped_neumann,locator=front_of_domain,face='z',dir=+1,canCorrect=.false. )
         call fs%add_bcond(name='back',  type=dirichlet,locator=back_of_domain,face='z',dir=-1,canCorrect=.false. )
         call fs%add_bcond(name='jet',   type=dirichlet,locator=jet_loc,face='z',dir=-1,canCorrect=.false. )
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Prepare and configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ! ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         !ps=fourier3d(cfg=cfg,name='Pressure',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
      end block create_and_initialize_flow_solver

      ! Create an LES model
      create_sgs: block
         call param_read(tag='Use SGS model',val=use_sgs,default=.false.)
         if (use_sgs) call param_read('SGS model type',sgs_type)
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.1_WP
      end block create_sgs

      ! Prepare initial velocity field
      initialize_velocity: block
         use incomp_class, only: bcond
         use random,   only: random_normal
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         use mathtools,only: Pi
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         ! Read in forcing, grid, and initial velocity field parameters
         ! Zero initial velocity field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         if (restarted) then
            call df%pullvar(name='U',var=fs%U)
            call df%pullvar(name='V',var=fs%V)
            call df%pullvar(name='W',var=fs%W)
            call df%pullvar(name='P',var=fs%P)
         else
            call fs%get_bcond('back', mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.1_WP 
            end do
           call fs%get_bcond('jet', mybc)
           do n=1,mybc%itr%no_
              i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
              fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.1_WP 
           end do
         end if
        
         meanU=0.0_WP; meanV=0.0_WP; meanW=0.0_WP
         call fs%get_mfr()
         call fs%correct_mfr()
         call fs%interp_vel(Ui,Vi,Wi)

         ! Compute turbulence stats
         ! call compute_stats()
      end block initialize_velocity
      
      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform, random_normal
         real(WP) :: dp
         integer :: i,np
         character(len=str_medium) :: timestamp
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle Stokes number', stk)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Check if a stochastic SGS model is used
         if (restarted) then
            call param_read('Restart from',timestamp,'r')
            ! Read the part file
            call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
         else
         ! Root process initializes np particles randomly
            if (lp%cfg%amRoot) then
               tau_eta = 1.0_WP !sqrt(visc/EPS0)
               dp = 1.0_WP !sqrt(18.0_WP*visc*stk*tau_eta/lp%rho)
               call lp%resize(np)
               do i=1,np
                  ! Give id
                  lp%p(i)%id=int(i,8)
                  ! Set the diameter
                  lp%p(i)%d=dp
                  ! Assign random position in the domain
                  lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
                  &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
                  &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
                  ! Give zero velocity
                  lp%p(i)%vel=0.0_WP
                  ! Give zero dt
                  lp%p(i)%dt=0.0_WP
                  ! temp
                  lp%p(i)%uf=0.0_WP
                  ! Locate the particle on the mesh
                  lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  ! Activate the particle
                  lp%p(i)%flag=0
               end do
            end if
            ! Distribute particles
            call lp%sync()
         end if
      end block initialize_lpt
      

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: np,i
         call param_read('Number of particles',np)
         pmesh=partmesh(nvar=1,nvec=2,name='lpt')
         pmesh%varname(1)="id"
         pmesh%vecname(1)="vel"
         pmesh%vecname(2)="fld"
         call lp%resize(np)
         call lp%update_partmesh(pmesh)
         do i = 1,lp%np_
            pmesh%var(1,i) = lp%p(i)%id
            pmesh%vec(:,1,i) = lp%p(i)%vel
            pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
         end do
      end block create_pmesh
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='HIT')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('visc',sgs%visc)
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
         call mfile%add_column(meanU,'Umean')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(meanV,'Vmean')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(meanW,'Wmean')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
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
   
   
   !> Time integrate our problem
   subroutine simulation_run
      use parallel,       only: parallel_time
      implicit none
      integer :: ii
      ! Perform time integration
      do while (.not.time%done())
         ! init wallclock
         ! wt_total%time_in=parallel_time()
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! wt_lpt%time_in=parallel_time()
         ! Advance particles by dt
         resU=fs%rho; resV=fs%visc
         if (stk.gt.0.0_WP) then
            call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV)
         else
            call lp%advance_tracer(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV)
         end if
         ! wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Turbulence modeling
         call fs%interp_vel(Ui,Vi,Wi)
         ! wt_sgs%time_in=parallel_time()
         if (use_sgs) then
            fs%visc=visc
            call fs%get_strainrate(SR=SR)
            call fs%get_gradu(gradu=gradu)
            resU=fs%rho
            call sgs%get_visc(type=sgs_type,rho=resU,gradu=gradu,dt=time%dt,SR=SR,Ui=Ui,Vi=Vi,Wi=Wi)
            where (sgs%visc.lt.-fs%visc)
               sgs%visc=-fs%visc
            end where
            fs%visc=fs%visc+sgs%visc
         end if
         ! wt_sgs%time=wt_sgs%time+parallel_time()-wt_sgs%time_in
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            ! wt_vel%time_in=parallel_time()
            
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
            
            ! wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in
            ! wt_force%time_in=parallel_time()
            ! wt_force%time=wt_force%time+parallel_time()-wt_force%time_in
            ! wt_vel%time_in=parallel_time()

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            
            
            
            dirichlet_velocity: block
               use incomp_class, only:bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k
               real(WP) :: ratio, tlin
               logical :: ramp_up,ramp_down,ramp_stop
               tlin=1.0_WP
               ratio = time%t / tlin
               ramp_up=.false.; ramp_down=.false.; ramp_stop=.false.
               if (ratio.lt.1.0_WP) ramp_up = .true.
               if (ratio.gt.1.0_WP) ramp_down = .true.
               if (ratio.gt.1.5_WP) ramp_stop = .true.
               call fs%get_bcond('back', mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.1_WP
               end do
               call fs%get_bcond('jet', mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP;fs%W(i,j,k)=0.1_WP
                  if (time%t.lt.1.0_WP) then
                     fs%W(i,j,k) = 0.75_WP
                  else
                     fs%W(i,j,k) = 0.1_WP
                  end if
                  ! elseif (ramp_down.and..not.ramp_stop) then
                  !    fs%W(i,j,k) = 0.1_WP
                  !    if (fs%W(i,j,k).lt.0.0_WP) fs%W(i,j,k) = 0.1_WP
                  ! else 
                  !    fs%W(i,j,k) = 0.1_WP
                  ! end if
               end do
            end block dirichlet_velocity
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            ! wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in
            
            ! Solve Poisson equation
            ! wt_pres%time_in=parallel_time()
            call fs%get_mfr()
            call fs%correct_mfr()
            call fs%get_div()
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
            ! wt_pres%time=wt_pres%time+parallel_time()-wt_pres%time_in
            
            ! Increment sub-iteration counter
            time%it=time%it+1
   
         end do
         
         ! wt_vel%time_in=parallel_time()
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! wt_vel%time=wt_vel%time+parallel_time()-wt_vel%time_in
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            do ii = 1,lp%np_
               pmesh%var(1,ii) = lp%p(ii)%id
               pmesh%vec(:,1,ii) = lp%p(ii)%vel
               pmesh%vec(:,2,ii) = lp%cfg%get_velocity(pos=lp%p(ii)%pos,i0=lp%p(ii)%ind(1),j0=lp%p(ii)%ind(2),k0=lp%p(ii)%ind(3),U=fs%U,V=fs%V,W=fs%W)
            end do
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()

         ! ! Monitor timing
         ! wt_total%time=parallel_time()-wt_total%time_in
         ! wt_vel%percent=wt_vel%time/wt_total%time*100.0_WP
         ! wt_pres%percent=wt_pres%time/wt_total%time*100.0_WP
         ! wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
         ! wt_sgs%percent=wt_sgs%time/wt_total%time*100.0_WP
         ! wt_stat%percent=wt_stat%time/wt_total%time*100.0_WP
         ! wt_force%percent=wt_force%time/wt_total%time*100.0_WP
         ! wt_rest%time=wt_total%time-wt_vel%time-wt_pres%time-wt_lpt%time-wt_sgs%time-wt_stat%time-wt_force%time
         ! wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
         ! call tfile%write()
         ! wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         ! wt_vel%time=0.0_WP;   wt_vel%percent=0.0_WP
         ! wt_pres%time=0.0_WP;  wt_pres%percent=0.0_WP
         ! wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         ! wt_sgs%time=0.0_WP;   wt_sgs%percent=0.0_WP
         ! wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
         ! wt_force%time=0.0_WP; wt_force%percent=0.0_WP
         ! wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              character(len=str_medium) :: timestamp
              ! Prefix for files
              write(timestamp,'(es12.5)') time%t
              ! Prepare a new directory
              if (fs%cfg%amRoot) call execute_command_line('mkdir -p restart')
              ! Populate df and write it
              call df%pushval(name='t' ,val=time%t     )
              call df%pushval(name='dt',val=time%dt    )
              call df%pushvar(name='U' ,var=fs%U       )
              call df%pushvar(name='V' ,var=fs%V       )
              call df%pushvar(name='W' ,var=fs%W       )
              call df%pushvar(name='P' ,var=fs%P       )
              call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
              ! Write particle file
              call lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
            end block save_restart
         end if
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
