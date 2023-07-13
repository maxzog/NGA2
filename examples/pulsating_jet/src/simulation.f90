!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
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
   type(incomp),      public :: fs
   type(ddadi),       public :: vs
   type(timetracker), public :: time
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu
   
   !> Fluid, forcing, and particle parameters
   real(WP) :: visc
   real(WP) :: stk,tau_eta
   logical  :: use_sgs
   integer  :: sgs_type

   !> Jet parameters
   real(WP) :: Djet
   
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
      
      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg,pfmg,gmres_smg
         use incomp_class,    only: clipped_neumann, dirichlet, slip
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Get jet diameter
         call param_read('Jet diameter',Djet,default=0.1_WP)
         ! Add BCs
         call fs%add_bcond(name='jet',   type=dirichlet,locator=jet_loc,face='x',dir=-1,canCorrect=.false. )
         call fs%add_bcond(name='coflow',type=dirichlet,locator=coflow_loc,face='x',dir=-1,canCorrect=.false. )
         call fs%add_bcond(name='right', type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true.)
         call fs%add_bcond(name='bottom',type=slip, locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top',   type=slip, locator=top_of_domain,face='y',dir=+1,canCorrect=.false.)
         call fs%add_bcond(name='front', type=slip, locator=front_of_domain,face='z',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='back',  type=slip, locator=back_of_domain,face='z',dir=+1,canCorrect=.false.)
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Prepare and configure pressure solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps, implicit_solver=vs)
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
         call fs%get_bcond('coflow', mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=0.05_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP 
         end do
         call fs%get_bcond('jet', mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=0.5_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP 
         end do

         call fs%get_mfr()
         call fs%correct_mfr()
         call fs%interp_vel(Ui,Vi,Wi)
      end block initialize_velocity
      
!!$      ! Initialize LPT solver
!!$      initialize_lpt: block
!!$         use random, only: random_uniform, random_normal
!!$         real(WP) :: dp
!!$         integer :: i,np
!!$         character(len=str_medium) :: timestamp
!!$         ! Create solver
!!$         lp=lpt(cfg=cfg,name='LPT')
!!$         ! Get drag model from the inpit
!!$         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
!!$         ! Get particle density from the input
!!$         call param_read('Particle density',lp%rho)
!!$         ! Get particle diameter from the input
!!$         call param_read('Particle Stokes number', stk)
!!$         ! Get number of particles
!!$         call param_read('Number of particles',np)
!!$         ! Root process initializes np particles randomly
!!$         if (lp%cfg%amRoot) then
!!$            tau_eta = 1.0_WP !sqrt(visc/EPS0)
!!$            dp = 1.0_WP !sqrt(18.0_WP*visc*stk*tau_eta/lp%rho)
!!$            call lp%resize(np)
!!$            do i=1,np
!!$               ! Give id
!!$               lp%p(i)%id=int(i,8)
!!$               ! Set the diameter
!!$               lp%p(i)%d=dp
!!$               ! Assign random position in the domain
!!$               lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
!!$                    &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
!!$                    &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
!!$               ! Give zero velocity
!!$               lp%p(i)%vel=0.0_WP
!!$               ! Give zero dt
!!$               lp%p(i)%dt=0.0_WP
!!$               ! temp
!!$               lp%p(i)%uf=0.0_WP
!!$               ! Locate the particle on the mesh
!!$               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
!!$               ! Activate the particle
!!$               lp%p(i)%flag=0
!!$            end do
!!$            ! Distribute particles
!!$            call lp%sync()
!!$         end if
!!$      end block initialize_lpt
      

!!$      ! Create partmesh object for Lagrangian particle output
!!$      create_pmesh: block
!!$         integer :: np,i
!!$         call param_read('Number of particles',np)
!!$         pmesh=partmesh(nvar=1,nvec=2,name='lpt')
!!$         pmesh%varname(1)="id"
!!$         pmesh%vecname(1)="vel"
!!$         pmesh%vecname(2)="fld"
!!$         call lp%resize(np)
!!$         call lp%update_partmesh(pmesh)
!!$         do i = 1,lp%np_
!!$            pmesh%var(1,i) = lp%p(i)%id
!!$            pmesh%vec(:,1,i) = lp%p(i)%vel
!!$            pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
!!$         end do
!!$      end block create_pmesh
      
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
         !call ens_out%add_particle('particles',pmesh)
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
         
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

!!$         ! Advance particles by dt
!!$         resU=fs%rho; resV=fs%visc
!!$         if (stk.gt.0.0_WP) then
!!$            call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV)
!!$         else
!!$            call lp%advance_tracer(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV)
!!$         end if
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Turbulence modeling
         call fs%interp_vel(Ui,Vi,Wi)
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

            ! Implicit velocity solve
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! From Monroe et al. 2021
            ! Q(t) = Q_0 * |exp{-t/tau} * sin(omega * t)|

!!$            !> Time-varying boundary condition
!!$            dirichlet_velocity: block
!!$               use incomp_class, only:bcond
!!$               type(bcond), pointer :: mybc
!!$               integer :: n,i,j,k
!!$               real(WP) :: ratio, tlin
!!$               logical :: ramp_up,ramp_down,ramp_stop
!!$               tlin=1.0_WP
!!$               ratio = time%t / tlin
!!$               ramp_up=.false.; ramp_down=.false.; ramp_stop=.false.
!!$               if (ratio.lt.1.0_WP) ramp_up = .true.
!!$               if (ratio.gt.1.0_WP) ramp_down = .true.
!!$               if (ratio.gt.1.5_WP) ramp_stop = .true.
!!$               call fs%get_bcond('jet', mybc)
!!$               do n=1,mybc%itr%no_
!!$                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
!!$                  fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP;fs%W(i,j,k)=0.1_WP
!!$                  if (time%t.lt.1.0_WP) then
!!$                     fs%W(i,j,k) = 0.75_WP
!!$                  else
!!$                     fs%W(i,j,k) = 0.1_WP
!!$                  end if
!!$                  ! elseif (ramp_down.and..not.ramp_stop) then
!!$                  !    fs%W(i,j,k) = 0.1_WP
!!$                  !    if (fs%W(i,j,k).lt.0.0_WP) fs%W(i,j,k) = 0.1_WP
!!$                  ! else 
!!$                  !    fs%W(i,j,k) = 0.1_WP
!!$                  ! end if
!!$               end do
!!$             end block dirichlet_velocity
             
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)

            ! Solve Poisson equation
            call fs%get_mfr()
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
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
         call fs%get_div()
         ! Output to ensight
         if (ens_evt%occurs()) then
!!$            call lp%update_partmesh(pmesh)
!!$            do ii = 1,lp%np_
!!$               pmesh%var(1,ii) = lp%p(ii)%id
!!$               pmesh%vec(:,1,ii) = lp%p(ii)%vel
!!$               pmesh%vec(:,2,ii) = lp%cfg%get_velocity(pos=lp%p(ii)%pos,i0=lp%p(ii)%ind(1),j0=lp%p(ii)%ind(2),k0=lp%p(ii)%ind(3),U=fs%U,V=fs%V,W=fs%W)
!!$            end do
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         !call lp%get_max()
         call mfile%write()
         call cflfile%write()
         
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
