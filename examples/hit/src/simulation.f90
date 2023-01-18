!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   !use pfft3d_class,      only: pfft3d
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp
   use lpt_class,         only: lpt
   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Single-phase incompressible flow solver, pressure and implicit solvers, and a time tracker
   !type(pfft3d),      public :: ps
   type(hypre_str),   public :: ps
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,hitfile,lptfile,sgsfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi,Lambda,EPSArr
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu

   !> Fluid and forcing parameters
   real(WP) :: visc
   real(WP) :: Urms0,KE0,KE,EPS,Re_L,Re_lambda,eta,Re_dom
   real(WP) :: Uvar,Vvar,Wvar,TKE,URMS,ell,sgsTKE
   real(WP) :: meanvisc,Lx,tau_eddy, tau ! tau_eddy is input, tau is calculated from tke and epsilon
   logical  :: linforce

   !> OU
   logical  :: use_crw

contains
   
   
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
         allocate(Lambda  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(EPSArr  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
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
         use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Prepare and configure pressure solver
         !ps=pfft3d(cfg=cfg,name='Pressure')
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         call param_read('Linear forcing',linforce)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)
      end block create_and_initialize_flow_solver


      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.17_WP
      end block create_sgs
      
      ! Prepare initial velocity field
      initialize_velocity: block
         use random,   only: random_normal
         use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
         use parallel, only: MPI_REAL_WP
         integer :: i,j,k,ierr
         real(WP) :: myKE
         ! Read in velocity rms and forcing time scale
         call param_read('Initial rms',Urms0)
         call param_read('Eddy turnover time',tau_eddy)
         ! Gaussian initial field
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Project to ensure divergence-free
         call fs%get_div()
         fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU/fs%rho
         fs%V=fs%V-time%dt*resV/fs%rho
         fs%W=fs%W-time%dt*resW/fs%rho
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         ! Calculate KE0
         myKE=0.0_WP
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  myKE=myKE+0.5_WP*(Ui(i,j,k)**2+Vi(i,j,k)**2+Wi(i,j,k)**2)!*fs%cfg%vol(i,j,k)
               end do
            end do
         end do
         call MPI_ALLREDUCE(myKE,KE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)!; KE = KE/fs%cfg%vol_total
         ! Read in Lx for calculation of KE0 - Assumes ell = 1/5 * Lx (Carrol & Blanquart 2013)
         call param_read('Lx', Lx)
         if (linforce) then
            KE0 = KE
         else
            KE0 = 0.0_WP
         end if

         Re_L = 0.0_WP
         Re_dom = 0.0_WP
         Re_lambda = 0.0_WP
         eta = 0.0_WP
         TKE = 0.0_WP
         meanvisc = 0.0_WP
         sgsTKE = 0.0_WP
         tau = 0.0_WP
         URMS = 0.0_WP
         ell = 0.0_WP
         EPS = 0.0_WP
         Lambda = 0.0_WP
         EPSArr = 0.0_WP

      end block initialize_velocity
      
      
      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform, random_normal
         real(WP) :: dp
         integer :: i,np
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter',dp)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Check if a stochastic SGS model is used
         call param_read('Use CRW', use_crw)
         ! Root process initializes np particles randomly
         if (lp%cfg%amRoot) then
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
               ! OU
               lp%p(i)%uf=0.0_WP
               if (use_crw) then
                  lp%p(i)%uf= [random_normal(m=0.0_WP,sd=0.1_WP),& 
                               random_normal(m=0.0_WP,sd=0.1_WP),&    
                               random_normal(m=0.0_WP,sd=0.1_WP)] 
               end if
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               lp%p(i)%a_crw=0.0_WP
               lp%p(i)%b_crw=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
      end block initialize_lpt
      

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: np,i
         call param_read('Number of particles',np)
         pmesh=partmesh(nvar=1,nvec=3,name='lpt')
         pmesh%varname(1)="id"
         pmesh%vecname(1)="vel"
         pmesh%vecname(2)="fld"
         pmesh%vecname(3)="uf"
         call lp%resize(np)
         call lp%update_partmesh(pmesh)
         do i = 1,lp%np_
            pmesh%var(1,i) = lp%p(i)%id
            pmesh%vec(:,1,i) = lp%p(i)%vel
            pmesh%vec(:,2,i) = lp%p(i)%uf + lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=Ui,V=Vi,W=Wi)
            pmesh%vec(:,3,i) = lp%p(i)%uf
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
         call ens_out%add_scalar('visc',fs%visc)
         call ens_out%add_scalar('lambda',Lambda)
         call ens_out%add_scalar('eps',EPSArr)
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
         call mfile%add_column(KE,'Kinetic energy')
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
         ! Create hit monitor
         hitfile=monitor(fs%cfg%amRoot,'hit')
         call hitfile%add_column(time%n,'Timestep number')
         call hitfile%add_column(time%t,'Time')
         call hitfile%add_column(Re_L,'Re_L')
         call hitfile%add_column(Re_dom,'Re_dom')
         call hitfile%add_column(Re_lambda,'Re_lambda')
         call hitfile%add_column(meanvisc,'mean visc') 
         call hitfile%add_column(eta,'eta')
         call hitfile%add_column(sgsTKE,'sgsTKE')
         call hitfile%add_column(TKE,'TKE')
         call hitfile%add_column(tau,'tau_eddy')
         call hitfile%add_column(URMS,'URMS')
         call hitfile%add_column(EPS,'Epsilon')
         call hitfile%add_column(ell,'Integral lengthscale')
         call hitfile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
         ! Create SGS monitor
         sgsfile=monitor(fs%cfg%amroot,'sgs')
         call sgsfile%add_column(time%n,'Timestep number')
         call sgsfile%add_column(time%t,'Time')
         call sgsfile%add_column(sgs%min_visc,'Min eddy visc')
         call sgsfile%add_column(sgs%max_visc,'Max eddy visc')
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Time integrate our problem
   subroutine simulation_run
      use sgsmodel_class, only: vreman, constant_smag, off
      implicit none
      integer :: ii
      ! Perform time integration
      do while (.not.time%done())
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         ! Advance particles by dt
         resU=fs%rho; resV=fs%visc
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV,use_crw=use_crw,sgs=sgs)
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         ! Reset here fluid properties
         fs%visc=visc
         ! Turbulence modeling
         call fs%get_strainrate(SR=SR)
         call fs%get_gradu(gradu=gradu)
         resU=fs%rho
         call sgs%get_visc(type=off,rho=resU,gradu=gradu)
         where (sgs%visc.lt.-fs%visc)
            sgs%visc=-fs%visc
         end where
         fs%visc=fs%visc+sgs%visc

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
            
            ! Add linear forcing term
            linear_forcing: block
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
               use parallel, only: MPI_REAL_WP
               real(WP) :: meanU,meanV,meanW,myKE
               real(WP) :: myUvar,myVvar,myWvar,mysgsTKE
               integer :: i,j,k,ierr
               ! Calculate KE
               call fs%interp_vel(Ui,Vi,Wi)
               
               ! Calculate mean velocity
               call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
               
               myUvar=0.0_WP
               mysgstke=0.0_WP
               myVvar=0.0_WP
               myWvar=0.0_WP
               myKE=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        myKE=myKE+0.5_WP*(Ui(i,j,k)**2+Vi(i,j,k)**2+Wi(i,j,k)**2)
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myKE,KE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)

               ! Add forcing term
               resU=resU+time%dt*KE0/KE*(fs%U-meanU)/(2.0_WP*tau_eddy)
               resV=resV+time%dt*KE0/KE*(fs%V-meanV)/(2.0_WP*tau_eddy)
               resW=resW+time%dt*KE0/KE*(fs%W-meanW)/(2.0_WP*tau_eddy)
            end block linear_forcing
            

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho
            
            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
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
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         sgsTKE = 0.0_WP 
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
            
         
         ! Compute turbulence statistics for monitor
         compute_stats: block 
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
               use parallel, only: MPI_REAL_WP
               real(WP) :: meanU,meanV,meanW,meanSR2
               real(WP) :: myUvar,myVvar,myWvar,mySR2
               real(WP) :: mysgstke,myvisc
               integer :: i,j,k,ierr
               
               call fs%interp_vel(Ui,Vi,Wi)
               call fs%get_strainrate(SR=SR)

               ! Calculate mean velocity
               call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
               call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total

               mySR2=0.0_WP
               myUvar=0.0_WP
               myVvar=0.0_WP
               myWvar=0.0_WP
               mysgstke=0.0_WP
               myvisc=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        mysgstke=mysgstke + fs%cfg%vol(i,j,k)*(sgs%visc(i,j,k)/0.067_WP/sgs%delta(i,j,k)/fs%rho)**2
                        myUvar=myUvar+(Ui(i,j,k)-meanU)**2*fs%cfg%vol(i,j,k)
                        myVvar=myVvar+(Wi(i,j,k)-meanW)**2*fs%cfg%vol(i,j,k)
                        myWvar=myWvar+(Wi(i,j,k)-meanV)**2*fs%cfg%vol(i,j,k)
                        mySR2=mySR2+(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)+2*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))*fs%cfg%vol(i,j,k)
                        myvisc=myvisc+fs%visc(i,j,k)*fs%cfg%vol(i,j,k)

                        ! Epsilon - on grid
                        EPSArr(i,j,k) = 2.0_WP*fs%visc(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)+2*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))/fs%rho

                        ! Lambda  - on grid
                        Lambda(i,j,k) = sqrt(15.0_WP*fs%visc(i,j,k)/EPSArr(i,j,k))
                     end do
                  end do
               end do

               call MPI_ALLREDUCE(myUvar,Uvar,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Uvar=Uvar/fs%cfg%vol_total
               call MPI_ALLREDUCE(myVvar,Vvar,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Vvar=Vvar/fs%cfg%vol_total
               call MPI_ALLREDUCE(myWvar,Wvar,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); Wvar=Wvar/fs%cfg%vol_total
               call MPI_ALLREDUCE(mysgstke,sgsTKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); sgsTKE=sgsTKE/fs%cfg%vol_total
               call MPI_ALLREDUCE(myvisc,meanvisc,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanvisc=meanvisc/fs%cfg%vol_total
               ! Calculate <s_ij s_ij>
               call MPI_ALLREDUCE(mySR2,meanSR2,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanSR2=meanSR2/fs%cfg%vol_total
               TKE = 0.5_WP*(Uvar+Vvar+Wvar)
               URMS = sqrt(2.0_WP/3.0_WP*(TKE+sgsTKE))
               Lambda = Lambda*URMS
               Re_dom = fs%Umax*Lx/meanvisc
               EPS = 2.0_WP*meanvisc*meanSR2/fs%rho
               Re_L = (TKE+sgsTKE)**2.0_WP/EPS/meanvisc 
               Re_lambda = sqrt(20.0_WP*Re_L/3.0_WP)
               eta = (meanvisc**3.0_WP/EPS)**0.25_WP
               ell = 1.5_WP**1.5_WP * URMS**3.0_WP / EPS
               tau  = (TKE+sgsTKE) / EPS

         end block compute_stats
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            do ii = 1,lp%np_
               pmesh%var(1,ii) = lp%p(ii)%id
               pmesh%vec(:,1,ii) = lp%p(ii)%vel
               pmesh%vec(:,2,ii) = lp%p(ii)%uf + lp%cfg%get_velocity(pos=lp%p(ii)%pos,i0=lp%p(ii)%ind(1),j0=lp%p(ii)%ind(2),k0=lp%p(ii)%ind(3),U=Ui,V=Vi,W=Wi)
            end do
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call hitfile%write()
         call lp%get_max()
         call lptfile%write()
         call sgsfile%write()
         
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR,gradu,Lambda,EPSArr)
      
   end subroutine simulation_final
   
   
   
   
   
end module simulation
