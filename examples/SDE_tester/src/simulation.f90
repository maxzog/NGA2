!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp
   use randomwalk_class,  only: lpt
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
   type(monitor) :: mfile,hitfile,lptfile,sgsfile,tfile,ssfile

   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi

   !> Fluid, forcing, and particle parameters
   real(WP) :: visc,meanU,meanV,meanW

   !> Wallclock time for monitoring
   type :: timer
      real(WP) :: time_in
      real(WP) :: time
      real(WP) :: percent
   end type timer
   type(timer) :: wt_total,wt_lpt,wt_rest,wt_stat
   
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
      end block allocate_work_arrays

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max iter',time%nmax)
         call param_read('Max time',time%tmax)
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

      ! Initialize timers
      initialize_timers: block
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP
         wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
      end block initialize_timers

      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
         ! Prepare and configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pfmg,nst=7)
         call fs%setup(pressure_solver=ps)
      end block create_and_initialize_flow_solver

      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.1_WP
         sgs%visc=0.005_WP
      end block create_sgs

      ! Prepare initial velocity field
      initialize_velocity: block
         fs%U = 0.0_WP
         fs%V = 0.0_WP
         fs%W = 0.0_WP

         ! Compute mean and remove it from the velocity field to obtain <U>=0
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
 
         fs%U = fs%U - meanU
         fs%V = fs%V - meanV
         fs%W = fs%W - meanW

         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)

         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
      end block initialize_velocity


      ! Initialize LPT solver
      initialize_lpt: block
         use random, only: random_uniform, random_normal
         integer :: i,np
         real(WP) :: dp
         character(len=str_medium) :: timestamp
         ! Create solver
         lp=lpt(cfg=cfg,name='CRW')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter from the input
         call param_read('Particle diameter', dp)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Check if spatially correlated
         call param_read('Correlation function', lp%corr_type)
         ! Check if a particles should be read in
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
               lp%p(i)%us=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Init Wiener increment
               lp%p(i)%dW=0.0_WP
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
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
            pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
            pmesh%vec(:,3,i) = lp%p(i)%us 
         end do
      end block create_pmesh

      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='SDE')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         time%cfl=0.9_WP
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
         call mfile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%Usvar,'Seen Uvar')
         call lptfile%add_column(lp%Vsvar,'Seen Vvar')
         call lptfile%add_column(lp%Wsvar,'Seen Wvar')
         call lptfile%add_column(lp%Usmin,'Seen Umin')
         call lptfile%add_column(lp%Usmax,'Seen Umax')
         call lptfile%add_column(lp%Vsmin,'Seen Vmin')
         call lptfile%add_column(lp%Vsmax,'Seen Vmax')
         call lptfile%add_column(lp%Wsmin,'Seen Wmin')
         call lptfile%add_column(lp%Wsmax,'Seen Wmax')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%write()
         ! Create timing monitor
         tfile=monitor(amroot=fs%cfg%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep number')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(wt_total%time,'Total [s]')
         call tfile%add_column(wt_lpt%time,'LPT [s]')
         call tfile%add_column(wt_lpt%percent,'LPT [%]')
         call tfile%add_column(wt_stat%time,'Stats [s]')
         call tfile%add_column(wt_stat%percent,'Stats [%]')
         call tfile%add_column(wt_rest%time,'Rest [s]')
         call tfile%add_column(wt_rest%percent,'Rest [%]')
         call tfile%write()
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
         wt_total%time_in=parallel_time() 
         
         call time%increment()
         
         !> ADVANCE PARTICLES
         wt_lpt%time_in=parallel_time()
         resU=fs%rho; resV=fs%visc-sgs%visc
         call lp%advance_scrw_tracer(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,sgs_visc=sgs%visc)
         wt_lpt%time=wt_lpt%time+parallel_time()-wt_lpt%time_in
         
         !> COMPUTE STATS
         wt_stat%time_in=parallel_time()
         call lp%get_max()
         if (mod(time%t,0.2_WP).lt.epsilon(1.0_WP)) call compute_pdfs(time%n)
         wt_stat%time=wt_stat%time+parallel_time()-wt_stat%time_in
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call lp%update_partmesh(pmesh)
            do ii = 1,lp%np_
               pmesh%var(1,ii) = lp%p(ii)%id
               pmesh%vec(:,1,ii) = lp%p(ii)%vel
               pmesh%vec(:,2,ii) = lp%cfg%get_velocity(pos=lp%p(ii)%pos,i0=lp%p(ii)%ind(1),j0=lp%p(ii)%ind(2),k0=lp%p(ii)%ind(3),U=fs%U,V=fs%V,W=fs%W)
               pmesh%vec(:,3,ii) = lp%p(ii)%us
            end do
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call lptfile%write()

         ! Monitor timing
         wt_total%time=parallel_time()-wt_total%time_in
         wt_lpt%percent=wt_lpt%time/wt_total%time*100.0_WP
         wt_stat%percent=wt_stat%time/wt_total%time*100.0_WP
         wt_rest%time=wt_total%time-wt_lpt%time-wt_stat%time
         wt_rest%percent=wt_rest%time/wt_total%time*100.0_WP
         call tfile%write()
         wt_total%time=0.0_WP; wt_total%percent=0.0_WP
         wt_lpt%time=0.0_WP;   wt_lpt%percent=0.0_WP
         wt_stat%time=0.0_WP;  wt_stat%percent=0.0_WP
         wt_rest%time=0.0_WP;  wt_rest%percent=0.0_WP

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
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi)
      
   end subroutine simulation_final

   !> Compute and output PDFs of fld vel and part vel
   subroutine compute_pdfs(step)
      use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM
      use parallel,only: MPI_REAL_WP
      character(5) :: stepstr
      character(len=8) :: fmt
      character(len=1024) :: filename
      real(WP), dimension(65) :: Uax,Vax
      real(WP), dimension(64) :: myUpdf,myVpdf,Updf,Vpdf
      real(WP) :: binw,Umean,Vmean,Uvar,Vvar,Vmin,Vmax,Umin,Umax
      integer, intent(in) :: step
      integer :: nbins,bin,i,j

      fmt = '(I5.5)'
      write(stepstr,fmt) step
      filename='./pdfs/UandV'//trim(stepstr)//'.dat'
      
      nbins = 64

      Updf=0.0_WP; Vpdf=0.0_WP

      Umin=-15.0_WP; Umax=15.0_WP
      Vmin=-15.0_WP; Vmax=15.0_WP

      binw=(Umax-Umin)/nbins

      do i=1,nbins+1
         Uax(i) = Umin + binw*i
         Vax(i) = Vmin + binw*i
      end do

      do i=1,lp%np_
         do j=1,3
            bin = CEILING(lp%p(i)%us(j)/binw, 1)
            if (bin.lt.0) then
               bin = (nbins-1.0_WP)/2.0_WP - ABS(bin)
            else
               bin = ABS(bin) + (nbins-1.0_WP)/2.0_WP
            end if
            Updf(bin) = Updf(bin) + 1.0_WP
            
            bin = CEILING(lp%p(i)%vel(j)/binw, 1)
            if (bin.lt.0) then
               bin = (nbins-1.0_WP)/2.0_WP - ABS(bin)
            else
               bin = ABS(bin) + (nbins-1.0_WP)/2.0_WP
            end if
            Vpdf(bin) = Vpdf(bin) + 1.0_WP
         end do
      end do

      Updf = Updf / sum(Updf) / binw
      Vpdf = Vpdf / sum(Vpdf) / binw

      open(1, file = filename, status='unknown')
      do i=1,nbins
         write(1,*) Updf(i), Vpdf(i)
      end do
      close(1)

      open(2, file = './pdfs/axes.dat', status='unknown')
      do i=1,nbins+1
         write(2,*) Uax(i), Vax(i)
      end do
      close(2)

   end subroutine compute_pdfs

end module simulation
