!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use sgsmodel_class,    only: sgsmodel
   use randomwalk_class,  only: lpt
   use datafile_class,    only: datafile
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single-phase incompressible flow solver and corresponding time tracker
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(lpt),         public :: lp
   type(sgsmodel),    public :: sgs
   type(datafile),    public :: df

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile,lptfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: S_
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:), allocatable :: uiuj
   real(WP), dimension(:,:,:,:), allocatable :: taudiv
   real(WP), dimension(:,:,:,:,:), allocatable :: gradu
   real(WP), dimension(:,:), allocatable :: avgUU
   real(WP), dimension(:), allocatable :: avgSS
   
   !> Statistics
   real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_,U2,U2_
   real(WP), dimension(:), allocatable :: c,c_,N_avg,N_inst,c_inst,c_inst_
   real(WP), dimension(:), allocatable :: pUavg,pVavg,pWavg,pU2,pV2,pW2
   real(WP), dimension(:), allocatable :: pUavg_,pVavg_,pWavg_,pU2_,pV2_,pW2_
   real(WP), dimension(:), allocatable :: counter, counter_
   
   !> Fluid viscosity
   real(WP) :: visc

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW

   !> Event for post-processing
   type(event) :: ppevt

   !> Misc
   logical :: restarted
   type(event) :: save_evt

contains

   !> Specialized subroutine that outputs particle concentration statistics
   subroutine postproc_particles()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      use mathtools, only: PI
      implicit none
      integer :: iunit,ierr,i,j,k
      character(len=str_medium) :: filename,timestamp

      c_inst=0.0_WP; c_inst_=0.0_WP
!      pUavg_=0.0_WP
!      pVavg_=0.0_WP
!      pWavg_=0.0_WP
!      pU2_=0.0_WP
!      pV2_=0.0_WP
!      pW2_=0.0_WP
!      counter_=0.0_WP

      ! Get expected density
      do i=1,fs%cfg%ny
         N_avg(i) = N_avg(i) + lp%np * fs%cfg%dy(i) / fs%cfg%yL * time%dt
         N_inst(i) = lp%np * fs%cfg%dy(i) / fs%cfg%yL * time%dt
      end do

      ! Accumulate particles
      do i=1,lp%np_
         k=lp%p(i)%ind(2)
         c_inst_(k)  = c_inst_(k) + time%dt
         c_(k)       = c_(k) + time%dt
         counter_(k) = counter_(k) + time%dt
         pUavg_(k)   = pUavg_(k) + time%dt * lp%p(i)%vel(1)
         pVavg_(k)   = pVavg_(k) + time%dt * lp%p(i)%vel(2)
         pWavg_(k)   = pWavg_(k) + time%dt * lp%p(i)%vel(3)
         pU2_(k)     = pU2_(k)   + time%dt * lp%p(i)%vel(1)**2
         pV2_(k)     = pV2_(k)   + time%dt * lp%p(i)%vel(2)**2
         pW2_(k)     = pW2_(k)   + time%dt * lp%p(i)%vel(3)**2
      end do

      ! All-reduce the data
      call MPI_ALLREDUCE(c_inst_,c_inst,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(c_,c,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pUavg_,pUavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pVavg_,pVavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pWavg_,pWavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pU2_,pU2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pV2_,pV2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(pW2_,pW2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(counter_,counter,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)

      do j=fs%cfg%jmin,fs%cfg%jmax
         if (counter(j).gt.0.0_WP) then
            pUavg(j)=pUavg(j)/counter(j)
            pVavg(j)=pVavg(j)/counter(j)
            pWavg(j)=pWavg(j)/counter(j)
            pU2(j)  =pU2(j)  /counter(j)
            pV2(j)  =pV2(j)  /counter(j)
            pW2(j)  =pW2(j)  /counter(j)
         else
            pUavg(j)=0.0_WP
            pVavg(j)=0.0_WP
            pWavg(j)=0.0_WP
            pU2(j)  =0.0_WP
            pV2(j)  =0.0_WP
            pW2(j)  =0.0_WP
         end if
         if (N_avg(j).gt.0.0_WP) then
            c(j) = c(j) / N_avg(j)
         else
            c(j) = 1.0_WP
         end if
         if (N_inst(j).gt.0.0_WP) then
            c_inst(j) = c_inst(j) / N_inst(j)
         else
            c_inst(j) = 1.0_WP
         end if
         
      end do
      
      ! If root, print it out
      if (fs%cfg%amRoot.and.ppevt%occurs()) then
         filename='./outs/Concentration_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'Height','c','c_inst','U','V','W','UU','VV','WW'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') & 
               & fs%cfg%ym(j),c(j),c_inst(j),pUavg(j),pVavg(j),pWavg(j),pU2(j)-pUavg(j)**2,pV2(j)-pVavg(j)**2,pW2(j)-pWavg(j)**2
         end do
         close(iunit)
      end if
   end subroutine postproc_particles

   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_vel()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      character(len=str_medium) :: filename,timestamp
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j) = vol_(j)+fs%cfg%vol(i,j,k)*time%dt
               Uavg_(j)=Uavg_(j)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)*time%dt
               U2_(j)=U2_(j)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)**2*time%dt
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE( vol_, vol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(U2_,U2,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      do j=fs%cfg%jmin,fs%cfg%jmax
         if (vol(j).gt.0.0_WP) then
            U2(j)=U2(j)/vol(j)
            Uavg(j)=Uavg(j)/vol(j)
         else
            Uavg(j)=0.0_WP
            U2(j)=0.0_WP
         end if
      end do
      ! If root, print it out
      if (fs%cfg%amRoot.and.ppevt%occurs()) then
         call execute_command_line('mkdir -p outs')
         filename='./outs/Uavg_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12)') 'Height','Uavg','U2'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5)') fs%cfg%ym(j),Uavg(j),U2(j)-Uavg(j)**2
         end do
         close(iunit)
      end if
     end subroutine postproc_vel


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
         allocate(SR(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(S_(  cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Allocate postproc arrays
         allocate(Uavg (cfg%jmin:cfg%jmax)); Uavg =0.0_WP
         allocate(Uavg_(cfg%jmin:cfg%jmax)); Uavg_=0.0_WP
         allocate(U2   (cfg%jmin:cfg%jmax)); U2   =0.0_WP
         allocate(U2_  (cfg%jmin:cfg%jmax)); U2_  =0.0_WP
         allocate(vol_ (cfg%jmin:cfg%jmax)); vol_ =0.0_WP
         allocate(vol  (cfg%jmin:cfg%jmax)); vol  =0.0_WP
         ! Allocate concentration arrays
         allocate(c      (1:cfg%ny)); c      =0.0_WP
         allocate(c_     (1:cfg%ny)); c_     =0.0_WP
         allocate(N_avg  (1:cfg%ny)); N_avg  =0.0_WP
         allocate(N_inst (1:cfg%ny)); N_inst =0.0_WP
         allocate(c_inst (1:cfg%ny)); c_inst =0.0_WP
         allocate(c_inst_(1:cfg%ny)); c_inst_=0.0_WP
         ! Allocate particle velocity stats
         allocate(pUavg (1:cfg%ny)); pUavg =0.0_WP
         allocate(pVavg (1:cfg%ny)); pVavg =0.0_WP
         allocate(pWavg (1:cfg%ny)); pWavg =0.0_WP
         allocate(pUavg_(1:cfg%ny)); pUavg_=0.0_WP
         allocate(pVavg_(1:cfg%ny)); pVavg_=0.0_WP
         allocate(pWavg_(1:cfg%ny)); pWavg_=0.0_WP
         allocate(pU2   (1:cfg%ny)); pU2   =0.0_WP
         allocate(pV2   (1:cfg%ny)); pV2   =0.0_WP
         allocate(pW2   (1:cfg%ny)); pW2   =0.0_WP
         allocate(pU2_  (1:cfg%ny)); pU2_  =0.0_WP
         allocate(pV2_  (1:cfg%ny)); pV2_  =0.0_WP
         allocate(pW2_  (1:cfg%ny)); pW2_  =0.0_WP
         allocate(counter (1:cfg%ny)); counter =0.0_WP
         allocate(counter_(1:cfg%ny)); counter_=0.0_WP
         ! Allocate arrays for Reynolds stress estimation
         allocate(uiuj(6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); uiuj=0.0_WP
         allocate(taudiv(3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); taudiv=0.0_WP
         allocate(avgUU(6,cfg%jmin:cfg%jmax)); avgUU=0.0_WP
         allocate(avgSS(  cfg%jmin:cfg%jmax)); avgSS=0.0_WP
         allocate(gradu(3,3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); gradu=0.0_WP
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

      ! Handle restart/saves here
      restart_and_save: block
        use string, only: str_medium
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
           df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=7)
           df%valname(1)='t'
           df%valname(2)='dt'
           df%varname(1)='U'
           df%varname(2)='V'
           df%varname(3)='W'
           df%varname(4)='P'
           df%varname(5)='LM'
           df%varname(6)='MM'
           df%varname(7)='visc'
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

      ! Create a single-phase flow solver without bconds
      create_and_initialize_flow_solver: block
        use incomp_class,   only: dirichlet, bcond
         use mathtools, only: twoPi
         use random, only: random_normal,random_uniform
         integer :: i,j,k,n
         real(WP) :: amp,vel
         real(WP), dimension(2) :: flucs
         type(bcond), pointer :: mybc
         ! Create flow solver
         fs=incomp(cfg=cfg,name='NS solver')
         ! Define boundary conditions
         call fs%add_bcond(name='bottom',type=dirichlet,locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top',type=dirichlet,locator=top_of_domain,face='y',dir=+1,canCorrect=.false. )
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',fs%rho)
        ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initialize velocity based on specified bulk
         call param_read('Ubulk',Ubulk)
         call param_read('Wbulk',Wbulk)
         where (fs%umask.eq.0) fs%U=Ubulk
         where (fs%wmask.eq.0) fs%W=Wbulk
         meanU=Ubulk
         meanW=Wbulk
         ! To facilitate transition
         call param_read('Perturbation',amp)
         vel=sqrt(Ubulk**2+Wbulk**2)
         if (restarted) then
            call df%pullvar(name='U',var=fs%U)
            call df%pullvar(name='V',var=fs%V)
            call df%pullvar(name='W',var=fs%W)
            call df%pullvar(name='P',var=fs%P)
         else
            do k=fs%cfg%kmino_,fs%cfg%kmaxo_
               do j=fs%cfg%jmino_,fs%cfg%jmaxo_
                  do i=fs%cfg%imino_,fs%cfg%imaxo_
                     !if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+amp*vel*cos(16.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)!*random_normal(m=0.0_WP,sd=0.5_WP)
                     !if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+amp*vel*cos(16.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)!*random_normal(m=0.0_WP,sd=0.5_WP)
                     if (fs%umask(i,j,k).eq.0)fs%U(i,j,k)=fs%U(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)*cos(8.0_WP*twoPi*fs%cfg%ym(j)/fs%cfg%yL)
                     if (fs%wmask(i,j,k).eq.0)fs%W(i,j,k)=fs%W(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  end do
               end do
            end do
         end if
         ! Set no-slip walls
         call fs%get_bcond('bottom',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=0.0_WP
         end do
         call fs%get_bcond('top',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%V(i,j,k)=0.0_WP
         end do
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      ! Initialize the LES model
      initialize_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         sgs%Cs_ref=0.1_WP
      end block initialize_sgs
      
      ! Revisit subgrid model to update arrays (dynamic) if this is a restart
      update_sgs: block
         if (restarted) then
            call df%pullvar(name='LM'   ,var=sgs%LM)
            call df%pullvar(name='MM'   ,var=sgs%MM)
            call df%pullvar(name='visc' ,var=sgs%visc)
         end if
      end block update_sgs

      ! Initialize LPT solver
      initialize_lpt: block
         use string, only: str_medium
         use random, only: random_uniform, random_normal
         real(WP) :: dp
         integer :: i,np
         character(len=str_medium) :: timestamp
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')
         ! Get drag model from the inpit
         ! call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
         ! Get particle density from the input
         call param_read('Particle density',lp%rho)
         ! Get particle diameter
         call param_read('Particle diameter',dp)
         ! Get number of particles
         call param_read('Number of particles',np)
         ! Randomly distribute particles w/o checking for overlap
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
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Init stochastic process at zero
               lp%p(i)%us=0.0_WP
               ! Init Wiener process to zero
               lp%p(i)%dW=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         lp%meanUn=0.0_WP
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
            pmesh%vec(:,2,i) = lp%cfg%get_velocity(pos=lp%p(i)%pos,i0=lp%p(i)%ind(1),j0=lp%p(i)%ind(2),k0=lp%p(i)%ind(3),U=fs%U,V=fs%V,W=fs%W)
            pmesh%vec(:,3,i) = lp%p(i)%us 
         end do
      end block create_pmesh

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('uu_diagonal',uiuj(1,:,:,:),uiuj(2,:,:,:),uiuj(3,:,:,:))
         call ens_out%add_vector('uu_offdiagonal',uiuj(4,:,:,:),uiuj(5,:,:,:),uiuj(6,:,:,:))
         call ens_out%add_vector('UUn_d',sgs%UUn(1,:,:,:),sgs%UUn(2,:,:,:),sgs%UUn(3,:,:,:)) 
         call ens_out%add_vector('UUn_o',sgs%UUn(4,:,:,:),sgs%UUn(5,:,:,:),sgs%UUn(6,:,:,:)) 
         call ens_out%add_scalar('dSS',sgs%dSS)
         call ens_out%add_scalar('viscosity',fs%visc)
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
         ! Create forcing monitor
         forcefile=monitor(fs%cfg%amRoot,'forcing')
         call forcefile%add_column(time%n,'Timestep number')
         call forcefile%add_column(time%t,'Time')
         call forcefile%add_column(meanU,'Bulk U')
         call forcefile%add_column(meanW,'Bulk W')
         call forcefile%write()
         ! Create LPT monitor
         call lp%get_max()
         lptfile=monitor(amroot=lp%cfg%amRoot,name='lpt')
         call lptfile%add_column(time%n,'Timestep number')
         call lptfile%add_column(time%t,'Time')
         call lptfile%add_column(lp%np,'Particle number')
         call lptfile%add_column(lp%ncol,'Collisions')
         call lptfile%add_column(lp%meanUn(1),'Deposition X vel')
         call lptfile%add_column(lp%meanUn(2),'Deposition Y vel')
         call lptfile%add_column(lp%meanUn(3),'Deposition Z vel')
         call lptfile%add_column(lp%Umin,'Particle Umin')
         call lptfile%add_column(lp%Umax,'Particle Umax')
         call lptfile%add_column(lp%Vmin,'Particle Vmin')
         call lptfile%add_column(lp%Vmax,'Particle Vmax')
         call lptfile%add_column(lp%Wmin,'Particle Wmin')
         call lptfile%add_column(lp%Wmax,'Particle Wmax')
         call lptfile%add_column(lp%dmin,'Particle dmin')
         call lptfile%add_column(lp%dmax,'Particle dmax')
         call lptfile%write()
      end block create_monitor

      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Perform the output
         call postproc_vel()
         call postproc_particles()
      end block create_postproc

   end subroutine simulation_init


   !> Time integrate our problem
   subroutine simulation_run
      use sgsmodel_class, only: dynamic_smag
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

         subgrid: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM,MPI_IN_PLACE
            use parallel, only: MPI_REAL_WP
            integer :: i,j,k,ierr
            avgUU = 0.0_WP; avgSS = 0.0_WP
            ! Get eddy viscosity
            call fs%interp_vel(Ui,Vi,Wi)
            call fs%get_strainrate(SR=SR)
            S_=sqrt(SR(1,:,:,:)**2+SR(2,:,:,:)**2+SR(3,:,:,:)**2+2.0_WP*(SR(4,:,:,:)**2+SR(5,:,:,:)**2+SR(6,:,:,:)**2))
            fs%visc=visc; resU=fs%rho; resV=fs%visc

            call sgs%get_visc(type=dynamic_smag,rho=resU,dt=time%dt,SR=SR,Ui=Ui,Vi=Vi,Wi=Wi)
            where (sgs%visc.lt.-fs%visc)
               sgs%visc=-fs%visc
            end where
            fs%visc=fs%visc + sgs%visc
            
            ! Average the Reynolds stress estimation in homogeneous directions
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     avgUU(:,j) = avgUU(:,j) + sgs%UUn(:,i,j,k) * fs%cfg%vol(i,j,k)
                     avgSS(j) = avgSS(j) + sgs%dSS(i,j,k) * fs%cfg%vol(i,j,k)
                  end do
               end do
            end do
            call MPI_ALLREDUCE(MPI_IN_PLACE,avgUU,6*fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,avgSS,  fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
            avgUU(1,:) = avgUU(1,:) / vol
            avgUU(2,:) = avgUU(2,:) / vol
            avgUU(3,:) = avgUU(3,:) / vol
            avgUU(4,:) = avgUU(4,:) / vol
            avgUU(5,:) = avgUU(5,:) / vol
            avgUU(6,:) = avgUU(6,:) / vol
            avgSS = avgSS / vol

            ! Compute the Reynolds stress
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     uiuj(:,i,j,k) = sgs%delta(i,j,k)**2 * S_(i,j,k)**2 * avgUU(:,j) / avgSS(j)
                  end do
               end do
            end do
         end block subgrid

         call fs%get_taurdiv(tau=uiuj,taudiv=taudiv)
         call fs%get_gradu(gradu=gradu)

         ! Advance particles
!!         call lp%advance_crw(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV,sgs_visc=sgs%visc)
         call lp%advance_crw_anisotropic(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=resU,visc=resV, &
                                 & sgs_visc=sgs%visc,gradu=gradu,taudiv=taudiv,uiuj=uiuj)

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

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

            ! Add body forcing
            forcing: block
               use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
               use parallel, only: MPI_REAL_WP
               integer :: i,j,k,ierr
               real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
               myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) then
                           myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*(2.0_WP*fs%U(i,j,k)-fs%Uold(i,j,k))
                           myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                        end if
                        if (fs%wmask(i,j,k).eq.0) then
                           myW   =myW   +fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)*(2.0_WP*fs%W(i,j,k)-fs%Wold(i,j,k))
                           myWvol=myWvol+fs%cfg%dx(i)*fs%cfg%dy(j)*fs%cfg%dzm(k)
                        end if
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               where (fs%umask.eq.0) resU=resU+fs%rho*(Ubulk-meanU)
               call MPI_ALLREDUCE(myWvol,Wvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myW   ,meanW,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanW=meanW/Wvol
               where (fs%wmask.eq.0) resW=resW+fs%rho*(Wbulk-meanW)
            end block forcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)

            ! Reset Dirichlet BCs
            dirichlet_velocity: block
              use incomp_class, only: bcond
              type(bcond), pointer :: mybc
              integer :: n,i,j,k
              call fs%get_bcond('bottom',mybc)
              do n=1,mybc%itr%no_
                 i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                 fs%V(i,j,k)=0.0_WP
              end do
              call fs%get_bcond('top',mybc)
              do n=1,mybc%itr%no_
                 i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                 fs%V(i,j,k)=0.0_WP
              end do
            end block dirichlet_velocity

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

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         output_ens: block
            integer :: ii
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
         end block output_ens

         ! Perform and output monitoring
         call fs%get_max()
         call lp%get_max()
         call mfile%write()
         call cflfile%write()
         call forcefile%write()
         call lptfile%write()
        
         ! Fluid phase post-processing 
         call postproc_vel()
         call postproc_particles()

         ! Finally, see if it's time to save restart files
         if (save_evt%occurs()) then
            save_restart: block
              use string, only: str_medium
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
              call df%pushvar(name='LM'   ,var=sgs%LM       )
              call df%pushvar(name='MM'   ,var=sgs%MM       )
              call df%pushvar(name='visc' ,var=sgs%visc     )
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
      deallocate(resU,resV,resW,Ui,Vi,Wi,SR)
      ! Deallocate work arrays
      deallocate(Uavg,Uavg_,U2,U2_,vol,vol_)

   end subroutine simulation_final





end module simulation
