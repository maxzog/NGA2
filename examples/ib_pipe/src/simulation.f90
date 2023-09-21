!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,D
!   use hypre_str_class,   only: hypre_str
   use fft3d_class,       only: fft3d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
!   use sgsmodel_class,    only: sgsmodel
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   ! use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use datafile_class,    only: datafile
   use string,            only: str_medium
   ! use lpt_class,         only: lpt
   implicit none
   private
   
   !> Get an an incompressible solver, pressure solver, and corresponding time tracker
   type(incomp),      public :: fs
   type(fft3d),       public :: ps
   !type(hypre_str),   public :: vs
   type(ddadi),       public :: vs
   ! type(lpt),         public :: lp
!   type(sgsmodel),    public :: sgs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt

   !> Provide a datafile and an event tracker for saving restarts
   type(event)    :: save_evt
   type(datafile) :: df
   ! type(partmesh) :: pmesh
   logical :: restarted
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile
   
   public :: simulation_init,simulation_run,simulation_final

   !> Timestep restrictions for LPT
   ! real(WP) :: lp_dt, lp_dt_max
   
   !> Work arrays
   !   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   ! real(WP), dimension(:,:,:), allocatable :: G
   real(WP) :: visc,mfr_target,mfr,bforce

   ! Stats
   integer :: nr
   real(WP), dimension (:), allocatable :: Uxr, Urr, Utr, Ux2, Ur2, Ut2,vol
   real(WP), dimension (:), allocatable :: Uxr_, Urr_, Utr_, Ux2_, Ur2_, Ut2_,vol_

   !> Event for calling post-processing script
   type(event) :: ppevt
   
contains
   
   !> Compute massflow rate
   function get_bodyforce_mfr(srcU) result(mfr)
      use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      real(WP), dimension(fs%cfg%imino_:,fs%cfg%jmino_:,fs%cfg%kmino_:), intent(in), optional :: srcU !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,ierr
      real(WP) :: vol,myRhoU,myUvol,Uvol,mfr
      myRhoU=0.0_WP; myUvol=0.0_WP
      if (present(srcU)) then
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*(fs%U(i,j,k)*fs%rho+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*fs%U(i,j,k)*fs%rho
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); mfr=mfr/Uvol
   end function get_bodyforce_mfr

   !> Subroutine for getting pipe velocity stats
   subroutine postproc_vel()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE, MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit, ierr, i, j, k, rind
      real(WP) :: rnorm, rmax, dr
      real(WP), dimension(3) :: r, vel, rhat, rtheta
      character(len=str_medium) :: filename, timestamp

      rmax = 0.5_WP * D
      dr = rmax / real(nr,WP)

      ! GET MEAN PROFILES
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               r = [0.0_WP, fs%cfg%ym(j), fs%cfg%zm(k)]    ! get radial position on axial slice
               rnorm = norm2(r)                          
               if (rnorm.gt.rmax) cycle                  ! cycle if r > d/2
               rind = floor(rnorm/dr) + 1  
               rhat = r / rnorm
               rtheta = cross(rhat, [0.0_WP, 1.0_WP, 0.0_WP])  ! assemble theta-wise unit vector
               rtheta = cross(rtheta, rhat) / norm2(cross(rtheta, rhat))
               vel = [Ui(i,j,k), Vi(i,j,k), Wi(i,j,k)]
               vol_(rind)=vol_(rind)+fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Uxr_(rind) = Uxr_(rind) + Ui(i,j,k)*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Urr_(rind) = Urr_(rind) + dot_product(vel, rhat)*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Utr_(rind) = Utr_(rind) + dot_product(vel, rtheta)*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Ux2_(rind) = Ux2_(rind) + Ui(i,j,k)**2*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Ur2_(rind) = Ur2_(rind) + dot_product(vel, rhat)**2*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
               Ut2_(rind) = Ut2_(rind) + dot_product(vel, rtheta)**2*fs%cfg%vol(i,j,k)*fs%cfg%VF(i,j,k)*time%dt
            end do
         end do
      end do

      call MPI_ALLREDUCE(vol_, vol, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Uxr_, Uxr, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Urr_, Urr, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Utr_, Utr, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Ux2_, Ux2, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Ur2_, Ur2, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)
      call MPI_ALLREDUCE(Ut2_, Ut2, nr, MPI_REAL_WP, MPI_SUM, fs%cfg%comm, ierr)

      do i=1,nr
         if (vol(i).gt.0.0_WP) then
            Uxr(i)=Uxr(i)/vol(i)
            Urr(i)=Urr(i)/vol(i)
            Utr(i)=Utr(i)/vol(i)
            Ux2(i)=Ux2(i)/vol(i)
            Ur2(i)=Ur2(i)/vol(i)
            Ut2(i)=Ut2(i)/vol(i)
         end if
      end do
       ! If root, print it out
      if (fs%cfg%amRoot.and.ppevt%occurs()) then
         filename='rstat_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12,3x,a12,3x,a12,3x,a12)') 'r','Ux','Ux_rms','Ut_rms','Ur_rms'
         do i=1,nr
            write(iunit,'(es12.5,3x,es12.5,3x,es12.5,3x,es12.5,3x,es12.5)') dr*real(i-1,WP),Uxr(i),Uxr(i)**2-Ux2(i),Utr(i)**2-Ut2(i),Urr(i)**2-Ur2(i)
         end do
         close(iunit)
      end if
      
      contains

         function cross(a, b) result(c)
            real(WP), dimension(3) :: c
            real(WP), dimension(3), intent(IN) :: a, b
            c(1) = a(2) * b(3) - a(3) * b(2)
            c(2) = a(3) * b(1) - a(1) * b(3)
            c(3) = a(1) * b(2) - a(2) * b(1)
         end function cross

   end subroutine postproc_vel   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max time',time%tmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      
      ! Create an incompressible flow solver without bconds
      create_flow_solver: block
      !   use hypre_str_class, only: pcg_pfmg
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Configure pressure solver
         ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         !ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         !ps%maxlevel=14
         !call param_read('Pressure iteration',ps%maxit)
         !call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         !vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg,nst=7)
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
!         allocate(SR(1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))   
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! allocate(G   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

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
      
      ! Initialize our velocity field
      initialize_velocity: block
         use mathtools, only: twoPi
         use random,    only: random_uniform
         integer :: i,j,k
         real(WP) :: Ubulk,amp
         ! Initial fields
         call param_read('Bulk velocity',Ubulk)
         fs%U=Ubulk; fs%V=0.0_WP; fs%W=0.0_WP; fs%P=0.0_WP
         ! For faster transition
         call param_read('Fluctuation amp',amp,default=0.0_WP)
         if (restarted) then
            call df%pullvar(name='U',var=fs%U)
            call df%pullvar(name='V',var=fs%V)
            call df%pullvar(name='W',var=fs%W)
            call df%pullvar(name='P',var=fs%P)
         else
            do k=fs%cfg%kmin_,fs%cfg%kmax_
               do j=fs%cfg%jmin_,fs%cfg%jmax_
                  do i=fs%cfg%imin_,fs%cfg%imax_
                     fs%U(i,j,k)=fs%U(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)*cos(8.0_WP*twoPi*fs%cfg%ym(j)/fs%cfg%yL)
                     fs%V(i,j,k)=fs%V(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                     fs%W(i,j,k)=fs%W(i,j,k)+Ubulk*random_uniform(lo=-0.5_WP*amp,hi=0.5_WP*amp)+amp*Ubulk*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  end do
               end do
            end do
         end if
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Compute cell-centered velocity
         call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         call fs%get_div()
         ! Get target MFR and zero bodyforce
         mfr=get_bodyforce_mfr()
         mfr_target=mfr
         bforce=0.0_WP
      end block initialize_velocity
      
      
      ! ! Initialize IBM fields
      ! initialize_ibm: block
      !    integer :: i,j,k
      !    do k=fs%cfg%kmino_,fs%cfg%kmaxo_
      !       do j=fs%cfg%jmino_,fs%cfg%jmaxo_
      !          do i=fs%cfg%imino_,fs%cfg%imaxo_
      !             G(i,j,k)=0.5_WP*D-sqrt(fs%cfg%ym(j)**2+fs%cfg%zm(k)**2)
      !          end do
      !       end do
      !    end do
      ! end block initialize_ibm

      ! ! Initialize LPT solver
      ! initialize_lpt: block
      !    use random, only: random_uniform, random_normal
      !    real(WP) :: dp
      !    integer :: i,np
      !    character(len=str_medium) :: timestamp
      !    ! Create solver
      !    lp=lpt(cfg=cfg,name='LPT')
      !    ! Get drag model from the inpit
      !    call param_read('Drag model',lp%drag_model,default='Schiller-Naumann')
      !    ! Get particle density from the input
      !    call param_read('Particle density',lp%rho)
      !    ! Get particle diameter from the input
      !    call param_read('Particle Stokes number', stk)
      !    ! Get number of particles
      !    call param_read('Number of particles',np)
      !    ! Check if a stochastic SGS model is used
      !    if (restarted) then
      !       call param_read('Restart from',timestamp,'r')
      !       ! Read the part file
      !       call lp%read(filename='restart/part_'//trim(adjustl(timestamp)))
      !    else
      !    ! Root process initializes np particles randomly
      !       if (lp%cfg%amRoot) then
      !          tau_eta = sqrt(visc/EPS0)
      !          dp = sqrt(18.0_WP*visc*stk*tau_eta/lp%rho)
      !          call lp%resize(np)
      !          do i=1,np
      !             ! Give id
      !             lp%p(i)%id=int(i,8)
      !             ! Set the diameter
      !             lp%p(i)%d=dp
      !             ! Assign random position in the domain
      !             lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
      !             &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
      !             &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
      !             ! Give zero velocity
      !             lp%p(i)%vel=0.0_WP
      !             ! Give zero dt
      !             lp%p(i)%dt=0.0_WP
      !             ! temp
      !             lp%p(i)%uf=0.0_WP
      !             ! Locate the particle on the mesh
      !             lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
      !             ! Activate the particle
      !             lp%p(i)%flag=0
      !          end do
      !       end if
      !       ! Distribute particles
      !       call lp%sync()
      !    end if
      ! end block initialize_lpt
      
      
!      ! Create an LES model
!      create_sgs: block
!         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
!      end block create_sgs


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='pipe')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         ! call ens_out%add_scalar('levelset',G)
         call ens_out%add_scalar('pressure',fs%P)
         ! call ens_out%add_scalar('visc_sgs',sgs%visc)
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
         call mfile%add_column(mfr,'MFR')
         call mfile%add_column(bforce,'Body force')
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

      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         nr=fs%cfg%ny/2
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Allocate
         allocate(vol(nr)); vol=0.0_WP
         allocate(Uxr(nr)); Uxr=0.0_WP
         allocate(Urr(nr)); Urr=0.0_WP
         allocate(Utr(nr)); Utr=0.0_WP
         allocate(Ux2(nr)); Ux2=0.0_WP
         allocate(Ur2(nr)); Ur2=0.0_WP
         allocate(Ut2(nr)); Ut2=0.0_WP
         allocate(vol_(nr)); vol_=0.0_WP
         allocate(Uxr_(nr)); Uxr_=0.0_WP
         allocate(Urr_(nr)); Urr_=0.0_WP
         allocate(Utr_(nr)); Utr_=0.0_WP
         allocate(Ux2_(nr)); Ux2_=0.0_WP
         allocate(Ur2_(nr)); Ur2_=0.0_WP
         allocate(Ut2_(nr)); Ut2_=0.0_WP
      end block create_postproc
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
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
         
!         ! Turbulence modeling
!         sgs_modeling: block
!               use sgsmodel_class, only: vreman
!               resU=fs%rho
!               call fs%get_gradu(gradU)
!               call fs%interp_vel(Ui,Vi,Wi)
!               call fs%get_strainrate(SR)
!               call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
!               fs%visc=visc+sgs%visc
!         end block sgs_modeling
         
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
            bodyforcing: block
               real(WP) :: mfr
               mfr=get_bodyforce_mfr(resU)
               bforce=(mfr_target-mfr)/time%dtmid
               resU=resU+time%dt*bforce
            end block bodyforcing
            
            ! Apply IB forcing to enforce BC at the pipe walls
            !ibforcing: block
            !   integer :: i,j,k
            !   do k=fs%cfg%kmino_,fs%cfg%kmaxo_
            !      do j=fs%cfg%jmino_,fs%cfg%jmaxo_
            !         do i=fs%cfg%imino_,fs%cfg%imaxo_
            !            resU(i,j,k)=resU(i,j,k)-(1.0_WP-get_VF(i,j,k,'U'))*fs%rho*fs%U(i,j,k)
            !            resV(i,j,k)=resV(i,j,k)-(1.0_WP-get_VF(i,j,k,'V'))*fs%rho*fs%V(i,j,k)
            !            resW(i,j,k)=resW(i,j,k)-(1.0_WP-get_VF(i,j,k,'W'))*fs%rho*fs%W(i,j,k)
            !         end do
            !      end do
            !   end do
            !end block ibforcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply IB forcing to enforce BC at the pipe walls
            ibforcing: block
              integer :: i,j,k
              do k=fs%cfg%kmin_,fs%cfg%kmax_
                 do j=fs%cfg%jmin_,fs%cfg%jmax_
                    do i=fs%cfg%imin_,fs%cfg%imax_
                       fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                       fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                       fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                    end do
                 end do
              end do
              call fs%cfg%sync(fs%U)
              call fs%cfg%sync(fs%V)
              call fs%cfg%sync(fs%W)
            end block ibforcing

           
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
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()

         ! Comput stats
         call postproc_vel()
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         
         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()

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
              ! call lp%write(filename='restart/part_'//trim(adjustl(timestamp)))
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
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,gradU)
      deallocate(Uxr, Urr, Utr, Ux2, Ur2, Ut2,vol,Uxr_, Urr_, Utr_, Ux2_, Ur2_, Ut2_,vol_)
      
   end subroutine simulation_final
   
end module simulation
