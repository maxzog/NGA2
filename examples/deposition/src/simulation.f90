!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use lowmach_class,     only: lowmach
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use datafile_class,    only: datafile
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single-phase lowmach flow solver and corresponding time tracker
   type(lowmach),     public :: fs
   type(lpt),         public :: lp
   type(timetracker), public :: time
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs

   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(partmesh) :: pmesh
   type(event)    :: ens_evt

   !> Restart data
   type(datafile) :: df
   logical :: restarted
   type(event) :: save_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,forcefile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: rho0,dRHOdt
   real(WP), dimension(:,:,:), allocatable :: srcUlp,srcVlp,srcWlp
   real(WP), dimension(:,:,:), allocatable :: tmp1,tmp2,tmp3
   real(WP), dimension(:,:,:,:), allocatable :: SR

   !> Timesteps for LPT
   real(WP) :: lp_dt,lp_dt_max

   !> Fluid viscosity
   real(WP) :: visc,rho,mfr,mfr_target,bforce

   !> Channel forcing
   real(WP) :: Ubulk,Wbulk
   real(WP) :: meanU,meanW

   !> Event for post-processing
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
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*(fs%rhoU(i,j,k)+srcU(i,j,k))
               end do
            end do
         end do
      else
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  vol=fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                  myUvol=myUvol+vol
                  myRhoU=myRhoU+vol*fs%rhoU(i,j,k)
               end do
            end do
         end do
      end if
      call MPI_ALLREDUCE(myUvol,Uvol,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(myRhoU,mfr ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); mfr=mfr/Uvol
   end function get_bodyforce_mfr


   !> Specialized subroutine that outputs the velocity distribution
   subroutine postproc_vel()
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      implicit none
      integer :: iunit,ierr,i,j,k
      real(WP), dimension(:), allocatable :: Uavg,Uavg_,vol,vol_
      character(len=str_medium) :: filename,timestamp
      ! Allocate vertical line storage
      allocate(Uavg (fs%cfg%jmin:fs%cfg%jmax)); Uavg =0.0_WP
      allocate(Uavg_(fs%cfg%jmin:fs%cfg%jmax)); Uavg_=0.0_WP
      allocate(vol_ (fs%cfg%jmin:fs%cfg%jmax)); vol_ =0.0_WP
      allocate(vol  (fs%cfg%jmin:fs%cfg%jmax)); vol  =0.0_WP
      ! Integrate all data over x and z
      do k=fs%cfg%kmin_,fs%cfg%kmax_
         do j=fs%cfg%jmin_,fs%cfg%jmax_
            do i=fs%cfg%imin_,fs%cfg%imax_
               vol_(j) = vol_(j)+fs%cfg%vol(i,j,k)
               Uavg_(j)=Uavg_(j)+fs%cfg%vol(i,j,k)*fs%U(i,j,k)
            end do
         end do
      end do
      ! All-reduce the data
      call MPI_ALLREDUCE( vol_, vol,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      call MPI_ALLREDUCE(Uavg_,Uavg,fs%cfg%ny,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
      do j=fs%cfg%jmin,fs%cfg%jmax
         if (vol(j).gt.0.0_WP) then
            Uavg(j)=Uavg(j)/vol(j)
         else
            Uavg(j)=0.0_WP
         end if
      end do
      ! If root, print it out
      if (fs%cfg%amRoot) then
         filename='Uavg_'
         write(timestamp,'(es12.5)') time%t
         open(newunit=iunit,file=trim(adjustl(filename))//trim(adjustl(timestamp)),form='formatted',status='replace',access='stream',iostat=ierr)
         write(iunit,'(a12,3x,a12)') 'Height','Uavg'
         do j=fs%cfg%jmin,fs%cfg%jmax
            write(iunit,'(es12.5,3x,es12.5)') fs%cfg%ym(j),Uavg(j)
         end do
         close(iunit)
      end if
      ! Deallocate work arrays
      deallocate(Uavg,Uavg_,vol,vol_)
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

      restart_and_save: block
         use string, only: str_medium
         character(len=str_medium) :: timestamp
         save_evt=event(time,'Restart output')
         call param_read('Restart output period', save_evt%tper)
         ! Check for restart
         call param_read(tag='Restart from',val=timestamp,short='r',default='')
         restarted=.false.; if(len_trim(timestamp).gt.0) restarted=.true.
         if (restarted) then
            call param_read('Restart from',timestamp,'r')
            df=datafile(pg=cfg,fdata='restart/data_'//trim(adjustl(timestamp)))
         else
            df=datafile(pg=cfg,filename=trim(cfg%name),nval=2,nvar=4)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%varname(1)='U'
            df%varname(2)='V'
            df%varname(3)='W'
            df%varname(4)='P'
         end if
      end block restart_and_save

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(dRHOdt (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho0   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcUlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcVlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcWlp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp1   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp2   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(tmp3   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR   (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Create a single-phase flow solver without bconds
      create_flow_solver: block
         use lowmach_class,   only: dirichlet, bcond
         use mathtools, only: twoPi
         integer :: i,j,k,n
         real(WP) :: amp,vel
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='NS solver')
         ! Define boundary conditions
         call fs%add_bcond(name='bottom',type=dirichlet,locator=bottom_of_domain,face='y',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='top',type=dirichlet,locator=top_of_domain,face='y',dir=+1,canCorrect=.true. )
         ! Assign constant viscosity
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Assign constant density
         call param_read('Density',rho); fs%rho=rho
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver

      ! Initialize our LPT solver
      initialize_lpt: block
        use random, only: random_lognormal,random_uniform
        use mathtools, only: Pi,twoPi
        use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE,MPI_INTEGER
        use parallel, only: MPI_REAL_WP
        real(WP) :: VFavg,Vol_,sumVolp,dp,up,amp
        integer :: i,j,k,ii,jj,kk,nn,ip,jp,kp,np,offset,ierr
        integer, dimension(:,:,:), allocatable :: npic      !< Number of particle in cell
        integer, dimension(:,:,:,:), allocatable :: ipic    !< Index of particle in cell
        logical :: overlap
        ! Create solver
        lp=lpt(cfg=cfg,name='LPT')
        ! Get mean volume fraction from input
        call param_read('Particle volume fraction',VFavg)
        ! Get drag model from input
        ! Get particle density from input
        call param_read('Particle density',lp%rho)
        ! Get particle diameter from input
        call param_read('Particle diameter',dp)
        ! Get particle velocity from input
        call param_read('Particle u velocity',up,default=0.0_WP)
        call param_read('Particle fluctuations',amp,default=0.0_WP)
        ! Maximum timestep size used for particles
        call param_read('Particle timestep size',lp_dt_max,default=huge(1.0_WP))
        lp_dt=lp_dt_max
        ! Set collision timescale
        call param_read('Collision timescale',lp%tau_col,default=15.0_WP*time%dt)
        ! Set coefficient of restitution
        call param_read('Coefficient of restitution',lp%e_n)
        call param_read('Friction coefficient',lp%mu_f,default=0.0_WP)
        ! Set gravity
        call param_read('Gravity',lp%gravity)
        ! Set filter scale to 3.5*dx
        lp%filter_width=3.5_WP*cfg%min_meshsize
        ! Initialize particles
        ! Get volume of domain belonging to this proc
        Vol_=0.0_WP
        do k=lp%cfg%kmin_,lp%cfg%kmax_
           do j=lp%cfg%jmin_,lp%cfg%jmax_
              do i=lp%cfg%imin_,fs%cfg%imax_
                 Vol_=Vol_+lp%cfg%dx(i)*lp%cfg%dy(j)*lp%cfg%dz(k)
              end do
           end do
        end do
        ! Get particle diameters
        np=ceiling(VFavg*Vol_/(pi*dp**3/6.0_WP))
        call lp%resize(np)
        ! Allocate particle in cell arrays
        allocate(npic(     lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); npic=0
        allocate(ipic(1:40,lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_)); ipic=0
        ! Distribute particles
        sumVolp=0.0_WP
        do i=1,np
           ! Set the diameter
           lp%p(i)%d=dp
           ! Give position (avoid overlap)
           overlap=.true.
           do while (overlap)
              lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin_),lp%cfg%x(lp%cfg%imax_+1)-dp),&
                   &       random_uniform(lp%cfg%y(lp%cfg%jmin_),lp%cfg%y(lp%cfg%jmax_+1)-dp),&
                   &       random_uniform(lp%cfg%z(lp%cfg%kmin_),lp%cfg%z(lp%cfg%kmax_+1)-dp)]
              lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
              overlap=.false.
              do kk=lp%p(i)%ind(3)-1,lp%p(i)%ind(3)+1
                 do jj=lp%p(i)%ind(2)-1,lp%p(i)%ind(2)+1
                    do ii=lp%p(i)%ind(1)-1,lp%p(i)%ind(1)+1
                       do nn=1,npic(ii,jj,kk)
                          j=ipic(nn,ii,jj,kk)
                          if (sqrt(sum((lp%p(i)%pos-lp%p(j)%pos)**2)).lt.0.5_WP*(lp%p(i)%d+lp%p(j)%d)) overlap=.true.
                       end do
                    end do
                 end do
              end do
           end do
             
           ! Activate the particle
           lp%p(i)%flag=0
           ip=lp%p(i)%ind(1); jp=lp%p(i)%ind(2); kp=lp%p(i)%ind(3)
           npic(ip,jp,kp)=npic(ip,jp,kp)+1
           ipic(npic(ip,jp,kp),ip,jp,kp)=i
           ! Give zero velocity
           lp%p(i)%vel=0.0_WP
           lp%p(i)%vel(1)=up
           if (amp.gt.0.0_WP) lp%p(i)%vel=lp%p(i)%vel+[random_uniform(-amp,amp),random_uniform(-amp,amp),random_uniform(-amp,amp)]

           ! Give zero collision force
           lp%p(i)%Acol=0.0_WP
           lp%p(i)%Tcol=0.0_WP
           ! Give zero dt
           lp%p(i)%dt=0.0_WP
           ! Sum up volume
           sumVolp=sumVolp+Pi/6.0_WP*lp%p(i)%d**3
        end do
        deallocate(npic,ipic)
        call lp%sync()
        ! Set ID
        offset=0
        do i=1,lp%cfg%rank
           offset=offset+lp%np_proc(i)
        end do
        do i=1,lp%np_
           lp%p(i)%id=int(i+offset,8)
        end do
        ! Get mean diameter and volume fraction
        call MPI_ALLREDUCE(sumVolp,VFavg,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); VFavg=VFavg/(lp%cfg%xL*lp%cfg%yL*lp%cfg%zL)
        if (lp%cfg%amRoot) then
           print*,"===== Particle Setup Description ====="
           print*,'Number of particles', lp%np
           print*,'Mean volume fraction',VFavg
        end if
        ! Get initial particle volume fraction
        call lp%update_VF()
      end block initialize_lpt

      initialize_flow: block
         use lowmach_class,   only: dirichlet, bcond
         use mathtools, only: Pi,twoPi
         type(bcond), pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: vel,amp
         ! Initialize velocity based on specified bulk
         call param_read('Ubulk',Ubulk)
         call param_read('Wbulk',Wbulk)
         ! Store fluid density
         fs%rho=rho; rho0=rho
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
                     if (fs%umask(i,j,k).eq.0) fs%U(i,j,k)=fs%U(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%zm(k)/fs%cfg%zL)
                     if (fs%wmask(i,j,k).eq.0) fs%W(i,j,k)=fs%W(i,j,k)+amp*vel*cos(8.0_WP*twoPi*fs%cfg%xm(i)/fs%cfg%xL)
                  end do
               end do
            end do
         end if
         ! Set no-slip walls
         call fs%get_bcond('bottom',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
         end do
         call fs%get_bcond('top',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
         end do
         ! Set density based on particle VF
         fs%rho=rho*(1.0_WP-lp%VF) 
         ! Get rho*U
         call fs%rho_multiply
         mfr=get_bodyforce_mfr()
         bforce=0.0_WP
         fs%rhoUold=fs%rhoU
         fs%rhoVold=fs%rhoV
         fs%rhoWold=fs%rhoW
         dRHOdt=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=dRHOdt)
         call fs%get_mfr()
         mfr_target=mfr
      end block initialize_flow

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='channel')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('viscosity',fs%visc)
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
      end block create_monitor


      ! Create a specialized post-processing file
      create_postproc: block
         ! Create event for data postprocessing
         ppevt=event(time=time,name='Postproc output')
         call param_read('Postproc output period',ppevt%tper)
         ! Perform the output
         if (ppevt%occurs()) call postproc_vel()
      end block create_postproc


   end subroutine simulation_init


   !> Time integrate our problem
   subroutine simulation_run
      implicit none
      real(WP) :: cfl
      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call lp%get_cfl(time%dt,cflc=time%cfl)
         call fs%get_cfl(time%dt,time%cfl); time%cfl=max(time%cfl,cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old velocity
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU 
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         if (lp%np.gt.0) then
            lpt: block
               real(WP) :: dt_done, mydt
               ! Get fluid stress
               call fs%get_div_stress(resU,resV,resW)
               resU=resU+bforce
               ! Zero-out LPT source terms
               srcUlp=0.0_WP; srcVlp=0.0_WP; srcWlp=0.0_WP
               ! Sub-iterate
               call lp%get_cfl(lp_dt,cflc=cfl,cfl=cfl)
               if (cfl.gt.0.0_WP) lp_dt=min(lp_dt*time%cflmax/cfl,lp_dt_max)
               dt_done=0.0_WP
               do while (dt_done.lt.time%dtmid)
                  ! Get timestep
                  mydt=min(lp_dt,time%dtmid-dt_done)
!                  call lp%collide(dt=mydt)
                  call lp%advance(dt=mydt,U=fs%U,V=fs%V,W=fs%W,rho=rho0,visc=fs%visc,stress_x=resU,stress_y=resV,stress_z=resW,&
                  &    srcU=tmp1,srcV=tmp2,srcW=tmp3)
                  srcUlp=srcUlp+tmp1
                  srcVlp=srcVlp+tmp2
                  srcWlp=srcWlp+tmp3
                  ! Increment
                  dt_done=dt_done+mydt
               end do
               fs%rho=rho*(1.0_WP-lp%VF)
               dRHOdt=(fs%RHO-fs%RHOold)/time%dtmid
               ! Compute PTKE and store source terms
               call lp%get_ptke(dt=time%dtmid,Ui=Ui,Vi=Vi,Wi=Wi,visc=fs%visc,rho=rho0,srcU=tmp1,srcV=tmp2,srcW=tmp3)
               srcUlp=srcUlp+tmp1
               srcVlp=srcVlp+tmp2
               srcWlp=srcWlp+tmp3
            end block lpt
         end if
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Build mid-time momentum
            fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rhoU-fs%rhoUold)+time%dtmid*resU
            resV=-2.0_WP*(fs%rhoV-fs%rhoVold)+time%dtmid*resV
            resW=-2.0_WP*(fs%rhoW-fs%rhoWold)+time%dtmid*resW

            ! Add momentum source from LPT
            add_lpt_source: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*srcUlp(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*srcUlp(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*srcUlp(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
            end block add_lpt_source

!            ! Add body forcing
!            forcing: block
!               mfr=get_bodyforce_mfr(resU)
!               bforce=(mfr_target-mfr)/time%dtmid
!               resU=resU+time%dtmid*bforce
!            end block forcing

!            ! Add body forcing
            forcing: block
               use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
               use parallel, only: MPI_REAL_WP
               integer :: i,j,k,ierr
               real(WP) :: myU,myUvol,Uvol
               myU=0.0_WP; myUvol=0.0_WP
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        if (fs%umask(i,j,k).eq.0) then
                           myU   =myU   +fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)*(2.0_WP*fs%rhoU(i,j,k)-fs%rhoUold(i,j,k))
                           myUvol=myUvol+fs%cfg%dxm(i)*fs%cfg%dy(j)*fs%cfg%dz(k)
                        end if
                     end do
                  end do
               end do
               call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr)
               call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); meanU=meanU/Uvol
               bforce=rho*Ubulk-meanU
               where (fs%umask.eq.0) resU=resU+bforce
            end block forcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions on the resulting fields
            call fs%apply_bcond(time%t,time%dt)
            call fs%rho_multiply()

            ! Reset Dirichlet BCs
            dirichlet_velocity: block
              use lowmach_class, only: bcond
              type(bcond), pointer :: mybc
              integer :: n,i,j,k
              call fs%get_bcond('bottom',mybc)
              do n=1,mybc%itr%no_
                 i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                 fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
              end do
              call fs%get_bcond('top',mybc)
              do n=1,mybc%itr%no_
                 i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                 fs%U(i,j,k)=0.0_WP; fs%V(i,j,k)=0.0_WP; fs%W(i,j,k)=0.0_WP
              end do
            end block dirichlet_velocity

            ! Solve Poisson equation
            call fs%correct_mfr(drhodt=dRHOdt)
            call fs%get_div(drhodt=dRHOdt)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid !TODO: Does this need *rho?
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide

            ! Increment sub-iteration counter
            time%it=time%it+1
         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=dRHOdt)
      
         mfr=get_bodyforce_mfr()

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i=1,lp%np_
                  pmesh%var(1,i)=lp%p(i)%d
                  pmesh%vec(:,1,i)=lp%p(i)%vel
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if 

         ! Perform and output monitoring
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call forcefile%write()

         ! Specialized post-processing
         if (ppevt%occurs()) call postproc_vel()

         ! Output restart files
         if (save_evt%occurs()) then
            save_restart: block
               use string, only: str_medium
               character(len=str_medium) :: timestamp
               ! File prefix
               write(timestamp,'(es12.5)') time%t
               ! Make save directory
               if (fs%cfg%amRoot) call execute_command_line('mkdir -p restart')
               ! Populate df and write
               call df%pushval(name='t' ,val=time%t)
               call df%pushval(name='dt',val=time%dt)
               call df%pushvar(name='U' ,var=fs%U)
               call df%pushvar(name='V' ,var=fs%V)
               call df%pushvar(name='W' ,var=fs%W)
               call df%pushvar(name='P' ,var=fs%P)
               call df%write(fdata='restart/data_'//trim(adjustl(timestamp)))
               ! Write particles
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

   end subroutine simulation_final





end module simulation
