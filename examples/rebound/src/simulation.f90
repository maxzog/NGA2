!> Various definitions and tools for running an NGA2 simulation
module simulation
  use precision,         only: WP
  use geometry,          only: cfg
  use lpt_class,         only: lpt
  use timetracker_class, only: timetracker
  use ensight_class,     only: ensight
  use partmesh_class,    only: partmesh
  use event_class,       only: event
  use monitor_class,     only: monitor
  use messager,          only: die
  implicit none
  private

  !> Only get a LPT solver and corresponding time tracker
  type(lpt),         public :: lp
  type(timetracker), public :: time

  !> Ensight postprocessing
  type(partmesh) :: pmesh
  type(ensight)  :: ens_out
  type(event)    :: ens_evt

  !> Simulation monitor file
  type(monitor) :: mfile,cflfile
  real(WP), dimension(3) :: pos,vel,angvel

  public :: simulation_init,simulation_run,simulation_final


  integer :: save_count

contains


   !> Get partoicle properties to all procs
   subroutine communicate_part
   use mpi_f08
   use parallel, only: MPI_REAL_WP
   implicit none
   integer :: i,ierr
   real(WP), dimension(3) :: buf
   pos=-huge(1.0_WP)
   vel=-huge(1.0_WP)
   angvel=-huge(1.0_WP)
   do i=1,lp%np_
      pos=lp%p(i)%pos
      vel=lp%p(i)%vel
      angvel=lp%p(i)%angvel
   end do
   call MPI_ALLREDUCE(pos,buf,3,MPI_REAL_WP,MPI_MAX,lp%cfg%comm,ierr); pos=buf
   call MPI_ALLREDUCE(vel,buf,3,MPI_REAL_WP,MPI_MAX,lp%cfg%comm,ierr); vel=buf
   call MPI_ALLREDUCE(angvel,buf,3,MPI_REAL_WP,MPI_MAX,lp%cfg%comm,ierr); angvel=buf
 end subroutine communicate_part


  !> Initialization of problem solver
  subroutine simulation_init
    use param, only: param_read
    implicit none


    ! Initialize time tracker with 1 subiterations
    initialize_timetracker: block
      time=timetracker(amRoot=cfg%amRoot)
      call param_read('Max timestep size',time%dtmax)
      call param_read('Max cfl number',time%cflmax)
      time%dt=time%dtmax
      time%itmax=1
    end block initialize_timetracker


    ! Initialize our LPT
    initialize_lpt: block
      use random, only: random_uniform
      use mathtools, only: PI
      real(WP) :: d,dx,dtheta,v0,y0
      integer :: i,np
      ! Create solver
      lp=lpt(cfg=cfg,name='LPT')
      ! Get particle density from the input
      call param_read('Particle density',lp%rho)
      ! Set gravity
      call param_read('Gravity', lp%gravity)
      ! Turn off drag
      lp%drag_model='none'
      ! Set collision timescale
      lp%tau_col=15.0_WP*time%dt
      ! Particle properties
      call param_read('Coefficient of restitution',lp%e_n)
      call param_read('Wall restitution',lp%e_w)
      call param_read('Friction coefficient',lp%mu_f)
      call param_read('Particle diameter',d)
      call param_read('Particle position',y0)
      pos=0.0_WP
      call param_read('Particle velocity',v0)
      vel=0.0_WP
      call param_read('Number of particles',np)
      call param_read('Elastic modulus',lp%E)
      call param_read('Shear modulus',lp%Eshear)
      call param_read('Free surface energy',lp%gamma)
      call param_read('Yield strength',lp%sigma_y)
      if (lp%cfg%amRoot) then
         dx=lp%cfg%xL/np
         dtheta=90.0_WP/np
         call lp%resize(np)
         do i=1,np
            ! print *, (i-1)*dtheta+1.0_WP
            print *, i*v0
            ! Initialize with one particle
            lp%p(i)%d=d
            lp%p(i)%col=0
            lp%p(i)%colId=0
            lp%p(i)%pos=[(i-1.0_WP)*dx+0.5_WP*dx-0.5_WP*lp%cfg%xL, y0, 0.0_WP]
            lp%p(i)%debug=0.0_WP
            lp%p(i)%id=int(i,8)
            lp%p(i)%delta_t=0.0_WP
            ! lp%p(i)%vel=v0*[COSD((i-1)*dtheta+1.0_WP), -SIND((i-1)*dtheta+1.0_WP), 0.0_WP]
            lp%p(i)%vel=i*v0*[0.0_WP, -1.0_WP, 0.0_WP]
            lp%p(i)%Acol=0.0_WP
            lp%p(i)%Tcol=0.0_WP
            lp%p(i)%dt=time%dt
            lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
            lp%p(i)%flag=0
         end do

         ! lp%p(1)%vel=v0*[ 0.0_WP, 0.0_WP, 0.0_WP]
         ! lp%p(2)%vel=v0*[ 1.0_WP, 0.0_WP, 0.0_WP]
         ! lp%p(3)%vel=v0*[-1.0_WP, 0.0_WP, 0.0_WP]
         ! lp%p(4)%vel=v0*[ 0.0_WP, 1.0_WP, 0.0_WP]
         ! lp%p(5)%vel=v0*[ 0.0_WP,-1.0_WP, 0.0_WP]

         ! lp%p(1)%pos=[ 0.0_WP, 0.05_WP, 0.0_WP]
         ! lp%p(2)%pos=[real(-5.5E-6, WP), 0.05_WP, 0.0_WP]
         ! lp%p(3)%pos=[real( 5.5E-6, WP), 0.05_WP, 0.0_WP]
         ! lp%p(4)%pos=[ 0.0_WP, 0.05_WP-real(5.5E-6, WP), 0.0_WP]
         ! lp%p(5)%pos=[ 0.0_WP, 0.05_WP+real(5.5E-6, WP), 0.0_WP]
      end if
      call lp%sync()
    end block initialize_lpt


    ! Create partmesh object for Lagrangian particle output
    create_pmesh: block
      integer :: i
      pmesh=partmesh(nvar=2,nvec=2,name='lpt')
      pmesh%varname(1)='id'
      pmesh%varname(2)='diameter'
      pmesh%vecname(1)='velocity'
      pmesh%vecname(2)='ang_vel'
      call lp%update_partmesh(pmesh)
      do i=1,lp%np_
         pmesh%var(1,i)=real(lp%p(i)%id,WP)
         pmesh%var(2,i)=lp%p(i)%d
         pmesh%vec(:,1,i)=lp%p(i)%vel
         pmesh%vec(:,2,i)=lp%p(i)%debug
      end do
    end block create_pmesh


    ! Add Ensight output
    create_ensight: block
      save_count=0
      ! Create Ensight output from cfg
      ens_out=ensight(cfg=lp%cfg,name='rebound')
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
      call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
      call lp%get_max()
      call communicate_part
      ! Create simulation monitor
      mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
      call mfile%add_column(time%n,'Timestep number')
      call mfile%add_column(time%t,'Time')
      call mfile%add_column(time%dt,'Timestep size')
      call mfile%add_column(lp%np,'Particle number')
      call mfile%add_column(lp%ncol,'Particle collisions')
      call mfile%add_column(lp%p(1)%Acol(1),'Particle X')
      call mfile%add_column(pos(2),'Particle Y')
      call mfile%add_column(pos(3),'Particle Z')
      call mfile%add_column(vel(1),'Particle U')
      call mfile%add_column(vel(2),'Particle V')
      call mfile%add_column(vel(3),'Particle W')
      call mfile%add_column(angvel(1),'Particle WX')
      call mfile%add_column(angvel(2),'Particle WY')
      call mfile%add_column(angvel(3),'Particle WZ')
      call mfile%write()
      ! Create CFL monitor
      cflfile=monitor(lp%cfg%amRoot,'cfl')
      call cflfile%add_column(time%n,'Timestep number')
      call cflfile%add_column(time%t,'Time')
      call cflfile%add_column(lp%CFLp_x,'Particle xCFL')
      call cflfile%add_column(lp%CFLp_y,'Particle yCFL')
      call cflfile%add_column(lp%CFLp_z,'Particle zCFL')
      call cflfile%add_column(lp%CFL_col,'Collision CFL')
      call cflfile%write()
    end block create_monitor

  end subroutine simulation_init


  !> Perform an NGA2 simulation
  subroutine simulation_run
    implicit none

    ! Perform time integration
    do while (.not.time%done())

       ! Increment time
       call lp%get_cfl(time%dt,cflc=time%cfl,cfl=time%cfl)
       call time%adjust_dt()
       call time%increment()

       ! Collide particles
       call lp%collide_bons(dt=time%dt)

       ! Advance particles by dt
       call lp%advance(dt=time%dt)

       ! Output to ensight
       if (ens_evt%occurs()) then
          update_pmesh: block
            integer :: i
            call lp%update_partmesh(pmesh)
            save_count=save_count + 1
            do i=1,lp%np_
               pmesh%var(1,i)=real(lp%p(i)%id,WP)
               pmesh%var(2,i)=lp%p(i)%d
               pmesh%vec(:,1,i)=lp%p(i)%vel
               pmesh%vec(:,2,i)=lp%p(i)%debug
            end do
          end block update_pmesh
          call ens_out%write_data(time%t)
       end if

       ! Perform and output monitoring
       call lp%get_max()
       call communicate_part
       call mfile%write()
       call cflfile%write()
        
       if (save_count.gt.1000) call die("Done!")
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

  end subroutine simulation_final

end module simulation
