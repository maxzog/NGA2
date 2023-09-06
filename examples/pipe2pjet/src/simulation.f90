!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,       only: WP
   use pjet_class,      only: pjet
   use pipe_class,      only: pipe
   implicit none
   private
   
   !> Pipe simulation
   type(pipe) :: turb
   logical :: isInPIPEGrp
   
   !> Pulsed jet simulation
   type(pjet) :: flow
   
   public :: simulation_init,simulation_run,simulation_final
   
contains
   
   
   !> Initialization of our simulation
   subroutine simulation_init
      use mpi_f08, only: MPI_Group
      implicit none
      type(MPI_Group) :: hit_group
      
      ! Initialize crossflow simulation
      call flow%init()
      
      ! Create an MPI group using leftmost processors only
      create_hit_group: block
         use parallel, only: group
         use mpi_f08,  only: MPI_Group_incl
         integer, dimension(:), allocatable :: ranks
         integer, dimension(3) :: coord
         integer :: ngrp,ierr,ny,nz
         ngrp=flow%cfg%npy*flow%cfg%npz
         allocate(ranks(ngrp))
         ngrp=0
         do nz=1,flow%cfg%npz
            do ny=1,flow%cfg%npy
               ngrp=ngrp+1
               coord=[0,ny-1,nz-1]
               call MPI_CART_RANK(flow%cfg%comm,coord,ranks(ngrp),ierr)
            end do
         end do
         call MPI_Group_incl(group,ngrp,ranks,hit_group,ierr)
         if (flow%cfg%iproc.eq.1) then
            isInPIPEGrp=.true.
         else
            isInPIPEGrp=.false.
         end if
       end block create_hit_group
      
      ! Prepare HIT simulation
      if (isInPIPEGrp) then
         prepare_hit: block
            real(WP) :: dt
            ! Initialize HIT
            call turb%init(group=hit_group)
            ! Run HIT until t/tau_eddy=20
            dt=0.15_WP*turb%cfg%min_meshsize/turb%Urms_tgt !< Estimate maximum stable dt
            do while (turb%time%t.lt.2.0_WP*turb%tau_tgt)
               call turb%step(dt)
            end do
            ! Reset time
            turb%time%t=0.0_WP
         end block prepare_hit
      end if      
      
   end subroutine simulation_init
   
   
   !> Run the simulation
   subroutine simulation_run
     implicit none
      
      ! Crossflow drives overall time integration
      do while (.not.flow%time%done())
         
         ! Advance crossflow simulation
         call flow%step()
         
         ! Advance HIT simulation and transfer velocity info
         if (isInPIPEGrp) then
            ! Advance HIT
            call turb%step(flow%time%dt)
            ! Transfer turbulent velocity from hit to flow
            apply_boundary_condition: block
               use incomp_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n,i,j,k,ihit
               real(WP) :: rescaling
               rescaling=1.0_WP!turb%ti/turb%Urms_tgt
               call flow%fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  ihit=i-flow%fs%cfg%imin+turb%fs%cfg%imax+1
                  flow%fs%U(i  ,j,k)=flow%Ubulk+turb%fs%U(ihit  ,j,k)*rescaling
                  flow%fs%V(i-1,j,k)=       turb%fs%V(ihit-1,j,k)*rescaling
                  flow%fs%W(i-1,j,k)=       turb%fs%W(ihit-1,j,k)*rescaling
               end do
            end block apply_boundary_condition
         end if
         
      end do
      
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Finalize crossflow simulation
      call flow%final()
      
      ! Finalize HIT simulation
      if (isInPIPEGrp) call turb%final()
      
   end subroutine simulation_final
   
   
end module simulation
