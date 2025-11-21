!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class, only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg
   
   real(WP), public :: Rcyl

   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      real(WP) :: B,T,P,C,E,R
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))
         
         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)/real(nx,WP)*Lx-0.4_WP*Lx
         end do
         do j=1,ny+1
            y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
         end do
         do k=1,nz+1
            z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
         end do
         
         ! General serial grid object
         grid=sgrid(coord=cartesian,no=1,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='blade')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use messager, only: die
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         if (partition(2).ne.1) call die("[create_cfg] IB translation not compatible with y-partitioning")
         ! Create partitioned grid
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      

      ! Create IB walls for this config
      create_walls: block
         use mathtools,      only: Pi, twoPi
         use ibconfig_class, only: sharp
         real(WP) :: theta,Xtheta,Ytheta,dTheta,dist,dmin,minTheta
         real(WP) :: dxmin,dymin,dXdTheta,dYdTheta,right,left,nx,ny,inside
         integer :: i,j,k,n,nTheta,nMin,nLeft,nRight,nIter,blade,nBlade,iMin,nOver
         real(WP) :: x0, y0, x1, x2, y1, y2, cx, cy, dx, dy, kCoeff
         real(WP) :: xTop,xBot,xCenter,yTop,yBot,yCenter
         real(WP) :: dt, d, sign, cutThresh
         real(WP), dimension(:,:), allocatable :: Xv, Yv, XvAll, YvAll, bladeCenters
         real(WP), dimension(:), allocatable :: phi, signs, distances 
         logical, dimension(:), allocatable :: flip
         real(WP), dimension(2) :: normal, x, v1, v2

         call param_read('Base shape coefficient', B, default=2.0_WP)
         call param_read('Thickness fraction', T, default=0.2_WP)
         call param_read('Taper exponent', P, default=1.0_WP)
         call param_read('Camber', C, default=0.05_WP)
         call param_read('Camber exponent', E, default=1.0_WP)
         call param_read('Reflex parameter', R, default=0.0_WP)
         call param_read('Cut threshold', cutThresh, default=0.95_WP)

         nTheta=800
         dTheta=twoPi/(real(nTheta, WP))

         nBlade=3

         ! Blade locations, rotation angles, and x-axis relfection bools
         allocate(bladeCenters(1:nBlade, 2)); bladeCenters=0.0_WP
         allocate(phi(1:nBlade)); phi=0.0_WP
         allocate(flip(1:nBlade)); flip=.false.
         allocate(signs(1:nBlade)); signs=1.0_WP
         allocate(distances(1:nBlade)); distances=huge(1.0_WP)

         allocate(XvAll(1:nBlade, 1:nTheta)); XvAll=0.0_WP
         allocate(YvAll(1:nBlade, 1:nTheta)); YvAll=0.0_WP

         phi = [Pi/1.5_WP, -Pi/6.0_WP, -Pi/6.0_WP]!, Pi/1.5_WP]

         bladeCenters(1, 1) = -0.8_WP
         bladeCenters(1, 2) = -0.6_WP

         !bladeCenters(4, 1) = -0.5_WP
         !bladeCenters(4, 2) = +0.33_WP

         bladeCenters(2, 1) = +0.3_WP
         bladeCenters(2, 2) = -0.2_WP

         bladeCenters(3, 1) = +0.3_WP
         bladeCenters(3, 2) = +0.6_WP

         flip(1) = .true.

         do blade=1,nBlade
            nOver=0
            do n=1,ntheta
               theta=n*dTheta
               x0=1.5_WP*getX(theta)
               y0=1.5_WP*getY(theta)

               if (x0.le.1.5_WP*cutThresh) nOver=nOver+1

               x1 = x0*cos(phi(blade)) - y0*sin(phi(blade))
               y1 = x0*sin(phi(blade)) + y0*cos(phi(blade))

               if (flip(blade)) then
                  XvAll(blade, n) = -1.0_WP*x1 + bladeCenters(blade, 1)
               else
                  XvAll(blade, n) = x1 + bladeCenters(blade, 1)
               end if

               y1 = y1 + bladeCenters(blade, 2)

               YvAll(blade, n) = y1
            end do
         end do

         allocate(Xv(1:nBlade, 1:nOver)); Xv=0.0_WP
         allocate(Yv(1:nBlade, 1:nOver)); Yv=0.0_WP

         do blade=1,nBlade
            i=0
            do n=1,ntheta
               theta=n*dTheta
               x0=1.5_WP*getX(theta)
               y0=1.5_WP*getY(theta)

               if (x0.ge.1.5_WP*cutThresh) cycle

               x1 = x0*cos(phi(blade)) - y0*sin(phi(blade))
               y1 = x0*sin(phi(blade)) + y0*cos(phi(blade))

               i=i+1
               if (flip(blade)) then
                  Xv(blade, i) = -1.0_WP*x1 + bladeCenters(blade, 1)
               else
                  Xv(blade, i) = x1 + bladeCenters(blade, 1)
               end if

               y1 = y1 + bladeCenters(blade, 2)

               Yv(blade, i) = y1
            end do
         end do

         
         ! Create IB field
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_

                  x = [cfg%xm(i), cfg%ym(j)]

                  distances=huge(1.0)
                  signs=1.0_WP

                  ! Loop over number of blades
                  do blade=1,nBlade

                     ! Set dmin
                     dmin = huge(1.0)

                     ! Minimum distance to any edge of the current blade
                     do n = 1, nOver 
                        v1 = [Xv(blade,n), Yv(blade,n)]
                        if (n < nOver) then
                           v2 = [Xv(blade,n+1), Yv(blade,n+1)]
                        else
                           v2 = [Xv(blade,1), Yv(blade,1)]
                        end if
                        d = point_segment_distance(x, v1, v2)
                        if (d < dmin) dmin = d
                     end do

                     ! Sign from inside test
                     if (is_inside_polygon(x, Xv(blade,:), Yv(blade,:), nOver)) then
                        signs(blade) = +1.0_WP
                     else
                        signs(blade) = -1.0_WP
                     end if

                     ! Store dmin
                     distances(blade)=dmin
                  end do

                  ! Find the closest blade
                  iMin=minloc(distances, dim=1)

                  ! Store signed distance
                  cfg%Gib(i,j,k) = signs(iMin)*distances(iMin)
               end do
            end do
         end do
         ! Get normal vector
         call cfg%calculate_normal()
         ! Get VF field
         call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
      end block create_walls

      contains 

         function getX(theta) result(X)
            implicit none
            real(WP) :: theta, X
            X = 0.5_WP + 0.5_WP*abs(cos(theta))**B/cos(theta)
         end function getX

         function getY(theta) result(Y)
            use mathtools, only: Pi, twoPi
            implicit none
            real(WP) :: theta, X, Y
            X=getX(theta)
            Y=0.5_WP*T*abs(sin(theta))**B/sin(theta)*(1.0_WP - X**P) + C*sin(X**E*Pi) + R*sin(X*twoPi)
         end function getY

        ! Compute the distance between a point x and a line segment v1--v2
        real(WP) function point_segment_distance(x, v1, v2)
          implicit none
          real(WP), intent(in) :: x(2), v1(2), v2(2)
          real(WP) :: seg(2), t, proj(2), dx(2)

          seg = v2 - v1
          if (dot_product(seg, seg) == 0.0) then
             ! Degenerate segment (point)
             point_segment_distance = sqrt(sum((x - v1)**2))
             return
          end if

          t = dot_product(x - v1, seg) / dot_product(seg, seg)
          t = max(0.0, min(1.0, t))

          proj = v1 + t * seg
          dx = x - proj
          point_segment_distance = sqrt(sum(dx**2))
        end function point_segment_distance

        ! Ray casting algorithm: returns .true. if point is inside polygon
        logical function is_inside_polygon(x, Xv, Yv, nvert)
          implicit none
          real(WP), intent(in) :: x(2)
          real(WP), intent(in) :: Xv(:), Yv(:)
          integer, intent(in) :: nvert
          integer :: i, j
          logical :: c
          real :: x1, y1, x2, y2, x0, y0, xcross

          x0 = x(1)
          y0 = x(2)
          c = .false.

          j = nvert
          do i = 1, nvert
             x1 = Xv(i)
             y1 = Yv(i)
             x2 = Xv(j)
             y2 = Yv(j)

             if ( ((y1 > y0) .neqv. (y2 > y0)) ) then
                xcross = x1 + (y0 - y1) * (x2 - x1) / (y2 - y1)
                if (xcross > x0) c = .not. c
             end if
             j = i
          end do

          is_inside_polygon = c
        end function is_inside_polygon

   end subroutine geometry_init
   
   
end module geometry
