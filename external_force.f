
      subroutine external_force(mspace,ntype,npairs,pairsv,nfilasv,mrijv
     +     ,mxijv,extforce)

c     mspace(1,i) = id -- label of each particle                   [out]      
c     mspace(2 to 4,i) = x-- coordinates of particles              [in/out]
c     mspace(5 to 7,i) = vx-- velocities of particles              [in/out]
c     mspace(8,i) = mass-- mass of particles                       [in]
c     mspace(9,i) = rho-- dnesities of particles                   [in/out]
c     mspace(10,i) = p-- pressure of particles                     [in/out]
c     mspace(11,i) = u-- internal energy of particles              [in/out]
c     mspace(12,i) = itype-- types of particles                    [in]
c     mspace(13,i) = hsml-- smoothing lengths of particles         [in/out]
c     mspace(14,i) = c-- sound velocity of particle                [out]
c     mspace(15,i) = s-- entropy of particles                      [out]
c     mspace(16,i) = e-- total energy of particles                 [out]
c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [out]
c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [out]
c     mspace(23,i) = du        : du = du/dt                        [out]
c     mspace(24,i) = ds        : ds = ds/dt                        [out]
c     mspace(25,i) = drho      : drho = drh,o/dt                   [out]

c     npairs    : maximum number of pairs interaction             [in]
c     nfilas    : number of interaction for each fluid particle   [in]
c     mrij      : matrix of ri-rj for all fluid particles for each interaction [in]
c     mxij      : matrix of xi-xj for all fluid particles for each interaction [in]

      implicit none
      include 'param.inc'

      integer i,j,npairs,ntype(2)
      double precision mspace(25,nmax)
      integer pairsv(npairs,ntype(1)),nfilasv(ntype(1))
      double precision mrijv(npairs,ntype(1)),mxijv(3,npairs,ntype(1))
      double precision extforce(3,ntype(1)),dd,p1,p2,f,fx,fy,fz,r0

      double precision rijv,xijv,yijv,zijv

      do i=1,ntype(1)
         extforce(1,i) = 0.0d0
         extforce(2,i) = 0.0d0
         extforce(3,i) = 0.0d0
      enddo
      
      do i=1,ntype(1)-ntype(2)
         extforce(3,i) = -g
      enddo

      dd = g * ht !/ 10000000000000. 
      p1 = 12
      p2 = 6
      r0 = h0/1.1!kappa0 * h0!h0
      
      do i=1,ntype(1)-ntype(2)
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         f = 0.0
         do j=1,nfilasv(i)
         if(mrijv(j,i).lt.r0)then
            f = dd*( (r0/mrijv(j,i))**p1 - (r0/mrijv(j,i))**p2  ) /
     +           mrijv(j,i)**2
            fx = fx + f*mxijv(1,j,i)
            fy = fy + f*mxijv(2,j,i)
            fz = fz + f*mxijv(3,j,i)
c            write(*,*)i,j,nfilasv(i),mrijv(j,i),fx,fy,fz
c            call radioij(mspace(2,i),mspace(3,i),mspace(4,i),
c     +           mspace(2,pairsv(j,i)),mspace(3,pairsv(j,i)),
c     +           mspace(4,pairsv(j,i)),rijv,xijv,yijv,zijv)
c     write(*,*)i,j,r0,mrijv(j,i),f,mxijv(1,j,i),mxijv(3,j,i),fx,fy,fz
            
            endif
         enddo
         extforce(1,i) = extforce(1,i) + fx
         extforce(2,i) = extforce(2,i) + fy
         extforce(3,i) = extforce(3,i) + fz
      enddo
      
      end
