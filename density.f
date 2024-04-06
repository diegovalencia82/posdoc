c     Aaproximación SPH para la densida

      subroutine density(mspace,ntype,npairs,pairs,nfilas,w)
      implicit none
      include 'param.inc'

      integer i,j
      integer ntype(2),npairs,pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision mspace(25,nmax),w(npairs,ntype(1)),rho,r,dx(dim),
     +     selfdens

      do i=1,dim
         dx(i)=0.e0
      enddo
      
      r = 0.
      do i=1,ntype(1)
         call Kernel(r,dx,mspace(13,i),selfdens,dx)
         mspace(9,i) = selfdens * mspace(8,i)
      enddo
      
      do i=1,ntype(1)
         rho = 0
         do j=1,nfilas(i)
            mspace(9,i) = mspace(9,i) + mspace(8,pairs(j,i))*w(j,i)
c            rho = mspace(9,pairs(j,i)) + mspace(8,pairs(j,i))*w(j,i)
         enddo
c     mspace(9,i) = rho
c         write(*,*)i,mspace(9,i),mspace(8,i)
      enddo
      
      end
