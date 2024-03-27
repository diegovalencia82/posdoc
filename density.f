c     Aaproximación SPH para la densida

      subroutine density(mspace,ntype,npairs,pairs,nfilas,w)
      implicit none
      include 'param.inc'

      integer i,j
      integer ntype(2),npairs,pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision mspace(25,nmax),w(npairs,ntype(1)),rho

      do i=1,ntype(1)
c         if(i.eq.150)write(*,*)mspace(2,i),mspace(4,i)
         rho = 0
         do j=1,nfilas(i)
            rho = rho + mspace(8,pairs(j,i))*w(j,i)
c         if(i.eq.150)write(*,*)mspace(2,pairs(j,i)),mspace(4,pairs(j,i))
         enddo
c         write(*,*)i,j,rho
         mspace(9,i) = rho
      enddo
      
      end
