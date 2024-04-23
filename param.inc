c-------------------------------------------------------
c     Including file for parameters and constants used
c     in the entire SPH software packages.
c-------------------------------------------------------

      double precision pi
      parameter ( pi = 3.14159265358979323846 )

c     dim : Dimension of the problem (1, 2 or 3)
      integer dim
      parameter ( dim = 2)

c     nmax	: Maximum number of particles
      integer nmax
      parameter ( nmax = 4538 )
c
      integer maxtimestep, initime, fintime, inioutfile
      parameter ( initime = 0, fintime = 26000 ) 
      parameter ( maxtimestep = 26001, inioutfile = 0000) 
      character infile*80
      parameter ( infile = 'snapshot_000' )

c     Time step integration
      double precision dt0
      parameter ( dt0 = 0.001/5. )
c      parameter ( dt0 = 0.001/2.5 )	

c     Definition of the support domain, take odd numbers
c     h0 is taken from the space of particles in the initial conditions
      double precision h0, kappa0
      parameter ( kappa0 = 2., h0 = 9.2592596990504743E-003)
c     Smoothing kernel function
c     skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
c         = 2, Gauss kernel (Gingold and Monaghan 1981)
c         = 3, Quintic kernel (Morris 1997)
      integer skf
      parameter(skf = 1)

c     ht 	    Maximum high of the block, it is equal to the initial conditions
c     The measurement must be in meters
      double precision ht
      parameter ( ht = 1.e-3*500)

c     g = 9.8. Gravity
      double precision g
      parameter ( g = 9.82d0 )

c     Artificial EOS.
c     Initial density of particles, rho0
      double precision gamma, beta, rho0
      parameter ( gamma=7.0, beta = 1.0, rho0 = 1000.0)
c      parameter ( gamma=7.0, beta = 5.0, rho0 = 916.0)

c     viscous force
      double precision mu, muc
      parameter ( mu=1.002e-3, muc=1000)

c      To save steps
       integer save_step
       parameter ( save_step = 10)