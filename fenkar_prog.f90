MODULE fenkar

	IMPLICIT NONE

	! dimension of the state
	! Note: when nx*ny*nz*nv >= (256*256*8*3) = 1572864: segfault
    	! extend stack size with `ulimit -s unlimited` at command line
	INTEGER, PARAMETER	:: nx = 200
	INTEGER, PARAMETER	:: ny = 200
	INTEGER, PARAMETER	:: nz = 50

CONTAINS

	FUNCTION H(x,k)RESULT(Hkx)
	! !DESCRIPTION:
	!	This function computes H(k*x) where H is a sigmoid of bias-tanh type
	IMPLICIT NONE
	REAL			:: x, k, Hkx
	Hkx = 0.5*(1.0 + tanh(k*x))
	END FUNCTION H

	SUBROUTINE init_state(x, p)
	!$ use omp_lib

	IMPLICIT NONE

	! !ARGUMENTS:
	REAL, INTENT(inout)	:: x(nx,ny,nz,3)    ! Model state
	REAL, INTENT(in)	:: p(13)
	INTEGER			:: ix,iy,iz,iv
	REAL			:: xx,yy,zz,th,pi

	pi = 4.0*atan(1.0)
	!$omp parallel do collapse(3) private(ix,iy,iz,xx,yy,zz,th)
	DO iz = 1,nz
		DO iy = 1,ny
			DO ix = 1,nx
				zz = real(iz)/real(nz)
				yy = real(iy)/real(ny) - 0.0
				xx = real(ix)/real(nx) - 0.0
				th = ATAN2(yy,xx)
				x(ix,iy,iz,1) = 0.00 + 1.0*H(xx-0.6,10.0)
				x(ix,iy,iz,2) = 1.00 - 1.0*H(yy-0.4,10.0)
				x(ix,iy,iz,3) = 1.00 - 1.0*H(zz-0.5,1.0)
				!x(ix,iy,iz,1) = 0.00 + 0.20*xx + 0.00*yy + 0.00*zz
				!x(ix,iy,iz,2) = 1.00 - 0.00*xx - 0.20*yy - 0.00*zz
				!x(ix,iy,iz,3) = 1.00 - 0.10*xx - 0.10*yy - 0.00*zz
				!x(ix,iy,iz,1)	= cos(0.5*th)**2.0
				!x(ix,iy,iz,2)	= 1.0 - 0.50*cos(0.5*th - (pi/3))**2.0
				!x(ix,iy,iz,3) 	= 1.0 - 0.25*cos(0.5*th - (pi/3) + (pi/3)*zz)**2.0
			END DO
		END DO
	END DO
	!$omp end parallel do

	END SUBROUTINE init_state

	SUBROUTINE read_state(x, FNAME)

	IMPLICIT NONE

	REAL			:: x(nx,ny,nz,3) ! fields
	CHARACTER(len=36)	:: FNAME

	PRINT *, 'Reading: ', FNAME
	OPEN(99, FILE=FNAME, FORM="unformatted")
	READ(99) x
	CLOSE(99)

	END SUBROUTINE read_state

	SUBROUTINE write_state(x, FNAME)

	IMPLICIT NONE

	REAL			:: x(nx,ny,nz,3) ! fields
	CHARACTER(len=36)	:: FNAME

	PRINT *, 'Writing: ', FNAME
	OPEN(99, FILE=FNAME, FORM="unformatted")
	WRITE(99) x
	CLOSE(99)

	END SUBROUTINE write_state

	SUBROUTINE update_state(x, dxdt, dt)
	!$ use omp_lib

	IMPLICIT NONE

	! !ARGUMENTS:
	REAL, INTENT(in)	:: dxdt(nx,ny,nz,3) ! Time derivative
	REAL, INTENT(inout)	:: x(nx,ny,nz,3)    ! Model state
	REAL, INTENT(in)	:: dt
	INTEGER			:: ix,iy,iz,iv

	! spatial loop(s) go here; parallelize with OpenMP?
	!$omp parallel do simd collapse(4)
	DO iz=1,nz
		DO iy=1,ny
			DO ix=1,nx
				DO iv=1,3
					x(ix,iy,iz,iv) = x(ix,iy,iz,iv) + dt * dxdt(ix,iy,iz,iv)
				END DO
			END DO
		END DO
	END DO
	!$omp end parallel do simd

	END SUBROUTINE update_state

	SUBROUTINE fenkar_dxdt(dxdt, x, p, t)	! Would be nice for (dxdt,x,p,t) signature for varying model params

	! !DESCRIPTION:
	! 	This subroutine computes the time derivative for the Fenton-Karma model

	!$ use omp_lib

	IMPLICIT NONE

	! !ARGUMENTS:
	REAL, INTENT(out)	:: dxdt(nx,ny,nz,3) 	! Time derivative
	REAL, INTENT(in)	:: x(nx,ny,nz,3)    	! Model state
	REAL, INTENT(in)	:: p(13)		! Model parameters
	REAL, INTENT(in)	:: t    		! Model time
	REAL			:: tsi, tv1m, tv2m, tvm, tvp, twm, twp
	REAL			:: td, to, tr, xk, uc, uv, ucsi, dx, D
	REAL			:: Isi, Iso, Ifi, Istim
	INTEGER			:: ix,iy,iz

	! assign model parameter values
	tsi=p(1) 	!30.0
	tv1m=p(2)	!1250.0
	tv2m=p(3)	!19.6
	tvp=p(4)	!3.33
	twm=p(5)	!41.0
	twp=p(6)	!870.0
	td=p(7)		!0.25
	to=p(8)		!12.5
	tr=p(9)		!33.0
	xk=p(10)	!10.0
	uc=p(11)	!0.13
	uv=p(12)	!0.04
	ucsi=p(13)	!0.85
	dx=0.02 	! cm
	D=0.0011/(dx*dx)! (cm^2/ms) /(cm^2) = 1/ms

	! spatial loop(s) go here; parallelize with OpenMP?
	!$omp parallel do collapse(3) private(ix,iy,iz,Iso,Isi,Ifi,tvm,Istim)
	DO iz=1,nz
		DO iy=1,ny
			DO ix=1,nx

				! Compute diffusion
				! NOTE: this uses a simple no-flux condition (homogeneous von Neumann)
				!	du(x=0,L,t)/dx = 0
				! with no fiber orientation
				! TODO: add (x,y) fiber orientation with depth z?
				! This impacts the boundary conditions on each z-level:
				!	(n.D.grad(u)) = [n_x,n_y,0]'*D*[u_x,u_y,u_z] = 0, for z=/=za,zb
				!	(n.D.grad(u)) = [0,0,n_z]'*D*[u_x,u_y,u_z] = u_z = 0, for z=za,zb
				! and impacts the diffusion stencil through the depth
				!IF (nx.gt.1) THEN
					IF ((ix.gt.1) .AND. (ix.lt.nx)) THEN 	! not at boundary
						dxdt(ix,iy,iz,1) = D*(x(ix-1,iy,iz,1) + x(ix+1,iy,iz,1) - 2*x(ix,iy,iz,1))
					ELSEIF (ix.eq.1) THEN			! at left boundary
						dxdt(ix,iy,iz,1) = D*(2*x(ix+1,iy,iz,1) - 2*x(ix,iy,iz,1))
					ELSEIF (ix.eq.nx) THEN			! at right boundary
						dxdt(ix,iy,iz,1) = D*(2*x(ix-1,iy,iz,1) - 2*x(ix,iy,iz,1))
					ENDIF
				!ENDIF
				!IF (ny.gt.1) THEN
					IF ((iy.gt.1) .AND. (iy.lt.ny)) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(x(ix,iy-1,iz,1) + x(ix,iy+1,iz,1) - 2*x(ix,iy,iz,1))
					ELSEIF (iy.eq.1) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(2*x(ix,iy+1,iz,1) - 2*x(ix,iy,iz,1))
					ELSEIF (iy.eq.ny) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(2*x(ix,iy-1,iz,1) - 2*x(ix,iy,iz,1))
					ENDIF
				!ENDIF
				!IF (nz.gt.1) THEN
					IF ((iz.gt.1) .AND. (iz.lt.nz)) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(x(ix,iy,iz-1,1) + x(ix,iy,iz+1,1) - 2*x(ix,iy,iz,1))
					ELSEIF (iz.eq.1) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(2*x(ix,iy,iz+1,1) - 2*x(ix,iy,iz,1))
					ELSEIF (iz.eq.nz) THEN
						dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) &
						+ D*(2*x(ix,iy,iz-1,1) - 2*x(ix,iy,iz,1))
					ENDIF
				!ENDIF

				! tvm
				tvm = tv1m*H(uv-x(ix,iy,iz,1),100.0) + tv2m*H(x(ix,iy,iz,1)-uv,100.0)

				! Istim
				dxdt(ix,iy,iz,1) = dxdt(ix,iy,iz,1) + 0.0 & ! 0.05*sin(t)**500.0
				! Isi
				- -x(ix,iy,iz,3)*H(x(ix,iy,iz,1)-ucsi,xk)/tsi &
				! Iso
				- x(ix,iy,iz,1)*H(uc-x(ix,iy,iz,1),100.0)/to &
				- H(x(ix,iy,iz,1)-uc,100.0)/tr &
				! Ifi
				- -x(ix,iy,iz,2)*H(x(ix,iy,iz,1)-uc,100.0)*(1.0-x(ix,iy,iz,1))*(x(ix,iy,iz,1)-uc)/td

				! update dx/dt
				!dxdt(ix,iy,iz,1) = Istim - Isi - Iso - Ifi
				dxdt(ix,iy,iz,2) = H(uc-x(ix,iy,iz,1),100.0)*(1.0-x(ix,iy,iz,2))/tvm &
				 - H(x(ix,iy,iz,1)-uc,100.0)*x(ix,iy,iz,2)/tvp
				dxdt(ix,iy,iz,3) = H(uc-x(ix,iy,iz,1),100.0)*(1.0-x(ix,iy,iz,3))/twm &
				 - H(x(ix,iy,iz,1)-uc,100.0)*x(ix,iy,iz,3)/twp

			END DO
		END DO
	END DO
	!$omp end parallel do

	END SUBROUTINE fenkar_dxdt

END MODULE fenkar

PROGRAM main

!$ use omp_lib
USE fenkar

IMPLICIT NONE

	REAL			:: dt, TT, t0	! time-step, time interval, initial time
	REAL			:: t		! current time
	INTEGER			:: it, Nt, wt	! iteration, time-steps, write interval
	REAL			:: dxdt(nx,ny,nz,3), x(nx,ny,nz,3) ! fields
	CHARACTER(len=36)	:: FNAME
	CHARACTER(len=36)	:: FINIT
	REAL			:: p(13)	! parameter vector
	INTEGER			:: ix, iy, iz	! iteration

	! set time step, time interval, and initial time
	dt = 0.01	! [ms]
	TT = 200.00
	t0 = 0.00

	! set number of time-steps
	Nt = Int(ceiling(TT/dt))	! Use Ceil because dt_true < dt_in is okay, dt_true > dt_in is not
	dt = TT/Nt			! update dt_in <- dt_true
	wt = Int(ceiling(2.0/dt))	! wt*dt == time interval

	! print out how large the system is and how long we're integrating
	PRINT *, 'Domain is: (', nx, ', ', ny, ', ', nz, '), Timespan: ', TT, ', timestep: ', dt

	! set initial conditions for x and t, set dxdt to zero
	dxdt(:,:,:,:) = 0.0
	p = (/ 30.0, 50.0, 19.6, 3.33, 41.0, 200.0, 0.25, 12.5, 33.0, 10.0, 0.13, 0.04, 0.85 /)
	t = t0

	!CALL init_state(x,p)
	!WRITE (FINIT, '(A13,I3,A1,I3,A1,I3,A1,A4,A4)') './out/fenkar_', nx, '_', ny, '_', nz, '_', 'init', '.dat' ! 6+7+4+1+4+1+4+1+4+4=36
	!PRINT *, FINIT
	!CALL read_state(x, FINIT)

	it=0
	WRITE (FNAME, '(A13,I3,A1,I3,A1,I3,A1,I4,A4)') './out/fenkar_', nx, '_', ny, '_', nz, '_', int(it/wt), '.dat' ! 6+7+4+1+4+1+4+1+4+4=36
	CALL write_state(x, FNAME)
	!PRINT *, FNAME
	!OPEN(99, FILE=FNAME, FORM="unformatted")
	!WRITE(99) x
	!CLOSE(99)

	!$ACC DATA COPY(x,dxdt,p,t)
	! temporal loop goes here; this is serial (parallel-in-time solvers are a research project)
	DO it=1,Nt

		! evaluate update
		CALL fenkar_dxdt(dxdt, x, p, t)		! compute time derivative

		! update state and time
		CALL update_state(x, dxdt, dt) 		!x = x + dxdt*dt
		t = it*dt 				!t + dt

		IF (mod(it,wt).EQ.0) THEN
			! print progress
			PRINT * , 't=', t!, 't+dt=', t+dt, ' (it-1)*dt=', (it-1)*dt

			WRITE (FNAME, '(A13,I3,A1,I3,A1,I3,A1,I4,A4)') './out/fenkar_', nx, '_', ny, '_', nz, '_', int(it/wt), '.dat' ! 6+7+4+1+4+1+4+1+4+4=36
			!PRINT *, FNAME
			!OPEN(99, FILE=FNAME, FORM="unformatted")
			!WRITE(99) x
			!CLOSE(99)
			CALL write_state(x, FNAME)
		END IF

	END DO
	!$ACC END DATA
	! ending clean up
	! ideally we would write the final state of x=(u,v,w) to a file

END PROGRAM main
