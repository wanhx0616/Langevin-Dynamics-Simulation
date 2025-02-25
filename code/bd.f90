!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CODE STARTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     HERE


PROGRAM bd
!
! This is the main driver for the brownian dynamics program.
!
USE Particles
USE Simulation_Control
USE Potential
IMPLICIT NONE
DOUBLE PRECISION :: ReqTRequested, StepTRequested, rescale2
DOUBLE PRECISION, DIMENSION(3) :: ReqBoxSize, StepBoxSize, rescale
INTEGER :: i, k

CALL Initialize


!Parameters of interactions

CALL Pair(1, 1, 1.d0, 3.d0)
CALL Pair(1, 2, 0.1d0, 3.d0)
CALL Pair(1, 3, 0.1d0, 3.d0)
CALL Pair(2, 2, 3.d0, 2.d0**(1.d0/6.d0))      
CALL Pair(2, 3, 3.d0, 2.d0**(1.d0/6.d0)) 
CALL Pair(3, 3, 3.d0, 2.d0**(1.d0/6.d0))

CALL Evolve_Sample(6000000)

    
CALL Terminate
END PROGRAM bd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  INITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ROUTINES


SUBROUTINE Initialize
!
! Initialization procedure
!
USE Simulation_Control
IMPLICIT NONE
INTEGER :: i, lis

        CALL Read_Sample

CALL Initial_Printout
END SUBROUTINE Initialize


SUBROUTINE Read_Sample
!
!  Reads the initial sample from files
!
use Particles
USE Potential
use Simulation_Control
USE ran_mod
implicit none
CHARACTER(LEN=10) :: strStep
CHARACTER(LEN=9) :: form
double precision, dimension(DIM) :: Mass_center, pij
DOUBLE PRECISION :: diameter1, diameter2, shift1, shift2, atheta, sumz, min1, dis
integer :: i, j, k, m2, totalz, x, y, z, m, a, b, L, wij, M1, rotation
double precision,allocatable::vert(:)

!Read the bond informtion of the soft nanoparticle
OPEN(UNIT=1, FILE='PARA', STATUS='OLD')
READ(1,'(3I7,3E23.15)') N, maxatomtype, BListLength, BoxSize

N = 10501
m2 = 676
memnum = N
nanonum = m2
N = N + nanonum + 1
maxatomtype = maxatomtype + 4
BoxSize(1) = 0.912304068485089E+02
BoxSize(2) = 0.912304068485089E+02
BoxSize(3) = 500
volume = PRODUCT(BoxSize)
density = N / volume

CALL Allocate_Particles
DO i = 1, BListLength
    READ(1,*) BondList(:,i)
END DO
CLOSE(UNIT=1)

!Read the information  of the membrane
OPEN(UNIT=1, FILE='out.xyz', STATUS='OLD')
READ(1,*)
READ(1,*)
DO i = 1, memnum
    READ(1,'(2X,3E23.15)') pos(:,i)
END DO
CLOSE(UNIT=1)

OPEN(UNIT=1, FILE='aout.xyz', STATUS='OLD')
READ(1,*)
READ(1,*)
DO i = 1, memnum
    READ(1,'(3E23.15)') axes(:,i)
END DO
CLOSE(UNIT=1)

OPEN(UNIT=1, FILE='uvw', STATUS='OLD')
DO i = 1, memnum
    READ(1,'(3E23.15)') vel(:,i)
END DO
DO i = 1, memnum
    READ(1,'(3E23.15)') acc(:,i)
END DO
CLOSE(UNIT=1)

OPEN(UNIT=1, FILE='auvw', STATUS='OLD')
DO i = 1, memnum
    READ(1,'(3E23.15)') angvel(:,i)
END DO
DO i = 1, memnum
    READ(1,'(3E23.15)') angacc(:,i)
END DO
CLOSE(UNIT=1)

!Shift the simulation box to position the membrane at the center of the box.
Mass_Center = sum(pos(:,1:memnum) , dim=2) / (N-m2-1)
sumz = 0.d0
totalz = 0
do i = 1, memnum
	if ((abs(pos(1,i)) < 10) .and. (abs(pos(2,i)) < 10)) then
		sumz = sumz + pos(3,i)
		totalz = totalz + 1
	END IF
END DO
sumz = sumz / totalz
pos(3,1:memnum) = pos(3,1:memnum) - sumz

ptype(1:memnum) = 1


!Read the informtion of the soft nanoparticle
OPEN(UNIT=1, FILE='outnp.xyz', STATUS='OLD')
READ(1,*)
READ(1,*)
DO i = 1, m2
    READ(1,'(2X,3E23.15)') pos(:,i+memnum)
END DO
CLOSE(UNIT=1)

OPEN(UNIT=1, FILE='uvwnp', STATUS='OLD')
DO i = 1, m2
    READ(1,'(3E23.15)') vel(:,i+memnum)
END DO
DO i = 1, m2
    READ(1,'(3E23.15)') acc(:,i+memnum)
END DO
CLOSE(UNIT=1)

ptype(memnum+1:N-1) = 2

Mass_Center = sum(pos(:,memnum+1:N-1) , dim=2) / nanonum
pos(:,N) = Mass_Center
ptype(N) = 4

!Shift the simulation box to position the soft particle above the center of the membrane.
do i = 1,m2
    pos(1,i+memnum) = pos(1,i+memnum) - pos(1,N)
    pos(2,i+memnum) = pos(2,i+memnum) - pos(2,N)
    pos(3,i+memnum) = pos(3,i+memnum) - pos(3,N) + sumz + 9
end do

Mass_Center = sum(pos(:,memnum+1:N-1) , dim=2) / nanonum
pos(:,N) = Mass_Center


!Identify the orientation of the particle according to the given initial angle.
atheta = 90.d0
atheta = atheta / 180.d0 * pi
diameter1 = 16.d0

 min1 = 10000
j = 0
do i=1,m2
    dis = (pos(1,i+memnum)- pos(1,N)-diameter1/2.d0*cos(atheta))**2 +&
           (pos(2,i+memnum)- pos(2,N))**2 +&
          (pos(3,i+memnum)-pos(3,N)+diameter1/2.d0*sin(atheta))**2
    if (dis < min1) then
        min1 = dis
        j = i
    end if
end do
ptype(memnum+j) = 3


step = 0
Mass_Center = sum(pos , dim=2) / (N-1)
do k=1,DIM
    pos(k,:) = pos(k,:) - Mass_Center(k)
end do


do i= 1, N-1
   xcm(:,i)= pos(:,i)
end do

END SUBROUTINE Read_Sample


SUBROUTINE Allocate_Particles
!
!   Allocate the arrays containing atomic information
!
use Particles
USE Potential
use Neighbor_List
USE Statistics
implicit none
allocate( ptype(N) );          ptype = 0
allocate( deltar(DIM,N) );     deltar = 0.d0
allocate( pos(DIM,N) );        pos = 0.d0
allocate( vel(DIM,N) );        vel = 0.d0
allocate( acc(DIM,N) );        acc = 0.d0
ALLOCATE(fcm(3,N));            fcm = 0.d0
allocate( axes(DIM,N) );       axes = 0.d0
allocate( angvel(DIM,N) );     angvel = 0.d0
allocate( angacc(DIM,N) );     angacc = 0.d0
allocate( ene_pot(N) );        ene_pot = 0.d0
ALLOCATE( ene_pot_exbond(N) ); ene_pot_exbond = 0.d0
allocate( ene_kin(N) );        ene_kin = 0.d0
 ALLOCATE(xcm(3,N));         xcm = 0.d0
BMaxListLength = MaxBondsPerAtom*N
allocate( BondList(3, BMaxListLength) ); BondList = 0
allocate( BondMarker1(10000) ); BondMarker1=0
allocate( BondMarker2(10000) ); BondMarker2=0
! Dynamic bonds
ALLOCATE(DBType(N)); DBType = 0
ALLOCATE(ndbond(N)); ndbond = 0
ALLOCATE(DBMap1(N)); DBMap1 = 0
ALLOCATE(DBMap2(N)); DBMap2 = 0
DBNeigMaxLength = 3 * N
DBNeigLength = 0
ALLOCATE(DBNeig1(DBNeigMaxLength)); DBNeig1 = 0
ALLOCATE(DBNeig2(DBNeigMaxLength)); DBNeig2 = 0
! also those related to the neighbor list:
MaxListLength = MaxPairsPerAtom*N
allocate( List(MaxListLength) )
allocate( Marker1(N) )
allocate( Marker2(N) )
allocate( Advance(N) )
allocate( DisplaceList(DIM,N) )
! cell list
allocate( CellList(N) )
! potential
ALLOCATE(epsil(maxatomtype,maxatomtype))
ALLOCATE(Rcutoffsq(maxatomtype,maxatomtype))
ALLOCATE(phicutoff(maxatomtype,maxatomtype))
ALLOCATE(PotShift(maxatomtype,maxatomtype)); PotShift = 0.d0
! Statistics
ALLOCATE(sta_data(8,iStaSave));     sta_data = 0
PS = 0.d0
G = DIM * N
END SUBROUTINE Allocate_Particles


SUBROUTINE Pair(i, j, ep, rc)
USE Particles
USE Potential
IMPLICIT NONE
INTEGER :: i, j
DOUBLE PRECISION :: ep, rc

IF ( i > maxatomtype .or. j > maxatomtype ) THEN
    PRINT *, 'Pair: Atom type is greater than Maxatomtype!'
    STOP
END IF

IF ( i .eq. 0 .or. j .eq. 0 ) THEN
    epsil     = ep
    Rcutoffsq = (rc+PotShift)**2
    phicutoff = 4.d0 * epsil * (1.d0 / Rcutoffsq**6 - 1.d0 / Rcutoffsq**3)
ELSE
    epsil(i,j)     = ep
    Rcutoffsq(i,j) = (rc+PotShift(i,j))**2
    phicutoff(i,j) = 4.d0 * epsil(i,j) * (1.d0 / Rcutoffsq(i,j)**6 - 1.d0 / Rcutoffsq(i,j)**3)
    IF ( i .ne. j ) THEN
            epsil(j,i) =     epsil(i,j)
        Rcutoffsq(j,i) = Rcutoffsq(i,j)
        phicutoff(j,i) = phicutoff(i,j)
    END IF
END IF
END SUBROUTINE Pair


SUBROUTINE Initial_Printout
!
!   Prints information on the run parameters on the standard output
!   Leading # are to directly use the output as a gunplot input.
!
use Particles
use Neighbor_List
use Simulation_Control
implicit none
print '(a,/,a,/,a)','#','# Brownian: a minimal molecular dynamics program','#'
print '(2a)','# ',title(1:len_trim(title))
print '(2a)','#  Input sample: ', argum
if ( argum == 'new' ) then
    print '(a)','#                                                              (new model)'
else
    print '(a)','#                                               (positions read from file)'
end if
print '(2a)','# Output sample: out-step'
print '(a,f7.4)', '# Time step:',deltat
print '(a,i6)','# Number of atoms:',N
print '(a,i6)','# Number of bonds:',BListLength
print '(a,3f12.6,a,f15.3)', &
      '# Box size:',BoxSize,', Volume:',volume
print '(a,f12.6)','# Constant T run with T = ',TRequested
print '(a,f8.4,a,i9)', &
      '# Skin:',skin,' , maximum neighbor list length:',MaxListLength

print '(a,/,a,/,a)', '#', &
'#     Step   Temperature     Kinetic      Potential   Total Energy    Pressure',&
'# --------  ------------  ------------  ------------  ------------  ------------'
END SUBROUTINE Initial_Printout


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TIME EVOLUTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     ROUTINES


SUBROUTINE Evolve_Sample(runsteps)
!
!   This is the main subroutine, controlling the time evolution of the system.
!
USE Particles
USE Neighbor_List
USE Simulation_Control
USE Statistics
USE ran_mod
IMPLICIT NONE
INTEGER :: currentstep, runsteps
LOGICAL, EXTERNAL :: MovedTooMuch
LOGICAL :: ListUpdateRequested
INTEGER :: i, j
DOUBLE PRECISION :: sumvel, rescale
DOUBLE PRECISION, DIMENSION(DIM) :: xi, Mass_Center
DOUBLE PRECISION, DIMENSION(3) :: boundary

ListUpdateRequested = .TRUE. ! always need an update when starting


currentstep = 0
DO WHILE (currentstep .lt. runsteps)
    currentstep = currentstep + 1
    step = step + 1
 fcm = 0.d0   
   !$omp parallel do private(i)
   
    DO i = 1,N-1
        IF (ptype(i) .eq. 1) THEN
            vel(:,i) = vel(:,i) + 0.5d0 * deltat * acc(:,i)
            deltar(:,i) = deltat * vel(:,i)
            pos(:,i) = pos(:,i) + deltar(:,i)           
                axes(:,i) = axes(:,i) + 0.5d0 * (deltat**2) * angacc(:,i) &
					+ deltat * angvel(:,i)
                angvel(:,i) = angvel(:,i) + 0.5d0 * deltat * angacc(:,i)
        ELSE IF (ptype(i) .eq. 2 .or. ptype(i) .eq. 3) THEN
		    vel(:,i) = vel(:,i) + 0.5d0 * deltat * acc(:,i)
            deltar(:,i) = deltat * vel(:,i)
            pos(:,i) = pos(:,i) + deltar(:,i)      
		    xcm(:,i) =  xcm(:,i) + deltar(:,i)	             		 			 		 		 			  				
        END IF
    END DO
	!$omp end parallel do
           Mass_Center = sum(xcm(:,memnum+1:N-1) , dim=2) / nanonum
					 pos(:,N) = Mass_Center	

    sumvel = 0
    do i=1,memnum
        sumvel = sumvel + dot_product(vel(:,i),vel(:,i)) 
    end do
    
    tension = (-sumvel  - virial) *0.5d0 / BoxSize(1) / BoxSize(2)

    CALL Refold_Positions
    
    IF ( NSigmaT .and. (step > 1)) THEN
               
		rescale = 1 + deltat / 3 / 1.d0 * (ten0 - tension)
		
        DO i = 1, 2
            BoxSize(i) = BoxSize(i)* rescale 
            pos(i,:) = pos(i,:)* rescale 
			xcm(i,:) = xcm(i,:) * rescale
        END DO
          

    END IF
    
    IF (ListUpdateRequested) THEN
        CALL Update_List
        ListUpdateRequested = .FALSE.
    END IF
    
    CALL Compute_Forces
    
    !$omp parallel do private(i)
    DO i = 1,N-1
        IF (ptype(i) .eq. 1) THEN
            vel(:,i) = vel(:,i) + 0.5d0 * deltat * acc(:,i)            
            angvel(:,i) = angvel(:,i) + 0.5d0 * deltat * angacc(:,i)
        ELSE  IF (ptype(i) .eq. 2 .or. ptype(i) .eq. 3)  THEN			   
			vel(:,i) = vel(:,i) + 0.5d0 * deltat * acc(:,i)			   
        END IF
    END DO
	!$omp end parallel do

    DisplaceList = DisplaceList + deltar
    ListUpdateRequested = MovedTooMuch(Skin)
    
    IF ( lppp .or. lpps ) THEN
        CALL Print_Statistics
    END IF
    
    IF ( MOD(step, iStaInt) == 0 ) THEN
        CALL data_output
    END IF
    
    IF ( (MOD(step, iSave) .eq. 0) .or. (step == 1)) THEN
        CALL Write_Sample
    END IF
	
	IF ( (MOD(step, iSave2) .eq. 0) .or. (step == 1)) THEN
        CALL Write_Sample2
    END IF    
END DO
END SUBROUTINE Evolve_Sample


SUBROUTINE Refold_Positions
!
!   Particles that left the box are refolded back into the box by
!   periodic boundary conditions
!
USE Particles
USE Simulation_Control
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(3) :: boundary
DOUBLE PRECISION, DIMENSION(3) :: numbox
INTEGER :: i
boundary = BoxSize / 2.d0

DO i = 1, N-1
     IF (ptype(i) .eq. 1) THEN
        axes(:,i)=axes(:,i)/sqrt(dot_product(axes(:,i), axes(:,i)))
   END IF

    WHERE ( pos(:,i) >  boundary )   pos(:,i) = pos(:,i) - BoxSize
    WHERE ( pos(:,i) < -boundary )   pos(:,i) = pos(:,i) + BoxSize
    IF ( ANY(ABS(pos(:,i)) > boundary) ) THEN
        PRINT '(a,I6,a,I10)', 'Missing particle ', i, 'at step ', step
        CALL Write_Sample
        !pause
        STOP
    END IF
END DO

   WHERE ( pos(:,N) >  boundary )   
       numbox = ceiling((pos(:,N) - boundary)/BoxSize) 
       pos(:,N) = pos(:,N) - BoxSize * numbox
    end where
    WHERE ( pos(:,N) < -boundary )   
	   numbox = floor((pos(:,N) + boundary)/BoxSize)
       pos(:,N) = pos(:,N) - BoxSize * numbox
	end where

END SUBROUTINE Refold_Positions


SUBROUTINE Compute_Forces
!
!  Compute forces on atoms from the position, this version makes use of
!  an existing neighbor list, and runs in O(N) time.
!
use Particles
use Simulation_Control
USE Potential
use Neighbor_List
use ran_mod
use Statistics
implicit none
DOUBLE PRECISION, DIMENSION(DIM) :: Sij1, Sij, rbij, dangphii, dangphij, dangphir
DOUBLE PRECISION :: Rsqij, R, Rshift, rm2, rm4, phi, dphi, rm6, rm12, RH, RH2, diameter0
DOUBLE PRECISION :: rc, rmin, rminsq, temp
DOUBLE PRECISION :: dotmuls1, dotmuls2, dotmuls3, dotmuls4, dotmuls5, dotmuls6
DOUBLE PRECISION :: angphi, lambdai, lambdaj, dforce
DOUBLE PRECISION, DIMENSION(DIM) :: dphini, dphinj, dphir
INTEGER :: i, j, L, pti, ptj, rotation

ene_pot = 0.d0
ene_adh = 0.d0
ene_pot_exbond = 0.d0
acc = 0.d0
angacc = 0.d0
virial = 0.d0
rmin = 2.d0**(1.d0 / 6.d0)
rminsq = rmin**2
dforce = 0.004d0
isneighbor = 0

diameter0 = 16.d0            
RH = diameter0/2.0           
RH2 = RH**2                


CALL Compute_Bonds


do i = 1,N-1
    if (ptype(i) == 3) then
        rotation = i
        exit
    end if
end do


DO i = 1, N-1
    pti = ptype(i)
    DO L = Marker1(i),Marker2(i)
        j = List(L)
        ptj = ptype(j)
        Sij = pos(:,i) - pos(:,j)
        WHERE ( abs(Sij) > 0.5d0 * BoxSize )
            Sij = Sij - SIGN(BoxSize, Sij)
        END WHERE
        Rsqij = dot_product(Sij,Sij)
		
		R = SQRT(Rsqij)
		
		If  ( (pti .eq. 1) .and. (ptj .eq. 1) ) THEN
		    rbij = Sij / R   
		    dotmuls1 = dot_product(axes(:,i),axes(:,j))
		    dotmuls2 = dot_product(axes(:,i), rbij(:))
		    dotmuls3 = dot_product(axes(:,j), rbij(:))
		    angphi = 1 + mu * (dotmuls1 - dotmuls2 * dotmuls3 &
			    + sin(theta0) * (dotmuls3-dotmuls2) - sin(theta0)**2-1)
		    dangphir(:) = mu * (-axes(:,i) * dotmuls3 - axes(:,j) * dotmuls2 &
			    + sin(theta0) * (axes(:,j) - axes(:,i)))

		    dangphii(:) = mu * (axes(:,j)-rbij(:)*dotmuls3-sin(theta0)*rbij(:))
		    dangphij(:) = mu * (axes(:,i)-rbij(:)*dotmuls2+sin(theta0)*rbij(:))
    		
            IF ( Rsqij .lt. rminsq ) THEN
                rm2  = rminsq/Rsqij
                rm4  = rm2**2
                phi  =  epsil(pti,ptj) * (rm4 - 2 * rm2)
                dphi = 4.d0 * epsil(pti,ptj) * (rm4 - rm2) / Rsqij
                ene_pot(i) = ene_pot(i) + 0.5d0 * (phi+(1-angphi)*epsil(pti, ptj))  
                ene_pot(j) = ene_pot(j) + 0.5d0 * (phi+(1-angphi)*epsil(pti, ptj))    
                virial = virial - dphi * Rsqij
                acc(:,i) = acc(:,i) + dphi * Sij
                acc(:,j) = acc(:,j) - dphi * Sij
    			
			    dphir = epsil(pti,ptj) * dangphir
			    dotmuls4 = dot_product(dphir(:), rbij(:))
			    dotmuls5 = dot_product(dphir(:), Sij(:))
			    dotmuls6 = dot_product(rbij(:), Sij(:))
			    virial = virial -  (dotmuls5 - dotmuls4 * dotmuls6) / R
			    acc(:,i) = acc(:,i) + (dphir - dotmuls4 * rbij) / R
                acc(:,j) = acc(:,j) - (dphir - dotmuls4 * rbij) / R

    			
			    dphini = epsil(pti,ptj) * dangphii
			    lambdai = -dot_product(dphini,axes(:,i))!-dot_product(angvel(:,i), angvel(:,i))
			    dphinj = epsil(pti,ptj) * dangphij
			    lambdaj = -dot_product(dphinj,axes(:,j))!-dot_product(angvel(:,j), angvel(:,j))
			    angacc(:,i)= angacc(:,i) + (dphini+lambdai*axes(:,i))    
			    angacc(:,j)= angacc(:,j) + (dphinj+lambdaj*axes(:,j))    			
   			    			
            ELSE IF ((Rsqij .gt. rminsq) .and. (Rsqij .lt. Rcutoffsq(pti,ptj))) THEN
                rc = sqrt(Rcutoffsq(pti,ptj))
			    temp = pi * (R-rmin) / 2.d0 / (rc-rmin)
                phi = -epsil(pti,ptj) * cos(temp)**(2.d0*xi)
			    dphi = -epsil(pti,ptj) * xi * (cos(temp)**(2.d0*xi-1))&
				    *sin(temp)* pi / (rc-rmin) / R
			    ene_pot(i) = ene_pot(i) + 0.5d0 * phi * angphi
                ene_pot(j) = ene_pot(j) + 0.5d0 * phi * angphi
                virial = virial - dphi * angphi * Rsqij
                acc(:,i) = acc(:,i) + dphi * angphi * Sij
                acc(:,j) = acc(:,j) - dphi * angphi * Sij
    			
			    dphir = - phi * dangphir
			    dotmuls4 = dot_product(dphir(:), rbij(:))
			    dotmuls5 = dot_product(dphir(:), Sij(:))
			    dotmuls6 = dot_product(rbij(:), Sij(:))
			    virial = virial -  (dotmuls5 - dotmuls4 * dotmuls6) / R
			    acc(:,i) = acc(:,i) + (dphir - dotmuls4 * rbij) / R
                acc(:,j) = acc(:,j) - (dphir - dotmuls4 * rbij) / R

    			
			    dphini = - phi * dangphii
			    lambdai = -dot_product(dphini,axes(:,i))!-dot_product(angvel(:,i), angvel(:,i))
			    dphinj = - phi * dangphij
			    lambdaj = -dot_product(dphinj,axes(:,j))!-dot_product(angvel(:,j), angvel(:,j))
			    angacc(:,i)= angacc(:,i) + (dphini+lambdai*axes(:,i)) 
			    angacc(:,j)= angacc(:,j) + (dphinj+lambdaj*axes(:,j))
            END IF
        ELSE IF ( ((pti .eq. 1) .and. (ptj .eq. 2)) .or. ((ptj .eq. 1) .and. (pti .eq. 2)) .or. &
		        ((pti .eq. 1) .and. (ptj .eq. 3)) .or. ((ptj .eq. 1) .and. (pti .eq. 3)) ) THEN
            IF ( Rsqij .lt. Rcutoffsq(pti,ptj) ) THEN
                R = SQRT(Rsqij)
                Rshift = R - PotShift(pti,ptj)
                IF ( Rshift < 0 ) THEN
                    print *, 'barrier'
                END IF
                rm2  = 1.d0/(Rshift**2)
                rm6  = rm2**3
                rm12 = rm6**2
                phi  =  4.d0 * epsil(pti,ptj) * (rm12 - rm6)
                dphi = 24.d0 * epsil(pti,ptj) * (2.d0 * rm12 - rm6) / Rshift / R
                ene_pot(i) = ene_pot(i) + 0.5d0 * phi
				ene_adh = ene_adh + phi
                ene_pot(j) = ene_pot(j) + 0.5d0 * phi
                acc(:,i) = acc(:,i) + dphi * Sij
                acc(:,j) = acc(:,j) - dphi * Sij
				IF (Rshift < 1.3d0) then
					isneighbor = isneighbor + 1
				END IF
            end if
		ELSE IF ( ((pti .eq. 2) .and. (ptj .eq. 2)) .or. ((pti .eq. 2) .and. (ptj .eq. 3)) .or. &
		         ((pti .eq. 3) .and. (ptj .eq. 2))) THEN
            IF ( Rsqij .lt. Rcutoffsq(pti,ptj) ) THEN   		
                R = SQRT(Rsqij)
                Rshift = R
                IF ( Rshift < 0 ) THEN
                    print *, 'barrier'
                END IF
                rm2  = 1.d0/(Rshift**2)
                rm6  = rm2**3
                rm12 = rm6**2
                phi  =  4.d0 * epsil(pti,ptj) * (rm12 - rm6)   
                dphi = 24.d0 * epsil(pti,ptj) * (2.d0 * rm12 - rm6) / Rshift / R
                ene_pot(i) = ene_pot(i) + 0.5d0 * phi				
                ene_pot(j) = ene_pot(j) + 0.5d0 * phi
                acc(:,i) = acc(:,i) + dphi * Sij
                acc(:,j) = acc(:,j) - dphi * Sij
			END IF
        END IF
    END DO
		
  ! Hertzian potential
 
  IF ((pti == 2) .or. (pti == 3))then
  
        Sij = pos(:,i) - pos(:,N)
        WHERE ( abs(Sij) > 0.5d0 * BoxSize )
            Sij = Sij - SIGN(BoxSize, Sij)
        END WHERE
        Rsqij = dot_product(Sij,Sij)		
		    		 			   
			   IF ( Rsqij .lt. RH2 ) THEN				   
                   R = SQRT(Rsqij)			   
				  phi  = U*((1- R/RH)**(5/2))
				  dphi = U*(5/2)*((1- R/RH)**(3/2))/ RH/ R
                  ene_pot(i) = ene_pot(i) + 0.5d0 * phi				  
                  ene_pot(N) = ene_pot(N) + 0.5d0 * phi
                  acc(:,i) = acc(:,i) + dphi * Sij
                  acc(:,N) = acc(:,N) - dphi * Sij

		        END IF
	END IF 		
	
	
	
	
    ! frictional forces and random forces
     if (pti == 1) then
        acc(:,i) = acc(:,i) - gamma * vel(:,i)
        acc(:,i) = acc(:,i) + ( (/ran(),ran(),ran()/) - 0.5d0 ) * &
                SQRT(24.d0*TRequested*gamma/deltat)
	end if 
    
    if (pti == 2 .or. pti == 3) then          
               Sij1 = pos(:,rotation) - pos(:,N)	
                WHERE ( abs(Sij1) > 0.5d0 * BoxSize )
                Sij1 = Sij1 - SIGN(BoxSize, Sij1)
                END WHERE		
        Sij1 = Sij1 / SQRT(dot_product(Sij1,Sij1))
        acc(:,i) = acc(:,i) + dforce * Sij1
		acc(:,i) = acc(:,i) - gamma * vel(:,i)
        acc(:,i) = acc(:,i) + ( (/ran(),ran(),ran()/) - 0.5d0 ) * &
                SQRT(24.d0*TRequested*gamma/deltat)
    end if
		
END DO

Do i = memnum+1, N-1      
	  acc(:,i) =  acc(:,i) + acc(:,N) / (N-1-memnum)
end do


virial = - virial
END SUBROUTINE Compute_Forces


SUBROUTINE Compute_Bonds
!
!   Spring force
!
USE Simulation_Control
use Particles
use Potential
implicit none
double precision, dimension(DIM) :: Sij
double precision :: Rsqij, R, Rshift, phi,dphi
INTEGER :: pti, ptj
integer :: i,j,k,L
	
do k = 1, BListLength
       i = BondList(2, k)
       j = BondList(3, k)
	 Sij = pos(:,i) - pos(:,j) 
       pti = ptype(i)
       ptj = ptype(j)
        where ( abs(Sij) > 0.5d0 * BoxSize )
            Sij = Sij - sign(BoxSize,Sij)
        end where
        Rsqij = dot_product(Sij,Sij)
        R = SQRT(Rsqij)
		IF ( R > FENE_Rmax ) THEN
           PRINT '(a,I6,a,I10)', 'Missing bond ', k, 'at step ', step
           STOP
         END IF		
        Rshift = R
        Rsqij = Rshift**2
        ! phi  = -0.5d0 * FENE_K * (FENE_Rmax**2) * log(1 - Rsqij/(FENE_Rmax**2))
        dphi = -FENE_K * Rshift / R / (1 - Rsqij/(FENE_Rmax**2))
        acc(:,i) = acc(:,i) + dphi*Sij
        acc(:,j) = acc(:,j) - dphi*Sij    
end do
END SUBROUTINE Compute_Bonds

SUBROUTINE Compute_Temperature(ene_kin_aver,temperature)
use Particles
implicit none
double precision :: ene_kin_aver, temperature
integer :: i
do i=1,N
    ene_kin(i) = 0.5d0*dot_product(vel(:,i),vel(:,i))   ! kin en of each atom
end do
ene_kin_aver = sum( ene_kin ) / N
temperature = 2.d0*ene_kin_aver/DIM
END SUBROUTINE Compute_Temperature


SUBROUTINE Update_List
!
!  Update the neighbor list, including all atom pairs that are within a
!  distance Range to each other.
!
use Particles
use Simulation_Control
use Potential
use Neighbor_List
implicit none
double precision, dimension(maxatomtype,maxatomtype) :: Rrange, Rrangesq
double precision, dimension(DIM) :: Sij
double precision :: Rsqij
INTEGER :: pti, ptj
integer :: i,j,k,L,m(3)


Rrange=0
do i = 1,maxatomtype
    do j = 1, maxatomtype
        Rrange(i,j) = SQRT(Rcutoffsq(i,j)) + Skin
    end do
end do
Rrangesq = Rrange**2
call new_clist(maxval(Rrange))

L = 1
do i = 1, N
    pti = ptype(i)
    Marker1(i) = L
    iCell = 1 + int( (pos(:,i) / BoxSize + 0.5d0) * nPart )
    where(iCell .gt. nPart) iCell = nPart

    do k=0,13
        select case(k)
            case(0)
                m = iCell
            case(1)
                m = iCell + (/ 1,0,0 /)
            case(2:4)
                m = iCell + (/k-3,1,0/)
            case(5:7)
                m = iCell + (/k-6,-1,1/)
            case(8:10)
                m = iCell + (/k-9,0,1/)
            case(11:13)
                m = iCell + (/k-12,1,1/)
        end select
        ! PBC
        if ( all(m == iCell) ) then
            j = CellList(i)
        else
            where (m > nPart) m = m - nPart
            where (m <= 0) m = m + nPart
            j = Head(m(1),m(2),m(3))
        end if
        ! Linked list
        do while (j .ne. 0)
            Sij = pos(:,i) - pos(:,j)
            where ( abs(Sij) > 0.5d0 * BoxSize )
                Sij = Sij - sign(BoxSize,Sij)
            end where
            Rsqij = dot_product(Sij,Sij)
            if ( Rsqij < Rrangesq(pti,ptype(j)) ) then
                if ( L > MaxListLength ) then
                    print *,'Update_List: FATAL: List too small'
                    !pause
                    stop
                end if
                List(L) = j
                L = L + 1
            end if
            j = CellList(j)
        end do
    end do
    Marker2(i) = L - 1
end do
ListLength = L - 1
DisplaceList = 0.d0
END SUBROUTINE Update_List


SUBROUTINE new_clist(Range)
use Particles
use Neighbor_List
implicit none
double precision :: Range
integer :: i
nPart = minval(int(BoxSize / Range))
if ( nPart .lt. 3 ) then
    print *,'Range is too large!'
    stop
end if
Head = 0
DO i = 1,N
    iCell = 1 + int( (pos(:,i) / BoxSize + 0.5d0) * nPart)
    where(iCell .gt. nPart) iCell = nPart
    CellList(i) = Head(iCell(1),iCell(2),iCell(3))
    Head(iCell(1),iCell(2),iCell(3)) = i
end do
END SUBROUTINE new_clist


LOGICAL FUNCTION MovedTooMuch(Skin)
!
!   Scans the DisplaceList vector, containing the displacement vectors of
!   particles since the last time the list was updated.  It looks for the
!   magnitude of the two largest displacements in the system.   These are
!   summed together. If the sum is larger than Skin, the function returns
!   .TRUE., otherwise .FALSE.
!
use Particles
use Neighbor_List
implicit none
double precision :: Skin
double precision :: Displ,Displ1,Displ2
integer :: i
Displ1 = 0.d0
Displ2 = 0.d0
do i=1,N
	
	
    Displ = sqrt( dot_product(DisplaceList(:,i),DisplaceList(:,i)) )
    if (Displ >= Displ1) then       ! 1st position taken
        Displ2 = Displ1             ! push current 1st into 2nd place
        Displ1 = Displ              ! and put this one into current 1st
    else if (Displ >= Displ2) then  ! 2nd position taken
        Displ2 = Displ
    end if
end do
MovedTooMuch = ( Displ1 + Displ2 > Skin )
END FUNCTION MovedTooMuch


SUBROUTINE Print_Statistics
!
!   Print the mean value, averaged during the run, of the statistical
!   quantities which have been accumulated
!
USE Particles
USE Simulation_Control
USE Statistics

IMPLICIT NONE

double precision :: test

IF ( iPrint <= 0 ) RETURN

CALL Compute_Temperature(ene_kin_aver,temperature)
test = SUM( ene_pot )
ene_pot_aver = SUM( ene_pot ) / N
ene_pot_aver_exbond = SUM( ene_pot_exbond ) / N
ene_tot_aver = ene_kin_aver + ene_pot_aver
pressure = density*temperature + virial/volume/DIM

IF (lpps) THEN
    PRINT '(1x,i10,5f14.6)',step,temperature, &
                       ene_kin_aver,ene_pot_aver,ene_tot_aver,pressure, tension
END IF

IF (lppp) THEN
    temperature_sum = temperature_sum + temperature
    ene_kin_sum     = ene_kin_sum     + ene_kin_aver
    ene_pot_sum     = ene_pot_sum     + ene_pot_aver
    pressure_sum    = pressure_sum    + pressure
    IF ( MOD(step,iPrint).eq.0 ) THEN
        PRINT '(1x,i10,5f14.6)',step, &
                            temperature_sum / iPrint , &
                            ene_kin_sum / iPrint , &
                            ene_pot_sum / iPrint , &
                            ( ene_kin_sum + ene_pot_sum ) / iPrint , &
                            pressure_sum / iPrint
        temperature_sum = 0
        ene_kin_sum = 0
        ene_pot_sum = 0
        pressure_sum = 0
    END IF
END IF
END SUBROUTINE Print_Statistics


SUBROUTINE data_output
USE Particles
USE Simulation_Control
USE Statistics
IMPLICIT NONE
INTEGER :: nb1
INTEGER :: MSD_np, MSD_nt
DOUBLE PRECISION :: MSD_particle, MSD_tether
DOUBLE PRECISION, DIMENSION(3) :: Sij
INTEGER :: i, j, k


! Save
do i = 1,N
    if (ptype(i) == 3) then
        j = i
        exit
    end if
end do

k = MOD(step/iStaInt, iStaSave)
if (k == 0) k = iStaSave
sta_data(1:3,k) = pos(:,N)
sta_data(4:6,k) = pos(:,j)
sta_data(7,k) = ene_adh
sta_data(8,k) = isneighbor

IF ( k == iStaSave ) THEN
    OPEN(UNIT=2,FILE='data',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*8,STATUS='OLD')
    DO i = 1, iStaSave
        WRITE(2,REC=i+step/iStaInt-iStaSave) sta_data(:,i)
    END DO
    CLOSE(UNIT=2)
END IF

END SUBROUTINE data_output


SUBROUTINE Write_Sample
!
!  Writes the final sample to file unit 2
!
use Particles
use Simulation_Control
implicit none
integer :: i
character(LEN=10) :: strStep
CHARACTER(LEN=4) :: form
DOUBLE PRECISION, DIMENSION(3) :: Sij
WRITE(strStep, '(I10.10)') step

OPEN(unit=2, file='out2-'//strStep//'.xyz', status='replace')
WRITE(form,'(a,I1,a)') '(I', INT(LOG10(10.*N)), ')'
WRITE(2,form) N
WRITE(2,'(a,I9,a,3E23.15)') 'Atoms. Timestep: ', step, ', BoxSize: ', BoxSize
DO i = 1, N
    WRITE(2,'(I1,1X,3E23.15)') ptype(i), pos(:,i)
END DO

CLOSE(UNIT=2)

END SUBROUTINE Write_Sample

SUBROUTINE Write_Sample2
!
!  Writes the final sample to file unit 2
!
use Particles
use Simulation_Control
use Potential
implicit none
integer :: i
character(LEN=10) :: strStep
CHARACTER(LEN=20) :: form
DOUBLE PRECISION, DIMENSION(3) :: Sij
WRITE(strStep, '(I10.10)') step

OPEN(unit=2, file='out-'//strStep//'.xyz', status='replace')
WRITE(form,'(a,I1,a)') '(I', INT(LOG10(10.*N)), ')'
WRITE(2,form) N
WRITE(2,'(a,I9,a,3E23.15)') 'Atoms. Timestep: ', step, ', BoxSize: ', BoxSize
DO i = 1, N
    WRITE(2,'(I1,1X,3E23.15)') ptype(i), pos(:,i)
END DO
CLOSE(UNIT=2)
!
OPEN(UNIT=2, FILE='uvw-'//strStep, STATUS='REPLACE')
DO i = 1, N
    WRITE(2,'(3E23.15)') vel(:,i)
END DO
DO i = 1, N
    WRITE(2,'(3E23.15)') acc(:,i)
END DO
CLOSE(UNIT=2)

OPEN(unit=2, file='aout-'//strStep//'.xyz', status='replace')
WRITE(form,'(a,I1,a)') '(I', INT(LOG10(10.*N)), ')'
WRITE(2,form) memnum
WRITE(2,'(a,I9,a,3E23.15)') 'Atoms. Timestep: ', step, ', BoxSize: ', BoxSize
DO i = 1, memnum
    WRITE(2,'(3E23.15)') axes(:,i)	
END DO
CLOSE(UNIT=2)


OPEN(UNIT=2, FILE='auvw-'//strStep, STATUS='REPLACE')
DO i = 1, memnum
    WRITE(2,'(3E23.15)') angvel(:,i)
END DO
DO i = 1, memnum
    WRITE(2,'(3E23.15)') angacc(:,i)
END DO
CLOSE(UNIT=2)

open(unit=2,file='PA',status='replace')
write(2,'(3I7,3E23.15)') N, maxatomtype, BListLength, BoxSize
do i=1,BListLength
    write(2,'(3I7)'), BondList(:,i)
end do
close(unit=2)
END SUBROUTINE Write_Sample2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TERMINATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    ROUTINES


SUBROUTINE Terminate
!
!   Termination procedure
!
use Particles
USE Potential
use Neighbor_List
USE Statistics
implicit none

deallocate( ptype )
deallocate( deltar )
deallocate( pos )
deallocate( vel )
deallocate( acc )
deallocate( axes )
deallocate( angvel )
deallocate( angacc )
deallocate( ene_pot )
DEALLOCATE( ene_pot_exbond )
deallocate( ene_kin )
!
deallocate( BondList )
deallocate( BondMarker1 )
deallocate( BondMarker2 )
!
DEALLOCATE( DBType )
DEALLOCATE( ndbond )
DEALLOCATE( DBMap1 )
DEALLOCATE( DBMap2 )
DEALLOCATE( DBNeig1 )
DEALLOCATE( DBNeig2 )
!
deallocate( List )
deallocate( Marker1 )
deallocate( Marker2 )
deallocate( Advance )
deallocate( DisplaceList )
!
deallocate( CellList )
!
DEALLOCATE( epsil )
DEALLOCATE( Rcutoffsq )
DEALLOCATE( phicutoff )
DEALLOCATE( PotShift )
!
DEALLOCATE( sta_data )
END SUBROUTINE Terminate

SUBROUTINE Crossmuls(vec1,vec2,crossmul1)
implicit none
real:: vec1(3),vec2(3)
real:: crossmul1(3)
crossmul1(1)=vec1(2)*vec2(3)-vec2(2)*vec1(3)
crossmul1(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
crossmul1(3)=vec1(1)*vec2(2)-vec2(1)*vec1(2)
return
END SUBROUTINE Crossmuls

