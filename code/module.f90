MODULE Particles
!
!   Contains all the data structures containing informations
!   about atom properties
!
    INTEGER, PARAMETER :: DIM = 3   ! change DIM to switch between 2D and 3D !
    LOGICAL :: VelAcc = .FALSE.     ! velocities and accelerations in file?
    INTEGER :: N, partnum, memnum, nanonum
    DOUBLE PRECISION :: density
    DOUBLE PRECISION, DIMENSION(3) :: BoxSize
!       the following arrays are allocated at run time
!       when the number of atoms N is known
    integer :: maxatomtype
    integer, dimension(:), allocatable :: ptype                 ! particle type
    double precision, dimension(:,:), allocatable :: deltar     ! displacement
    double precision, dimension(:,:), allocatable :: pos        ! positions
    double precision, dimension(:,:), allocatable :: vel        ! velocities
    double precision, dimension(:,:), allocatable :: acc        ! accelerations
    double precision, dimension(:), allocatable :: ene_pot      ! potential energies
	double precision  :: ene_adh      ! potential energies
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ene_pot_exbond ! potential energies exclude bond potential
    double precision, dimension(:), allocatable :: ene_kin      ! kinetic energies
	
	double precision, dimension(:,:), allocatable :: axes        ! positions
    double precision, dimension(:,:), allocatable :: angvel        ! velocities
    double precision, dimension(:,:), allocatable :: angacc        ! accelerations
	
	!integer, dimension(:), allocatable :: body        ! body
!       scalars derived from the above quantities
    double precision :: volume      ! these are constants
    double precision :: virial=0.d0      ! virial term, used to compute the presure
!       bond list
    double precision :: PS=0.d0
    integer, parameter :: MaxBondsPerAtom = 6
    integer, dimension(:), allocatable :: BondMarker1,BondMarker2
	integer, dimension(:,:), allocatable :: BondList
    integer :: BMaxListLength,BListLength
    ! Dynamic bonds
    INTEGER, DIMENSION(:), ALLOCATABLE :: DBType, ndbond
    INTEGER, DIMENSION(:), ALLOCATABLE :: DBMap1, DBMap2
    INTEGER :: DBMapLength
    INTEGER, DIMENSION(:), ALLOCATABLE :: DBNeig1, DBNeig2
    INTEGER :: DBNeigMaxLength, DBNeigLength    
	! rigid body
    ! INTEGER :: nbody
    ! INTEGER, DIMENSION(:), ALLOCATABLE :: body
    ! DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: masstotal
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: displace
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: inertial
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xcm
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vcm
     DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fcm
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: quat
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: torque
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: omega
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: angmom
    ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: conjqm
    ! DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rotation_matrix 
    !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dforce
	
    INTEGER :: G
	LOGICAL :: NSigmaT = .TRUE.
    double precision, parameter :: Q = 4500.d0
END MODULE Particles


module Simulation_Control
!
!   Contains the parameters controlling the simulation, which will be
!   supplied as input by the user.
!
    character*80     :: title   = 'BROWNIAN DYNAMICS'   ! Arbitrary 'title' string
    CHARACTER*20     :: argum
    integer          :: iPrint  = 1000             ! Number of time setps to print information
    integer          :: iSave   = 5000        ! Number of time steps to save sample
	integer          :: iSave2  = 100000
    double precision :: deltat  = 0.005d0                ! Time steps (reduced units)
    double precision :: gamma   = 1.d0
    double precision :: TRequested = 0.23d0             ! Desired temperature
    double precision :: Skin = 0.5d0                  ! extra range for neighbor list
    double precision :: ten0 = 0.d0                ! extra range for neighbor list
    logical          :: lpps = .FALSE.                  ! print per step
    logical          :: lppp = .TRUE.                   ! print per peroid
    ! global variables
    integer          :: step = 0
end module Simulation_Control


module Potential
    ! pairwise interactions (LJ and WCA)
    double precision, dimension(:,:), allocatable :: epsil
    double precision, dimension(:,:), allocatable :: Rcutoffsq
    double precision, dimension(:,:), allocatable :: phicutoff
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PotShift
    double precision, parameter :: U = 130.d0
    ! bond interactions (FENE)
    double precision, parameter :: FENE_K = 1.0d0
    double precision, parameter :: FENE_Rmax = 1.5d0
	double precision, parameter :: mu = 3.0d0, xi = 4.d0, theta0 = 0.d0, pi = 3.1415926d0
	
end module Potential


module Statistics
!
!   Contains statistical quantities accumulated during the run.
!   All quantities are resetted to zero.
!
    DOUBLE PRECISION :: ene_kin_aver, ene_pot_aver, ene_tot_aver, ene_pot_aver_exbond
    DOUBLE PRECISION :: temperature, pressure, tension
    double precision :: temperature_sum = 0.d0
    double precision :: ene_kin_sum     = 0.d0
    double precision :: ene_pot_sum     = 0.d0
    double precision :: pressure_sum    = 0.d0
    INTEGER :: iStaInt = 1
    INTEGER :: iStaSave = 100000
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sta_data
    INTEGER :: MSD_int = 100000
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MSD_ri0
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MSD_dr
	INTEGER :: isneighbor = 0
end module Statistics


module Neighbor_List
    integer, parameter :: MaxPairsPerAtom = 200
    integer, dimension(:), allocatable :: Advance, Marker1, Marker2, List
    integer :: MaxListLength, ListLength
    double precision, dimension(:,:), allocatable :: DisplaceList
!   cell list
    integer :: nPart
    integer, dimension(3) :: iCell
    integer, dimension(200,200,200) :: Head
    integer, dimension(:), allocatable :: CellList
end module Neighbor_List


MODULE ran_mod
! module contains two FUNCTIONs
! ran RETURN a uniform random number between 0-1
! normal RETURN a normal distribution
CONTAINS
    FUNCTION ran()  ! RETURN random number between 0-1
        IMPLICIT NONE
        INTEGER,SAVE :: flag = 0
        DOUBLE PRECISION :: ran
        
        IF (flag == 0) THEN
            CALL init_random_seed(); flag = 1
        END IF
        CALL random_number(ran) ! built in fortran 90 random number FUNCTION
        RETURN
    END FUNCTION ran
    
    FUNCTION normal() ! RETURNs a normal distribution
        IMPLICIT NONE
        INTEGER, SAVE :: flag = 0
        DOUBLE PRECISION :: normal
        DOUBLE PRECISION :: fac, gsave, rsq, r1, r2
        SAVE gsave
        
        IF (flag.eq.0) THEN
            rsq = 2.0d0
            DO WHILE (rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do
                r1 = 2.0d0 * ran() - 1.0d0
                r2 = 2.0d0 * ran() - 1.0d0
                rsq = r1 * r1 + r2 * r2
            END DO
            fac = sqrt(-2.0d0*log(rsq)/rsq)
            gsave = r1 * fac
            normal = r2 * fac
            flag = 1
        ELSE
            normal = gsave
            flag = 0
        ENDIF
        
        RETURN
    END FUNCTION normal
    
    SUBROUTINE init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        CALL SYSTEM_CLOCK(COUNT = clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
    END SUBROUTINE    
END MODULE ran_mod
