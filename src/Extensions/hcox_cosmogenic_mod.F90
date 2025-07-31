!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_cosmogenic_mod.F90
!
! !DESCRIPTION: Defines the HEMCO extension for cosmogenic production
!\\
! More detail
!\\
! !INTERFACE:
!
MODULE HCOX_Cosmogenic_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State   ! Derived type for HEMCO state
  USE HCOX_State_Mod, ONLY : Ext_State   ! Derived type for External state

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HcoX_Cosmogenic_Run
  PUBLIC  :: HcoX_Cosmogenic_Init
  PUBLIC  :: HcoX_Cosmogenic_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
!  PRIVATE :: Init_7Be_Emissions
!
! !REMARKS:
!  References:
!  ============================================================================
!
! !REVISION HISTORY:
!  19 Jul 2025 - Lee T. Murray - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst

   ! Emissions indices etc.
   INTEGER               :: Instance
   INTEGER               :: ExtNr         ! Main Extension number
   INTEGER               :: IDT14CO       

   ! For tracking emissions
   REAL(hp), POINTER     :: Emiss14CO(:,:,:)
   
   TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER  :: AllInst => NULL()

  ! 14C cosmogenic production table
  REAL(hp), ALLOCATABLE :: C14Prod_Table(:,:,:)
  REAL(hp), ALLOCATABLE :: C14Prod_phi(:)      
  REAL(hp), ALLOCATABLE :: C14Prod_gcr(:)      
  REAL(hp), ALLOCATABLE :: C14Prod_depth(:)  
!
! !DEFINED PARAMETERS:
!
  ! To convert kg to atoms
  REAL*8,  PARAMETER     :: XNUMOL_14CO  = ( 6.022140857d23 / 29.9981566086d-3 )
  ! To convert radians/degrees
  REAL*8,  PARAMETER     :: r2d = 180d0 / 3.14159d0
  REAL*8,  PARAMETER     :: d2r = 3.14159d0 / 180d0
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Cosmogenic_run
!
! !DESCRIPTION: Subroutine HcoX\_Cosmogenic\_Run computes cosmogenic production
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Cosmogenic_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,    ONLY : HCO_EvalFld
    USE HCO_FluxArr_Mod, ONLY : HCO_EmisAdd
    USE HCO_CLOCK_MOD,   ONLY : HcoClock_Get
    USE NETCDF    
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState    ! Options for extension
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  19 Jul 2025 - Lee T. Murray - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scalars
    INTEGER            :: I,        J,          L,          N
    INTEGER            :: HcoID,    DOY,     Year,      Month
    INTEGER            :: YearI
    INTEGER            :: VARID,   NCID,    NCID2,      NCID3
    REAL(hp)           :: DTSRCE
    REAL(hp)           :: ADD_14CO, C14_TMP, P_TMP, LAT_TMP
    REAL(hp)           :: DEPTH, GCR, PHI, POLE_LON, POLE_LAT
    REAL(hp)           :: Val2D(1,1), Val3D(1,1,1)
    CHARACTER(LEN=255) :: MSG, LOC
    
    ! Pointers
    TYPE(MyInst), POINTER :: Inst
    REAL(hp),     POINTER :: Arr2D(:,:  )
    REAL(hp),     POINTER :: Arr3D(:,:,:)

    !=======================================================================
    ! HCOX_Cosmogenic_RUN begins here!
    !=======================================================================
    LOC = 'HCOX_Cosmogenic_RUN (HCOX_COSMOGENIC_MOD.F90)'

    ! Return if extension not turned on
    IF ( ExtState%Cosmogenic <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Set error flag
    !ERR = .FALSE.

    ! Get instance
    Inst   => NULL()
    CALL InstGet ( ExtState%Cosmogenic, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find Cosmogenic instance Nr. ', ExtState%Cosmogenic
       CALL HCO_ERROR(MSG,RC)
       RETURN
    ENDIF

    ! Emission timestep [s]
    DTSRCE = HcoState%TS_EMIS

    ! Nullify
    Arr2D => NULL()
    Arr3D => NULL()

    !=================================================================
    ! Update Solar Modulation Potential and Magnetic Pole
    !=================================================================
      
    ! Get todays date
    CALL HcoClock_Get( HcoState%Clock, cDOY = DOY, cYYYY = Year, cMM = Month, RC=RC )

    ! Solar modulation potential from 
    ! References:                 
    ! [1] Usoskin, I.G., K. Alanko-Huotari, G.A. Kovaltsov, and K. Mursula, 
    ! Heliospheric modulation of cosmic rays: Monthly reconstruction for 1951-2004, 
    ! J. Geophys. Res., 110, A12108, 2005 (doi:10.1029/2005JA011250).      
    ! [2] Usoskin, I. G., G. A. Bazilevskaya, and G. A. Kovaltsov, 
    ! Solar modulation parameter for cosmic rays since 1936 reconstructed from ground-based neutron monitors
    ! and ionization chambers, J. Geophys. Res., 116, 2011 (doi:10.1029/2010JA016105).                
    ! [3] McCracken, K.G., and J. Beer, 
    ! Long-term changes in the cosmic ray intensity at Earth, 
    ! 1428-2005, J. Geophys. Res., 112, A10101, 2007 (doi:10.1029/2006JA012117).        	
    
    !call check( nf90_open( 'solar_mod_pot.1936-2022.nc', nf90_nowrite, ncid2) )
    !call check( nf90_inq_varid( ncid2, "PHI", varid) )
    !call check( nf90_get_var( ncid2, varid, Val2D, (/Month,Year-1936+1/), (/1,1/) ) )
    !PHI = Val2D(1,1)
    !call check( nf90_close(ncid2) )
 
    IF ( YEAR < 1850 ) THEN
      YEARI = 1850
    ELSE
      YEARI = YEAR
    ENDIF

    call check( nf90_open( 'solar_mod_pot.1850-2023.nc', nf90_nowrite, ncid2) )
    call check( nf90_inq_varid( ncid2, "PHI", varid) )
    call check( nf90_get_var( ncid2, varid, Val2D, (/Month,YearI-1850+1/), (/1,1/) ) )
    PHI = Val2D(1,1)
    call check( nf90_close(ncid2) )

    IF ( YEAR < 1700 ) THEN
       YEARI = 1700
    ELSE
       YEARI = YEAR
    ENDIF

    ! Get magnetic north pole reconstruction/forecast from NOAA
    ! https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy
    call check( nf90_open( 'mag_pole.1700-2024.nc', nf90_nowrite, ncid3) )
    call check( nf90_inq_varid( ncid3, "MAG_POLE", varid) )
    call check( nf90_get_var( ncid3, varid, Val3D, (/YearI-1700+1,DOY,1/), (/1,1,1/) ) )
    POLE_LON = Val3D(1,1,1)
    call check( nf90_get_var( ncid3, varid, Val3D, (/YearI-1700+1,DOY,2/), (/1,1,1/) ) )
    POLE_LAT = Val3D(1,1,1)
    call check( nf90_close(ncid3) )
    !print*,'Magnetic Pole and PHI:',Year,DOY,POLE_LON,POLE_LAT,PHI

    !=======================================================================
    ! Compute 14CO emissions [kg/m2/s]
    !=======================================================================
    IF ( Inst%IDT14CO > 0 ) THEN
!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I, J, L, LAT_TMP, P_TMP, C14_TMP, ADD_14CO         ) &
!$OMP PRIVATE( DEPTH, GCR                                         ) &
!$OMP SCHEDULE( DYNAMIC )
       DO L = 1, HcoState%Nz
       DO J = 1, HcoState%Ny
       DO I = 1, HcoState%Nx

          ! Reset
          ADD_14CO = 0d0

          ! Determine midpoint pressure [Pa]
          P_TMP = ( HcoState%Grid%PEDGE%Val(I,J,L) + &
                    HcoState%Grid%PEDGE%Val(I,J,L+1)    ) / 2.0_hp 

          ! Determine atmospheric depth from TOA in [g cm-2]
          DEPTH = 1d-1 * ( P_TMP / HcoState%Phys%g0 )
          
          ! Determine geomagnetic cutoff rigidity in [GV]
          GCR   = GEOMAGNETIC_CUTOFF_RIGIDITY( HcoState%Grid%XMID%Val(I,J), &
                                               HcoState%Grid%YMID%Val(I,J), &
                                               POLE_LON, POLE_LAT )

          ! Get 14CO production from Poluianov et al., 2016 [atoms/g/s]
          ! Assuming 93% yield of 14C + O2 -> 14CO [Mak et al., 1994]
          C14_TMP = 0.93 * get_14C_prod( PHI, GCR, DEPTH )

          ! Calculate mass in layer 
          P_TMP = ( HcoState%Grid%PEDGE%Val(I,J,L) - &
                    HcoState%Grid%PEDGE%Val(I,J,L+1) ) / HcoState%Phys%g0 ! kg/m2
          P_TMP = P_TMP * HcoState%Grid%AREA_M2%Val(I,J) * 1e3            ! kg/m2 -> g
          
          ! Convert production to atoms/s
          C14_TMP = C14_TMP * P_TMP
          
          ! ADD_C14 = [atoms/s] / [atom/kg] / [m2] = 14CO emissions [kg/m2/s]
          ADD_14CO  = ( C14_TMP / XNUMOL_14CO  ) / HcoState%Grid%AREA_M2%Val(I,J)

          ! Save emissions into an array for use below
          Inst%Emiss14CO(I,J,L) = ADD_14CO
          
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO

       !------------------------------------------------------------------------
       ! Add 14CO emissions to HEMCO data structure & diagnostics
       !------------------------------------------------------------------------

       ! Add emissions
       IF ( Inst%IDT14CO > 0 ) THEN
          Arr3D => Inst%Emiss14CO(:,:,:)
          CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDT14CO, &
                            RC,       ExtNr=Inst%ExtNr )
          Arr3D => NULL()
          IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( &
                             'HCO_EmisAdd error: Emiss14CO', RC )
             RETURN
          ENDIF
       ENDIF

    ENDIF !IDT14CO > 0

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Nullify pointers
    Inst    => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Cosmogenic_Run
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Cosmogenic_Init
!
! !DESCRIPTION: Subroutine HcoX\_Cosmogenic\_Init initializes the HEMCO
! GC\_Rn-Pb-Be extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Cosmogenic_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
    USE NETCDF
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! Hemco state
    INTEGER,          INTENT(INOUT)  :: RC

! !REVISION HISTORY:
!  07 Jul 2014 - R. Yantosca - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, nSpc, ExtNr
    INTEGER                        :: nDepth,  nGCR,  nPhi,  ncid,  varid
    CHARACTER(LEN=255)             :: MSG, LOC

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    ! Pointers
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_Cosmogenic_INIT begins here!
    !=======================================================================
    LOC = 'HCOX_COSMOGENIC_INIT (HCOX_COSMOGENIC.F90)'

    ! Get the main extension number
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%Cosmogenic, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( 'Cannot create Cosmogenic instance', RC )
       RETURN
    ENDIF
    ! Also fill the extension numbers in the Instance object
    Inst%ExtNr      = ExtNr
    
    ! Set HEMCO species IDs
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'Could not set HEMCO species IDs', RC )
       RETURN
    ENDIF

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN

       ! Write the name of the extension regardless of the verbose setting
       msg = 'Using HEMCO extension: Cosmogenic (radionuclide emissions)'
       IF ( HCO_IsVerb( HcoState%Config%Err ) ) THEN
          CALL HCO_Msg( HcoState%Config%Err, sep1='-' ) ! with separator
       ELSE
          CALL HCO_Msg( msg, verb=.TRUE.              ) ! w/o separator
       ENDIF

       ! Write all other messages as debug printout only
       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
    ENDIF

    ! Set up tracer and HEMCO indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'C14O', '14CO' )
             Inst%IDT14CO     = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! ERROR: No tracer defined
    IF ( Inst%IDT14CO <= 0 ) THEN       
       CALL HCO_ERROR( 'Cannot use Cosmogenic extension: no valid species!', RC )
    ENDIF

    ! Activate met fields required by this extension
    ExtState%FRCLND%DoUse  = .TRUE.
    ExtState%T2M%DoUse     = .TRUE.
    ExtState%AIR%DoUse     = .TRUE.
    ExtState%TropLev%DoUse = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    IF ( Inst%IDT14CO > 0 ) THEN

       ! Open the file. 
       call check( nf90_open( 'C14_production.nc', nf90_nowrite, ncid) ) 
       
       ! Get the dimensions
       call check( nf90_inq_dimid( ncid, "phi", varid ) )
       call check( nf90_inquire_dimension(ncid, varid, len = nPhi ) )
       call check( nf90_inq_dimid( ncid, "Pc", varid ) )
       call check( nf90_inquire_dimension(ncid, varid, len = nGCR ) )
       call check( nf90_inq_dimid( ncid, "Z", varid ) )
       call check( nf90_inquire_dimension(ncid, varid, len = nDepth ) )
       
       ! Array for Production Table
       ALLOCATE( C14Prod_Table( nPhi, nGCR, nDepth ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate C14Prod_table', RC )
          RETURN
       ENDIF
       call check( nf90_inq_varid(ncid, "Q", varid) )
       call check( nf90_get_var(ncid, varid, C14Prod_Table ) )
       
       ! Array for Solar Modulation Potential
       ALLOCATE( C14Prod_phi( nPhi ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate C14Prod_phi', RC )
          RETURN
       ENDIF
       call check( nf90_inq_varid(ncid, "phi", varid) )
       call check( nf90_get_var(ncid, varid, C14Prod_phi ) )
       C14Prod_phi = C14Prod_phi * 1e3 ! Convert from GV to MV
       WRITE(6,*) nPhi, 'Solar Modulation Potential Values Found'
       WRITE(6,*) C14Prod_phi
       WRITE(6,*) ''
       
       ! Array for Geomagnetic Cutoff Rigidity
       ALLOCATE( C14Prod_gcr( nGCR ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate C14Prod_gcr', RC )
          RETURN
       ENDIF
       call check( nf90_inq_varid(ncid, "Pc", varid) )
       call check( nf90_get_var(ncid, varid, C14Prod_gcr ) )
       WRITE(6,*) nGCR, "Geomagnetic Cutoff Rigidity Values Found"
       WRITE(6,*) C14Prod_gcr
       WRITE(6,*) ''
       
       ! Array for Atmospheric Depth (Pressure)
       ALLOCATE( C14Prod_depth( nDepth ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate C14Prod_depth', RC )
          RETURN
       ENDIF
       call check( nf90_inq_varid(ncid, "Z", varid) )
       call check( nf90_get_var(ncid, varid, C14Prod_depth ) )
       WRITE(6,*) nDepth, 'Atmospheric Depth Values Found'
       WRITE(6,*) C14Prod_depth
       WRITE(6,*) ''
       
       ! Close the production table
       call check( nf90_close(ncid) )
       
       ALLOCATE( Inst%Emiss14CO( HcoState%Nx, HcoState%NY, HcoState%NZ ), STAT=RC )
       IF ( RC /= 0 ) THEN
          CALL HCO_ERROR ( 'Cannot allocate Emiss14CO', RC )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    ! Nullify pointers
    Inst    => NULL()

    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Cosmogenic_Init
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Cosmogenic_Final
!
! !DESCRIPTION: Subroutine HcoX\_Cosmogenic\_Final finalizes the HEMCO
!  extension for the GEOS-Chem Rn-Pb-Be specialty simulation.  All module
!  arrays will be deallocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Cosmogenic_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_GC_RNPBBE_FINAL begins here!
    !=======================================================================

    CALL InstRemove ( ExtState%Cosmogenic )

  END SUBROUTINE HCOX_Cosmogenic_Final
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       !---------------------------------------------------------------------
       ! Deallocate fields of Inst before popping Inst off the list
       ! in order to avoid memory leaks
       !---------------------------------------------------------------------

       IF ( ASSOCIATED( Inst%Emiss14CO ) ) THEN
          DEALLOCATE( Inst%Emiss14CO )
       ENDIF
       Inst%Emiss14CO => NULL()

       !---------------------------------------------------------------------
       ! Pop off instance from list
       !---------------------------------------------------------------------
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)

    ENDIF

    ! Free pointers before exiting
    PrevInst => NULL()
    Inst     => NULL()

  END SUBROUTINE InstRemove
!EOC
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function get_14C_prod( phi, gcr, depth )

    ! Arguments
    real*8 :: get_14C_prod ! Output, 14 C production rate in g cm-2 s-1
    real*8 :: phi          ! Input, Solar modulation parameter
    real*8 :: gcr          ! Input, geomagnetic cutoff rigidity
    real*8 :: depth        ! Input, atmospheric depth [ g cm-2 ]

    integer :: i, j, k

    ! What is the closest atmospheric depth?
    k = minloc( abs( C14Prod_depth - depth ), 1 )

    ! What is the closest geomagnetic cutoff rigidity?
    j = minloc(abs( C14Prod_gcr - gcr ), 1 )

    ! What is the closet solar modulation potential?
    i = minloc(abs( C14Prod_phi - phi ), 1 )

    ! Return the closest C-14 Production rate
    get_14C_prod = C14Prod_Table( i, j, k )

    ! Use smarter interpolater later

  end function get_14C_prod

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function geomagnetic_cutoff_rigidity( lon, lat, pole_lon, pole_lat )

    ! http://www.cosmicrays.org/muon-cutoff-rigidity.php

    ! Arguments
    real*8 :: geomagnetic_cutoff_rigidity ! Output [GV]
    real*8 :: lat, lon           ! Input, geograph coords [Degrees]
    real*8 :: pole_lat, pole_lon ! Input, magnetic pole location [Degrees]

    ! Internals
    real*8 :: gmlat

    ! Get geomagnetic latitude
    gmlat = geomag_lat( lon, lat, pole_lon, pole_lat )

    ! Use dipole approximation [Smart and Shea, 2005]
    geomagnetic_cutoff_rigidity = ( 14.5d0 * ( cos( gmlat * d2r ) )**4d0 )

  end function geomagnetic_cutoff_rigidity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
  function geomag_lat( lon, lat, pole_lon, pole_lat )
    ! Returns the geomagnetic latitude for a given lat, lon and year
    ! i.e., the latitude of the sphere with the central axis
    !       being the Earth's geomagnetic dipole
    
    ! Arguments
    real*8 :: geomag_lat ! Output, geomagnetic lat   [Degrees]
    real*8 :: lat, lon   ! Input, geographic coords  [Degrees]
    real*8 :: pole_lat, pole_lon ! Input, magnectic north pole [Degrees]
    
    ! Determine the Geomagnetic Latitude relative to mag. North [Radians]
    geomag_lat = asin( sin(lat*d2r) * sin(pole_lat*d2r) + &
        cos(lat*d2r) * cos(pole_lat*d2r) * cos(lon*d2r-pole_lon*d2r) )
    
    ! Convert to Degrees
    geomag_lat = geomag_lat * r2d
    
  end function geomag_lat
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
  
END MODULE HCOX_Cosmogenic_Mod
