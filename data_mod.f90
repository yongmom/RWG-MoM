!----------------------------------------------------------------------!
!    This file defines the common variables which will be used in many !
! subroutines.                                                         !
!               Designed by WangYong in November 1th 2012              !
!                     Update by Wang Yong 24.02.2016                   !
!----------------------------------------------------------------------!
module MoM_VARIABLES
    integer,parameter :: DBL = kind(1.0D0)
    real(DBL),parameter :: pi = 3.1415926535898D0  ! parameter pi
    complex(DBL),parameter :: cj = ( 0.0D0, 1.0D0 )

    !------------- Set the parameters of Gauss integration ------------!
    ! the number of Gauss sampling point ( default Value: 4 )
    integer :: NGauss = 4
    !------------------------------------------------------------------!
    real(DBL) Freq0             ! the frequency to be analysised
    real(DBL) FreqMin, FreqMax  ! start and end of sweeping frequency
    integer FreqNum             ! the number of sweep frequency
    integer SwpType             ! the sampling method of sweeping
    integer SymSet(3)           ! setting for symmetry
    integer NunJ, NunM, NUN     ! number of unknowmn

    ! the name and ProjPath of the model file
    character(len=500) ProjPath, ProjFile, meshFile, PostFile
    integer nmed
    complex(DBL),allocatable :: eps(:), miu(:)
    real(DBL),allocatable :: sgm(:)

end module MoM_VARIABLES
!****************************** The End *******************************!

!****************** the parameters of the model mesh ******************!
module mod_Mesh
    use MoM_VARIABLES, only : DBL
    ! the number of vertices, wires triangles and junction
    integer VertNum, WireNum, TriaNum, NmomTria, NmomWire
    integer NDomain, NBound      ! number of wire domain and boundary
    integer,allocatable :: DomaMedi(:)    ! medium of each domain
    integer,allocatable :: BoudMedi(:,:)  ! medium across boundary
    integer,allocatable :: DomaWire(:,:)  ! medium wire located in
    integer,allocatable :: BoudTria(:,:)  ! triange in each boundary
    integer,allocatable :: WireVert(:,:)  ! index of wire vertex
    integer,allocatable :: TriaVert(:,:)  ! index of triangle vertex
    ! the coordinate of triangle vertices
    real(DBL),allocatable :: Vertex_r(:,:)
    real(DBL),allocatable :: radius(:)    ! radius of each wire
    real(DBL),allocatable :: TriaArea(:)  ! area of each triangle
end module mod_Mesh
!****************************** The End *******************************!

!----------------------------------------------------------------------!
!   This module defines the variables which will be used in scattering !
! problems.                                                            !
module mod_Scatter
    use MoM_VARIABLES, only : DBL
    !---------------- the parameters of the plane wave ----------------!
    integer ScatNum              ! number of plane wave
    integer,allocatable :: InciNum(:,:)
    ! min, max, polarization, axial ratio of incident
    real(DBL),allocatable :: PW(:,:)
    ! polarization( 1:right-hand; -1:left-hand, 0:linear )
    integer,allocatable :: LR_rot(:)
    !------------------------------------------------------------------!
end module mod_Scatter
!****************************** The End *******************************!

!----------------------------------------------------------------------!
!   This module defines the variables which will be used in WirePort   !
! problems.                                                            !
module mod_WirePort
    use MoM_VARIABLES, only : DBL
    !---------------- the parameters of the excitation ----------------!
    integer WGapNum              ! number of wire gap
    integer,allocatable :: WireGap(:,:)
    real(DBL),allocatable :: GapLoc(:,:)
    real(DBL),allocatable :: GapEin(:,:)
    !------------------------------------------------------------------!
end module mod_WirePort
!****************************** The End *******************************!

!----------------------------------------------------------------------!
!   This module defines the variables which will be used in LumpPort   !
! problems.                                                            !
module mod_LumpPort
    use MoM_VARIABLES, only : DBL
    !---------------- the parameters of the excitation ----------------!
    integer LumpNum              ! number of lump port
    real(DBL),allocatable :: LumpV(:,:)
    real(DBL),allocatable :: LumpE(:,:)
    real(DBL),allocatable :: LumpH(:,:)
    !------------------------------------------------------------------!
end module mod_LumpPort
!****************************** The End *******************************!

!----------------------------------------------------------------------!
!     This module defines the variables which will be used in Waveport !
! problems.                                                            !
module mod_Waveport
    use MoM_VARIABLES, only : DBL
    !---------------- the parameters of the  excitation ---------------!
    ! number of the waveport. if not given, the default value is 0
    integer PortNum
    ! the kinds of waveports : 
    ! 1 -- Rectangle ; 2 -- Circle ; 3 -- Coaxial cable
    integer,allocatable :: PortKind(:)
    integer :: ModeNum(3) = (/5, 5, 1/)
    ! the width and height of each port
    real(DBL),allocatable :: width(:), heigh(:)
    ! the relative origin of each port
    real(DBL),allocatable :: origin(:,:)
    ! the directions of the transmission associated with each port
    real(DBL),allocatable :: a_in(:,:)
    ! the directions of electric field associated with each port
    real(DBL),allocatable :: e_in(:,:)
    ! the directions of magnetic field associated with each port
    real(DBL),allocatable :: h_in(:,:)
    !------------------------------------------------------------------!
end module mod_Waveport
!****************************** The End *******************************!

!---------------- Parameters of post-processing setting ---------------!
module MoM_POST
    use MoM_VARIABLES, only : DBL
    integer VoltNum, VoltGrp
    integer SwpFreq
    ! set the far field plane
    integer Far_Num
    integer NearNum
    integer isCurrt
    integer,allocatable :: NumPh(:),NumTh(:)
    real(DBL),allocatable :: MinPh(:), MaxPh(:)
    real(DBL),allocatable :: MinTh(:), MaxTh(:)
    complex(DBL),allocatable :: Excit(:,:)

end module MoM_POST
!****************************** The End *******************************!