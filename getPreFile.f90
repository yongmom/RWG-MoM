!*********************** output the project file **********************!
subroutine getPreFile( )
    use MoM_VARIABLES
    use mod_Mesh
    use MoM_POST
    use mod_Scatter
    use mod_WirePort
    use mod_LumpPort
    use mod_Waveport
    implicit none
    integer m, n, I, J
    real(DBL),allocatable :: Rpack(:)
    integer,allocatable :: Ipack(:)

    ! open the project file
    open( unit=10, file=trim(ProjFile)//".tgm", Form='binary' )

    !#------- This is a project file for RWG-MoM -------#!
    write(10) FreqNum
    write(10) SwpType
    write(10) SwpFreq
    write(10) NGauss
    write(10) SymSet(:)
    write(10) nmed
    write(10) ScatNum
    write(10) WgapNum
    write(10) LumpNum
    write(10) PortNum
    write(10) Far_Num
    write(10) NearNum
    write(10) isCurrt
    write(10) VoltNum
    write(10) VoltGrp
    write(10) VertNum
    write(10) WireNum
    write(10) TriaNum
    write(10) NDomain
    write(10) NBound
    write(10) NmomTria
    write(10) NmomWire

    write(10) Freq0
    write(10) FreqMin
    write(10) FreqMax
    write(10) eps(:)
    write(10) miu(:)
    write(10) sgm(:)

    !#----------- Information of excitations -----------#!
    m = 6 * ScatNum + 9 * LumpNum + 11 * PortNum
    n = 3 * ScatNum + 2 * WGapNum + PortNum
    allocate( Rpack( m ) )
    allocate( Ipack( n ) )
    
    do m = 1, ScatNum
        I = 3 * m - 2
        J = 6 * m - 5
        Ipack( I:I+1 ) = InciNum(:,m)
        Ipack( I + 2 ) = LR_rot(m)
        Rpack( J:J+5 ) = PW(:,m)
    end do

    I = 3 * ScatNum + 1
    do m = 1, WgapNum
        Ipack( I:I+1 ) = WireGap(:,m)
        I = I + 2
    end do

    I = 6 * ScatNum + 1
    do m = 1, LumpNum
        Rpack( I:I+2 ) = LumpV(:,m)
        I = I + 3
        Rpack( I:I+2 ) = LumpE(:,m)
        I = I + 3
        Rpack( I:I+2 ) = LumpH(:,m)
        I = I + 3
    end do

    I = 3 * ScatNum + 2 * WgapNum
    J = 6 * ScatNum + 9 * LumpNum
    do m = 1, PortNum
        Ipack( I + m ) = PortKind(m)
        Rpack( J + 1 ) = width(m)
        Rpack( J + 2 ) = heigh(m)
        Rpack( J+3:J+5 ) = origin(:,m)
        Rpack( J+6:J+8 ) = a_in(:,m)
        Rpack( J+9:J+11) = e_in(:,m)
        J = J + 11
    end do
    ! write  the file
    write(10) Ipack(:)
    write(10) Rpack(:)

    !#-------- Parameters for the post process ---------#!
    write(10) NumPh(:),NumTh(:)
    write(10) MinPh(:),MaxPh(:)
    write(10) MinTh(:),MaxTh(:)

    do m = 1, VoltGrp
        write(10) Excit( :, m )
    end do
    !#------------ Model mesh of the target ------------#!
    write(10) Vertex_r(:,:)

    if( WireNum > 0 ) then
        write(10) WireVert(:,:)
		write(10) DomaWire(:,:)
        write(10) radius(:)
        write(10) DomaMedi(:)
    end if

    if( TriaNum > 0 ) then
        write(10) TriaVert(:,1:TriaNum)
		write(10) BoudTria(:,1:NBound )
        write(10) BoudMedi(:,1:NBound)
    end if

    close( 10 )
    !---------------------------- All done ----------------------------!

end subroutine getPreFile
!******************************* The End ******************************!