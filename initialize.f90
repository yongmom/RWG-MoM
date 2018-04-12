!----------------------------------------------------------------------!
!This source file includes the subroutines which are used to initialize!
!data for the program.                                                 !
!              Designed by Wang Yong in November 1th 2012              !
!----------------------------------------------------------------------!

!*********************** initial the parameters ***********************!
subroutine readConf( L )
    use mod_Mesh
    use MoM_POST
    use MoM_VARIABLES
    use mod_Scatter
    use mod_WirePort
    use mod_LumpPort
    use mod_Waveport
    implicit none
    integer L
    character String*500, Intro*10, Values*490
    logical alive, param(5)
    integer ios, I, J
    real(DBL) mag, pha, epsr, tanE, miur, tanU

    !--------------Check whether the configure file exist--------------!
    inquire( file = trim(ProjFile), exist = alive )
    if( .NOT. alive ) then
        write(*,*)"->Can't find the configure file, Please check it !"
        write(L,*)"Can't find the configure file, Please check it !"
        stop
    end if

    !------------------- initialize the parameters --------------------!
    FreqNum = 0
    ScatNum = 0
    WGapNum = 0
    LumpNum = 0
    PortNum = 0
    SwpFreq = 0
    Far_Num = 0
    NearNum = 0
    isCurrt = 0
    !-------------------- Read the configure file ---------------------!
    open( unit = 10, file = trim(ProjFile) )
    param(:) = .FALSE.
    do
        read( 10,"(A500)",IOSTAT = ios ) String
        ! judge whether it has reached the end of the file
        if( ios /= 0 ) exit

        Intro(:) = String(1:10)
        Values(:) = String(11:490)
        ! set the number of frequencies
        if(Intro == "Set.FMesh:") then
            meshFile = trim(adjustl(Values))
            param( 1 ) = .TRUE.
        
        ! set the analysis frequency
        else if(Intro == "Set.Freq0:") then
            read( Values, * ) Freq0
            param( 2 ) = .TRUE.

        ! set for sweeping the frequency
        else if(Intro == "Set.Sweep:") then
            read( Values, * ) FreqMin, FreqMax, FreqNum, SwpType

        ! whether calculate field of each frequency
        else if( Intro == "Set.SavFr:" ) then
            read( Values, * ) SwpFreq

        ! set the number of Gauss point
        else if( Intro == "Set.IntAc:" ) then
            read(Values,*) NGauss

        ! set the symmetric planes
        else if( Intro == "Set.Symme:" ) then
            read(Values,*) SymSet(:)
            param( 3 ) = .TRUE.

        ! set the material
        else if( Intro == "Set.Mediu:" ) then
            read(Values,*) nmed
            allocate( eps(nmed) )
            allocate( miu(nmed) )
            allocate( sgm(nmed) )
            do I = 1, nmed 
                read (10,*) J, epsr, tanE, miur, tanU
                if( J == 0 ) then
                    sgm(I) = 0.0D0
                else
                    sgm(I) = tanE
                    tanE = 0.0D0
                end if
                tanE = epsr * tanE
                tanU = miur * tanU
                eps(I) = epsr - cj * tanE
                miu(I) = miur - cj * tanU
            end do
            param( 4 ) = .TRUE.

        ! set the parameters of plane wave
        else if( Intro == "Set.PWave:" ) then
            read(Values,*) ScatNum
            allocate( InciNum(2,ScatNum) )
            allocate( PW(6,ScatNum) )
            allocate( LR_rot(ScatNum) )
            do I = 1, ScatNum
                read(10,*) InciNum(:,I), PW(:,I), LR_rot(I)
            end do

        ! set the delta voltage on wire gap
        else if( Intro == "Set.WGaps:" ) then
            read( Values,* ) WgapNum
            allocate( GapLoc( 3,WgapNum ) )
            allocate( GapEin( 3,WgapNum ) )
            allocate( WireGap(2,WgapNum ) )
            do I = 1, WgapNum
                read(10,*) GapLoc(:,I), GapEin(:,I)
            end do

        ! set the parameters of lump port
        else if( Intro == "Set.Lumps:" ) then
            read( Values,* ) LumpNum
            allocate( LumpV( 3,LumpNum ) )
            allocate( LumpE( 3,LumpNum ) )
            allocate( LumpH( 3,LumpNum ) )
            do I = 1, LumpNum
                read(10,*) LumpV(:,I)
                read(10,*) LumpE(:,I)
                read(10,*) LumpH(:,I)
            end do
            
        ! set the parameters of waveport
        else if( Intro == "Set.Ports:" ) then
            read( Values,* ) PortNum
            ! apply the memory space
            allocate( PortKind( PortNum ) )
            allocate( origin( 3,PortNum ) )
            allocate( width( PortNum ) )
            allocate( heigh( PortNum ) )
            allocate( a_in( 3,PortNum ) )
            allocate( e_in( 3,PortNum ) )
            allocate( h_in( 3,PortNum ) )
            ! input the information of wave ports
            do I = 1, PortNum
                read(10,*) PortKind(I)
                if( 1 == PortKind(I) ) then
                    ! rectangular waveport
                    read(10,*) width(I), heigh(I)
                    read(10,*) origin(:,I)
                    read(10,*) a_in(:,I)
                    read(10,*) e_in(:,I)
                    mag = sum( a_in(:,I) * a_in(:,I) )
                    a_in(:,I) = a_in(:,I) / sqrt(mag)
                    mag = sum( e_in(:,I) * e_in(:,I) )
                    e_in(:,I) = e_in(:,I) / sqrt(mag)
                    ! get the direction of magnetic field
                    call cross_mult( a_in(:,I), e_in(:,I), h_in(:,I) )
                    mag = heigh(I) / 2.0D0
                    pha = width(I) / 2.0D0
                    origin(:,I) = origin(:,I) - e_in(:,I) * mag
                    origin(:,I) = origin(:,I) + h_in(:,I) * pha
                else if( 2 == PortKind( I ) ) then
                    ! circle waveport
                    read(10,*) width(I)
                    read(10,*) origin(:,I)
                    read(10,*) a_in(:,I)
                    read(10,*) e_in(:,I)
                    heigh(I) = 0.0D0
                    mag = sum( a_in(:,i) * a_in(:,I) )
                    a_in(:,I) = a_in(:,I) / sqrt(mag)
                    mag = sum( e_in(:,i) * e_in(:,I) )
                    e_in(:,I) = e_in(:,I) / sqrt(mag)
                else
                    ! coaxial cable
                    read(10,*) width(I), heigh(I)
                    read(10,*) origin(:,I)
                    read(10,*) a_in(:,I)
                    mag = sum( a_in(:,i) * a_in(:,I) )
                    a_in(:,I) = a_in(:,I) / sqrt(mag)
                    e_in(:,I) = 0.0D0
                end if
            end do

        ! Set the far field setting
        else if( Intro == "Set.NFarF:" ) then
            read( Values, * ) Far_Num
            allocate( NumPh( Far_Num ) )
            allocate( NumTh( Far_Num ) )
            allocate( MinPh( Far_Num ) )
            allocate( MaxPh( Far_Num ) )
            allocate( MinTh( Far_Num ) )
            allocate( MaxTh( Far_Num ) )
            ! read the far field setting
            do I = 1, Far_Num
                read(10,*) NumPh(I),NumTh(I),MinPh(I),MaxPh(I),MinTh(I),MaxTh(I)
            end do

        else if( Intro == "Set.NNear:" ) then
            read( Values, * ) NearNum
            
        else if( Intro == "Set.isCur:" ) then
            read( Values, * ) isCurrt
            
        ! set the voltage of excitations
        else if( Intro == "Set.Excit:" ) then
            read(Values,*) VoltNum, VoltGrp
            allocate( Excit( VoltNum, VoltGrp ) )
            do J = 1, VoltGrp
                do I = 1, VoltNum
                    read(10,*) mag, pha
                    pha = pha * pi/180.0D0
                    Excit(I,J) = mag * exp( cj*pha )
                end do
            end do
            param( 5 ) = .TRUE.

        end if
    end do
    close( unit=10, STATUS='KEEP' )
    !-------------------- Read configure file end ---------------------!

    write(*,*) "  The project is : ", trim(ProjFile)
    write(L,*) "The project is : ", trim(ProjFile)

    !------------------- Check the input information ------------------!
    if( 0 == (ScatNum + WgapNum + LumpNum + PortNum) ) then
        write(*,*) " Warn : You did not set any excitation ! "
        write(L,*) "Warn : You did not set any excitation ! "
    end if

    if( .NOT. param(1) ) then
        write(*,*) " Error : Please set the mesh file !"
        write(L,*) "Error : Please set the mesh file !"
        stop
    end if
    if( .NOT. param(2) ) then
        write(*,*) " Error : Please set the frequency !"
        write(L,*) "Error : Please set the frequency !"
        stop
    end if
    if( .NOT. param(3) ) then
        write(*,*) " Error : Please set the symmetry !"
        write(L,*) "Error : Please set the symmetry !"
        stop
    end if
    if( .NOT. param(4) ) then
        write(*,*) " Error : Please set the material !"
        write(L,*) "Error : Please set the material !"
        stop
    end if
    if( .NOT. param(5) ) then
        write(*,*) " Warn : You didn't set the voltage !"
        write(L,*) "Warn : You didn't set the voltage !"
    end if
    !------------------------------------------------------------------!

    !------------ print the information of the excitations ------------!
    if( ScatNum > 0 ) then
        ! print the incident field information
        I = sum( InciNum(1,:) * InciNum(2,:) )
        write(*,*) "->The number of plane waves is  :",I
        write(*,*) "---------------------------------------------------"

        write(L,*) "The number of plane waves is  :",I
        write(L,*) "---------------------------------------------------"
    end if
    
    if( WgapNum > 0 ) then
        ! print the incident field information
        write(*,*) "->The number of delta voltage is :", WgapNum
        write(*,*) "---------------------------------------------------"
        write(L,*) "The number of delta voltage is :", WgapNum
        write(L,*) "---------------------------------------------------"
    end if

    if( LumpNum > 0 ) then
        ! print the incident field information
        write(*,*) "->The number of lump ports is :", LumpNum
        write(*,*) "---------------------------------------------------"
        write(L,*) "The number of lump ports is :", LumpNum
        write(L,*) "---------------------------------------------------"
    end if

    if( PortNum > 0 ) then
        ! print the incident field information
        write(*,*) "->The number of wave ports is :", PortNum
        write(*,*) "---------------------------------------------------"
        write(L,*) "The number of wave ports is :", PortNum
        write(L,*) "---------------------------------------------------"
    end if
    !------------------------------------------------------------------!
    I = sum( InciNum(1,:) * InciNum(2,:) )
    I = I + WgapNum + LumpNum + PortNum

    if( I /= VoltNum ) then
        write(*,*) "Error: The number of excitation is not matched !"
        stop
    end if

    ! modify the file name
    I = index( ProjFile, '.', .TRUE. )
    if( I > 1 ) ProjFile = ProjFile(1:I-1)

end subroutine readConf
!******************************* The End ******************************!

!*********************** initial the model data ***********************!
subroutine initModel( log )
    use mod_Mesh
    use MoM_VARIABLES
    use mod_WirePort
    use mod_LumpPort
    use mod_Waveport
    implicit none
    integer log
    logical alive          ! the flag of whether the model file exist
    integer,allocatable :: L_TriaNum(:)
    integer,allocatable :: P_TriaNum(:)
    integer m, Idx, V(2)
    real(DBL) norm(3), cent(3), Lm, dist, angle

    !--scan the model file and statistic the number of the triangles --!
    write(*,*) "Now check the information of the model..."

    !--------------- Check whether the model file exist ---------------!
    inquire( file = trim(meshFile)//'.gdm', exist = alive )
    if( .NOT. alive ) then
        write(*,*) "->Error : There is no such model file !"
        write(log,*)"Error : There is no such model file !"
        stop
    end if

    !---------------------- Read the model file -----------------------!
    call readMesh( log )
    ! initialize the mom triangle and wire
    NmomTria = 0
    NmomWire = 0
    NunJ  = 0
    NunM  = 0
    if( TriaNum > 0 ) call define_RWG( )   ! define RWG on triangle
    if( WireNum > 0 ) call define_Wire( )  ! define wire basis function
    ! get the overall unknowns of the matrix equation
    if( PortNum > 0 ) then
        m = sum( ModeNum( PortKind(:) ) )
        NunM = NunM + m
    end if
    NUN  = NunJ + NunM

    !-------- display the number of the vertices and triangles --------!
    write(*,*)"->Number of the vertices     : ",VertNum
    write(*,*)"->Number of the segments     : ",WireNum
    write(*,*)"->Number of the triangles    : ",TriaNum
    write(*,*)"->Number of electric current : ",NunJ
    write(*,*)"->Number of magnetic current : ",NunM
    write(*,*)"->Number of MoM unknowns     : ",NUN

    write(log,*)"Number of the vertices     : ",VertNum
    write(log,*)"Number of the segments     : ",WireNum
    write(log,*)"Number of the triangles    : ",TriaNum
    write(log,*)"Number of electric current : ",NunJ
    write(log,*)"Number of magnetic current : ",NunM
    write(log,*)"Number of MoM unknowns     : ",NUN

    !----------------- find the wire with delta gap ------------------!
    if( WgapNum > 0 ) then
        do Idx = 1, WgapNum
            WireGap(:,Idx) = 0     ! initialize the wire gap
            do m = 1, WireNum
                V(:) = WireVert(:,m)
                norm(:) = Vertex_r(:,V(2)) - Vertex_r(:,V(1))
                ! get the length of wire m
                Lm = sum( norm * norm )
                Lm = sqrt( Lm )
                
                cent(:) = GapLoc(:,Idx) - Vertex_r(:,V(1))
                dist = sqrt( sum( cent * cent ) )
                cent(:) = Vertex_r(:,V(2)) - GapLoc(:,Idx)
                dist = dist + sqrt( sum( cent * cent ) )
                
                if( abs( dist - Lm ) < 1.0e-8 ) then
                    ! port Idx is in wire m
                    WireGap(1,Idx) = m
                    angle = sum( norm * GapEin(:,Idx) )
                    if( angle >= 0.0 ) then
                        WireGap(2,Idx) = 1
                    else
                        WireGap(2,Idx) = -1
                    end if
                    ! find port Idx, exit the cycle
                    exit
                end if
            end do
            
            if( 0 == WireGap(1,Idx) ) then
                write(*,*) "Error : Can't find wire port", Idx
                write(log,*) "Error : Can't find wire port", Idx
                stop
            end if
        end do
        ! release the space of wire port
        deallocate( GapLoc, GapEin )
    end if
    !------------------------------------------------------------------!
    
    !----------------- find the triangle on waveport ------------------!
    if( LumpNum > 0 ) then
        allocate( L_TriaNum( LumpNum ) )
        ! initialize the number of patches in each port
        L_TriaNum(:) = 0
        do m = 1, TriaNum
            ! get the center of patch m
            cent(:) = Vertex_r( :, TriaVert(1,m) )
            cent(:) = cent(:) + Vertex_r( :, TriaVert(2,m) )
            cent(:) = cent(:) + Vertex_r( :, TriaVert(3,m) )
            cent(:) = cent(:)/3.0D0
            ! check whether the patch is in port and in which port
            call IsInLump( cent, Idx )
            if( 0 /= Idx ) then
                L_TriaNum( Idx ) = L_TriaNum( Idx ) + 1
            end if
        end do
        ! check the setting of waveport
        do m = 1, LumpNum
            if( 0 == L_TriaNum(m) ) then
                write(*,*) "->Can't find triangle on lump port:", m
                write(log,*) "Can't find triangle on lump port:", m
                stop
            end if
        end do
        deallocate( L_TriaNum )
    end if
    !------------------------------------------------------------------!
    
    !----------------- find the triangle on waveport ------------------!
    if( PortNum > 0 ) then
        allocate( P_TriaNum( PortNum ) )
        ! initialize the number of patches in each port
        P_TriaNum(:) = 0
        do m = 1, TriaNum
            ! get the center of patch m
            cent(:) = Vertex_r( :, TriaVert(1,m) )
            cent(:) = cent(:) + Vertex_r( :, TriaVert(2,m) )
            cent(:) = cent(:) + Vertex_r( :, TriaVert(3,m) )
            cent(:) = cent(:)/3.0D0
            ! check whether the patch is in port and in which port
            call IsInPort( cent, Idx )
            if( 0 /= Idx ) then
                P_TriaNum( Idx ) = P_TriaNum( Idx ) + 1
            end if
        end do
        ! check the setting of waveport
        do m = 1, PortNum
            if( 0 == P_TriaNum(m) ) then
                write(*,*) "->Can't find triangle on waveport:", m
                write(log,*) "Can't find triangle on waveport:", m
               ! stop
            end if
        end do
        deallocate( P_TriaNum )
    end if
    !------------------------------------------------------------------!

end subroutine initModel
!******************************* The End ******************************!

!**************** Check whether the patch in Lump ports ***************!
! Input : cent(3)--The coordinate of the triangle center               !
! Output: Indx----The index of the port in which the patch is located  !
!         if the patch is not in ports Indx = 0                        !
!----------------------------------------------------------------------!
subroutine IsInLump( cent, Indx )
    use MoM_VARIABLES, only : FreqMax
    use mod_LumpPort
    implicit none
    real(DBL) cent(3)           ! the index of the patch
    integer Indx
    integer P                   ! the index of the port
    ! the 3 vertices of the port P
    real(DBL) V1(3), V2(3), V3(3), norm(3)
    real(DBL) r0(3), absR

    ! set the patch not in port as default
    Indx = 0
    ! cycle on each port
    do P = 1, LumpNum
        !-- check whether the point locates in the wave ports --!
        ! get the rectangular waveport's vertices
        V1(:) = LumpV( :, P )
        V2(:) = V1(:) + LumpH( :, P )
        V3(:) = V1(:) + LumpE( :, P )
        call cross_mult( LumpE(:,P), LumpH(:,P), norm )

        r0(:) = cent(:) - LumpV( :, P )
        ! get the distance between cent and port P
        absR = abs( sum( norm * r0 ))
        absR = absR * FreqMax /0.03D0
        if( absR > 0.1 ) cycle
            
        ! point cent is in the plane of port
        call lineProj( V1, V2, cent, r0, absR )
        ! if cent is between V1 and V2
        if( absR > 0.0 .AND. absR < 1.0 ) then
            call lineProj( V1, V3, cent, r0, absR )
            ! if cent is between V2 and V3
            if( absR > 0 .AND. absR < 1 ) then
                ! the patch is in port P
                Indx = P
                return
            end if
        end if
    end do

end subroutine IsInLump
!******************************* The End ******************************!

!****************** Check whether the patch in ports ******************!
! Input : cent(3)--The coordinate of the triangle center               !
! Output: Indx----The index of the port in which the patch is located  !
!         if the patch is not in ports Indx = 0                        !
!----------------------------------------------------------------------!
subroutine IsInPort( cent, Indx )
    use MoM_VARIABLES, only : Freq0
    use mod_Waveport
    implicit none
    real(DBL) cent(3)            ! the index of the patch
    integer Indx
    integer P                    ! the index of the port
    ! the 3 vertices of the port P
    real(DBL) V1(3), V2(3), V3(3)
    real(DBL) r0(3), absR

    ! set the patch not in port as default
    Indx = 0
    ! cycle on each port
    do P = 1, PortNum
        !-- check whether the point locates in the wave ports --!
        r0(:) = cent(:) - origin( :, P )
        ! get the distance between cent and port P
        absR = abs(sum( a_in(:,P) * r0 ))
        absR = absR * Freq0 /0.03
        
        if( absR > 0.1 ) cycle
            
        if( 1 == PortKind(P) ) then
            ! get the rectangular waveport's vertices
            V1(:) = origin( :, P )
            V2(:) = V1(:) - width( P ) * h_in( :, P )
            V3(:) = V1(:) + heigh( P ) * e_in( :, P )
            ! point cent is in the plane of port
            call lineProj( V1, V2, cent, r0, absR )
            ! if cent is between V1 and V2
            if( absR > 0 .AND. absR < 1 ) then
                call lineProj( V1, V3, cent, r0, absR )
                ! if cent is between V2 and V3
                if( absR > 0 .AND. absR < 1 ) then
                    ! the patch is in port P
                    Indx = P
                    return
                end if
            end if
        else if( 2 == PortKind(P) ) then
            absR = sum( r0 * r0 )
            absR = sqrt( absR )
            if( absR < width( P ) ) then
                ! the patch is in port P
                Indx = p
                return
            end if
        else
            absR = sum( r0 * r0 )
            absR = sqrt( absR )
            if( absR > width(P) .AND. absR < heigh(P) ) then
                ! the patch is in port p
                Indx = p
                return
            end if
        end if
    end do
end subroutine IsInPort
!******************************* The End ******************************!

!----------------------------------------------------------------------!
! This subroutine is designed to calculate the projection of point r in!
! the given line.                                                      !
!----------------------------------------------------------------------!
subroutine lineProj( V1, V2, r, r_proj, k )
    ! Input : V1(3)---The first vertex of the line                     !
    !       : V2(3)---The second vertex of the line                    !
    !       : r(3)----The point to be projected                        !
    ! Output: r_proj(3)---The projection of r in the line              !
    !       : k---Is the project point in or outside the line segment  !
    !------------------------------------------------------------------!
    use MoM_VARIABLES, only : DBL
    implicit none
    real(DBL) V1(3), V2(3)    ! the vertices of the line
    real(DBL) r(3), r_proj(3) ! r_proj is the projection of r
    real(DBL) k
    real(DBL) line(3)         ! the vector of the line

    line = V2 - V1            ! the vector along the line
    k = sum( (r-V1)*line )/sum( line*line )
    r_proj = V1 + k*line
end subroutine
!***************************** The End ********************************!

!********* Determine the cross multiplication of two vectors *********!
!                                c = a X b
subroutine cross_mult( a, b, c )
    use MoM_VARIABLES, only : DBL
    implicit none
    real(DBL) a(3), b(3), c(3)
    real(DBL) temp(3)

    temp(1) = a(2)*b(3) - a(3)*b(2)
    temp(2) = a(3)*b(1) - a(1)*b(3)
    temp(3) = a(1)*b(2) - a(2)*b(1)
    c(:) = temp(:)
end subroutine
!***************************** The End *******************************!
