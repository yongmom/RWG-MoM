!************** Read the model file from gdm file format **************!
subroutine readMesh( log )
    use mod_Mesh
    use MoM_VARIABLES  
    implicit none
    integer log
    integer,allocatable :: WireMedi(:)
    integer,allocatable :: MediTemp(:,:)
    integer,allocatable :: BoudTemp(:)
    integer,allocatable :: newBoud(:)
    integer,allocatable :: Temp(:,:)
    integer m, n, I, J, V(3)
    real(DBL) rad

    !------------------------- Read mesh file -------------------------!
    ! open the original model file
    open( unit = 1, file = trim(meshFile)//".gdm" )

    read( 1, * ) VertNum, WireNum, TriaNum, n
    read( 1, * ) NBound

    if( VertNum < 2 ) then
        write(*,*) "Error : There is at least 2 nodes."
        write(log,*) "Error : There is at least 2 nodes."
        stop
    end if
    if( WireNum < 1 .AND. TriaNum < 1 ) then
        write(*,*) "Error : There is no wire and triangle."
        write(log,*) "Error : There is no wire and triangle."
        stop
    end if
    
    if( TriaNum > 0 .AND. NBound <= 0 ) then
        write(*,*) "Error : There is at least one medium surface."
        write(log,*) "Error : There is at least one medium surface."
        stop
    end if

    if( NBound > 0 ) then
        if( TriaNum <= 0 ) then
            ! read and skip the boundary
            do m = 1, NBound
                read( 1, * ) I, J
            end do
            ! if there is no triangle, delete the boundary
            NBound = 0
        else
            allocate( MediTemp( 2,NBound ) )
            allocate( Temp( 2,NBound ) )
            allocate( newBoud(NBound) )
            ! read the Mediums across each boundry
            do m = 1, NBound
                read( 1, * ) I, J
                if( I * J < 0 ) then
                    write(*,*) "Error : Wrong medium setting on boundry", m
                    write(log,*)"Error : Wrong medium setting on boundry",m
                    stop
                end if
                if( I < J ) then
                    Temp(1,m) = I
                    Temp(2,m) = J
                else
                    Temp(1,m) = J
                    Temp(2,m) = I
                end if

                if( Temp(2,m) == 0 ) then
                    write(*,*) "Error : Wrong medium setting on boundry", m
                    write(log,*)"Error : Wrong medium setting on boundry",m
                    stop
                end if
            
                if( abs(I) > nmed .OR. abs(J) > nmed ) then
                    write(*,*) "Error : Medium on boundry", m, 'exceeds'
                    write(log,*)"Error : Medium on boundry",m, 'exceeds'
                    stop
                end if
        
                newBoud( m ) = m
            end do

            ! remove the duplicate boundary
            n = 0
            do m = 1, NBound
                if( Temp(1,m) == Temp(2,m) ) then
                    newBoud( m ) = 0
                else
                    J = 0
                    do I = 1, m-1
                        if( all( Temp(:,m) == Temp(:,I) ) ) then
                            newBoud( m ) = newBoud( I )
                            J = 1
                            exit
                        end if
                    end do
                    if( 0 == J ) then
                        n = n + 1
                        newBoud( m ) = n
                        MediTemp(:,n) = Temp(:,m)
                    end if
                end if
            end do
            NBound = n
        
            deallocate( Temp )
        end if        
    end if
    
    ! apply memory for the mesh
    allocate( Vertex_r( 3,VertNum ) ) ! the coordinate of vertex
    allocate( WireVert( 2,WireNum ) ) ! the index of vertex
    allocate( TriaVert( 3,TriaNum ) ) ! the index of vertex
    allocate( BoudTemp( TriaNum ) )   ! boundry of triangles
    allocate( WireMedi( WireNum ) )   ! which domain each wire in
    allocate( radius( WireNum ) )     ! radius of each wire

    ! begin to read the mesh
    do m = 1, VertNum
        read( 1, * ) Vertex_r(:,m)
    end do

    do m = 1, WireNum
        read( 1, * ) WireVert(:,m), radius(m), WireMedi(m)
        if( radius(m) <= 0 ) then
            write(*,*) "Error : Wrong radius of wire no.", m
            write(log,*)"Error : Wrong radius of wire no.",m
            stop
        end if
        if( WireMedi(m) <= 0 .OR. WireMedi(m) > nmed ) then
            write(*,*) "Error : Wrong medium of wire no.", m
            write(log,*)"Error : Wrong medium of wire no.",m
            stop
        end if
    end do

    do m = 1, TriaNum
        read( 1, * ) TriaVert(:,m), n
        ! mark the boundary of triangle
        BoudTemp(m) = newBoud( n )
    end do
    ! close the model file
    close( 1 )
    !------------------------- Read file end --------------------------!

    !---------------- Arrange the medium for triangle -----------------!
    if( NBound > 0 ) then
        allocate( BoudMedi( 2,NBound ) )
        allocate( BoudTria( 2,NBound ) )     ! boundry of triangles
        
        BoudMedi( :,1:NBound ) = MediTemp( :,1:NBound )
        ! re-sort the triangles according to the boundary
        I = 0
        do m = 1, NBound
            BoudTria( 1, m ) = I + 1
            do n = I+1, TriaNum
                if( BoudTemp(n) == m ) then
                    I = I + 1
                    V(:) = TriaVert(:,I)
                    TriaVert(:,I) = TriaVert(:,n)
                    TriaVert(:,n) = V(:)
                    BoudTemp(n) = BoudTemp(I)
                end if
            end do
            BoudTria( 2, m ) = I
        end do
        TriaNum = BoudTria( 2, NBound )
        
        deallocate( BoudTemp )
        deallocate( MediTemp )
    end if
    !------------------------------------------------------------------!
    
    !------------------- Arrange the medium for wire ------------------!
    if( WireNum > 0 ) then
        ! count the number of wire domain
        allocate( Temp(WireNum,1) )
        NDomain = 1
        Temp(1,1) = WireMedi(1)
        do I = 2, WireNum
            m = 0
            do J = 1, NDomain
                if( WireMedi(I) == Temp(J,1) ) then
                    m = 1
                    exit
                end if
            end do
            if( m == 0 ) then
                NDomain = NDomain + 1
                Temp( NDomain, 1 ) = WireMedi(I)
            end if
        end do
        
        ! count the number of wire in each domain
        allocate( DomaMedi( NDomain ) )
        allocate( DomaWire( 2,NDomain ) )
        DomaMedi(:) = Temp( 1:NDomain, 1 )
        deallocate( Temp )

        m = 0
        do I = 1, NDomain
            DomaWire( 1, I ) = m + 1
            do J = m+1, WireNum
                if( DomaMedi( I ) == WireMedi(J) ) then
                    m = m + 1
                    V(1:2) = WireVert( :, m )
                    WireVert( :, m ) = WireVert( :, J )
                    WireVert( :, J ) = V(1:2)
                    rad = radius( m )
                    radius( m ) = radius( J )
                    radius( J ) = rad
                    WireMedi(J) = WireMedi(m)
                end if
            end do
            DomaWire( 2, I ) = m
        end do
    end if
    !------------------------------------------------------------------!
end subroutine readMesh
!******************************* The End ******************************!

!**************** Read the model file of msh file type ****************!
subroutine read_msh( )
    use mod_Mesh
    use MoM_VARIABLES
    implicit none
    integer ios       ! the flag of the model file status
    integer(kind=4) m, n

    !--------------------------- Open files ---------------------------!
    ! open the original model file
    open( unit = 1, file = trim(ProjFile) )
    VertNum = 0
    WireNum = 0
    TriaNum = 0
    read( 1, * )  ! MESH dimension 3 ElemType Triangle Nnode 3
    read( 1, * )  ! Coordinates     
    do
        read( 1, *, IOSTAT = ios ) m
        if( ios /= 0 ) exit    ! end coordinates     
        VertNum = VertNum + 1  ! number of the vertex
    end do
      
    read( 1, * )  ! ' '
    read( 1, * )  ! Elements

    do
        read( 1, *, IOSTAT = ios ) m
        if( ios /= 0 ) exit    ! end elements
        TriaNum = TriaNum + 1  ! number of the triangle
    end do

    read( 1, * )  ! MESH dimension 3 ElemType Linear Nnode 2
    read( 1, * )  ! Coordinates
    read( 1, * )  ! End Coordinates
    read( 1, * )  ! ' '
    read( 1, * )  ! Elements

    do
        read( 1, *, IOSTAT = ios ) m
        if( ios /= 0 ) exit    ! end elements
        WireNum = WireNum + 1  ! number of the triangle
    end do
    close( 1 )                 ! close the model file

    ! begin to read the mesh file
    m = VertNum + 6 * WireNum
    n = TriaNum + 4 * WireNum
    allocate( Vertex_r( 3, m ) )    ! the coordinate of vertex
    allocate( TriaVert(3, n ) )     ! the index of vertex
    allocate( WireVert(2, WireNum ) )

    ! load the model datum
    open( unit = 1, file = trim(ProjFile) )

    read( 1, * )  ! MESH dimension 3 ElemType Triangle Nnode 3
    read( 1, * )  ! Coordinates     
    do m = 1, VertNum
        read( 1, * ) n, Vertex_r(:,m)  
    end do
    
    read( 1, * )  ! End coordinates
    read( 1, * )  ! ' '
    read( 1, * )  ! Elements
    do m = 1, TriaNum
        read( 1, * ) n, TriaVert(:,m)
    end do

    read( 1, * )  ! End Elements
    read( 1, * )  ! MESH dimension 3 ElemType Linear Nnode 2
    read( 1, * )  ! Coordinates
    read( 1, * )  ! End Coordinates
    read( 1, * )  ! ' '
    read( 1, * )  ! Elements

    do m = 1, WireNum
        read( 1, * ) n, WireVert(:,m)
    end do

    close( 1 )    ! close the model file
end subroutine read_msh
!******************************* The End ******************************!

!**************** Read the model file of nas file type ****************!
subroutine read_nas( )
    use mod_Mesh
    use MoM_VARIABLES
    implicit none
    integer ios       ! the flag of the model file status
    integer(kind=4) m, n
    character OneLine*80
    character XCoord*16, YCoord*16, ZCoord*16

    !--------------------------- Open files ---------------------------!
    ! open and check the original model file
    open( unit = 1, file = trim( ProjFile ) )
    VertNum = 0
    TriaNum = 0
  
    do
        read( 1, "(A80)", IOSTAT=ios ) OneLine
        if( ios /= 0 ) exit    ! end of file

        if( "GRID*" == OneLine(1:5) ) then
            VertNum = VertNum + 1
        end if
        if( "CTRIA3" == OneLine(1:6) ) then
            TriaNum = TriaNum + 1
        end if

        if( "ENDDATA" == OneLine(1:7) ) exit

    end do
    close( 1 )
    !------------------------------------------------------------------!

    allocate( Vertex_r( 3, VertNum ) )  ! the coordinate of vertex
    allocate( TriaVert( 3, TriaNum ) )  ! the index of vertex
    ! load the model datum
    open( unit = 1, file = trim(ProjFile) )
    m = 0
    n = 0
    do
        read( 1, "(A80)", IOSTAT=ios ) OneLine
        if( ios /= 0 ) exit      ! end coordinates
        if( "ENDDATA" == OneLine(1:7) ) exit

        if( "GRID*" == OneLine(1:5) ) then
            XCoord(:) = OneLine(41:56)
            YCoord(:) = OneLine(57:72)
            read( 1, "(A80)" ) OneLine
            ZCoord(:) = OneLine(9:24)
            m = m + 1
            read( XCoord , * ) Vertex_r( 1, m )
            read( YCoord , * ) Vertex_r( 2, m )
            read( ZCoord , * ) Vertex_r( 3, m )
        end if

        if( "CTRIA3" == OneLine(1:6) ) then
            n = n + 1
            XCoord(:) = OneLine(25:32)
            YCoord(:) = OneLine(33:40)
            ZCoord(:) = OneLine(41:48)
            read( XCoord,* ) TriaVert( 1, n )
            read( YCoord,* ) TriaVert( 2, n )
            read( ZCoord,* ) TriaVert( 3, n )
        end if
    end do
    close( 1 )
end subroutine read_nas
!******************************* The End ******************************!

!**************** Read the model file of raw file type ****************!
subroutine read_raw( )
    use mod_Mesh
    use MoM_VARIABLES
    implicit none
    logical ok
    character(len=250) object
    integer node, FileEnd    ! the flag of the model file status
    integer i
    integer(kind=4) m, n
    real(DBL) RawValue( 9 )
    real(DBL) dist
    real(DBL),allocatable::tri_xyz( : ,: )

    !--------------------------- Open files ---------------------------!
    ! open the original model file
    open( unit = 1, file = trim(ProjFile) )
    !---- Scan the model file and count the number of the triangles ---!
    TriaNum = 0            ! initial the number of the triangles
    do
        read( 1,"(A250)",iostat = FileEnd ) object
        if( FileEnd /= 0 ) exit		! if reach the end of the file
        node=len_trim(object(30:32))+len_trim(object(32:34))
        if( node /= 0 ) then
            TriaNum = TriaNum + 1
        end if
    end do
    close( 1 )   ! close the model file

    !---------------------- load the model datum ----------------------!
    allocate( tri_xyz( 3 , 3*TriaNum ) )   ! the coordinate of vertex
    allocate( TriaVert( 3, TriaNum ) ) ! the index of vertex

    open( unit = 1, file = trim(ProjFile) )
    VertNum = 0
    n = 0
    do
        read( 1,"(A250)",iostat = FileEnd ) object
        if( FileEnd /= 0 ) exit		! if reach the end of the file
        node=len_trim(object(30:32))+len_trim(object(32:34))
        if( node /= 0 ) then
            n = n + 1
            read( object, * ) RawValue( : )
            if( 0 == VertNum ) then
                tri_xyz( :, 1 ) = RawValue( 1:3 )
                tri_xyz( :, 2 ) = RawValue( 4:6 )
                tri_xyz( :, 3 ) = RawValue( 7:9 )
                VertNum = 3
                TriaVert( :, 1 ) = (/1,2,3/)
                cycle
            end if
            do i = 1, 3
                ok = .TRUE.
                do m = 1, VertNum
                    dist=sum( abs( RawValue(3*i-2:3*i)-tri_xyz(:,m) ) )
                    ! find the coincident vertex
                    if( dist < 1.0e-10 ) then
                        ok = .false.		! it's not a new vertex
                        TriaVert( i, n ) = m
                        exit				! break the inner cycle
                    end if
                end do
                if( ok ) then				! it's a new vertex
                    ! add the number of vertice
                    VertNum = VertNum + 1
                    TriaVert( i, n ) = VertNum
                    tri_xyz( :, VertNum ) = RawValue( 3*i-2 : 3*i )
                end if
            end do
        end if
    end do
    allocate( Vertex_r( 3, VertNum ) )
    Vertex_r(:,:) = tri_xyz(:,1:VertNum )
    ! release the storage space which will not be used
    deallocate( tri_xyz )
    close( 1 )
    
end subroutine read_raw
!******************************* The End ******************************!

!**************** Read the model file of stl file type ****************!
subroutine read_stl( )
    use mod_Mesh
    use MoM_VARIABLES
    implicit none
    logical ok
    character(len=50) object
    integer FileEnd    ! the flag of the model file status
    integer I
    integer(kind=4) m, n
    real(DBL) node( 3, 3 )
    real(DBL) dist
    real(DBL),allocatable::tri_xyz( : ,: )

    !--------------------------- Open files ---------------------------!
    ! open the original model file
    open( unit = 1, file = trim(ProjFile) )
    !---- Scan the model file and count the number of the triangles ---!
    TriaNum = 0            ! initial the number of the triangles
    read( 1, * )
    do
        read( 1,*, iostat = FileEnd ) object
        if( FileEnd /= 0 ) exit		! if reach the end of the file
        if( trim(object) == "endfacet" ) then
            TriaNum = TriaNum + 1
        end if
    end do
    close( 1 )   ! close the model file

    !---------------------- load the model datum ----------------------!
    allocate( tri_xyz( 3 , 3*TriaNum ) )   ! the coordinate of vertex
    allocate( TriaVert( 3, TriaNum ) ) ! the index of vertex

    open( unit = 1, file = trim(ProjFile) )
    VertNum = 0
    read( 1, * )
    do n = 1, TriaNum
        read( 1, * )
        read( 1, * )
        read( 1, * ) object, node( :,1 )
        read( 1, * ) object, node( :,2 )
        read( 1, * ) object, node( :,3 )

        if( 0 == VertNum ) then
            tri_xyz( :, 1 ) = node( :,1 )
            tri_xyz( :, 2 ) = node( :,2 )
            tri_xyz( :, 3 ) = node( :,3 )
            TriaVert( :, 1 ) = (/1,2,3/)
            VertNum = 3
        else
            do I = 1, 3
                ok = .TRUE.
                do m = 1, VertNum
                    dist = sum( abs( node(:,I)-tri_xyz(:,m) ) )
                    ! find the coincident vertex
                    if( dist < 1.0e-10 ) then
                        ok = .FALSE.		! it's not a new vertex
                        TriaVert( I, n ) = m
                        exit				! break the inner cycle
                    end if
                end do
                if( ok ) then				! it's a new vertex
                    ! add the number of vertice
                    VertNum = VertNum + 1
                    TriaVert( I, n ) = VertNum
                    tri_xyz( :, VertNum ) = node( :,I )
                end if
            end do
        end if

        read( 1, * )
        read( 1, * )

    end do
    allocate( Vertex_r( 3, VertNum ) )
    Vertex_r(:,:) = tri_xyz(:,1:VertNum)
    ! release the storage space which will not be used
    deallocate( tri_xyz )
    close( 1 )
    
end subroutine read_stl
!******************************* The End ******************************!