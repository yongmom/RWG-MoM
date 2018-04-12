!----------------------------------------------------------------------!
! This source file is written to define the RWG basis function on each !
! triangle for general RWG-MoM program. It is the key part of MoM.     !
!                Designed by Wang Yong in Oct. 15th 2016               !
!----------------------------------------------------------------------!

!******** find the common edges and define RWG basis functions ********!
subroutine define_RWG( )
    use mod_Mesh
    !use mod_RWG
    use MoM_VARIABLES
    implicit none
    integer m, n, I, J, K, P, T          ! loop control variable
    integer EdgeNum, iTria, iNode, Idx, flag
    integer newBound, extra, M1, V(2)
    !----temporary variables of triangle information----!
    integer,allocatable :: indx(:)
    integer,allocatable :: allEdge(:,:)
    integer,allocatable :: TriaToBound(:,:)
    integer,allocatable :: BoudTriaTmp(:,:)
    integer,allocatable :: LocToGlob(:)
    integer,allocatable :: MediTmp(:,:)
    integer,allocatable :: JTemp(:,:)

    !------------------- define the basis functions -------------------!
    n = 2 * NBound
    allocate( MediTmp( 2,n ) )
    allocate( BoudTriaTmp(2,n ) )

    ! copy the boundary triangle index and medium
    MediTmp(  :, 1:NBound ) = BoudMedi(:,1:NBound)
    BoudTriaTmp(:,1:NBound) = BoudTria(:,1:NBound)
    ! check how many microstrip boundary
    newBound = NBound
    do m = 1, NBound
        if( BoudMedi(2,m) < 0 ) then
            ! modify the medium of microstrip
            MediTmp(1,m) = 0
            MediTmp(2,m) = -BoudMedi(2,m)
            ! copy boundary for microstrip
            newBound = newBound + 1
            MediTmp(1,newBound) = 0
            MediTmp(2,newBound) = -BoudMedi(1,m)
            ! copy the triangle
            BoudTriaTmp( :,newBound ) = BoudTriaTmp( :,m )
        end if
    end do
    
    ! map the triangle to boundary
    allocate( TriaToBound(2,TriaNum) )
    allocate( LocToGlob( newBound ) )
    TriaToBound(:,:) = 0
    LocToGlob(:) = 0     ! microstrip tria. local index to global index

    ! map the boundary information to each tria.
    do I = 1, NBound
        m = BoudTriaTmp(1,I)
        n = BoudTriaTmp(2,I)
        TriaToBound(1,m:n) = I
    end do
    ! mark for the copied microstrip boundary
    K = TriaNum + 1
    do I = NBound+1, newBound
        m = BoudTriaTmp(1,I)
        n = BoudTriaTmp(2,I)
        TriaToBound(2,m:n) = I
        LocToGlob(I) = K - m
        BoudTriaTmp(1,I) = K
        BoudTriaTmp(2,I) = K + n - m
        K = BoudTriaTmp(2,I) + 1
    end do
    
    EdgeNum = 3 * TriaNum               ! number of all triangle edge
    allocate( allEdge( 2, EdgeNum ) )   ! the vector of common edge
    allocate( indx( EdgeNum ) )
    ! fetch all of the triangles edges
    m = 0
    do n = 1, TriaNum
        m = m + 1
        ! the first edge is between point 2 and point 3
        allEdge(1,m) = TriaVert(2,n)
        allEdge(2,m) = TriaVert(3,n)

        m = m + 1
        ! the second edge is between point 3 and point 1
        allEdge(1,m) = TriaVert(3,n)
        allEdge(2,m) = TriaVert(1,n)

        m = m + 1
        ! the third edge is between point 1 and point 2
        allEdge(1,m) = TriaVert(1,n)
        allEdge(2,m) = TriaVert(2,n)
    end do
    ! make sure each edge's first point index less than the second point 
    do m = 1, EdgeNum
        if( allEdge(1,m) > allEdge(2,m) ) then
            ! exchange the vertex
            n = allEdge(1,m)
            allEdge(1,m) = allEdge(2,m)
            allEdge(2,m) = n
        end if
    end do

    ! sort all the triangle edges
    call EdgeSort( EdgeNum, allEdge, indx )
    ! sort triangle in each patch junction
    call LocalSort( MediTmp, TriaNum, TriaToBound, allEdge, indx )

    ! initialize the number of mom triangle
    NmomTria = BoudTriaTmp(2,newBound)
    allocate( JTemp( 3, NmomTria ) )

    ! initialize the mom tirangle
    do I = 1, NmomTria
        JTemp( :, I ) = 0
    end do
    !----------- check how many triangles need to be copied -----------!
    extra = 0
    I = 1
    do while( I < EdgeNum )
        ! define mom tirangle and basis functions
        if( indx(I) > 0 ) then
            ! it's not a junction, only 2 triangle share a common edge
            iTria = ceiling( indx(I) / 3.0 )
            iNode = mod( indx(I), 3 )
            if( 0 == iNode ) iNode = 3

            if( TriaToBound( 2,iTria ) == 0 ) then
                ! deal with general triangle pair
                NunJ = NunJ + 1
                JTemp( iNode, iTria ) = NunJ
                Idx = TriaToBound( 1, iTria )
                if( MediTmp(1,Idx) > 0 ) then
                    NunM = NunM + 1
                end if

                iTria = ceiling( indx(I+1) / 3.0 )
                iNode = mod( indx(I+1), 3 )
                if( 0 == iNode ) iNode = 3

                JTemp( iNode, iTria ) = NunJ
            else
                ! deal with microstrip triangles
                do m = 1, 2
                    iTria = ceiling( indx(I) / 3.0 )
                    iNode = mod( indx(I), 3 )
                    if( 0 == iNode ) iNode = 3

                    NunJ = NunJ + 1
                    Idx = TriaToBound( m, iTria )
                    K = iTria + LocToGlob( Idx )
                    JTemp( iNode, K ) = NunJ
                    
                    iTria = ceiling( indx(I+1) / 3.0 )
                    iNode = mod( indx(I+1), 3 )
                    if( 0 == iNode ) iNode = 3

                    K = iTria + LocToGlob( Idx )
                    JTemp( iNode, K ) = NunJ
                end do
            end if
            ! move down the index
            I = I + 2
        else if( indx(I) == 0 ) then
            ! it is not a common edge
            I = I + 1
        else
            ! define mom triangle and basis function in junction
            indx(I) = -indx(I)
            ! check how many triangles attached to the junction
            J = I + 1
            do while( J < EdgeNum )
                if( any(allEdge(:,J) /= allEdge(:,J+1)) ) exit
                J = J + 1
            end do

            iTria = ceiling( indx(I) / 3.0 )
            Idx = TriaToBound( 1, iTria )
            M1 = MediTmp( 1,Idx )

            if( M1 > 0 ) then
                ! deal with dielectric junction
                NunJ = NunJ + 1
                NunM = NunM + 1
            else
                ! deal with composite junction
                flag = 1
                do m = I, J
                    n = m + 1
                    if( m == J ) then
                        n = I
                        if( 1 == flag ) exit
                    end if

                    ! here skip closed metallic
                    K = ceiling( indx(m) / 3.0 )
                    P = ceiling( indx(n) / 3.0 )
                    T = TriaToBound( 2, K ) + TriaToBound( 2, P )
                    K = TriaToBound( 1, K )
                    P = TriaToBound( 1, P )
                    T = T + MediTmp( 1, K ) + MediTmp( 1, P )
                    K = MediTmp( 2, K )
                    P = MediTmp( 2, P )
                    if( 0 == T .AND. K /= P ) cycle

                    ! deal with the first triangle
                    iTria = ceiling( indx(m) / 3.0 )
                    iNode = mod( indx(m), 3 )
                    if( 0 == iNode ) iNode = 3
                    
                    Idx = TriaToBound( 2, iTria )
                    if( Idx /= 0 ) then
                        flag = 0
                        ! a new electric current
                        NunJ = NunJ + 1
                        ! deal with microstrip
                        call SortMicro( indx,MediTmp, TriaToBound,I,J,m,V )
                        Idx = V(2)
                        K = iTria + LocToGlob( Idx )
                        JTemp( iNode, K ) = NunJ
                        M1 = MediTmp( 2, Idx )
                    else
                        Idx = TriaToBound( 1, iTria )
                        if( MediTmp(1,Idx) == 0 ) then
                            ! a new electric current
                            NunJ = NunJ + 1
                            K = JTemp( iNode, iTria )
                            if( 0 == K ) then
                                JTemp( iNode, iTria ) = NunJ
                            else
                                ! need to be copied
                                extra = extra + 1
                            end if
                            M1 = MediTmp( 2, Idx )
                        end if
                    end if

                    ! deal with the second triangle
                    iTria = ceiling( indx(n) / 3.0 )
                    iNode = mod( indx(n), 3 )
                    if( 0 == iNode ) iNode = 3

                    Idx = TriaToBound( 2, iTria )
                    if( 0 == Idx ) then
                        ! it is not a microstrip triangle
                        Idx = TriaToBound( 1, iTria )
                    else
                        ! it is a microstrip triangle
                        call SortMicro( indx,MediTmp,TriaToBound,I,J,n,V )
                        Idx = V(1)
                    end if

                    if( MediTmp(1,Idx) == 0 ) then
                        ! it is a metallic triangle
                        P = iTria + LocToGlob( Idx )
                        K = JTemp( iNode, P )
                        if( 0 == K ) then
                            JTemp( iNode, P ) = NunJ
                        else
                            ! need to be copied
                            extra = extra + 1
                        end if
                    else
                        ! it is a dielectric triangle
                        flag = 0
                        JTemp( iNode, iTria ) = NunJ
                        if( M1 == MediTmp(1,Idx) ) then
                            M1 = MediTmp(2,Idx)
                        else
                            M1 = MediTmp(1,Idx)
                        end if
                    end if
                end do
            end if
            ! move down the index
            I = J + 1
        end if
    end do

    NmomTria = NmomTria + extra

    ! release the temporary memory
    deallocate( indx, allEdge )
    deallocate( TriaToBound )
    deallocate( LocToGlob )
    deallocate( BoudTriaTmp, JTemp )
    !---------------------- check processing done ---------------------!
    
end subroutine define_RWG
!******************************* The End ******************************!

!*********************** Sort the triangle edge ***********************!
!   This routine performs an in-memory sort of the first N elements of !
! array datum, returning into array indx the indices of elements of    !
! datum arranged in ascending order.                                   !
!----------------------------------------------------------------------!
! EdgeSort use a hybrid QuickSort algorithm. In particular, the "pivot !
! key" for dividing each subsequence is chosen to be the Median of the !
! first, last, and middle values of the subsequence; and the QuickSort !
! is cut off when a subsequence has 9 or fewer elements,and a straight !
! insertion sort of the entire array is done at the end. The result is !
! comparable to a pure insertion sort for very short arrays, and very  !
! fast for very large array(of order 12 micro-sec/element on the 3081K !
! for arrays of 10K elements). It is also not subject to the poor      !
! performance of the pure QuickSort on partially ordered datum.        !
!                   Created:  15 Jul 1986  Len Moss                    !
!                  Modified:  22 May 2014  Wang Yong                   !
!----------------------------------------------------------------------!
subroutine EdgeSort( N, datum, indx )
! edge in each triangle                                                !
! Input : N --The number of the element of datum                       !
!         datum--The element which needed to be sorted                 !
! Output: Indx----The index of the common edge in the first triangle   !
!----------------------------------------------------------------------!
    implicit none
    integer N
    integer datum( 2, N )
    integer indx( N )
    ! local variable
    integer,parameter::Cutoff = 9
    integer LSTK(31), RSTK(31), ISTK
    integer L, R, I, J, P
    integer indxP, indxT
    integer datP(2), datT(2)
    logical FLAG
    logical,external :: EdgeGT

    ! deal with the error input
    if( N <= 0 ) return
    ! Make initial guess for indx
    do I = 1, N
        indx( I ) = I
    end do

    !------------------------------------------------------------------!
    ! if the array is short, use the straight insertion sort directly. !
    if( N <= Cutoff ) then
        do I = 2, N
            FLAG = EdgeGT( datum( :, I-1 ), datum( :, I ) )
            if( FLAG ) then
                indxP = indx( I )
                datP = datum( :, I )
                indx( I ) = indx( I-1 )
                datum( :, I ) = datum( :, I-1 )
                do P = I-1, 2, -1
                    FLAG = EdgeGT( datum( :, P-1 ), datP )
                    if( .NOT. FLAG ) exit
                    indx(P) = indx(P-1)
                    datum( :, P ) = datum( :, P-1 )
                end do
                indx( P ) = indxP
                datum( :, P ) = datP
            end if
        end do
        ! the sort is completed successfully
        return
    end if
    !------------------------------------------------------------------!

    ! QuickSort
    ! Select the Median of the first, last, and middle elements as !
    ! the "pivot key" (in Knuth's notation, "K"). Also modified to !
    ! leave datum in place and produce an indx array.              ! 
    ! Q1: Initialize
    ISTK = 0
    L = 1
    R = N
 
    do while( .TRUE. )
 
        ! Q2: Sort the subsequence datum[L]..datum[R].
        ! At this point, datum[l] <= datum[m] <= datum[r] for all l<L
        ! r > R, and L <= m <= R.  (First time through, there is no
        ! datum for l < L or r > R.)
        I = L
        J = R
 
        ! Q2.5: Select pivot key
        ! Let the pivot, P, be the midpoint of this subsequence,
        ! P=(L+R)/2; then rearrange indx(L), indx(P), and indx(R)
        ! so the corresponding datum values are in increasing order.
        ! The pivot key, datP, is then datum[P].
        P = ( L + R )/2
        indxP = indx( P )
        datP = datum( :, P )
 
        if( EdgeGT( datum( :, L ), datP) ) then
            indx( P ) = indx( L )
            indx( L ) = indxP
            indxP = indx( P )
            datum( :, P )=datum( :, L )
            datum( :, L )=datP
            datP = datum( :, P )
        end if
 
        if( EdgeGT(datP , datum( :, R ) ) ) then
            if( EdgeGT(datum(:,L) , datum(:,R) ) ) then
                indx( P ) = indx( L )
                indx( L ) = indx( R )
                datum( :, P ) = datum( :, L )
                datum( :, L ) = datum( :, R )
            else
                indx( P ) = indx( R )
                datum( :, P ) = datum( :, R )
            end if
            indx( R ) = indxP
            datum( :, R )= datP
            indxP = indx( P )
            datP = datum( :, P )
        end if
 
        ! Now we swap values between the right and left sides and/or
        ! move datP until all smaller values are on the left and all
        ! larger values are on the right.  Neither the left or right
        ! side will be internally ordered yet; however, datP will be
        ! in its final position.
 
        do while( .TRUE. )
            ! Q3: Search for datum on left >= datP
            !
            ! At this point, datum[L] <= datP.  We can therefore start 
            ! scanning up from L, looking for a value >= datP (this 
            ! scan is guaranteed to terminate since we initially placed 
            ! datP near the middle of the subsequence).
            I = I + 1
            do while( EdgeGT( datP, datum( :, I ) ) )
                I = I + 1
            end do

            ! Q4: Search for datum on right <= datP
            !
            ! At this point, datum[R] >= datP.  We can therefore start
            ! scanning down from R, looking for a value <= datP (this
            ! scan is guaranteed to terminate since we initially placed
            ! datP near the middle of the subsequence).
            J = J - 1
            do while( EdgeGT(datum( :, J ),datP) )
                J = J -1
            end do

            ! Q5: Have the two scans collided?
            if( I < J ) then
                ! Q6: No, interchange datum[I] <-> datum[J] and continue
                indxT = indx( I )
                indx( I ) = indx( J )
                indx( J ) = indxT
                datT = datum( :, I )
                datum( :, I ) = datum( :, J )
                datum( :, J ) = datT
            else
                ! Q7: Yes, select next subsequence to sort. At this point, 
                ! I >= J and datum[l] <= datum[I] == datP <= datum[r],for
                ! all L <= l < I and J < r <= R.  if both subsequences are
                ! more than M elements long, push the longer one on the
                ! stack and go back to QuickSort the shorter; if only one
                ! is more than M elements long, go back and QuickSort it;
                ! otherwise, pop a subsequence off the stack and QuickSort it.
                if( R-J >= I-L .AND. I-L > Cutoff ) then
                    ISTK = ISTK + 1
                    LSTK(ISTK) = J + 1
                    RSTK(ISTK) = R
                    R=I-1
                else if( I-L > R-J .AND. R-J > Cutoff ) then
                    ISTK = ISTK + 1
                    LSTK(ISTK) = L
                    RSTK(ISTK) = I - 1
                    L = J + 1
                else if( R-J > Cutoff ) then
                    L = J+1
                else if( I-L > Cutoff ) then
                    R = I-1
                else
                ! Q8: Pop the stack, or terminate QuickSort if empty
                    if( ISTK < 1 ) then
                        !--------------------------------------------------!
                        !              Straight Insertion sort             !
                        !--------------------------------------------------!
                        do I = 2, N
                            FLAG = EdgeGT( datum( :, I-1 ), datum( :, I ) )
                            if( FLAG ) then
                                indxP = indx( I )
                                datP = datum( :, I )
                                indx( I ) = indx( I-1 )
                                datum( :, I ) = datum( :, I-1 )
                                do P = I-1, 2, -1
                                    FLAG = EdgeGT( datum( :, P-1 ), datP )
                                    if( .NOT. FLAG ) exit
                                    indx(P) = indx(P-1)
                                    datum( :, P ) = datum( :, P-1 )
                                end do
                                indx( P ) = indxP
                                datum( :, P ) = datP
                            end if
                        end do
                        !--------------------------------------------------!
                        ! the sort is completed successfully
                        return
                    end if
                    L = LSTK(ISTK)
                    R = RSTK(ISTK)
                    ISTK = ISTK-1
                end if
                ! jump out the cycle
                exit
            end if
        end do
    end do
end subroutine EdgeSort
!******************************* The End ******************************!

!****************** Check whether A is great than B *******************!
logical function EdgeGT( A, B )
    implicit none
    integer A(2), B(2)

    if( A(1) > B(1) ) then
        EdgeGT= .TRUE.
    else if( A(1) == B(1) ) then
        if( A(2) > B(2) ) then
            EdgeGT= .TRUE.
        else 
            EdgeGT= .FALSE.
        end if
    else 
        EdgeGT = .FALSE.
    end if
end function EdgeGT
!******************************* The End ******************************!

!******************** Sort the junction triangles *********************!
subroutine LocalSort( Medi, Num, TtoB, allEdge, indx )
use MoM_VARIABLES, only : DBL
    use mod_mesh, only : Vertex_r, TriaVert
    implicit none
    integer Num, Medi(2,*), TtoB(2,Num), allEdge(2,3*Num), indx(3*Num)
    integer EdgeNum, iNode, iTria, NT, m, n, I, J, K
    logical isCommon
    real(DBL) high(3), norm(3), edge(3), dist
    real,allocatable :: dot(:)

    ! find each junction
    EdgeNum = 3 * Num
    I = 1
    do while( I < EdgeNum )
        if( all( allEdge(:,I) == allEdge(:,I+1) ) ) then
            ! it is a common edge
            J = I + 1
            ! check how many triangles share a common edge
            do while( J < EdgeNum )
                if( any(allEdge(:,J) /= allEdge(:,J+1)) ) exit
                J = J + 1
            end do

            NT = J - I + 1      ! number of joint triangles
            
            if( NT > 2 ) then
                ! make sure the first triangle is a metal, if have
                K = 0
                do m = I, J
                    iTria = ceiling( indx(m) / 3.0 )

                    n = TtoB( 1, iTria )
                    if( 0 == Medi(1,n) ) then
                        n = indx( I )
                        indx( I ) = indx( m )
                        indx( m ) = n
                        K = K + 1
                    end if
                end do
                ! if only one metallic triangle
                if( 1 == K ) then
                    iTria = ceiling( indx(I) / 3.0 )
                    n = TtoB( 2, iTria )
                    ! it's not a microstrip
                    if( 0 == n ) then
                        ! kick out the alone open metallic
                        indx( I ) = 0
                        I = I + 1
                        NT = NT - 1
                    end if
                end if
            end if

            if( NT > 2 ) then
                ! more than 2 triangles joint
                allocate( dot(NT) )

                iTria = ceiling( indx(I) / 3.0 )
                iNode = mod( indx(I), 3 )
                if( 0 == iNode ) iNode = 3

                m = allEdge( 1, I )
                n = allEdge( 2, I )

                ! get the direction of common edge
                norm(:) = Vertex_r( :, n ) - Vertex_r( :,m )
                dist = sum( norm * norm )
                dist = sqrt( dist )
                norm(:) = norm(:) / dist
                
                ! get the direction of altitudes
                n = TriaVert( iNode, iTria )
                edge(:) = Vertex_r( :, n ) - Vertex_r( :,m )
                dist = sum( edge * norm )
                high(:) = edge(:) - dist * norm(:)
                dist = sum( high(:) * high(:) )
                dist = sqrt( dist )
                high(:) = high(:) / dist
                
                do K = I+1, J
                    iTria = ceiling( indx(K) / 3.0 )
                    iNode = mod( indx(K), 3 )
                    if( 0 == iNode ) iNode = 3

                    ! get the direction of altitudes
                    n = TriaVert( iNode, iTria )
                    edge(:) = Vertex_r(:,n) - Vertex_r(:,m)
                    dist = sum( edge * norm )
                    edge(:) = edge(:) - dist * norm(:)
                    dist = sum( edge * edge )
                    dist = sqrt( dist )
                    edge(:) = edge(:) / dist
                    dot(K-I) = sum( high * edge )
                    call cross_mult( high, edge, edge )
                    dist = sum( edge * norm )
                    ! deal with angle great than 180 degree
                    if( dist < 0.0 ) dot(K-I) = -2 - dot(K-I)
                end do
                ! sort the triangle in local junction
                do m = I+1, J
                    do n = m+1, J
                        if( dot(m-I) < dot(n-I) ) then
                            dist = dot(m-I)
                            dot(m-I) = dot(n-I)
                            dot(n-I) = dist
                            K = indx(m)
                            indx(m) = indx(n)
                            indx(n) = K
                        end if
                    end do
                end do

                ! mark for junction with more than 2 triangle
                indx(I) = -indx(I)
                ! free memory
                deallocate( dot )
            else
                ! mark for junction with microstrip triangle
                m = ceiling( indx(I) / 3.0 )
                n = ceiling( indx(J) / 3.0 )
                if( TtoB(2,m) == 0 .AND. TtoB(2,n) /= 0 ) then
                    K = indx(I)
                    indx(I) = indx(J)
                    indx(J) = K
                end if
                if( TtoB(2,m) /= TtoB(2,n) ) then
                    indx(I) = -indx(I)
                end if
            end if

            ! update the index of edge
            I = J
        end if
        ! move to next edge
        I = I + 1
    end do

    ! erase non-common-edge
    if( .NOT. all( allEdge(:,1) == allEdge(:,2) ) ) indx(1) = 0
    do I = 2, EdgeNum - 1
        isCommon = all( allEdge(:,I) == allEdge(:,I-1) )
        isCommon = isCommon .OR. all( allEdge(:,I) == allEdge(:,I+1) )
        if( .NOT. isCommon ) then
            indx(I) = 0
        end if
    end do
    I = EdgeNum
    if( .NOT. all( allEdge(:,I) == allEdge(:,I-1) ) ) indx(I) = 0

end subroutine LocalSort
!******************************* The End ******************************!

!******************** Sort the microstrip boundary ********************!
subroutine SortMicro( indx, Medi, TtoB, I, J, m, V )
    implicit none
    integer indx(*), Medi(2,*), TtoB(2,*), I, J, m, V(2)
    integer iTria, M1, M2, M3, M4, K, n
    logical isCommon

    ! re-mark the boundary order of junction for microstrip
    iTria = ceiling( abs(indx(m)) / 3.0 )
    V(:) = TtoB( :, iTria )

    M1 = Medi( 2, V(1) )
    M2 = Medi( 2, V(2) )
    ! mark no need to exchange the microstrip boundary
    isCommon = .FALSE.

    n = m - 1
    if( m == I ) n = J
    K = ceiling( abs(indx(n)) / 3.0 )
    n = TtoB( 2, K )
    if( n /= 0 ) then
        M3 = Medi( 2, TtoB(1,K) )
        M4 = Medi( 2, n )
    else
        M3 = Medi( 1, TtoB(1,K) )
        M4 = Medi( 2, TtoB(1,K) )
    end if
    if( M1 /= M3 .AND. M1 /= M4 ) then
        isCommon = .TRUE.
    end if
                
    n = m + 1
    if( m == J ) n = I
    K = ceiling( abs(indx(n)) / 3.0 )
    n = TtoB( 2, K )
    if( n /= 0 ) then
        M3 = Medi( 2, TtoB(1,K) )
        M4 = Medi( 2, n )
    else
        M3 = Medi( 1, TtoB(1,K) )
        M4 = Medi( 2, TtoB(1,K) )
    end if
    if( M2 /= M3 .AND. M2 /= M4 ) then
        isCommon = .TRUE.
    end if
                
    if( isCommon ) then
        V(1) = TtoB( 2, iTria )
        V(2) = TtoB( 1, iTria )
    end if

end subroutine SortMicro
!******************************* The End ******************************!
