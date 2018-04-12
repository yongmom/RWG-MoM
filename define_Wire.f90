!----------------------------------------------------------------------!
!  This source file is written to define the triangular basis function !
! on each wire for RWG-MoM program. It is the key part of MoM.         !
!                Designed by Wang Yong in Nov. 5th 2016                !
!----------------------------------------------------------------------!

!***** find the joined point and define triangular basis function *****!
subroutine define_Wire( )
    use mod_Mesh
    use MoM_VARIABLES
    implicit none
    integer I, J, m, n ! loop control variable
    integer extra, NunW
    !----temporary variables of triangle information----!
    integer,allocatable :: AllNode(:)
    integer,allocatable :: indx(:)
    !---------------------------------------------------!

    ! the number of joined point is no more than 2 times of wire
    NunW = 2 * WireNum
    allocate( allNode( NunW ) )       ! all wire nodes
    allocate( indx( NunW ) )

    !-------------------- find out the common node -------------------!
    do m = 1, WireNum
        n = 2 * m
        allNode(n-1) = WireVert( 1, m )
        allNode( n ) = WireVert( 2, m )
    end do

    call NodeSort( NunW, allNode, indx )

    ! initialize the extra wires
    extra = 0
    NunW  = 0
    I = 1
    do while( I < 2*WireNum )
        if( allNode(I) == allNode(I+1) ) then
            J = I + 1
            ! check how many wires attached to the junction
            do while( J < 2*WireNum )
                if( allNode(J) /= allNode(J+1) ) exit
                J = J + 1
            end do
            
            extra = extra + J - I - 1
            NunW = NunW + J - I
            
            I = J + 1
        else
            I = I + 1
        end if
    end do
    ! get the number of mom wire
    NmomWire = WireNum + extra

    NunJ = NunJ + NunW
    
    deallocate( allNode, indx )
    !------------------------------------------------------------------!
    
end subroutine define_Wire
!******************************* The End ******************************!

!*********************** Sort the triangle edge ***********************!
!   This routine performs an in-memory sort of the first N elements of !
! array datum, returning into array indx the indices of elements of    !
! datum arranged in ascending order.                                   !
!----------------------------------------------------------------------!
!                   Created:  15 Jul 1986  Len Moss                    !
!                  Modified:  22 May 2014  Wang Yong                   !
!----------------------------------------------------------------------!
subroutine NodeSort( N, datum, indx )
! edge in each triangle                                                !
! Input : N --The number of the element of datum                       !
!         datum--The element which needed to be sorted                 !
! Output: Indx----The index of node                                    !
!----------------------------------------------------------------------!
    implicit none
    integer N
    integer datum( N )
    integer indx( N )
    ! local variable
    integer,parameter::Cutoff = 9
    integer LSTK(31), RSTK(31), ISTK
    integer L, R, I, J, P
    integer indxP, indxT
    integer datP, datT
    logical FLAG

    ! Make initial guess for indx
    do I = 1, N
        indx( I ) = I
    end do

    !------------------------------------------------------------------!
    ! if the array is short, use the straight insertion sort directly. !
    if ( N <= Cutoff ) then
        !-------------------------------------------------------!
        !                 Straight Insertion sort               !
        !-------------------------------------------------------!
        do I = 2, N
            FLAG = datum(I-1) > datum(I)
            if( FLAG ) then
                indxP = indx( I )
                datP = datum( I )
                indx( I ) = indx( I-1 )
                datum( I ) = datum( I-1 )
                do P = I-1, 2, -1
                    FLAG = datum(P-1) > datP
                    if( .NOT. FLAG ) exit
                    indx(P) = indx(P-1)
                    datum( P ) = datum( P-1 )
                end do
                indx( P ) = indxP
                datum( P ) = datP
            end if
        end do
        ! the sort is completed successfully and leave back.
        return

    end if
    !------------------------------------------------------------------!

    ! QuickSort
    ! Select the median of the first, last, and middle elements as !
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
        ! Let the pivot P be the midpoint of this subsequence,
        ! P=(L+R)/2; then rearrange indx(L) indx(P) and indx(R). so
        ! the corresponding datum values are in increasing order.
        ! The pivot key, datP, is then datum[P].
        P = ( L + R )/2
        indxP = indx( P )
        datP = datum( P )
 
        if( datP < datum(L) ) then
            indx( P ) = indx( L )
            indx( L ) = indxP
            indxP = indx( P )
            datum( P )=datum( L )
            datum( L )=datP
            datP = datum( P )
        end if
 
        if( datP > datum(R) ) then
            if( datum(L) > datum(R) ) then
                indx( P ) = indx( L )
                indx( L ) = indx( R )
                datum( P ) = datum( L )
                datum( L ) = datum( R )
            else
                indx( P ) = indx( R )
                datum( P ) = datum( R )
            end if
            indx( R ) = indxP
            datum( R )= datP
            indxP = indx( P )
            datP = datum( P )
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
            do while( datP > datum(I) )
                I = I + 1
            end do

            ! Q4: Search for datum on right <= datP
            !
            ! At this point, datum[R] >= datP.  We can therefore start
            ! scanning down from R, looking for a value <= datP (this
            ! scan is guaranteed to terminate since we initially placed
            ! datP near the middle of the subsequence).
            J = J - 1
            do while( datP < datum(J) )
                J = J -1
            end do

            ! Q5: Have the two scans collided?
            if( I < J ) then
                ! Q6: No, interchange datum[I] <-> datum[J] and continue
                indxT = indx( I )
                indx( I ) = indx( J )
                indx( J ) = indxT
                datT = datum( I )
                datum( I ) = datum( J )
                datum( J ) = datT
            else
                ! Q7: Yes, select next subsequence to sort. At this point, 
                ! I >= J and datum[L] <= datum[I] == datP <= datum[R],for
                ! all L <= l < I and J < r <= R.  if both subsequences are
                ! more than M elements long, push the longer one on the
                ! stack and go back to QuickSort the shorter; if only one
                ! is more than M elements long, go back and QuickSort it;
                ! otherwise, pop a subsequence off the stack and QuickSort it.
                if( R+L >= I+J .AND. I-L > Cutoff ) then
                    ISTK = ISTK + 1
                    LSTK(ISTK) = J + 1
                    RSTK(ISTK) = R
                    R = I - 1
                else if( R+L < I+J .AND. R-J > Cutoff ) then
                    ISTK = ISTK + 1
                    LSTK(ISTK) = L
                    RSTK(ISTK) = I - 1
                    L = J + 1
                else if( R-J > Cutoff ) then
                    L = J + 1
                else if( I-L > Cutoff ) then
                    R = I - 1
                else
                ! Q8: Pop the stack, or terminate QuickSort if empty
                    if( ISTK < 1 ) then
                        !--------------------------------------------------!
                        !             Straight Insertion sort              !
                        !--------------------------------------------------!
                        do I = 2, N
                            FLAG = datum(I-1) > datum(I)
                            if( FLAG ) then
                                indxP = indx( I )
                                datP = datum( I )
                                indx( I ) = indx( I-1 )
                                datum( I ) = datum( I-1 )
                                do P = I-1, 2, -1
                                    FLAG = datum(P-1) > datP
                                    if( .NOT. FLAG ) exit
                                    indx(P) = indx(P-1)
                                    datum( P ) = datum( P-1 )
                                end do
                                indx( P ) = indxP
                                datum( P ) = datP
                            end if
                        end do
                        !--------------------------------------------------!
                        ! the sort is completed successfully and leave back.
                        return

                    end if

                    L = LSTK( ISTK )
                    R = RSTK( ISTK )
                    ISTK = ISTK - 1
                end if
                ! jump out the cycle
                exit
            end if
        end do
    end do
end subroutine
!******************************* The End ******************************!