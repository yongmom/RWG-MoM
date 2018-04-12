!----------------------------------------------------------------------!
!     This program generates the project file for RWG-MoM program      !
!   The necessary input is ***.myc and ***.gdm. ***.myc provids the    !
! electromagnetics parameters,  while  ***.gdm  containing the mesh    !
! information of the target.                                           !
!                        Author:         Wang Yong                     !
!                        Date  :       Jun. 20th 2014                  !
!            Compile platform  : Intel Fortran on CentOS 6.2 Linux     !
!                   Compaq Visual Fortran passed also                  !
!----------------------------------------------------------------------!
program Pre_Main
    use MoM_VARIABLES
    implicit none
    integer,parameter :: L = 25               ! port of the log file

    !----------------------------- Begin ------------------------------!
    call checkArg( ProjPath, ProjFile )

    open( unit=L, file = trim(ProjPath)//"log.out" )
    ! print the introducing information on the screen.
    write(*,*)"####################################################"
    write(*,*)"#  Generate the project file for RWG-MoM program   #"
    write(*,*)"#              Used as pre-processing              #"
    write(*,*)"#       Designed by Wang Yong in 2015.08.27        #"
    write(*,*)"####################################################"
    ! output the introducing information in the log file.
    write(L,*)"####################################################"
    write(L,*)"#  Generate the project file for RWG-MoM program   #"
    write(L,*)"#              Used as pre-processing              #"
    write(L,*)"#       Designed by Wang Yong in 2015.08.27        #"
    write(L,*)"####################################################"

    !----------------------------- Initial ----------------------------!
    call readConf(L)            ! initialize the parameter
    call initModel(L)           ! initialize the model data
    !------------------------------------------------------------------!

    !------------------ Generate the LU file for MoM ------------------!
    write(*,*) 'Now generate the project file...'
    call getPreFile( )          ! initialize the parameter
    !write(*,*) 'Now generate the post-setting file...'
    !call PostMsh( )
    !------------------------------------------------------------------!

    ! print the concluding information
    ! print the concluding information on the screen.
    write(*,*)" !------------------------------------------------!"
    write(*,*)" ! Completed successfully! The project file *.tgm !"
    write(*,*)" ! is stored in file folder : ./                  !"
    write(*,*)" !------------------------------------------------!"
    ! output the concluding information in the log file.
    write(L,*)"*-------------------------------------------------*"
    write(L,*)"* Completed successfully ! The project file *.tgm *"
    write(L,*)"* is stored in file folder : ./                   *"
    write(L,*)"*-------------------------------------------------*"
    close( L )

end program Pre_Main
!****************************** The End *******************************!