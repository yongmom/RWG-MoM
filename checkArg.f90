!----------------------------------------------------------------------!
!              Last update by Wang Yong in Apr. 24 2016                !
!----------------------------------------------------------------------!
    
!*** Check the argument of command line and get the path of project ***!
subroutine checkArg( ProjPath, ProjFile )
    implicit none
    character(len=500) ProjPath, ProjFile
    integer m, n
    logical alive

    call getArg( 1, ProjFile )
    
    m = index( ProjFile, '/', .TRUE. )
    if( 0 == m ) m = index( ProjFile, '\', .TRUE. )
    if( 0 == m ) then
        ProjPath = ''
    else
        n = len(trim(ProjFile))
        ProjPath = ProjFile(1:m)
    end if

    inquire( directory = trim(ProjPath)//'post', exist = alive )
    if( .NOT. alive ) call system( 'mkdir '//trim(ProjPath)//'post' )
    
end subroutine checkArg
!****************************** The End ******************************!
