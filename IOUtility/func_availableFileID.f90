
integer function availableFileID()
implicit none
integer fileID
logical statusOfOpened
!================================================
fileID=11  ! initial file id start from 11
do
    inquire(fileID,opened=statusOfOpened)
    if ( statusOfOpened .eqv. .false. ) then
        availableFileID=fileID
        exit
    endif
    fileID = fileID+1
enddo
!================================================
return
stop
end function availableFileID
