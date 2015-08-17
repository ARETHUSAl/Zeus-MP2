!-----------------------------------------------------------------------
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
program zmp_pp
!
! ZEUS-MP Post-processor: Main Program.
!
! PURPOSE:  This program asks the user what post-processing taskes
!           to perform, and calls the appropriate routines.
!
! LAST MODIFIED: RAF, 9/12/96.
!
!.......................................................................
!
      use real_prec
      use param
!
      implicit NONE
!
      integer :: i, ltask
      character :: task*12
!.......................................................................
!
! Write greeting.
!
       write(*,"(/'ZMP_PP: Welcome to the ZEUS-MP Post-Processor!')")
       ltask = len( task )
!
! Loop back here for next task.
!
   10  continue
!
! Determine task.
!
       do 20 i=1,ltask
         task(i:i) = ' '
   20  continue
       write(*,"(/'What task would you like performed?')")
       write(*,"( 'Type one of the following:'/)")
       write(*,"( '   h5splice   (concatenate HDF5 files one by one)')")
       write(*,"( '   auto_h5    (automated h5splice)')")
       write(*,"( '   quit       (exit ZMP_PP)'/)")
       read(*,"(a)") task
!
! QUIT
!
       if (task(1:4) .eq. 'quit') then
         stop
       else if (task(1:8) .eq. 'h5splice') then
!
! SPLICE HDF5 DUMPs -- HDF5 Scientific Data Sets "hdf<id><coords>.<ndump>"
!
         call h5splice
         goto 10
       else if (task(1:7) .eq. 'auto_h5') then
!
! AUTO-SPLICE HDF5 DUMPs -- HDF5 Scientific Data Sets "hdf<id><coords>.<ndump>"
!
         call auto_h5
         goto 10
       else
!
! WHAT?
!
         write(*,"(/'ZMP_PP: Please retype your task.')")
         goto 10
       endif
!
end program zmp_pp
