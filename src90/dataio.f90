!=======================================================================
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!               University of Illinois at Urbana-Champaign
!
      subroutine dataio(iswres,iswhdf,iswtsl,iswhst,iswusr)
!
!  PURPOSE:  Controls data I/O for restart, HDF, history, and ASCII
!  dumps.
!
!  INPUT ARGUMENTS: iswres,iswhdf,iswhst=switches for restart,hdf, and
!    history dumps.  Values of 1 ensure dumps will be made.
!
!  OUTPUT ARGUMENTS: none
!
!
!  LAST MODIFIED:  by John Hayes; 9/9/97
!  LAST MODIFIED:  by efh 04/15/99 including call to tslice
!  LAST MODIFIED:  by John Hayes, 05/26/2005; removed call to tslice
!-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use root
      use mpiyes
      use mpipar
!-----------------------------------------------------------------------
      implicit NONE
!
      integer  :: iswres,iswhdf,iswhst,iswusr,iswshl,iswtsl
      integer  :: incr,strtoi
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================
!
!FH   if (dtdump .gt. 0.0 .and. tused .ge. (tdump+dtdump)) then
      if (dtdump .gt. 0.0 .and. time .ge. (tdump+dtdump)) then
        tdump = tdump + dtdump
        iswres = 1
      endif
      if (dthdf  .gt. 0.0 .and. time  .ge. (thdf+dthdf)) then
        thdf = thdf + dthdf
        iswhdf = 1
      endif
      if (dttsl  .gt. 0.0 .and. time  .ge. (ttsl+dttsl)) then
        ttsl = ttsl + dttsl
        iswtsl = 1
      endif
      if (dthist .gt. 0.0 .and. time  .ge. (thist+dthist)) then
        thist = thist + dthist
        iswhst = 1
      endif
      if (dtusr .gt. 0.0 .and. time  .ge. (tusr+dtusr)) then
        tusr = tusr + dtusr
        iswusr = 1
      endif
!
!  restart dump
!
      if(xrestart) then
      if (iswres .eq. 1) then
        call msave(resfile)
        incr = strtoi(resfile,13,16) + 1
        write(resfile,"(a3,a2,3i2.2,'.',i4.4)") 'res',id,coords(1) &
                                            ,coords(2),coords(3),incr
        iswres=0
      endif
      endif ! xrestart
!
!  HDF dump
!
      if(xhdf) then
       if (iswhdf .eq. 1) then
        call hdfall(hdffile)
        incr = strtoi(hdffile,13,16) + 1
        write(hdffile,"(a3,a2,3i2.2,'.',i4.4)") 'hdf',id,coords(1) &
                                            ,coords(2),coords(3),incr
!
!  history dump
!
      if(xhst) then
        if (iswhst .eq. 1) then
           call printd 
           iswhst=0
         endif 
      endif
!
! If an array of output times t_out has been defined, set the value
! of (thdf + dthdf) to the next desired output time.
!
        if (incr .le. nbuff-8) then
          if (t_out(incr).ne.0.0) then
            thdf = t_out(incr) - dthdf
          endif
        endif
        iswhdf=0
       endif
      endif ! xhdf
!
!  ascii dump
!
      if(xascii) then
      if (iswusr .eq. 1) then
        call textdmp2
        incr = strtoi(usrfile,13,16) + 1
        write(usrfile,"(a3,a2,3i2.2,'.',i4.4)") 'usr',id,coords(1) &
                                            ,coords(2),coords(3),incr
!
! If an array of output times t_out has been defined, set the value
! of (tusr + dtusr) to the next desired output time.
!
        if (incr .le. nbuff-8) then
          if (t_out(incr).ne.0.0) then
            tusr = t_out(incr) - dtusr
          endif
        endif
        iswusr=0
      endif
      endif ! xascii
!
      return
      end
