subroutine ran1(idum,x) 
    implicit none
    integer, parameter          :: k4b=selected_int_kind(9)
    integer(k4b), intent(inout) :: idum
    real(8)                     :: x 
    integer(k4b), parameter     :: ia=16807,      &        
                                   im=2147483647, & 
                                   iq=127773,     & 
                                   ir=2836
    real(8), save               :: am
    integer(k4b), save          :: ix=-1,iy=-1,k
 
    if (idum <= 0 .or. iy < 0) then
        am=nearest(1.0,-1.0)/im
        iy=ior(ieor(888889999,abs(idum)),1)
        ix=ieor(777755555,abs(idum))
        idum=abs(idum)+1
    end if
 
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))
    k=iy/iq
    iy=ia*(iy-k*iq)-ir*k
    if (iy<0) iy=iy+im
    x=am*ior(iand(im,ieor(ix,iy)),1)
end subroutine ran1
