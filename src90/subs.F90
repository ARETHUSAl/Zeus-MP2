module subs
  implicit none
  contains 
    subroutine localize_xy(out, in)!, source_proc, dest_pz)
!
!  Localizes global 4D data first along the y, then along the x-direction to
!  the destination processor. The global data is supposed to include the outer
!  ghost layers. The returned data will include inner ghost layers.
!  Inner ghost layers are cut away during the combination of the data.
!
!  23-Apr-2012/Bourdin.KIS: adapted from non-torus-type localize_xy
!
      real, dimension(:,:,:,:), intent(out) :: out
      real, dimension(:,:,:,:), intent(in), optional :: in
      !integer, intent(in), optional :: source_proc, dest_pz
!
      !if (present (source_proc) .or. present (dest_pz)) continue
      if (present (in)) out = in
!
    end subroutine localize_xy
!
    subroutine backskip_to_time(lun,lroot)
!
!  Skips over possible persistent data blocks from end of snapshot to time record.
!
!  9-mar-15/MR: coded
!
      integer,           intent(in) :: lun
      logical, optional, intent(in) :: lroot

      integer :: i,id
      integer, parameter :: id_block_PERSISTENT = 2000 

      backspace(lun)
      read(lun) id

      if (id==id_block_PERSISTENT) then
        backspace(lun)

        do
          do i=1,3; backspace(lun); enddo
          read(lun) id
          if (id==id_block_PERSISTENT) exit
        enddo
        backspace(lun)
      endif

      backspace(lun)
      if (loptest(lroot)) backspace(lun)

    endsubroutine backskip_to_time

    logical function loptest(lopt,ldef)
!
!  returns value of optional logical parameter opt if present,
!  otherwise the default value ldef, if present, .false. if not
!
!  20-aug-13/MR: coded
!  26-aug-13/MR: optional default value ldef added
!
      logical, optional, intent(in) :: lopt, ldef

      if (present(lopt)) then
        loptest=lopt
      else if (present(ldef)) then
        loptest=ldef
      else
        loptest=.false.
      endif

    endfunction loptest

!***********************************************************************
    subroutine cross(a,b,c)
      use real_prec
!
!  Cross product, c = a x b, for simple 3-d vectors (independent of position).
!
      real(rl), dimension (3) :: a,b,c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
!
    endsubroutine cross
!***********************************************************************
    subroutine dot(a,b,c)
      use real_prec
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  11-mar-04/wolf: coded
!
      real(rl), dimension (:) :: a,b
      real(rl) :: c
!
      intent(in) :: a,b
      intent(out) :: c
!
      c = dot_product(a,b)
!
    endsubroutine dot
!***********************************************************************
    subroutine dot2(a,b)
      use real_prec
!
!  Dot product, c=a.b, of two simple 3-d arrays.
!
!  11-mar-04/wolf: coded
!
      real(rl), dimension (:) :: a
      real(rl) :: b
!
      intent(in) :: a
      intent(out) :: b
!
      b = dot_product(a,a)
!
    endsubroutine dot2
!***********************************************************************
end module subs
