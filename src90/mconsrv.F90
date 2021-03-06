#include "rtchem.def"
!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////               M C O N S R V               \\\\\\\\\\
!
!=======================================================================
!
      subroutine mconsrv
!
!
!  written by: Daniel Whalen 02.15.01
!
!  upgraded 09.02.06 by DJW to include HD chemistry and 3-body reactions
!
!  ported to ZEUS-MP 2.1 by DJW 11.29.06
!
!  PURPOSE:    ensures that the multiple species densities add up to the
!              total density d after the transport and chemistry steps;
!              also enforces charge conservation
!              
!  LOCAL VARIABLES:
!
!  EXTERNALS:  none
!
!  The advection and chemistry substeps compute delta(rho) for the total
!  density and individual species densities but don't guarantee that the
!  changes in the individual densites sum to equal the change in the 
!  total density.  MCONSRV first calculates the sum of the hydrogen 
!  densities and compares it to fh*(total density), where f is the 
!  hydrogen fraction of the baryon density.  The difference between the
!  two is the error and it is added to the larger of HI and HII to 
!  enforce H mass conservation.  MCONSRV then compares the sum of the 
!  He densities and compares this to (1-fh)*(total density), the He 
!  fraction of the total baryon density.  The error is lumped in with
!  the largest of the He densities.  In each case we add the error to 
!  the largest of the densities because it can be the same size as the
!  smallest densities and would swamp them if it were added to them.  
!  We call MCONSRV after the chemistry substep.
!
!
!-----------------------------------------------------------------------
!
      use real_prec
      use param
      use cons
      use chem
      use field
      use grid
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
!
      implicit NONE
       integer  ::  i, j, k       
       real(rl) ::  err,dh
!-----------------------------------------------------------------------
!
       dh = 1.0d-6
       do 30 k=ks,ke
       do 20 j=js,je
       do 10 i=is,ie
!      None of the hydrogen densities can exceed fh * (total baryon
!      density).  If HI or HII has somehow become greater than the
!      total hydrogen baryon density, knock it back down to just below
!      this density.  Do the same if any of the He densities exceed
!      (1-fh)*(total density).  Also, n_e,max = n_H + 2*n_He; this leads
!      to 
!
!                m_H * n_e,max = rho_H + 0.5 * (4*m_H*n_He)
! 
!                       de_max = rho_H + 0.5 * rho_He
!
!                       de_max = fh * d + 0.5 * (1-fh) * d 
!
!                              = 0.5 * (1+fh) * d
!
!      so if de exceeds 0.5 * (1+fh) * d we decrease it to just below
!      0.5 * (1+fh) * d
!
!      We're not in danger of having any of the H2 or HD species 
!      approaching the fh density fraction under normal circumstances,
!      but we should make provision for all H to be converted to H2
!      if 3-body production becomes dominant
#ifdef H
       if (abun(i,j,k,2)/d(i,j,k) .ge.         fh ) &
 &
         abun(i,j,k,2) =        fh *d(i,j,k)-tiny
       if (abun(i,j,k,1)/d(i,j,k) .ge.         fh ) &
 &
         abun(i,j,k,1) =        fh *d(i,j,k)-tiny
#endif /* H */
#if defined H || defined He
       if (abun(i,j,k,3)/d(i,j,k) .ge. 0.5*(1.+fh)) &
 &
         abun(i,j,k,3) =0.5*(1.+fh)*d(i,j,k)-tiny
#endif /* H or He */
#ifdef He
       if (abun(i,j,k,4)/d(i,j,k) .ge.      1.-fh ) &
 &
         abun(i,j,k,4) =    (1.-fh)*d(i,j,k)-tiny
       if (abun(i,j,k,5)/d(i,j,k) .ge.      1.-fh ) &
 &
         abun(i,j,k,5) =    (1.-fh)*d(i,j,k)-tiny
       if (abun(i,j,k,6)/d(i,j,k) .ge.      1.-fh ) &
 &
         abun(i,j,k,6) =    (1.-fh)*d(i,j,k)-tiny
#endif /* He */
#ifdef H2
       if (abun(i,j,k,8)/d(i,j,k) .ge.          fh ) &
 &
         abun(i,j,k,8) =        fh *d(i,j,k)-tiny
#endif /* H2 */
#ifdef HD
       if (abun(i,j,k,11)/d(i,j,k) .ge.         dh ) &
 &
         abun(i,j,k,11)=        dh *d(i,j,k)-tiny
       if (abun(i,j,k,10)/d(i,j,k) .ge.         dh ) &
 &
         abun(i,j,k,10)=        dh *d(i,j,k)-tiny
       if (abun(i,j,k,12)/d(i,j,k) .ge. 2.*thd* dh ) &
 &
         abun(i,j,k,12)= 2.*thd*dh *d(i,j,k)-tiny
#endif /* HD */
#ifdef H
       if (abun(i,j,k,2) .lt. tiny) abun(i,j,k,2)   = tiny
       if (abun(i,j,k,1) .lt. tiny) abun(i,j,k,1)   = tiny
#endif /* H */
#if defined H || defined He
       if (abun(i,j,k,3) .lt. tiny) abun(i,j,k,3 )  = tiny
#endif /* H or He */
#ifdef He
       if (abun(i,j,k,4) .lt. tiny) abun(i,j,k,4)   = tiny
       if (abun(i,j,k,5) .lt. tiny) abun(i,j,k,5)   = tiny
       if (abun(i,j,k,6) .lt. tiny) abun(i,j,k,6)   = tiny
#endif /* He */
#ifdef H2
       if (abun(i,j,k,7) .lt. tiny) abun(i,j,k,7)   = tiny
       if (abun(i,j,k,8) .lt. tiny) abun(i,j,k,8)   = tiny
       if (abun(i,j,k,9) .lt. tiny) abun(i,j,k,9)   = tiny
#endif /* H2 */
#ifdef HD
       if (abun(i,j,k,10) .lt. tiny) abun(i,j,k,10) = tiny
       if (abun(i,j,k,11) .lt. tiny) abun(i,j,k,11) = tiny
       if (abun(i,j,k,12) .lt. tiny) abun(i,j,k,12) = tiny
#endif /* HD */
!
!       Now enforce density conservation
!
#ifdef H
        err = fh*d(i,j,k)   - abun(i,j,k,2) - abun(i,j,k,1)
#ifdef H2 &
            - abun(i,j,k,7) - abun(i,j,k,8) - abun(i,j,k,9)
#endif /* H2 */
        if    (abun(i,j,k,1) .ge. abun(i,j,k,2) 
#ifdef H2 &
       .and.   abun(i,j,k,1) .ge. abun(i,j,k,8)
#endif /* H2 */ &
                                          ) then
          abun(i,j,k,1)  =   abun(i,j,k,1) + err
        endif
        if    (abun(i,j,k,2) .ge. abun(i,j,k,1)
#ifdef H2 &
       .and.   abun(i,j,k,2) .ge. abun(i,j,k,8)
#endif /* H2 */ &
                                          ) then
          abun(i,j,k,2)  =   abun(i,j,k,2) + err
        endif
#ifdef H2
        if (abun(i,j,k,8) .ge. abun(i,j,k,1)  .and. &
            abun(i,j,k,8) .ge. abun(i,j,k,2)) then
          abun(i,j,k,8)  =   abun(i,j,k,8) + err
        endif
#endif /* H2 */
#endif /* H */
#ifdef He
        err = (1.0-fh) * d(i,j,k) - abun(i,j,k,4) - abun(i,j,k,5) &
                                                  - abun(i,j,k,6)
        if(abun(i,j,k,4) .ge. abun(i,j,k,6)) &
           abun(i,j,k,4) = abun(i,j,k,4) + err
        if(abun(i,j,k,6) .ge. abun(i,j,k,4)) &
           abun(i,j,k,6) = abun(i,j,k,6) + err
#endif /* He */
#ifdef HD
        err = dh*d(i,j,k) - abun(i,j,k,11) - abun(i,j,k,10) &
                          -   2. *  thd    * abun(i,j,k,12)
        if (abun(i,j,k,10) .ge. abun(i,j,k,11)  .and. &
            abun(i,j,k,10 ).ge. abun(i,j,k,12)) then
           abun(i,j,k,10) = abun(i,j,k,10) + err
        else if (abun(i,j,k,11) .ge. abun(i,j,k,10)  .and.
                 abun(i,j,k,11) .ge. abun(i,j,k,12)) then
           abun(i,j,k,11) = abun(i,j,k,11) + err
        else if (abun(i,j,k,12) .ge. abun(i,j,k,10)  .and. &
                 abun(i,j,k,12) .ge. abun(i,j,k,11)) then
           abun(i,j,k,12) = abun(i,j,k,12) + err
        endif
#endif /* HD */
!
!       Enforce charge conservation
!
        abun(i,j,k,3) = 
#ifdef H &
                        abun(i,j,k,2) 
#endif /* H */
#ifdef He &
                      + abun(i,j,k,5)/4. + abun(i,j,k,6)/2.
#endif /* He */
#ifdef H2 &
                      - abun(i,j,k,7)    + abun(i,j,k,9)/2.
#endif /* H2 */
#ifdef HD &
                      + abun(i,j,k,11)/2.
#endif /* HD */
!
!       If any of the densities went negative make it positive
!
#ifdef H
        abun(i,j,k,2) = dabs(abun(i,j,k,2))
        abun(i,j,k,1) = dabs(abun(i,j,k,1))
#endif /* H */
#if defined H || defined He
        abun(i,j,k,3) = dabs(abun(i,j,k,3))
#endif /* H or He */
#ifdef He
        abun(i,j,k,4) = dabs(abun(i,j,k,4))
        abun(i,j,k,5) = dabs(abun(i,j,k,5))
        abun(i,j,k,6) = dabs(abun(i,j,k,6))
#endif /* He */
#ifdef H2
        abun(i,j,k,7) = dabs(abun(i,j,k,7))
        abun(i,j,k,8) = dabs(abun(i,j,k,8))
        abun(i,j,k,9) = dabs(abun(i,j,k,9))
#endif /* H2 */
#ifdef HD
        abun(i,j,k,10) = dabs(abun(i,j,k,10))
        abun(i,j,k,11) = dabs(abun(i,j,k,11))
        abun(i,j,k,12) = dabs(abun(i,j,k,12))
#endif /* HD */
!
!       If any of the densities go below 'tiny' set it to 'tiny'
!
#ifdef H
        if (abun(i,j,k,1) .lt. tiny) abun(i,j,k,1) = tiny
        if (abun(i,j,k,2) .lt. tiny) abun(i,j,k,2) = tiny
#endif /* H */
#if defined H || defined He
        if (abun(i,j,k,3) .lt. tiny) abun(i,j,k,3) = tiny
#endif /* H or He */
#ifdef He
        if (abun(i,j,k,4) .lt. tiny) abun(i,j,k,4) = tiny
        if (abun(i,j,k,5) .lt. tiny) abun(i,j,k,5) = tiny
        if (abun(i,j,k,6) .lt. tiny) abun(i,j,k,6) = tiny
#endif /* He */
#ifdef H2
        if (abun(i,j,k,7) .lt. tiny) abun(i,j,k,7) = tiny
        if (abun(i,j,k,8) .lt. tiny) abun(i,j,k,8) = tiny
        if (abun(i,j,k,9) .lt. tiny) abun(i,j,k,9) = tiny
#endif /* H2 */
#ifdef HD
        if (abun(i,j,k,10) .lt. tiny) abun(i,j,k,10) = tiny
        if (abun(i,j,k,11) .lt. tiny) abun(i,j,k,11) = tiny
        if (abun(i,j,k,12) .lt. tiny) abun(i,j,k,12) = tiny
#endif /* HD */
10     continue
20     continue
30     continue
       return
       end
!
!=======================================================================
!
!    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
!    //////////               M C O N S R V               \\\\\\\\\\
!
!=======================================================================
!
