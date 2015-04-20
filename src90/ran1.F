      function ran1(idum) 
      use real_prec
      implicit none                                           
c                                                                        
c    dac:zeus3d.ran1 <------------------- returns uniform random deviate 
c    from whp:numerical recipes                           december, 1992 
c                                                                        
c    written by: David Clarke                                            
c    modified 1:                                                         
c                                                                        
c    PURPOSE:  Returns a uniformly distrbuted random deviate between 0.0   
c    and 1.0 (see Numerical Recipes, 1st edition for FORTRAN, page 196).   
c                                                                        
c    INPUT VARIABLES:                                                      
c     idum     .ge. 0 => get next random number from sequence             
c              .lt. 0 => initialise sequence depending on absolute value  
c                        of "idum", and gets next random number           
c                                                                        
c    OUTPUT VARIABLES:                                                     
c      ran1     deviate with uniform distribution between 0.0 and 1.0      
c                                                                        
c    EXTERNALS: [NONE]                                                     
c                                                                        
c----------------------------------------------------------------------- 
c                                                                                                                                              
      real(rl) :: ran1
      integer, parameter   ::  imax=97                                        
      integer, parameter   ::  m1=259200, ia1=7141, ic1=54733   
      integer, parameter   ::  m2=134456, ia2=8121, ic2=28411 
      integer, parameter   ::  m3=243000, ia3=4561, ic3=51349    
      real(rl), parameter  ::  rm1=1.0/m1,  rm2=1.0/m2                              
      integer  :: idum, iff, ix1, ix2, ix3, i                                                   
      real(rl) :: r(imax)                                        
c                                                                       
c     Data statements                                                   
c                                                                       
      data          iff      / 0 /                                      
      save          iff     , ix1     , ix2     , ix3     , r           
c                                                                       
c---------------------------------------------------------------------- 
c                                                                       
c     If sequence is to be initialised, start here ...                  
c                                                                       
      if ( (idum .lt. 0) .or. (iff .eq. 0) ) then                       
        iff = 1                                                         
        ix1 = mod ( ic1 - idum   , m1 )                                 
        ix1 = mod ( ia1*ix1 + ic1, m1 )                                 
        ix2 = mod ( ix1          , m2 )                                 
        ix1 = mod ( ia1*ix1 + ic1, m1 )                                 
        ix3 = mod ( ix1          , m3 )                                 
        do 10 i=1,imax                                                  
          ix1  = mod ( ia1*ix1 + ic1, m1 )                              
          ix2  = mod ( ia2*ix2 + ic2, m2 )                              
          r(i) = (  real(ix1) +  real(ix2) * rm2 ) * rm1                
10      continue                                                        
        idum = 1                                                        
      endif                                                             
c                                                                       
c     Otherwise, start here ....                                        
c                                                                       
      ix1  = mod ( ia1*ix1 + ic1, m1 )                                  
      ix2  = mod ( ia2*ix2 + ic2, m2 )                                  
      ix3  = mod ( ia3*ix3 + ic3, m3 )                                  
      i    = 1 + ( imax * ix3 ) / m3                                    
      ran1 = r(i)                                                       
      r(i) = (  real(ix1) +  real(ix2) * rm2 ) * rm1                    
c                                                                       
      return                                                            
      end
