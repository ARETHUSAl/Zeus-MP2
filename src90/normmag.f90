!=======================================================================
!
!    \\\\\\\\\\      B E G I N   F U N C T I O N        //////////
!    //////////              N O R M M A G              \\\\\\\\\!
!=======================================================================
! normalizes magnetic - field to b_rms = rms
!
      function normmag(rms)
      use real_prec
      use param
      use field
      use bndry
      use grid
      use root
      use scratch
      use mpiyes
      use mpipar
      implicit none
      real(rl)  :: normmag
      real(rl)  :: norm, rms, av
      integer   :: i , j , k , ip, jp, kp
      norm=0.
      do 30 k=ks, ke
        kp = k+1
        do 20 j=js, je
          jp = j+1
          do 10 i=is, ie
            ip = i+1
            norm = norm + &
                   (b1(i,j,k)+b1(ip,j ,k ))**2. + &
                   (b2(i,j,k)+b2(i ,jp,k ))**2. + &
                   (b3(i,j,k)+b3(i ,j ,kp))**2.
10       continue
20      continue
30     continue
      norm = sqrt(norm*0.25/real(ijkn)**3/rms**2.)
      b1 = b1/norm
      b2 = b2/norm
      b3 = b3/norm
      av  = 0.0
      rms = 0.0
      do 60 k=ks, ke
        kp = k+1
        do 50 j=js, je
          jp = j+1
          do 40 i=is, ie
            ip = i+1
            av  = av  + ( b1(i,j,k) + b1(ip,j ,k ) &
                      +   b2(i,j,k) + b2(i ,jp,k ) &
                      +   b3(i,j,k) + b3(i, j ,kp) )*0.5
            rms = rms + ((b1(i,j,k) + b1(ip,j ,k ))**2 &
                      +  (b2(i,j,k) + b2(i ,jp,k ))**2 &
                      +  (b3(i,j,k) + b3(i ,j ,kp))**2)*0.25
40        continue
50      continue
60    continue
      av  = av/(3.0*real(ijkn)**3.0)
      rms = sqrt(rms/(real(ijkn)**3.0))
      if (myid_w .eq. 0) then
        write(6,*) 'NORMMAG  :'
        write(6,*) 'NORMMAG  : Norm     = ',norm
        write(6,*) 'NORMMAG  : Average  = ',av
        write(6,*) 'NORMMAG  : B_rms    = ',rms
        write(6,*) 'NORMMAG  :'
      endif
      normmag = norm
      return
      end                 
