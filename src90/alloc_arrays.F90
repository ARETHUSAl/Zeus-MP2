!=======================================================================
!
!    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
!    //////////         A L L O C _ A R R A Y S           \\\\\\\\\\
!
!                            Developed by
!                Laboratory of Computational Astrophysics
!                 University of California at San Diego
!
!=======================================================================
      subroutine alloc_arrays
!
      use param
      use config
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use field
      use grid
      use bndry
      use lor_scr
      use scratch
      use radiation
      use opac
      use chem
!
      implicit none
!
!=======================================================================
!     determine local array extents
!=======================================================================
!
      in   = izones + 5
      jn   = jzones + 5
      kn   = kzones + 5
      ijkn = maxijk + 5
!
!=======================================================================
!     FIELD arrays
!=======================================================================
!
      allocate(d (in,jn,kn))
      allocate(e (in,jn,kn))
      allocate(v1(in,jn,kn))
      allocate(v2(in,jn,kn))
      allocate(v3(in,jn,kn))
      allocate(p (in,jn,kn))
      allocate(tt(in,jn,kn))
      if(xtotnrg) then
       allocate(q11(in,jn,kn))
       allocate(q22(in,jn,kn))
       allocate(q33(in,jn,kn))
      endif
      if(nspec .gt. 1) then 
       allocate(abun(in,jn,kn,nspec))
      else
       allocate(abun(in,jn,kn,1))
      endif
!
      if(xmhd) then
       allocate(b1(in,jn,kn))
       allocate(b2(in,jn,kn))
       allocate(b3(in,jn,kn))
!       allocate(ar(in,jn,kn))
!       allocate(ai(in,jn,kn))
      endif
!
      if(xhse) then
       allocate(gp1(in,jn,kn))
       allocate(gp2(in,jn,kn))
       allocate(gp3(in,jn,kn))
      endif
      if(lrad .ne. 0) then
       allocate(er   (in,jn,kn))
       allocate(en   (in,jn,kn))
       allocate(ern  (in,jn,kn))
       allocate(de   (in,jn,kn))
       allocate(der  (in,jn,kn))
       allocate(pn   (in,jn,kn))
       allocate(dr1  (in,jn,kn))
       allocate(dr2  (in,jn,kn))
       allocate(dr3  (in,jn,kn))
       allocate(dvl11(in,jn,kn))
       allocate(dvl22(in,jn,kn))
       allocate(dvl12(in,jn,kn))
       allocate(dvl21(in,jn,kn))
       allocate(divvl(in,jn,kn))
       allocate(dvl33(in,jn,kn))
       allocate(bb    (in,jn,kn)) 
       allocate(kap   (in,jn,kn))
       allocate(sig   (in,jn,kn)) 
       allocate(bbn   (in,jn,kn))
       allocate(kapn  (in,jn,kn)) 
       allocate(dpde  (in,jn,kn))
       allocate(dbbde (in,jn,kn))
       allocate(dkapde(in,jn,kn)) 
       allocate(dkapdt(in,jn,kn))
       allocate(kapr  (in,jn,kn)) 
       allocate(kem   (in,jn,kn))
       allocate(dkemde(in,jn,kn)) 
       allocate(dkemdt(in,jn,kn))
       allocate(kemn  (in,jn,kn))
      endif
!
      if(xgrav .or. xgrvfft) then
       allocate(gp   (in,jn,kn))
       allocate(ggp  (in,jn,kn))
       allocate(oldgp(in,jn,kn))
       if(xsphgrv) allocate(intm(in))
      endif
!
!=======================================================================
!     GRID arrays
!=======================================================================
!
      allocate(x1a    (in))
      allocate(x1ai   (in))
      allocate(dx1a   (in))
      allocate(dx1ai  (in))
      allocate(vol1a  (in))
      allocate(dvl1a  (in))
      allocate(dvl1ai (in))
      allocate(g2a    (in))
      allocate(g2ai   (in))
      allocate(g31a   (in))
      allocate(g31ai  (in))
      allocate(dg2ad1 (in))
      allocate(dg31ad1(in))
!
      allocate(x2a    (jn))
      allocate(x2ai   (jn))
      allocate(dx2a   (jn))
      allocate(dx2ai  (jn))
      allocate(vol2a  (jn))
      allocate(dvl2a  (jn))
      allocate(dvl2ai (jn))
      allocate(g32a   (jn))
      allocate(g32ai  (jn))
      allocate(dg32ad2(jn))
      allocate(g4a    (jn))
!
      allocate(x3a    (kn))
      allocate(x3ai   (kn))
      allocate(dx3a   (kn))
      allocate(dx3ai  (kn))
      allocate(vol3a  (kn))
      allocate(dvl3a  (kn))
      allocate(dvl3ai (kn))
!
      allocate(x1b    (in))
      allocate(x1bi   (in))
      allocate(dx1b   (in))
      allocate(dx1bi  (in))
      allocate(vol1b  (in))
      allocate(dvl1b  (in))
      allocate(dvl1bi (in))
      allocate(g2b    (in))
      allocate(g2bi   (in))
      allocate(g31b   (in))
      allocate(g31bi  (in))
      allocate(dg2bd1 (in))
      allocate(dg31bd1(in))
!
      allocate(x2b    (jn))
      allocate(x2bi   (jn))
      allocate(dx2b   (jn))
      allocate(dx2bi  (jn))
      allocate(vol2b  (jn))
      allocate(dvl2b  (jn))
      allocate(dvl2bi (jn))
      allocate(g32b   (jn))
      allocate(g32bi  (jn))
      allocate(dg32bd2(jn))
      allocate(g4b    (jn))
!
      allocate(x3b    (kn))
      allocate(x3bi   (kn))
      allocate(dx3b   (kn))
      allocate(dx3bi  (kn))
      allocate(vol3b  (kn))
      allocate(dvl3b  (kn))
      allocate(dvl3bi (kn))
!
      allocate(vg1(in))
      allocate(vg2(jn))
      allocate(vg3(kn))
!
       allocate(x1ah    (in))
       allocate(dx1ah   (in))
       allocate(dvl1ah  (in))
       allocate(g2ah    (in))
       allocate(g31ah   (in))
       allocate(x1an    (in))
       allocate(dx1an   (in))
       allocate(dvl1an  (in))
       allocate(g2an    (in))
       allocate(g31an   (in))
       allocate(x1ahi   (in))
       allocate(dx1ahi  (in))
       allocate(dvl1ahi (in))
       allocate(g2ahi   (in))
       allocate(g31ahi  (in))
       allocate(x1ani   (in))
       allocate(dx1ani  (in))
       allocate(dvl1ani (in))
       allocate(g2ani   (in))
       allocate(g31ani  (in))
!
       allocate(x1bh    (in))
       allocate(dx1bh   (in))
       allocate(dvl1bh  (in))
       allocate(g2bh    (in))
       allocate(g31bh   (in))
       allocate(x1bn    (in))
       allocate(dx1bn   (in))
       allocate(dvl1bn  (in))
       allocate(g2bn    (in))
       allocate(g31bn   (in))
       allocate(x1bhi   (in))
       allocate(dx1bhi  (in))
       allocate(dvl1bhi (in))
       allocate(g2bhi   (in))
       allocate(g31bhi  (in))
       allocate(x1bni   (in))
       allocate(dx1bni  (in))
       allocate(dvl1bni (in))
       allocate(g2bni   (in))
       allocate(g31bni  (in))
!
       allocate(x2ah    (jn))
       allocate(dx2ah   (jn))
       allocate(dvl2ah  (jn))
       allocate(g32ah   (jn))
       allocate(g4ah    (jn))
       allocate(x2an    (jn))
       allocate(dx2an   (jn))
       allocate(dvl2an  (jn))
       allocate(g32an   (jn))
       allocate(g4an    (jn))
       allocate(x2ahi   (jn))
       allocate(dx2ahi  (jn))
       allocate(dvl2ahi (jn))
       allocate(g32ahi  (jn))
       allocate(x2ani   (jn))
       allocate(dx2ani  (jn))
       allocate(dvl2ani (jn))
       allocate(g32ani  (jn))
!
       allocate(x2bh    (jn))
       allocate(dx2bh   (jn))
       allocate(dvl2bh  (jn))
       allocate(g32bh   (jn))
       allocate(g4bh    (jn))
       allocate(x2bn    (jn))
       allocate(dx2bn   (jn))
       allocate(dvl2bn  (jn))
       allocate(g32bn   (jn))
       allocate(g4bn    (jn))
       allocate(x2bhi   (jn))
       allocate(dx2bhi  (jn))
       allocate(dvl2bhi (jn))
       allocate(g32bhi  (jn))
       allocate(x2bni   (jn))
       allocate(dx2bni  (jn))
       allocate(dvl2bni (jn))
       allocate(g32bni  (jn))
!
       allocate(x3ah    (kn))
       allocate(dx3ah   (kn))
       allocate(dvl3ah  (kn))
       allocate(x3an    (kn))
       allocate(dx3an   (kn))
       allocate(dvl3an  (kn))
       allocate(x3ahi   (kn))
       allocate(dx3ahi  (kn))
       allocate(dvl3ahi (kn))
       allocate(x3ani   (kn))
       allocate(dx3ani  (kn))
       allocate(dvl3ani (kn))
!
       allocate(x3bh    (kn))
       allocate(dx3bh   (kn))
       allocate(dvl3bh  (kn))
       allocate(x3bn    (kn))
       allocate(dx3bn   (kn))
       allocate(dvl3bn  (kn))
       allocate(x3bhi   (kn))
       allocate(dx3bhi  (kn))
       allocate(dvl3bhi (kn))
       allocate(x3bni   (kn))
       allocate(dx3bni  (kn))
       allocate(dvl3bni (kn))
!
!=======================================================================
!     BNDRY arrays
!=======================================================================
!
      allocate(niib(jn,kn))
      allocate(noib(jn,kn))
      allocate(nijb(in,kn))
      allocate(nojb(in,kn))
      allocate(nikb(in,jn))
      allocate(nokb(in,jn))
!
       allocate(niib2 (jn,kn))
       allocate(niib3 (jn,kn))
       allocate(niib23(jn,kn))
       allocate(noib2 (jn,kn))
       allocate(noib3 (jn,kn))
       allocate(noib23(jn,kn))
       allocate(nijb3 (in,kn))
       allocate(nijb1 (in,kn))
       allocate(nijb31(in,kn))
       allocate(nojb3 (in,kn))
       allocate(nojb1 (in,kn))
       allocate(nojb31(in,kn))
       allocate(nikb1 (in,jn))
       allocate(nikb2 (in,jn))
       allocate(nikb12(in,jn))
       allocate(nokb1 (in,jn))
       allocate(nokb2 (in,jn))
       allocate(nokb12(in,jn))
!
      if(lrad .ne. 0) then
       allocate(liib(jn,kn))
       allocate(loib(jn,kn))
       allocate(lijb(in,kn))
       allocate(lojb(in,kn))
       allocate(likb(in,jn))
       allocate(lokb(in,jn))
      endif
!
      allocate(diib(jn,kn,3))
      allocate(doib(jn,kn,3))
      allocate(dijb(in,kn,3))
      allocate(dojb(in,kn,3))
      allocate(dikb(in,jn,3))
      allocate(dokb(in,jn,3))
!
      allocate(eiib(jn,kn,2))
      allocate(eoib(jn,kn,2))
      allocate(eijb(in,kn,2))
      allocate(eojb(in,kn,2))
      allocate(eikb(in,jn,2))
      allocate(eokb(in,jn,2))
!
      allocate(v1iib(jn,kn,2))
      allocate(v1oib(jn,kn,2))
      allocate(v1ijb(in,kn,2))
      allocate(v1ojb(in,kn,2))
      allocate(v1ikb(in,jn,2))
      allocate(v1okb(in,jn,2))
      allocate(v2iib(jn,kn,2))
      allocate(v2oib(jn,kn,2))
      allocate(v2ijb(in,kn,2))
      allocate(v2ojb(in,kn,2))
      allocate(v2ikb(in,jn,2))
      allocate(v2okb(in,jn,2))
      allocate(v3iib(jn,kn,2))
      allocate(v3oib(jn,kn,2))
      allocate(v3ijb(in,kn,2))
      allocate(v3ojb(in,kn,2))
      allocate(v3ikb(in,jn,2))
      allocate(v3okb(in,jn,2))
!
      if(nspec .gt. 1) then
       allocate(abiib(jn,kn,3,nspec))
       allocate(aboib(jn,kn,3,nspec))
       allocate(abijb(in,kn,3,nspec))
       allocate(abojb(in,kn,3,nspec))
       allocate(abikb(in,jn,3,nspec))
       allocate(abokb(in,jn,3,nspec))
      endif ! nspec
!
      if(xmhd) then
       allocate(b1iib(jn,kn,2))
       allocate(b1oib(jn,kn,2))
       allocate(b1ijb(in,kn,2))
       allocate(b1ojb(in,kn,2))
       allocate(b1ikb(in,jn,2))
       allocate(b1okb(in,jn,2))
       allocate(b2iib(jn,kn,2))
       allocate(b2oib(jn,kn,2))
       allocate(b2ijb(in,kn,2))
       allocate(b2ojb(in,kn,2))
       allocate(b2ikb(in,jn,2))
       allocate(b2okb(in,jn,2))
       allocate(b3iib(jn,kn,2))
       allocate(b3oib(jn,kn,2))
       allocate(b3ijb(in,kn,2))
       allocate(b3ojb(in,kn,2))
       allocate(b3ikb(in,jn,2))
       allocate(b3okb(in,jn,2))
!
       allocate(emf1iib(jn,kn,3))
       allocate(emf1oib(jn,kn,3))
       allocate(emf1ijb(in,kn,3))
       allocate(emf1ojb(in,kn,3))
       allocate(emf1ikb(in,jn,3))
       allocate(emf1okb(in,jn,3))
       allocate(emf2iib(jn,kn,3))
       allocate(emf2oib(jn,kn,3))
       allocate(emf2ijb(in,kn,3))
       allocate(emf2ojb(in,kn,3))
       allocate(emf2ikb(in,jn,3))
       allocate(emf2okb(in,jn,3))
       allocate(emf3iib(jn,kn,3))
       allocate(emf3oib(jn,kn,3))
       allocate(emf3ijb(in,kn,3))
       allocate(emf3ojb(in,kn,3))
       allocate(emf3ikb(in,jn,3))
       allocate(emf3okb(in,jn,3))
      endif ! xmhd
!
      if(lrad .ne. 0) then
       allocate(eriib(jn,kn,2))
       allocate(eroib(jn,kn,2))
       allocate(erijb(in,kn,2))
       allocate(erojb(in,kn,2))
       allocate(erikb(in,jn,2))
       allocate(erokb(in,jn,2))
      endif ! lrad
!
      if(xgrav .or. xgrvfft) then
       allocate(gpiib(jn,kn,2))
       allocate(gpoib(jn,kn,2))
       allocate(gpijb(in,kn,2))
       allocate(gpojb(in,kn,2))
       allocate(gpikb(in,jn,2))
       allocate(gpokb(in,jn,2))
      endif ! xgrav
!
!=======================================================================
!     SCRATCH arrays
!=======================================================================
!
      allocate(w1da(ijkn))
      allocate(w1db(ijkn))
      allocate(w1dc(ijkn))
      allocate(w1dd(ijkn))
      allocate(w1de(ijkn))
      allocate(w1df(ijkn))
      allocate(w1dg(ijkn))
      allocate(w1dh(ijkn))
      allocate(w1di(ijkn))
      allocate(w1dj(ijkn))
      allocate(w1dk(ijkn))
      allocate(w1dl(ijkn))
      allocate(w1dm(ijkn))
      allocate(w1dn(ijkn))
      allocate(w1do(ijkn))
      allocate(w1dp(ijkn))
      allocate(w1dq(ijkn))
      allocate(w1dr(ijkn))
      allocate(w1ds(ijkn))
      allocate(w1dt(ijkn))
      allocate(w1du(ijkn))
!
      allocate(w3da(in,jn,kn))
      allocate(w3db(in,jn,kn))
      allocate(w3dc(in,jn,kn))
      allocate(w3dd(in,jn,kn))
      allocate(w3de(in,jn,kn))
      allocate(w3df(in,jn,kn))
      allocate(w3dg(in,jn,kn))
      allocate(w3dh(in,jn,kn))
      allocate(w3di(in,jn,kn))
      allocate(w3dj(in,jn,kn))
      if(nspec .gt. 1) allocate(w4da(in,jn,kn,nspec))
!
!=======================================================================
!     LOR_SCR arrays
!=======================================================================
!
      if(xmhd) then
       allocate(srd1(in,jn,kn))
       allocate(srd2(in,jn,kn))
       allocate(srd3(in,jn,kn))
      endif
!
!
!=======================================================================
!     CHEM/RT arrays
!=======================================================================
!
      allocate(DMa      (nratec))
      allocate(DM       (in, jn))
      if (nspec .gt. 1) then
       allocate(k1a     (nratec))
       allocate(k2a     (nratec))
       allocate(cehia   (nratec))
       allocate(cihia   (nratec))
       allocate(rehiia  (nratec))
       allocate(brema   (nratec))
       allocate(compt   (nratec))
       allocate(k1      (in, jn))
       allocate(k2      (in, jn))
       allocate(ceHI    (in, jn))
       allocate(ciHI    (in, jn))
       allocate(reHII   (in, jn))
       allocate(brem    (in, jn))
       allocate(cmpt    (in, jn))
       if (iRT .gt. 0) then
        allocate(k24    (in, jn))
        allocate(k24mv  (in, jn))
        allocate(piHI   (in, jn))
        allocate(nu      (nnu1+nnu2))
        allocate(sigma24 (nnu1+nnu2))
        allocate(fcentral(nnu1+nnu2))
        allocate(f_outer (nnu1+nnu2))
        allocate(f(in,jn, nnu1+nnu2))
        allocate(dv_ioniz(in,jn,kn ))
       endif
       allocate(tttt    (in, jn))
       allocate(indixe  (in, jn))
       allocate(logtem  (in, jn))
       allocate(dedot   (in, jn))
       allocate(edot  (in,jn,kn))
       allocate(tgas  (in,jn,kn))
      endif
      if (nspec .gt. 3) then
       allocate(k3a     (nratec))
       allocate(k4a     (nratec))
       allocate(k5a     (nratec))
       allocate(k6a     (nratec))
       allocate(ceheia  (nratec))
       allocate(ciheia  (nratec))
       allocate(ciheisa (nratec))
       allocate(ceheiia (nratec))
       allocate(ciheiia (nratec))
       allocate(reheii1a(nratec))
       allocate(reheii2a(nratec))
       allocate(reheiiia(nratec))
       allocate(k3      (in, jn))
       allocate(k4      (in, jn))
       allocate(k5      (in, jn))
       allocate(k6      (in, jn))
       allocate(ceHeI   (in, jn))
       allocate(ciHeI   (in, jn))
       allocate(ciHeIS  (in, jn))
       allocate(ceHeII  (in, jn))
       allocate(ciHeII  (in, jn))
       allocate(reHeII1 (in, jn))
       allocate(reHeII2 (in, jn))
       allocate(reHeIII (in, jn))
       if (iRT .gt. 0) then
         allocate(k25   (in, jn))
         allocate(k25mv (in, jn))
         allocate(k26   (in, jn))
         allocate(k26mv (in, jn))
         allocate(piHeI (in ,jn))
         allocate(piHeII(in ,jn))
         allocate(sigma25(nnu1+nnu2))
         allocate(sigma26(nnu1+nnu2))
       endif
      endif
      if (nspec .gt. 6) then
       allocate(k7a     (nratec))
       allocate(k8a     (nratec))
       allocate(k9a     (nratec))
       allocate(k10a    (nratec))
       allocate(k11a    (nratec))
       allocate(k12a    (nratec))
       allocate(k13a    (nratec))
       allocate(k14a    (nratec))
       allocate(k15a    (nratec))
       allocate(k16a    (nratec))
       allocate(k17a    (nratec))
       allocate(k18a    (nratec))
       allocate(k19a    (nratec))
       allocate(k22a    (nratec))
       allocate(vibha   (nratec))
       allocate(hyd01ka (nratec))
       allocate(h2k01a  (nratec))
       allocate(rotha   (nratec))
       allocate(rotla   (nratec))
       allocate(gpldla  (nratec))
       allocate(gphdla  (nratec))
       allocate(gaHIa   (nratec))
       allocate(gaH2a   (nratec))
       allocate(gaHea   (nratec))
       allocate(gaHpa   (nratec))
       allocate(gaela   (nratec))
       allocate(k7      (in, jn))
       allocate(k8      (in, jn))
       allocate(k9      (in, jn))
       allocate(k10     (in, jn))
       allocate(k11     (in, jn))
       allocate(k12     (in, jn))
       allocate(k13     (in, jn))
       allocate(k14     (in, jn))
       allocate(k15     (in, jn))
       allocate(k16     (in, jn))
       allocate(k17     (in, jn))
       allocate(k18     (in, jn))
       allocate(k19     (in, jn))
       allocate(k22     (in, jn))
       allocate(hyd01k  (in, jn))
       allocate(h2k01   (in, jn))
       allocate(roth    (in, jn))
       allocate(rotl    (in, jn))
       allocate(vibh    (in, jn))
       allocate(gpldl   (in, jn))
       allocate(gphdl   (in, jn))
       allocate(gaHI    (in, jn))
       allocate(gaH2    (in, jn))
       allocate(gaHe    (in, jn))
       allocate(gaHp    (in, jn))
       allocate(gael    (in, jn))
       allocate(galdl   (in, jn))
       if (iRT .gt. 0) then
         allocate(k27       (in, jn))
         allocate(k28       (in, jn))
         allocate(k29       (in, jn))
         allocate(k30       (in, jn))
         allocate(k31       (in, jn))
         allocate(sigma27(nnu1+nnu2))
         allocate(sigma28(nnu1+nnu2))
         allocate(sigma29(nnu1+nnu2))
         allocate(sigma30(nnu1+nnu2))
         allocate(sigma31(nnu1+nnu2))
       endif
      endif
      if (nspec .gt. 9) then
       allocate(k50a    (nratec))
       allocate(k51a    (nratec))
       allocate(k52a    (nratec))
       allocate(k53a    (nratec))
       allocate(k54a    (nratec))
       allocate(k55a    (nratec))
       allocate(k56a    (nratec))
       allocate(hdltea  (nratec))
       allocate(hdlowa  (nratec))
       allocate(k50     (in, jn))
       allocate(k51     (in, jn))
       allocate(k52     (in, jn))
       allocate(k53     (in, jn))
       allocate(k54     (in, jn))
       allocate(k55     (in, jn))
       allocate(k56     (in, jn))
       allocate(hdlte   (in, jn))
       allocate(hdlow   (in, jn))
      endif
      return
      end
