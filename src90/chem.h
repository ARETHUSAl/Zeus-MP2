      integer     ibrco, iceco       , icico       , ireco       ,
     .            ipiht, ih2co       , itmcool     , itmrate     ,
     .            iter , indixe(in,jn)             , i_start     ,
     .            icmpt, iDMco       , icycle      , iRT         ,
     .            ifreq, ispct       , iOTS        , nphbins     ,
     .            iLW  , iHM         , iPWA        , iHDco

      REAL temstart    , temend      , dlogtem     , tgas(in,jn,kn),
     .     logtem(in,jn)             , dedot(in,jn), rdshft        ,
     .     t_cutoff    , t_shutoff   , nphdot      , ephot         ,
     .     omega       , r_sep       , bbnm2       , t_star        ,
     .     hnumax      , alpha       , fh          , qsonm         ,
     .     nu(nnu1+nnu2), delnu      , nLWdot      , L_star        ,  
     .     xc          , jnu         , bbnm1       , nphdot1
#ifdef H
     .,    k1a (nratec), k2a (nratec), k4a (nratec)
#endif /* H */
#ifdef He
     .,    k3a (nratec), k5a (nratec), k6a (nratec) 
#endif /* He */
#ifdef H2
     .,    k7a (nratec), k8a (nratec), k9a (nratec), k10a(nratec), 
     .     k11a(nratec), k12a(nratec), k13a(nratec), k14a(nratec),
     .     k15a(nratec), k16a(nratec), k17a(nratec), k18a(nratec),
     .     k19a(nratec), k22a(nratec) 
#endif /* H2 */
#ifdef HD
     .,    k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     .     k54a(nratec), k55a(nratec), k56a(nratec) 
#endif /* HD */

      REAL brema   (nratec), compt   (nratec)
#ifdef H
     .,    cehia   (nratec), cihia   (nratec), rehiia  (nratec),
     .     DMa     (nratec), edot    (in,jn,kn)  
#endif /* H */
#ifdef He
     .,    ceheia  (nratec), ciheia  (nratec), ciheisa (nratec),
     .     ceheiia (nratec), ciheiia (nratec), reheii1a(nratec),
     .     reheii2a(nratec), reheiiia(nratec)
#endif /* He */
#ifdef H2
     .,    vibha   (nratec), hyd01ka (nratec), h2k01a  (nratec),
     .     rotha   (nratec), rotla   (nratec),
     .     gpldla  (nratec), gphdla  (nratec), gphdl1
#endif /* H2 */
#ifdef HD
     .,    hdltea  (nratec), hdlowa  (nratec)
#endif /* HD */

      common/coolblkc/ ibrco   , iceco   , icico   , ireco    , ipiht  , 
     .       ih2co   , i_start , tgas    , temstart, temend   , dlogtem, 
     .       logtem  , indixe  , dedot   , icmpt   , rdshft   , iDMco  ,
     .       icycle  , t_cutoff, iRT     , ifreq   , ispct    , ephot  , 
     .       omega   , r_sep   , bbnm1   , t_star  , hnumax   , fh     ,
     .       nphdot  , alpha   , iOTS    , qsonm   , t_shutoff, delnu  ,
     .       nu      , iPWA    , iLW     , iHM     , nLWdot   , L_star , 
     .       xc      , jnu     , bbnm2   , nphdot1 
#ifdef H
     .,      k1a     , k2a     , k4a     
#endif /* H */
#ifdef He
     .,      k3a     , k5a     , k6a       
#endif /* He */
#ifdef H2
     .,      k7a     , k8a     , k9a     , k10a     , k11a    , k12a   ,
     .       k13a    , k14a    , k15a    , k16a     , k17a    , k18a   , 
     .       k19a    , k22a      
#endif /* H2 */
#ifdef HD
     .,      k50a    , k51a    , k52a    , k53a     , k54a    , k55a   , 
     .       k56a          
#endif /* HD */
#ifdef H
     .,      cehia   , cihia   , rehiia  , DMa
#endif /* H */
#ifdef He
     .,      ceheia  , ciheia  , ciheisa ,
     .       ceheiia , ciheiia , reheii1a, reheii2a , reheiiia,
     .       edot 
#endif /* He */
#ifdef H2
     .,      vibha   , hyd01ka , h2k01a  , rotha    , rotla   ,
     .       gpldla  , gphdla  , gphdl1
#endif /* H2 */
#ifdef HD
     .,      hdltea  , hdlowa  
#endif /* HD */
     .,      brema   , itmcool , itmrate , iter     , compt

      REAL tttt   (in,jn)
#ifdef H
     .,    k1     (in,jn), k2     (in,jn)
#ifdef RT   
     .,    k24    (in,jn)
#endif /* RT */  
#endif /* H */
#ifdef He
     .,    k3     (in,jn), k4     (in,jn), k5     (in,jn), 
     .     k6     (in,jn)
#ifdef RT   
     .,    k25    (in,jn), k26    (in,jn)
#endif /* RT */   
#endif /* He */
#ifdef H2
     .,    k7     (in,jn), k8     (in,jn), k9     (in,jn),
     .     k10    (in,jn), k11    (in,jn), k12    (in,jn),
     .     k13    (in,jn), k14    (in,jn), k15    (in,jn),
     .     k16    (in,jn), k17    (in,jn), k18    (in,jn),
     .     k19    (in,jn), k22    (in,jn)
#ifdef RT   
     .,    k27    (in,jn), k28    (in,jn), k29    (in,jn),
     .     k30    (in,jn), k31    (in,jn)
#endif /* RT */   
#endif /* H2 */
#ifdef HD
     .,    k50    (in,jn), k51    (in,jn), k52    (in,jn), 
     .     k53    (in,jn), k54    (in,jn), k55    (in,jn), 
     .     k56    (in,jn) 
#endif /* HD */
#ifdef H
     .,    ceHI   (in,jn), ciHI   (in,jn), reHII  (in,jn) 
     .,    DM     (in,jn)
#ifdef RT   
     .,    piHI   (in,jn), sigma24(nnu1+nnu2)
#endif /* RT */   
#endif /* H */
#ifdef He
     .,    ceHeI  (in,jn), ciHeI  (in,jn), ciHeIS (in,jn),  
     .     ceHeII (in,jn), ciHeII (in,jn), reHeII1(in,jn), 
     .     reHeII2(in,jn), reHeIII(in,jn)
#ifdef RT   
     .,    piHeI  (in,jn), piHeII (in,jn), sigma25(nnu1+nnu2), 
     .     sigma26(nnu1+nnu2) 
#endif /* RT */   
#endif /* He */
#ifdef H2
     .,    hyd01k (in,jn), h2k01  (in,jn), roth   (in,jn),
     .     rotl   (in,jn), vibh   (in,jn),
     .     gpldl  (in,jn), gphdl  (in,jn)
#ifdef RT 
     .,    sigma27(nnu1+nnu2), sigma28(nnu1+nnu2), sigma29(nnu1+nnu2), 
     .     sigma30(nnu1+nnu2), sigma31(nnu1+nnu2), HMnm          (61)
#endif /* RT */   
#endif /* H2 */
#ifdef HD
     .,     hdlte  (in,jn), hdlow  (in,jn)
#endif /* HD */

     .,    brem   (in,jn), cmpt   (in,jn)
#ifdef RT 
     .,    fcentral(nnu1+nnu2), f(in,jn,nnu1+nnu2), f_outer(nnu1+nnu2)
c     .,    rad_press(in,jn,kn) 
#endif /* RT */   

      common/rate/    tttt     
#ifdef H
     .,      k1     , k2     
#ifdef RT 
     .,      k24    
#endif /* RT */   
#endif /* H */
#ifdef He
     .,      k3     , k4     , k5     , k6     
#ifdef RT 
     .,      k25    , k26    
#endif /* RT */   
#endif /* He */
#ifdef H2
     .,      k7     , k8     , k9     , k10    , k11    , k12    ,
     .       k13    , k14    , k15    , k16    , k17    , k18    ,
     .       k19    , k22    , HMnm
#ifdef RT 
     .,      k27    , k28    , k29    , k30    , k31    
#endif /* RT */   
#endif /* H2 */
#ifdef HD
     .,      k50    , k51    , k52    , k53    , k54    , k55    ,
     .       k56
#endif /* HD */
#ifdef H
     .,      ceHI   , ciHI   , reHII  , DM
#ifdef RT 
     .,      piHI   , sigma24
#endif /* RT */   
#endif /* H */
#ifdef He
     .,      ceHeI  , ciHeI  , ciHeIS , ceHeII ,
     .       ciHeII , reHeII1, reHeII2, reHeIII
#ifdef RT 
     .,      piHeI  , piHeII , sigma25, sigma26 
#endif /* RT */   
#endif /* He */
#ifdef H2
     .,      vibh   , hyd01k , h2k01  , roth   , rotl,   
     .       gpldl  , gphdl  
#ifdef RT 
     .,      sigma27, sigma28, sigma29, sigma30, sigma31 
#endif /* RT */   
#endif /* H2 */
#ifdef HD
     .,      hdlte  , hdlow  
#endif /* HD */
     .,      brem   , cmpt
#ifdef RT 
     .,      fcentral, f, f_outer !, rad_press  
#endif /* RT */   
