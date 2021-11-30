***
        SUBROUTINE WRITESPP(jp,tphys,evolve_type,
     &                      mass,kstar,aj,tm,mc,rad,m0,lum,
     &                      teff,rc,menv,renv,ospin,B_0,
     &                      bacc,tacc,epoch,bhspin)
        IMPLICIT NONE
*
* Write results to spp array.
*
*     Author : Katie Breivik
*     Date :   20 Nov 2021
*
        REAL*8 mass,evolve_type,tphys,aj,tm
        REAL*8 mc,rad,m0,lum,teff,rc,menv,renv,ospin
        REAL*8 B_0,bacc,tacc,epoch,bhspin
        INTEGER jp,kstar
        REAL scm(50000,16),spp(25,20)
        COMMON /SINGLE/ scm,spp

        jp = MIN(900,jp + 1)
        spp(jp,1) = tphys
        spp(jp,2) = mass
        spp(jp,3) = float(kstar)
        spp(jp,4) = evolve_type
        spp(jp,5) = aj
        spp(jp,6) = tm
        spp(jp,7) = mc
        spp(jp,8) = rad
        spp(jp,9) = m0
        spp(jp,10) = lum
        spp(jp,11) = teff
        spp(jp,12) = rc
        spp(jp,13) = menv
        spp(jp,14) = renv
        spp(jp,15) = ospin
        spp(jp,16) = B_0
        spp(jp,17) = bacc
        spp(jp,18) = tacc
        spp(jp,19) = epoch
        spp(jp,20) = bhspin
        END

***
        SUBROUTINE WRITESCM(ip,tphys,kstar,m0,mass,lum,rad,teff,
     &                      mc,rc,menv,renv,epoch,deltam,
     &                      ospin,B_0,SN)
        IMPLICIT NONE
*
* Write results to bcm array.
*
*     Author : Scott Coughlin
*     Date :   12th March 2019
*
        REAL*8 tphys,m0,mass,lum,rad,teff,mc,rc
        REAL*8 menv,renv,epoch,ospin,B_0,deltam
        INTEGER kstar,SN
        INTEGER ip
        REAL scm(50000,16),spp(25,20)
        COMMON /SINGLE/ scm,spp

        ip = ip + 1
        scm(ip,1) = tphys
        scm(ip,2) = float(kstar)
        scm(ip,3) = m0
        scm(ip,4) = mass
        scm(ip,5) = lum
        scm(ip,6) = rad
        scm(ip,7) = teff
        scm(ip,8) = mc
        scm(ip,9) = rc
        scm(ip,10) = menv
        scm(ip,11) = renv
        scm(ip,12) = epoch
        scm(ip,13) = ospin
        scm(ip,14) = deltam
        scm(ip,15) = B_0
        scm(ip,16) = float(SN)

        END