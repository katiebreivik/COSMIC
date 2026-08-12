      SUBROUTINE hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  bhspin,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      integer kw,id,metisse_id
*
      real*8 mass,aj,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20)
      real*8 bhspin
      real*8 r,lum,mc,rc,menv,renv,k2,mcx

      if (using_METISSE.eq.1) then
* persistent per-star track pool: remap the binary-component id (1
* or 2) to its pool slot, set once per evolv2 call in TRACKIDMAP.
* mc_he/mc_co below stay indexed by the original (unmapped) id --
* they are COSMIC's own fixed 2-element per-component arrays, not
* METISSE's track pool.
          metisse_id = id
          if (using_cmc.eq.1) metisse_id = track_id(id)
          CALL METISSE_hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  mcx,metisse_id)
*
* KB: Assign mc_he and mc_co from METISSE output.
* For kw<=5: mc = McHe (core_mass = He core), mcx = McCO
* For kw=6:  mc = McCO (core_mass switches to CO core), mcx = McCO
* For kw>=7: mc = McCO (He star CO core), mcx = McCO
*
           if(kw.le.1)then
              mc_he(id) = mc
              mc_co(id) = 0.d0
           elseif(kw.le.3)then
              mc_he(id) = mc
              mc_co(id) = 0.d0
           elseif(kw.eq.4)then
              mc_he(id) = mc-mcx
              mc_co(id) = mcx
           elseif(kw.eq.5)then
              mc_co(id) = mcx
              mc_he(id) = mc - mcx
           elseif(kw.eq.6)then
*
* KB NOTE: METISSE tracks may still have nonzero He shell mass
* (t% pars% McHe - McCO) during TPAGB. Setting mc_he=0 here
* follows the SSE convention for now.
              mc_co(id) = mc
              mc_he(id) = 0.d0
           elseif(kw.ge.7.and.kw.le.9)then
              mc_co(id) = mc
              mc_he(id) = mt - mc
           endif
           ! get_bhspin is defined in assign_commons_cosmic.f90
           if (kw==14) CALL get_bhspin(bhspin,metisse_id)
          
      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_hrdiag'
          CALL SSE_hrdiag(mass,aj,mt,tm,tn,tscls,lums,GB,zpars,
     &                  r,lum,kw,mc,rc,menv,renv,k2,
     &                  bhspin,id)
      endif

      END
