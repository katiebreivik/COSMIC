      SUBROUTINE star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,id)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      
      real*8 mass,mt,tm,tn,tscls(20),lums(10),GB(10),zpars(20),dtm
      integer kw ,id
      integer metisse_id

      if (using_METISSE.eq.1) then
* persistent per-star track pool: remap the binary-component id (1
* or 2) to its pool slot, set once per evolv2 call in TRACKIDMAP
          metisse_id = id
          if (using_cmc.eq.1) metisse_id = track_id(id)
          !WRITE(*,*) 'Calling METISSE_star'
          CALL METISSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars,dtm,
     &                       metisse_id)

      elseif (using_SSE.eq.1) then
          !WRITE(*,*) 'Calling SSE_star'
          CALL SSE_star(kw,mass,mt,tm,tn,tscls,lums,GB,zpars)
      endif

      END
