MODULE agrif_dom_update

   USE dom_oce
   USE agrif_parameters
   USE agrif_profiles
   USE agrif_recompute_scales
 
   IMPLICIT none
   PRIVATE

   PUBLIC agrif_update_all

CONTAINS 

#if defined key_agrif

   SUBROUTINE agrif_update_all
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE agrif_update_all  ***
      !!----------------------------------------------------------------------  
      !
      INTEGER :: ind1, ind2

      IF( Agrif_Root() ) return
      
      IF ( .NOT.ln_vert_remap ) THEN
         CALL agrif_update_variable(bottom_level_id,procname = update_bottom_level)
         Agrif_UseSpecialValueInUpdate = .FALSE.
         Agrif_SpecialValueFineGrid    = 0._wp         
         CALL agrif_update_variable(e3t_copy_id, procname = update_e3t_z) 
         !
      ELSE
         Agrif_UseSpecialValueInUpdate = .FALSE.
         Agrif_SpecialValueFineGrid    = 0._wp         
         CALL agrif_update_variable(e3t_id, procname = update_e3t_z_cons) 

         ! jc: extend update zone outside dynamical interface within sponge zone:
         ! Use max operator this time to account for cases for which Agrif_Rho > nbghostcells
         ind1 = CEILING(REAL(max(nbghostcells_x_w-1, nbghostcells_x_e-1), wp) / Agrif_Rhox() )
         ind2 = CEILING(REAL(max(nbghostcells_y_s-1, nbghostcells_y_n-1), wp) / Agrif_Rhoy() )
         CALL agrif_update_variable(e3t_copy_id, locupdate1=(/-ind1,0/), &
                             &                   locupdate2=(/-ind2,0/),procname = update_e3t_z_cons)
      ENDIF
      Agrif_UseSpecialValueInUpdate = .FALSE.
      !
      ! Update vertical scale factors at U, V and F-points:
      CALL Agrif_ChildGrid_To_ParentGrid()
      CALL agrif_recompute_scalefactors
      CALL Agrif_ParentGrid_To_ChildGrid()
      !    
   END SUBROUTINE agrif_update_all

   SUBROUTINE update_bottom_level( ptab, i1, i2, j1, j2, before)
      !!----------------------------------------------------------------------
      !!       ***  ROUTINE update_bottom_level  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL, DIMENSION(i1:i2,j1:j2)    , INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      !!---------------------------------------------------------------------- 
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = mbkt(i1:i2,j1:j2)*ssmask(i1:i2,j1:j2)
      ELSE
         mbkt(i1:i2,j1:j2) = nint(ptab(i1:i2,j1:j2))
         
         WHERE ( mbkt(i1:i2,j1:j2) .EQ. 0 )
            ssmask(i1:i2,j1:j2) = 0._wp
            mbkt(i1:i2,j1:j2)   = 1 
         ELSEWHERE
            ssmask(i1:i2,j1:j2) = 1._wp
         END WHERE 
      ENDIF
      !
   END SUBROUTINE update_bottom_level

   SUBROUTINE update_e3t_z( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** update_e3t_z ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj, jk
      !!---------------------------------------------
      !
      IF (before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  IF ( (ssmask(ji,jj) /=0._wp).AND.(mbkt(ji,jj).GE.jk) ) THEN
                     tabres(ji,jj,jk) = e3t_0(ji,jj,jk)
                  ELSE
                     tabres(ji,jj,jk) = 0._wp
                  ENDIF 
               END DO
            END DO
         END DO
      ELSE
         DO jk=1,jpk
            DO jj=j1,j2
               DO ji=i1,i2
                  IF ( ( mbkt(ji,jj).GE.jk ).AND.(ssmask(ji,jj)==1._wp) ) THEN
                     e3t_0(ji,jj,jk) = MAX(tabres(ji,jj,jk),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
                 !    e3t_0(ji,jj,jk) = tabres(ji,jj,jk)
                  ELSE
                     e3t_0(ji,jj,jk) = e3t_1d(jk)
                  ENDIF
               END DO
            END DO
         END DO
         !
      ENDIF
      ! 
   END SUBROUTINE update_e3t_z

   SUBROUTINE update_e3t_z_cons( tabres, i1, i2, j1, j2, k1, k2, before )
      !!---------------------------------------------
      !!           *** update_e3t_z_cons ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji, jj, jk, ik
      REAL(wp) :: zhmin, zdepth, zdepwp, ze3tp, ze3wp
      !!---------------------------------------------
      !
      IF (before) THEN
         DO jk = k1, k2-1
            DO jj = j1, j2
               DO ji = i1, i2
                   IF ( (ssmask(ji,jj) /=0._wp).AND.( mbkt(ji,jj) .GE. jk ) ) THEN
                      tabres(ji,jj,jk) = e3t_0(ji,jj,jk)
                   ELSE
                      tabres(ji,jj,jk) = 0._wp
                   endif
               END DO
            END DO
         END DO
         tabres(i1:i2,j1:j2,k2) = ssmask(i1:i2,j1:j2) ! To get fractional area
      ELSE
         IF( rn_hmin < 0._wp ) THEN   
            ik = - INT( rn_hmin )
         ELSE                          
            ik = MINLOC( gdepw_1d, mask = gdepw_1d > rn_hmin, dim = 1 )
         ENDIF
         zhmin = gdepw_1d(ik+1)

         ! Compute child bathymetry:
         bathy(i1:i2,j1:j2) = 0._wp
         DO jk=k1,k2-1   
            bathy(i1:i2,j1:j2) = bathy(i1:i2,j1:j2) + tabres(i1:i2,j1:j2,jk)
         END DO
         WHERE( bathy(i1:i2,j1:j2) == 0._wp )   ;   mbathy(i1:i2,j1:j2) = 0       
         ELSE WHERE                             ;   mbathy(i1:i2,j1:j2) = jpkm1  
         END WHERE

         DO jk = jpkm1, 1, -1
            zdepth = gdepw_1d(jk) ! + MIN( e3zps_min, e3t_1d(jk)*e3zps_rat )
            WHERE( 0._wp < bathy(i1:i2,j1:j2) .AND. bathy(i1:i2,j1:j2) <= zdepth )   mbathy(i1:i2,j1:j2) = jk-1
         END DO

         ! Scale factors and depth at T- and W-points
         DO jk = 1, jpk  
            gdept_0(i1:i2,j1:j2,jk) = gdept_1d(jk)
            gdepw_0(i1:i2,j1:j2,jk) = gdepw_1d(jk)
            e3t_0  (i1:i2,j1:j2,jk) = e3t_1d  (jk)
            e3w_0  (i1:i2,j1:j2,jk) = e3w_1d  (jk)
         END DO
         ! Scale factors and depth at T- and W-points
         DO jj = j1, j2
            DO ji = i1, i2 
               ik = mbathy(ji,jj)
               IF( ik > 0 ) THEN               ! ocean point only
                  ! max ocean level case
                  IF( ik == jpkm1 ) THEN
                     zdepwp = bathy(ji,jj)
                     ze3tp  = bathy(ji,jj) - gdepw_1d(ik)
                     ze3wp = 0.5_wp * e3w_1d(ik) * ( 1._wp + ( ze3tp/e3t_1d(ik) ) )
                     e3t_0(ji,jj,ik  ) = ze3tp
                     e3t_0(ji,jj,ik+1) = ze3tp
                     e3w_0(ji,jj,ik  ) = ze3wp
                     e3w_0(ji,jj,ik+1) = ze3tp
                     gdepw_0(ji,jj,ik+1) = zdepwp
                     gdept_0(ji,jj,ik  ) = gdept_1d(ik-1) + ze3wp
                     gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + ze3tp
                     !
                  ELSE                         ! standard case
                     IF( bathy(ji,jj) <= gdepw_1d(ik+1) ) THEN 
                        gdepw_0(ji,jj,ik+1) = bathy(ji,jj)
                     ELSE                                       
                        gdepw_0(ji,jj,ik+1) = gdepw_1d(ik+1)
                     ENDIF
                     gdept_0(ji,jj,ik) = gdepw_1d(ik) + ( gdepw_0(ji,jj,ik+1) - gdepw_1d(ik) )           &
                        &                             * ((gdept_1d(     ik  ) - gdepw_1d(ik) )           &
                        &                             / ( gdepw_1d(     ik+1) - gdepw_1d(ik) ))
                     e3t_0  (ji,jj,ik) = e3t_1d  (ik) * ( gdepw_0 (ji,jj,ik+1) - gdepw_1d(ik))           &
                        &                             / ( gdepw_1d(      ik+1) - gdepw_1d(ik)) 
                     e3w_0(ji,jj,ik) =                                                                   & 
                        &      0.5_wp * (gdepw_0(ji,jj,ik+1) + gdepw_1d(ik+1) - 2._wp * gdepw_1d(ik) )   &
                        &             * ( e3w_1d(ik) / ( gdepw_1d(ik+1) - gdepw_1d(ik) ) )
                     !       ... on ik+1
                     e3w_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                     e3t_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                     gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + e3t_0(ji,jj,ik)
                  ENDIF
               ENDIF
            END DO
         END DO
         !
         DO jj=j1,j2
            DO ji=i1,i2
               bathy(ji,jj) = SUM(e3t_0(ji,jj,1:mbkt(ji,jj) ) ) 
            END DO
         END DO
         !
         WHERE ( ( mbathy(i1:i2,j1:j2) .EQ. 0 )  & 
           & .OR.(tabres(i1:i2,j1:j2,k2)<0.5_wp) &
           & .OR.(bathy(i1:i2,j1:j2)<zhmin) )
            ssmask(i1:i2,j1:j2) = 0._wp
            mbathy(i1:i2,j1:j2) = 0 
         ELSEWHERE
            ssmask(i1:i2,j1:j2) = 1._wp
         END WHERE
         mbkt(i1:i2,j1:j2) = MAX( mbathy(i1:i2,j1:j2), 1 )
      ENDIF
      ! 
   END SUBROUTINE update_e3t_z_cons
      
#else
   SUBROUTINE agrif_update_all
   END SUBROUTINE agrif_update_all
#endif

END MODULE agrif_dom_update
