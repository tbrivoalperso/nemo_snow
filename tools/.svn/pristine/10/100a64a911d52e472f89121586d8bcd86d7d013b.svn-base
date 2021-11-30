MODULE agrif_connect

   USE dom_oce
   USE agrif_parameters
   USE agrif_profiles

   IMPLICIT NONE
   PRIVATE

   PUBLIC agrif_boundary_connections, agrif_bathymetry_connect 

CONTAINS

#if defined key_agrif

   SUBROUTINE agrif_boundary_connections
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE agrif_boundary_connections  ***
      !!----------------------------------------------------------------------  
      IF( Agrif_Root() ) return

      CALL agrif_connection()
      !
!      CALL Agrif_Bc_variable(bottom_level_id, procname = connect_bottom_level)
      ! 
!      CALL Agrif_Bc_variable(e3t_copy_id, procname = connect_e3t_copy)

      ALLOCATE(e3t_interp_done(jpi,jpj))
      e3t_interp_done(:,:) = .FALSE. 
      ! set extrapolation on for interpolation near the coastline:
      Agrif_UseSpecialValue = .TRUE.
      Agrif_SpecialValue = 0._wp
      CALL Agrif_Bc_variable(e3t_connect_id, procname = connect_e3t_connect)
      ! Override in ghost zone by nearest value:
      Agrif_UseSpecialValue = .FALSE.
      e3t_interp_done(:,:) = .FALSE.
      CALL Agrif_Bc_variable(e3t_copy_id,    procname = connect_e3t_connect)
      Agrif_UseSpecialValue = .FALSE.
      DEALLOCATE(e3t_interp_done)
      !
   END SUBROUTINE agrif_boundary_connections

   SUBROUTINE agrif_bathymetry_connect
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE agrif_bathymetry_connect  ***
      !!----------------------------------------------------------------------  
      IF( Agrif_Root() ) return

      CALL agrif_connection()
      !
      ALLOCATE(e3t_interp_done(jpi,jpj))
      e3t_interp_done(:,:) = .FALSE. 
      ! set extrapolation on for interpolation near the coastline:
      Agrif_UseSpecialValue = .TRUE.
      Agrif_SpecialValue = 0._wp
      CALL Agrif_Bc_variable(e3t_connect_id, procname = connect_bathy_connect)
      ! Override in ghost zone by nearest value:
      Agrif_UseSpecialValue = .FALSE.
      e3t_interp_done(:,:) = .FALSE.
      CALL Agrif_Bc_variable(e3t_copy_id,    procname = connect_bathy_connect)
      Agrif_UseSpecialValue = .FALSE.
      DEALLOCATE(e3t_interp_done)
      !
   END SUBROUTINE agrif_bathymetry_connect

   SUBROUTINE connect_e3t_copy( ptab, i1, i2, j1, j2, k1, k2, before, nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE connect_e3t_copy  ***
      !!----------------------------------------------------------------------  
      INTEGER                               , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) ::   ptab
      LOGICAL                               , INTENT(in   ) ::   before
      INTEGER                               , INTENT(in   ) ::   nb , ndir
      !
      !!---------------------------------------------------------------------- 
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2,k1:k2) = e3t_0(i1:i2,j1:j2,k1:k2)
      ELSE
         e3t_0(i1:i2,j1:j2,1:jpk) = ptab(i1:i2,j1:j2,1:jpk)
      ENDIF
      !
   END SUBROUTINE connect_e3t_copy
   
   SUBROUTINE connect_bottom_level( ptab, i1, i2, j1, j2, before, nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE connect_bottom_level  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      INTEGER                         , INTENT(in   ) ::   nb , ndir
      !
      !!---------------------------------------------------------------------- 
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = mbkt(i1:i2,j1:j2)*ssmask(i1:i2,j1:j2)
      ELSE
         mbkt(i1:i2,j1:j2) = nint(ptab(i1:i2,j1:j2))
         WHERE (mbkt(i1:i2,j1:j2)==0)
           ssmask(i1:i2,j1:j2) = 0.
         ELSEWHERE
           ssmask(i1:i2,j1:j2) = 1.
         END WHERE  
      ENDIF
      !
   END SUBROUTINE connect_bottom_level
   
   SUBROUTINE connect_e3t_connect( ptab, i1, i2, j1, j2, k1, k2, before, nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE connect_e3t_connect  ***
      !!----------------------------------------------------------------------  
      INTEGER                               , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) ::   ptab
      LOGICAL                               , INTENT(in   ) ::   before
      INTEGER                               , INTENT(in   ) ::   nb , ndir
      !
      !!---------------------------------------------------------------------- 
      INTEGER :: ji, jj, jk, ik 
      REAL(wp), DIMENSION(i1:i2,j1:j2) :: bathy_local, bathy_interp
      REAL(wp) :: zdepth, zdepwp, zmax, ze3tp, ze3wp, zhmin 
      !
      IF( before) THEN
         DO jk=k1, k2
            DO jj=j1,j2
               DO ji=i1,i2
                  IF( mbkt(ji,jj) .GE. jk ) THEN
                     ptab(ji,jj,jk) = e3t_0(ji,jj,jk)
                  ELSE
                     ptab(ji,jj,jk) = 0.
                  ENDIF
               END DO
            END DO
         END DO
         !
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj,k2) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj)
            END DO
         END DO
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               bathy_local (ji,jj) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj)
               bathy_interp (ji,jj) = ptab(ji,jj,k2)
               ! keep child masking in transition zone:
               IF ((ztabramp(ji,jj)/=1._wp).AND.(bathy_local(ji,jj)==0._wp)) bathy_interp(ji,jj)=0._wp
        ! Connected bathymetry
               IF( .NOT.e3t_interp_done(ji,jj) ) THEN
                  bathy_local(ji,jj)=(1.-ztabramp(ji,jj))*bathy_local(ji,jj)+ztabramp(ji,jj)*bathy_interp(ji,jj)
               ENDIF
            END DO
         END DO

        ! Update mbkt and ssmask
         IF( rn_hmin < 0._wp ) THEN
            ik = - INT( rn_hmin )
         ELSE
            ik = MINLOC( gdepw_1d, mask = gdepw_1d > rn_hmin, dim = 1 )
         ENDIF
         zhmin = gdepw_1d(ik+1)

         zmax = gdepw_1d(jpk) + e3t_1d(jpk)
         bathy_local(:,:) = MAX(MIN(zmax,bathy_local(:,:)),0._wp)
         WHERE( bathy_local(i1:i2,j1:j2) == 0._wp)
            mbathy(i1:i2,j1:j2) = 0
         ELSE WHERE 
            mbathy(i1:i2,j1:j2) = jpkm1 
            bathy_local(i1:i2,j1:j2) = MAX(  zhmin , bathy_local(i1:i2,j1:j2)  ) 
         END WHERE

         DO jk=jpkm1,1,-1
           zdepth = gdepw_1d(jk) + MIN(e3zps_min,e3t_1d(jk)*e3zps_rat)
           WHERE( 0._wp < bathy_local(:,:) .AND. bathy_local(:,:) <= zdepth ) mbathy(i1:i2,j1:j2) = jk-1
         ENDDO

         WHERE (mbathy(i1:i2,j1:j2) == 0); ssmask(i1:i2,j1:j2) = 0
         ELSE WHERE                      ; ssmask(i1:i2,j1:j2) = 1.
         END WHERE
         
         mbkt(i1:i2,j1:j2) = MAX( mbathy(i1:i2,j1:j2), 1 )
         !
         DO jj = j1, j2
            DO ji = i1, i2
               IF( .NOT.e3t_interp_done(ji,jj) ) THEN ! the connection has not yet been done
                  DO jk = 1, jpk    
                     gdept_0(ji,jj,jk) = gdept_1d(jk)
                     gdepw_0(ji,jj,jk) = gdepw_1d(jk)
                     e3t_0  (ji,jj,jk) = e3t_1d  (jk)
                     e3w_0  (ji,jj,jk) = e3w_1d  (jk)
                  END DO 
                  !
                  ik = mbathy(ji,jj)
                  IF( ik > 0 ) THEN               ! ocean point only
                     ! max ocean level case
                     IF( ik == jpkm1 ) THEN
                        zdepwp = bathy_local(ji,jj)
                        ze3tp  = bathy_local(ji,jj) - gdepw_1d(ik)
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
                        IF( bathy_local(ji,jj) <= gdepw_1d(ik+1) ) THEN
                           gdepw_0(ji,jj,ik+1) = bathy_local(ji,jj)
                        ELSE
                           gdepw_0(ji,jj,ik+1) = gdepw_1d(ik+1)
                        ENDIF
                        gdept_0(ji,jj,ik) = gdepw_1d(ik) + ( gdepw_0(ji,jj,ik+1) - gdepw_1d(ik) ) &
                              &                * ((gdept_1d(     ik  ) - gdepw_1d(ik) )           &
                              &                / ( gdepw_1d(     ik+1) - gdepw_1d(ik) ))
                        e3t_0  (ji,jj,ik) = e3t_1d  (ik) * ( gdepw_0 (ji,jj,ik+1) - gdepw_1d(ik)) &
                              &                / ( gdepw_1d(      ik+1) - gdepw_1d(ik))
                        e3w_0(ji,jj,ik) = &
                              & 0.5_wp * (gdepw_0(ji,jj,ik+1) + gdepw_1d(ik+1) - 2._wp * gdepw_1d(ik) )   &
                              &        * ( e3w_1d(ik) / ( gdepw_1d(ik+1) - gdepw_1d(ik) ) )
                        !       ... on ik+1
                        e3w_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                        e3t_0  (ji,jj,ik+1) = e3t_0  (ji,jj,ik)
                        gdept_0(ji,jj,ik+1) = gdept_0(ji,jj,ik) + e3t_0(ji,jj,ik)
                     ENDIF
                  ENDIF
               ENDIF
               e3t_interp_done(ji,jj) = .TRUE.
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE connect_e3t_connect

   SUBROUTINE connect_bathy_connect( ptab, i1, i2, j1, j2, k1, k2, before, nb,ndir)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE connect_e3t_connect  ***
      !!----------------------------------------------------------------------  
      INTEGER                               , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) ::   ptab
      LOGICAL                               , INTENT(in   ) ::   before
      INTEGER                               , INTENT(in   ) ::   nb , ndir
      !
      !!---------------------------------------------------------------------- 
      INTEGER :: ji, jj, jk
      !
      IF( before) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  IF( mbkt(ji,jj) .GE. jk ) THEN
                     ptab(ji,jj,jk) = e3t_0(ji,jj,jk)
                  ELSE
                     ptab(ji,jj,jk) = 0._wp
                  ENDIF
               END DO
            END DO
         END DO
         !
         DO jj=j1,j2
            DO ji=i1,i2
               ptab(ji,jj,k2) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj)
            END DO
         END DO
      ELSE
         DO jj=j1,j2
            DO ji=i1,i2
               ! keep child masking in transition zone:
               IF ((ztabramp(ji,jj)/=1._wp).AND.(bathy(ji,jj)==0._wp)) ptab(ji,jj,k2) = 0._wp
               ! Connected bathymetry
               IF( .NOT.e3t_interp_done(ji,jj) ) THEN
                  bathy(ji,jj)=(1._wp-ztabramp(ji,jj))*bathy(ji,jj)+ztabramp(ji,jj)*ptab(ji,jj,k2)
                  e3t_interp_done(ji,jj) = .TRUE.
               ENDIF
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE connect_bathy_connect
   
   SUBROUTINE agrif_connection
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE  Agrif_connection ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, ind1, ind2
      INTEGER  ::   ispongearea, istart
      REAL(wp) ::   z1_spongearea
      !!----------------------------------------------------------------------
      !
      ! Define ramp from boundaries towards domain interior at T-points
      ! Store it in ztabramp

      ALLOCATE(ztabramp(jpi,jpj))
      ispongearea = 1 + npt_connect * Agrif_iRhox()
      istart = npt_copy * Agrif_iRhox()
      z1_spongearea = 1._wp / REAL( ispongearea, wp )
      
      ztabramp(:,:) = 0._wp

      ! --- West --- !
      IF( lk_west ) THEN
         ind1 = nn_hls + nbghostcells + istart
         ind2 = ind1 + ispongearea 
         DO ji = mi0(ind1), mi1(ind2)   
            DO jj = 1, jpj               
               ztabramp(ji,jj) = REAL(ind2 - mig(ji), wp) * z1_spongearea
            END DO
         ENDDO
            ! ghost cells:
            ind1 = 1
            ind2 = nn_hls + nbghostcells + istart  ! halo + land + nbghostcells
            DO ji = mi0(ind1), mi1(ind2)   
               DO jj = 1, jpj               
                  ztabramp(ji,jj) = 1._wp
               END DO
            END DO
      ENDIF

      ! --- East --- !
      IF( lk_east ) THEN
         ind2 = jpiglo -  (nn_hls + nbghostcells -1 ) - istart
         ind1 = ind2 -ispongearea       
         DO ji = mi0(ind1), mi1(ind2)
            DO jj = 1, jpj
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( mig(ji) - ind1 ) * z1_spongearea )
            ENDDO
         ENDDO
            ! ghost cells:
            ind1 = jpiglo -  (nn_hls + nbghostcells - 1 ) - istart   ! halo + land + nbghostcells - 1
            ind2 = jpiglo - 1
            DO ji = mi0(ind1), mi1(ind2)
               DO jj = 1, jpj
                  ztabramp(ji,jj) = 1._wp
               END DO
            END DO
      ENDIF

      ispongearea = 1 + npt_connect * Agrif_iRhoy()
      istart = npt_copy * Agrif_iRhoy()
      z1_spongearea = 1._wp / REAL( ispongearea, wp )

      ! --- South --- !
      IF( lk_south ) THEN
         ind1 = nn_hls + nbghostcells + istart
         ind2 = ind1 + ispongearea 
         DO jj = mj0(ind1), mj1(ind2) 
            DO ji = 1, jpi
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( ind2 - mjg(jj) ) * z1_spongearea  )
            END DO
         ENDDO
            ! ghost cells:
            ind1 = 1
            ind2 = nn_hls + nbghostcells + istart                 ! halo + land + nbghostcells
            DO jj = mj0(ind1), mj1(ind2) 
               DO ji = 1, jpi
                  ztabramp(ji,jj) = 1._wp
               END DO
            END DO
      ENDIF

      ! --- North --- !
      IF( lk_north ) THEN
         ind2 = jpjglo - (nn_hls + nbghostcells - 1) - istart
         ind1 = ind2 -ispongearea
         DO jj = mj0(ind1), mj1(ind2)
            DO ji = 1, jpi
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( mjg(jj) - ind1 ) * z1_spongearea )
            END DO
         ENDDO
            ! ghost cells:
            ind1 = jpjglo - (nn_hls + nbghostcells - 1) - istart      ! halo + land + nbghostcells - 1
            ind2 = jpjglo
            DO jj = mj0(ind1), mj1(ind2)
               DO ji = 1, jpi
                  ztabramp(ji,jj) = 1._wp
               END DO
            END DO
      ENDIF
      !
   END SUBROUTINE agrif_connection

#else
   SUBROUTINE agrif_boundary_connections
   END SUBROUTINE agrif_boundary_connections
   SUBROUTINE agrif_bathymetry_connect 
   END SUBROUTINE agrif_bathymetry_connect 
#endif

END MODULE agrif_connect
