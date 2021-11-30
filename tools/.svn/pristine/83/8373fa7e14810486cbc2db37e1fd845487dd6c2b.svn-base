MODULE agrif_parameters
	
   USE par_kind

   PUBLIC 

#if defined key_agrif
        LOGICAL :: ln_remove_closedseas=.FALSE.
        LOGICAL :: ln_vert_remap=.FALSE. ! =T is using volume conserving update
	INTEGER :: npt_copy      ! area (in coarse grid points) with piecewise
                                 ! constant bathymetry inside child zoom: should equal the sponge length
	INTEGER :: npt_connect   ! area (in coarse grid points) of coarse/child
                                 ! bathymetry blending
	REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   ztabramp
	LOGICAL,  PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   e3t_interp_done

#endif

END MODULE agrif_parameters
