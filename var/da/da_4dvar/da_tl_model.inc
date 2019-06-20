subroutine da_tl_model ()

   !-----------------------------------------------------------------------
   ! Purpose: Acrivate tangent linear  model
   !    Author: Xin Zhang,  MMM/ESSL/NCAR,  Date: 05/19/2010
   !-----------------------------------------------------------------------

   implicit none

   integer :: time_step, fid, ierr, i, j, alarmid, gwd_opt
   character (len=80)     :: inpname 

   if (trace_use_dull) call da_trace_entry("da_tl_model")

#ifdef _OPENMP
   call setfeenv()
#endif

#if (WRF_CHEM == 1)
   IF ( use_nonchemobs ) THEN
#endif
      ! Set up perturbation output at each slots
      model_config_rec%auxhist18_outname = "tl_d<domain>_<date>"
      call nl_set_auxhist18_interval_s ( head_grid%id, var4d_bin ) 
      call nl_set_auxhist18_begin_s ( head_grid%id, 0 )
      call nl_set_auxhist18_end_h ( head_grid%id, 999 )
      call nl_set_io_form_auxhist18 ( head_grid%id, 2 )
      call nl_set_frames_per_auxhist18 ( head_grid%id, 1 )

      ! Set up basic states reading
      model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
      call nl_get_time_step ( head_grid%id, time_step )
      call nl_set_auxinput16_interval_s (head_grid%id, time_step )
      call nl_set_io_form_auxinput16 ( head_grid%id, 2 )
      call nl_set_frames_per_auxinput16 ( head_grid%id, 1 )
      call nl_get_gwd_opt ( head_grid%id, gwd_opt )
#if (WRF_CHEM == 1)
   ENDIF
#endif

   ! Set all perturbation to ZERO

   head_grid%g_scalar = 0.0
   head_grid%g_tke_2 = 0.0
   head_grid%g_al = 0.0
   head_grid%g_ustm = 0.0
   head_grid%g_hfx = 0.0
   head_grid%g_qfx = 0.0
   head_grid%g_h_diabatic = 0.0

   head_grid%g_u_1 = 0.0
   head_grid%g_v_1 = 0.0
   head_grid%g_w_1 = 0.0
   head_grid%g_t_1 = 0.0
   head_grid%g_ph_1 = 0.0
   head_grid%g_mu_1 = 0.0
   head_grid%g_tke_1 = 0.0

   head_grid%g_ru = 0.0
   head_grid%g_rv = 0.0
   head_grid%g_rw = 0.0
   head_grid%g_ww = 0.0
   head_grid%g_alt = 0.0
   head_grid%g_php = 0.0
   head_grid%g_mut = 0.0
   head_grid%g_muu = 0.0
   head_grid%g_muv = 0.0
   head_grid%g_p_hyd_w = 0.0
   head_grid%g_p_hyd = 0.0
   head_grid%g_z_at_w = 0.0
   head_grid%g_z = 0.0
   head_grid%g_zx = 0.0
   head_grid%g_zy = 0.0
   head_grid%g_rdz = 0.0
   head_grid%g_rdzw = 0.0
   head_grid%g_div = 0.0
   head_grid%g_nba_mij = 0.0
   head_grid%g_nba_rij = 0.0
   head_grid%g_defor11 = 0.0
   head_grid%g_defor22 = 0.0
   head_grid%g_defor33 = 0.0
   head_grid%g_defor12 = 0.0
   head_grid%g_defor13 = 0.0
   head_grid%g_defor23 = 0.0
   head_grid%g_bn2 = 0.0
   head_grid%g_xkmh = 0.0
   head_grid%g_xkmv = 0.0
   head_grid%g_xkhh = 0.0
   head_grid%g_xkhv = 0.0
   head_grid%g_ru_tend = 0.0
   head_grid%g_rv_tend = 0.0
   head_grid%g_u_save = 0.0
   head_grid%g_v_save = 0.0
   head_grid%g_t_save = 0.0
   head_grid%g_rthften = 0.0
   head_grid%g_rqvften = 0.0
   head_grid%g_rthcuten = 0.0
   head_grid%g_rqvcuten = 0.0
   head_grid%g_rthraten = 0.0
   head_grid%g_muts = 0.0
   head_grid%g_mudf = 0.0
   head_grid%g_ru_m = 0.0
   head_grid%g_rv_m = 0.0
   head_grid%g_ww_m = 0.0

   IF ( gwd_opt .EQ. 1 ) THEN
      head_grid%g_dtaux3d = 0.0
      head_grid%g_dtauy3d = 0.0
      head_grid%g_dusfcg = 0.0
      head_grid%g_dvsfcg = 0.0
   ENDIF

   head_grid%g_rainnc  = 0.0
   head_grid%g_rainncv = 0.0
   head_grid%g_rainc  = 0.0
   head_grid%g_raincv = 0.0
   head_grid%g_pratec = 0.0

#if (WRF_CHEM == 1)
   if (cv_options.eq.0 .and. cv_options_chem.gt.0) then
      ! Added these 6-30-2017, need to remove if chem ICs/BCs are control variables
      head_grid%g_chem = 0.0
      head_grid%g_chem_bxs = 0.0; head_grid%g_chem_bxe = 0.0
      head_grid%g_chem_bys = 0.0; head_grid%g_chem_bye = 0.0
      head_grid%g_chem_btxs = 0.0; head_grid%g_chem_btxe = 0.0
      head_grid%g_chem_btys = 0.0; head_grid%g_chem_btye = 0.0


      ! Likewise remove these if simultaneous physical/chemical DA
      head_grid%g_u_2 = 0.0
      head_grid%g_v_2 = 0.0
      head_grid%g_w_2 = 0.0
      head_grid%g_t_2 = 0.0
      head_grid%g_ph_2 = 0.0
      head_grid%g_p = 0.0
      head_grid%g_mu_2 = 0.0
      head_grid%g_moist = 0.0
      head_grid%g_tracer = 0.0

      head_grid%g_u_bxs = 0.0; head_grid%g_u_bxe = 0.0
      head_grid%g_u_bys = 0.0; head_grid%g_u_bye = 0.0
      head_grid%g_v_bxs = 0.0; head_grid%g_v_bxe = 0.0
      head_grid%g_v_bys = 0.0; head_grid%g_v_bye = 0.0
      head_grid%g_w_bxs = 0.0; head_grid%g_w_bxe = 0.0
      head_grid%g_w_bys = 0.0; head_grid%g_w_bye = 0.0
      head_grid%g_ph_bxs = 0.0; head_grid%g_ph_bxe = 0.0
      head_grid%g_ph_bys = 0.0; head_grid%g_ph_bye = 0.0
      head_grid%g_t_bxs = 0.0; head_grid%g_t_bxe = 0.0
      head_grid%g_t_bys = 0.0; head_grid%g_t_bye = 0.0
      head_grid%g_mu_bxs = 0.0; head_grid%g_mu_bxe = 0.0
      head_grid%g_mu_bys = 0.0; head_grid%g_mu_bye = 0.0
      head_grid%g_moist_bxs = 0.0; head_grid%g_moist_bxe = 0.0
      head_grid%g_moist_bys = 0.0; head_grid%g_moist_bye = 0.0
      head_grid%g_scalar_bxs = 0.0; head_grid%g_scalar_bxe = 0.0
      head_grid%g_scalar_bys = 0.0; head_grid%g_scalar_bye = 0.0
      head_grid%g_tracer_bxs = 0.0; head_grid%g_tracer_bxe = 0.0
      head_grid%g_tracer_bys = 0.0; head_grid%g_tracer_bye = 0.0

      head_grid%g_u_btxs = 0.0; head_grid%g_u_btxe = 0.0
      head_grid%g_u_btys = 0.0; head_grid%g_u_btye = 0.0
      head_grid%g_v_btxs = 0.0; head_grid%g_v_btxe = 0.0
      head_grid%g_v_btys = 0.0; head_grid%g_v_btye = 0.0
      head_grid%g_w_btxs = 0.0; head_grid%g_w_btxe = 0.0
      head_grid%g_w_btys = 0.0; head_grid%g_w_btye = 0.0
      head_grid%g_ph_btxs = 0.0; head_grid%g_ph_btxe = 0.0
      head_grid%g_ph_btys = 0.0; head_grid%g_ph_btye = 0.0
      head_grid%g_t_btxs = 0.0; head_grid%g_t_btxe = 0.0
      head_grid%g_t_btys = 0.0; head_grid%g_t_btye = 0.0
      head_grid%g_mu_btxs = 0.0; head_grid%g_mu_btxe = 0.0
      head_grid%g_mu_btys = 0.0; head_grid%g_mu_btye = 0.0
      head_grid%g_moist_btxs = 0.0; head_grid%g_moist_btxe = 0.0
      head_grid%g_moist_btys = 0.0; head_grid%g_moist_btye = 0.0
      head_grid%g_scalar_btxs = 0.0; head_grid%g_scalar_btxe = 0.0
      head_grid%g_scalar_btys = 0.0; head_grid%g_scalar_btye = 0.0
      head_grid%g_tracer_btxs = 0.0; head_grid%g_tracer_btxe = 0.0
      head_grid%g_tracer_btys = 0.0; head_grid%g_tracer_btye = 0.0
   end if
#endif
   ! Set the physics options to what they are in TL and AD
#if (WRF_CHEM == 1)
   IF ( use_chemobs) THEN
      if ( original_sf_surface_physics .ne. 1 .AND. &
           original_sf_surface_physics .ne. 7 ) then
         call nl_set_ra_lw_physics (head_grid%id, 0)
         call nl_set_ra_sw_physics (head_grid%id, 0)
      end if
      if ( original_sf_sfclay_physics .gt. 0 .and. &
           original_sf_sfclay_physics .ne. 7 .and. &
           original_sf_sfclay_physics .ne. 91 ) &
         call nl_set_sf_sfclay_physics (head_grid%id, 0)
      if ( original_bl_pbl_physics .gt. 0 .and. original_bl_pbl_physics .ne. 7) &
         call nl_set_bl_pbl_physics (head_grid%id, 98)
   ELSE
#endif
      call nl_set_ra_lw_physics (head_grid%id, 0)
      call nl_set_ra_sw_physics (head_grid%id, 0)
      call nl_set_sf_sfclay_physics (head_grid%id, 0)
      if (original_bl_pbl_physics .gt. 0) &
         call nl_set_bl_pbl_physics (head_grid%id, 98)
#if (WRF_CHEM == 1)
   END IF
#endif
   if (original_mp_physics .NE. 0 .and. original_mp_physics .NE. 98 .and. original_mp_physics .NE. 99) &
      call nl_set_mp_physics (head_grid%id, original_mp_physics_ad) 
   if (original_cu_physics .gt. 0) then
      call nl_set_cu_physics (head_grid%id, 98) 
      head_grid%cudt = 0
   endif
   !call nl_set_mp_zero_out (head_grid%id, 0)
   !call nl_set_sf_surface_physics (head_grid%id, 1)
   call nl_set_ifsnow (head_grid%id, 1)
   call nl_set_icloud (head_grid%id, 0)
   !call nl_set_isfflx (head_grid%id, 0)

   call nl_set_var4d_run (head_grid%id, .true.)

   ! WRF model time-stepping.  Calls integrate().
   IF ( checkpoint_interval .LE. 0 ) THEN
      call wrf_run_tl
   ELSE
      call da_set_run_hours ()

      io_form_checkpt = io_form_xtraj

      call wrf_run_tl_checkpt

      call da_reset_run_hours ()
   ENDIF

   ! restore the physics options to the originals
   call nl_set_mp_physics (head_grid%id, original_mp_physics) 
   call nl_set_ra_lw_physics (head_grid%id, original_ra_lw_physics)
   call nl_set_ra_sw_physics (head_grid%id, original_ra_sw_physics)
   call nl_set_sf_sfclay_physics (head_grid%id, original_sf_sfclay_physics)
   call nl_set_bl_pbl_physics (head_grid%id, original_bl_pbl_physics)
   call nl_set_cu_physics (head_grid%id, original_cu_physics)
   head_grid%cudt = original_cudt
   !call nl_set_mp_zero_out (head_grid%id, original_mp_zero_out)
   !call nl_set_sf_surface_physics (head_grid%id, original_sf_surface_physics)
   call nl_set_ifsnow (head_grid%id, original_ifsnow)
   call nl_set_icloud (head_grid%id, original_icloud)
   !call nl_set_isfflx (head_grid%id, original_isfflx)

   ! Turn off auxinput2 reading and auxhist18 output
   call nl_set_auxinput16_interval_s (head_grid%id, 0 )
   call nl_set_auxhist18_interval_s ( head_grid%id, 0 )

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   call Setup_Timekeeping ( head_grid )

   if (trace_use_dull) call da_trace_exit("da_tl_model")

end subroutine da_tl_model

