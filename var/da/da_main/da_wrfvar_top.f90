module da_wrfvar_top

   !-----------------------------------------------------------------------
   ! Purpose: 
   !-----------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type,model_config_rec, &
      model_to_grid_config_rec, get_config_as_buffer,set_config_as_buffer, &
      initial_config
   use module_domain, only : domain,alloc_and_configure_domain, head_grid, &
      program_name, domain_clock_get, domain_clock_set, x_type, dealloc_space_domain, &
      domain_destroy
   use module_driver_constants, only : max_comms
   use module_symbols_util, only : wrfu_finalize, wrfu_initialize, &
      wrfu_cal_gregorian
   use module_io_domain, only : close_dataset
   use module_timing
#ifdef VAR4D
   use da_4dvar, only : da_nl_model, model_grid, u6_2, v6_2, w6_2, t6_2, ph6_2, p6, &
#if (WRF_CHEM == 1)
      da_ad_model, da_init_model_input, da_init_model_output, &
#endif
#if defined(LAPACK) || (WRF_CHEM == 1)
      clean_4dvar, &
#endif
      mu6_2, psfc6, moist6, kj_swap, da_finalize_model, da_model_lbc_off
   use da_wrfvar_io, only : da_med_initialdata_output_lbc
#endif

#if defined(RTTOV) || defined(CRTM)
   use module_radiance, only : satinfo
#endif
#ifdef RTTOV
   use module_radiance, only : coefs, opts, sensor_descriptor
#endif

   use module_state_description, only : num_moist, num_a_moist, num_g_moist, &
      num_scalar, num_a_scalar, num_g_scalar, &
      !num_chem, PARAM_FIRST_SCALAR, num_tracer 
      num_tracer, num_chem, &
#if (WRF_CHEM == 1)
      num_a_tracer, num_g_tracer, &
      num_a_chem, num_g_chem, &
      num_emis_ant, num_scaleant, num_a_scaleant, num_g_scaleant, &
      num_ebu, num_ebu_in, num_scalebb, num_a_scalebb, num_g_scalebb, &
      num_ext_coef, & !num_a_ext_coef, num_g_ext_coef, &
      num_bscat_coef, & !num_a_bscat_coef, num_g_bscat_coef, &
      num_asym_par, & !num_a_asym_par, num_g_asym_par, &
      !num_eghg_bio,num_emis_dust,num_emis_seas,num_emis_vol
      num_chem_surf, num_chem_acft, &
      num_surf_hx, num_g_surf_hx, num_a_surf_hx, &
      num_acft_hx, num_g_acft_hx, num_a_acft_hx, &
#endif
      PARAM_FIRST_SCALAR
#if (WRF_CHEM == 1)
   use da_chem, only:  da_retrieve_chem_hx, da_write_obs_chem, da_calculate_chem_forcing_ad, &
      da_read_obs_chem_again
   use da_chem_tools, only: da_hdgn, da_setup_osse_chem

!   use da_control, only: use_synopobs, use_shipsobs, use_metarobs, use_soundobs, &
!              use_mtgirsobs, use_tamdarobs, use_bogusobs, use_pilotobs, &
!              use_airepobs, use_geoamvobs, use_polaramvobs, use_buoyobs, &
!              use_profilerobs, use_satemobs, use_gpspwobs, use_gpsztdobs, &
!              use_gpsrefobs, use_ssmiretrievalobs, use_ssmitbobs, &
!              use_ssmt1obs, use_ssmt2obs, use_qscatobs, use_hirs2obs, &
!              use_hirs3obs, use_hirs4obs, use_mhsobs, use_msuobs, &
!              use_amsuaobs, use_amsubobs, use_mwtsobs, use_mwhsobs, &
!              use_atmsobs, use_airsobs, use_eos_amsuaobs, use_hsbobs, use_airsretobs, &
!              use_rainobs, use_iasiobs, use_radarobs, use_radar_rv, use_radar_rf
#endif
   use module_tiles, only : set_tiles


#ifdef DM_PARALLEL
   use module_dm, only : local_communicator, local_communicator_x, &
      local_communicator_y, ntasks_x, ntasks_y, data_order_xyz, mytask, &
      ntasks, data_order_xy,wrf_dm_initialize
   use module_comm_dm, only : halo_radar_xa_w_sub,halo_ssmi_xa_sub, &
      halo_sfc_xa_sub, halo_xa_sub, halo_psichi_uv_adj_sub, halo_bal_eqn_adj_sub, &
      halo_psichi_uv_sub, halo_init_sub, halo_psichi_uv_adj_sub, halo_2d_work_sub,&
      halo_wpec_sub, halo_wpec_adj_sub, halo_xa_all_sub, halo_xb_all_sub
#endif

   ! too many namelist options to list
   use da_control
   use da_define_structures, only : y_type, j_type, iv_type, be_type, &
      xbx_type,da_deallocate_background_errors,da_initialize_cv, &
      da_zero_vp_type,da_allocate_y,da_deallocate_observations, &
#if (WRF_CHEM == 1)
      da_zero_xch_type, da_allocate_y_chem, &
#endif
#if defined(LAPACK)
      da_random_seed, &
#endif
      da_deallocate_y, da_zero_x, hessian_type
   use da_minimisation, only : da_get_innov_vector,da_minimise_cg, &
      da_minimise_lz, da_write_diagnostics, da_calculate_residual, &
      da_kmat_mul, &
#ifdef VAR4D
#if (WRF_CHEM == 1)
      da_calculate_aminusb, &
#endif
      da_dgn, da_transform_xatowrftl_lbc_wrapper, &
#endif
      da_calculate_grady, da_sensitivity, da_lanczos_io, da_calculate_j

   use da_linear_ops, only: da_output_increments, da_cv_io

#if defined(LAPACK)
   use da_randomisation, only: &
      da_randomise_svd, da_block_lanczos, &
#if (WRF_CHEM == 1)
      da_evaluate_increment, da_evaluate_hessian, &
#endif
      da_gen_omega, randvecdir

   use f95_precision, only: WP => DP
#endif
   use da_obs, only : da_transform_xtoy_adj 
   use da_obs_io, only : da_write_filtered_obs, da_write_obs, da_final_write_obs , &
#if (WRF_CHEM == 1)
      da_read_obs_chem, &
#endif
      da_write_obs_etkf, da_write_modified_filtered_obs
   use da_par_util, only : da_system,da_copy_tile_dims,da_copy_dims, da_cv_indices
   use da_physics, only : da_uvprho_to_w_lin
#if defined (CRTM) || defined (RTTOV)
   use da_radiance, only : da_deallocate_radiance
   use da_radiance1, only : num_tovs_before, tovs_recv_pe,tovs_copy_count, &
      tovs_send_pe,tovs_send_count,tovs_recv_start, num_tovs_after, &
      tovs_send_start, da_oi_stats_rad, da_write_oa_rad_ascii, da_setup_satcv
   use da_varbc, only : da_varbc_init,da_varbc_update
#endif
   use da_reporting, only : message, da_warning, da_error, da_message
   use da_setup_structures, only : da_setup_obs_structures, &
      da_setup_background_errors,da_setup_flow_predictors, &
      da_setup_cv, da_scale_background_errors, da_scale_background_errors_cv3
   use da_test, only : da_check, da_check_gradient
   use da_tools_serial, only : da_get_unit, da_free_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace, da_trace_report
   use da_transfer_model, only : da_transfer_xatoanalysis,da_setup_firstguess, &
#if (WRF_CHEM == 1)
       da_transfer_wrftoxb_chem, &
#endif
       da_transfer_headtomodel, da_transfer_wrftltoxa_adj, &
       da_transfer_xatowrftl, &
       da_transfer_wrftoxb
   use da_vtox_transforms, only : da_transform_vtox, da_transform_xtoxa, &
      da_transform_xtoxa_adj
   use da_wrfvar_io, only : da_med_initialdata_input, da_update_firstguess
   use da_tools, only : da_set_randomcv, da_get_julian_time

   use da_tools, only : map_info,map_info_ens,proj_merc, proj_ps,proj_lc,proj_latlon, &
      da_llxy_default,da_llxy_wrf,da_xyll,da_diff_seconds,da_map_set, &
      da_set_boundary_xb,da_togrid

#ifdef CRTM
   use module_radiance, only : crtm_destroy
   use da_crtm, only : channelinfo, sensor_descriptor
#endif

   use da_airep, only : da_oi_stats_airep
   use da_airsr , only : da_oi_stats_airsr
   use da_bogus, only : da_oi_stats_bogus
   use da_buoy , only : da_oi_stats_buoy
   use da_geoamv, only : da_oi_stats_geoamv
   use da_gpspw, only : da_oi_stats_gpspw
   use da_gpsref, only : da_oi_stats_gpsref
   use da_metar, only : da_oi_stats_metar
   use da_pilot, only : da_oi_stats_pilot
   use da_polaramv, only : da_oi_stats_polaramv
   use da_profiler, only : da_oi_stats_profiler
   use da_qscat, only : da_oi_stats_qscat
   use da_mtgirs, only : da_oi_stats_mtgirs
   use da_radar, only : da_oi_stats_radar, da_write_oa_radar_ascii
   use da_satem, only : da_oi_stats_satem
   use da_ships, only : da_oi_stats_ships
   use da_sound, only : da_oi_stats_sound, da_oi_stats_sonde_sfc
   use da_ssmi, only : da_oi_stats_ssmt1, da_oi_stats_ssmt2, da_oi_stats_ssmi_tb, da_oi_stats_ssmi_rv
   use da_synop, only : da_oi_stats_synop  
   use da_rain, only : da_oi_stats_rain

   use da_wrf_interfaces

   use da_netcdf_interface, only : da_get_var_2d_real_cdf

   implicit none

   integer :: loop, levels_to_process

   type (domain) , pointer :: keep_grid, grid_ptr, null_domain
   type (domain) , pointer :: another_grid, parent_grid, ensemble_grid, input_grid
   type (grid_config_rec_type), save :: config_flags
   integer                 :: number_at_same_level
   integer                 :: time_step_begin_restart

   integer :: domain_id , fid , oid , idum1 , idum2

#ifdef DM_PARALLEL
   integer                 :: nbytes
   integer, parameter      :: configbuflen = 4* CONFIG_BUF_LEN
   integer                 :: configbuf( configbuflen )
#endif

   character (len=80)      :: rstname

#ifdef RTTOV
#include "rttov_dealloc_coefs.interface"
#endif

contains

#include "da_wrfvar_init1.inc"
#include "da_wrfvar_init2.inc"
#include "da_wrfvar_run.inc"
#include "da_wrfvar_interface.inc"
#include "da_wrfvar_finalize.inc"
#include "da_solve.inc"

end module da_wrfvar_top
