module da_4dvar

use da_tracing, only : da_trace_entry, da_trace_exit
use da_reporting, only : da_error, message, da_message
use da_control, only : comm, var4d_bin, var4d_lbc, trace_use_dull, num_fgat_time, multi_inc, &
                       run_hours, run_days, adtl_run_hours, evalj, &
                       rootproc, &
#if (WRF_CHEM == 1)
                       init_osse_chem, calc_hx_only, &
                       cv_options, cv_options_chem, &
                       use_chemobs, use_nonchemobs, use_offlineobs, &
                       init_chem_scale, &
#endif
                       checkpoint_interval, write_checkpoints, disable_traj, cycle_interval

#ifdef VAR4D
use da_module_timing

#if (WRF_CHEM == 1)
   use, intrinsic :: iso_c_binding,                       &
                     ONLY: c_int32_t, C_CHAR, C_NULL_CHAR
#endif

use module_streams, only : MAX_WRF_ALARMS
use module_io_domain, only : open_r_dataset,close_dataset,input_auxhist16
use module_wrf_top, only : domain, head_grid, config_flags, &
             wrf_init, wrf_run, wrf_run_tl, wrf_run_ad, wrf_finalize, &
             wrf_run_tl_checkpt, wrf_run_ad_checkpt, wrf_run_cycle, &
             Setup_Timekeeping, gradient_out, io_form_checkpt
use mediation_pertmod_io, only : xtraj_io_initialize, adtl_initialize, &
             save_ad_forcing, read_ad_forcing, read_nl_xtraj, save_tl_pert, &
             read_tl_pert, swap_ad_forcing, read_nl_xtraj_checkpt
use module_configure, only :  model_to_grid_config_rec, grid_config_rec_type, model_config_rec
use module_domain, only : wrfu_timeinterval, domain_clock_get, domain_clock_set
use module_utility
use module_big_step_utilities_em, only : calc_mu_uv
use g_module_big_step_utilities_em, only : g_calc_mu_uv
use a_module_big_step_utilities_em, only : a_calc_mu_uv

#ifdef DM_PARALLEL
use module_dm, only : local_communicator
#endif

#if (WRF_CHEM == 1)
   use da_chem, only: da_copy_chem_to_model
#endif

type (domain), pointer :: model_grid
type (grid_config_rec_type) :: model_config_flags


character*256 :: timestr
character*20  :: xtrajprefix="xtraj_for_obs"
integer :: io_form_xtraj = 2 !Change this if alternative I/O MECHANISM is preferred
integer :: nl_called = 0

! Define some variables to save the NL physical option
integer :: original_mp_physics, original_ra_lw_physics, original_ra_sw_physics, &
           original_sf_sfclay_physics, original_bl_pbl_physics, original_cu_physics, &
           original_ifsnow, original_icloud, original_cudt, original_mp_physics_ad, &
           original_sf_surface_physics, original_interval_seconds, original_restart_interval

REAL , DIMENSION(:,:,:) , ALLOCATABLE  :: ubdy3dtemp1 , vbdy3dtemp1 , tbdy3dtemp1 , pbdy3dtemp1 , qbdy3dtemp1
REAL , DIMENSION(:,:,:) , ALLOCATABLE  :: ubdy3dtemp2 , vbdy3dtemp2 , tbdy3dtemp2 , pbdy3dtemp2 , qbdy3dtemp2
REAL , DIMENSION(:,:,:) , ALLOCATABLE  :: mbdy2dtemp1,  mbdy2dtemp2 , wbdy3dtemp1 , wbdy3dtemp2

REAL , DIMENSION(:,:,:) , ALLOCATABLE  :: u6_2, v6_2, w6_2, t6_2, ph6_2, p6
REAL , DIMENSION(:,:,:,:) , ALLOCATABLE  :: moist6
REAL , DIMENSION(:,:) , ALLOCATABLE  :: mu6_2, psfc6

integer :: end_year0, end_month0, end_day0, end_hour0

contains

#include "da_nl_model.inc"
#include "da_tl_model.inc"
#include "da_ad_model.inc"
#include "da_finalize_model.inc"
#include "da_4dvar_io.inc"
#include "da_4dvar_lbc.inc"
#include "da_set_run_hours.inc"
#include "da_reset_run_hours.inc"
#if (WRF_CHEM == 1)
#include "da_init_model_output.inc"
#include "da_init_model_input.inc"
#include "da_add_var_to_firstguess.inc"
#endif
#endif

end module da_4dvar
