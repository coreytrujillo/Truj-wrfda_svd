subroutine da_calculate_gradj(it, iter, cv_size, cv_size_jb, cv_size_je, cv_size_jp, &
                              cv_size_jl, xbx, be, iv, cv, y, grad, grid, config_flags, re )

   !---------------------------------------------------------------------------
   ! Purpose: Calculates the gradient of the cost function w/r to cv
   !
   ! Called from da_minimise_cg (or da_minimise_lz)
   !
   ! History: 12/12/08 - Creation from da_calculate_j (Tom Auligne)
   !
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(in)                :: it     ! external iteration #.
   integer, intent(in)                :: iter   ! internal iteration #.
   integer, intent(in)                :: cv_size    ! Total cv size.
   integer, intent(in)                :: cv_size_jb, cv_size_je, cv_size_jp, cv_size_jl
   type (xbx_type),intent(inout)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)         :: be     ! background error structure.
   type (iv_type), intent(inout)      :: iv     ! innovation vector (o-b).
   real*8, intent(in)                 :: cv     (1:cv_size)   ! control variables.
   type (y_type), intent(inout)       :: y
   real*8, intent(out)                :: grad(cv_size)        ! gradient of cost function
   type (y_type), optional, intent(inout) :: re     ! residual vector (o-a).
   
   type(domain), intent(inout)  :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   integer          :: je_start, je_end             ! Start/end indices of Je.
   integer          :: jl_start, jl_end             ! Start/end indices of Jl.
   real             :: jo_partial                   ! jo for this processor
   type (y_type)    :: jo_grad_y                    ! Grad_y(jo)
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   integer          :: mz(13)
#else
   integer          :: mz(7)
#endif
#else
#ifdef CLOUD_CV
   integer          :: mz(15)
#else
   integer          :: mz(9)
#endif
#endif
   real             :: grad_jo(cv_size)
   real             :: grad_jb(cv_size)
   real             :: grad_je(cv_size)
   real             :: grad_jd(cv_size)
   real             :: grad_jp(cv_size)
   real             :: grad_js(cv_size)
   real             :: grad_jl(cv_size)
   real             :: gnorm_j, gnorm_jo, gnorm_jb, gnorm_je, gnorm_jd, gnorm_jp, gnorm_js, gnorm_jl
   logical          :: jcdf_flag

   real            :: cv_dummy(cv_size)

   ! Variables for VarBC background constraint
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                           :: inst, ichan, npred, ipred, id
   real                              :: bgerr
   integer                           :: n

   if (trace_use) call da_trace_entry("da_calculate_gradj")

   !-------------------------------------------------------------------------
   ! [0.0] initialization:
   !-------------------------------------------------------------------------
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /)
#endif
#else
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#endif
#endif
   je_start   = cv_size_jb + 1
   je_end     = cv_size_jb + cv_size_je
   jp_start   = cv_size_jb + cv_size_je + 1
   jp_end     = cv_size_jb + cv_size_je + cv_size_jp
   jl_start   = cv_size_jb + cv_size_je + cv_size_jp + 1
   jl_end     = cv_size_jb + cv_size_je + cv_size_jp + cv_size_jl
   
   grad_jo = 0.0
   grad_jb = 0.0
   grad_je = 0.0
   grad_jd = 0.0
   grad_jp = 0.0
   grad_js = 0.0
   grad_jl = 0.0

   jcdf_flag = .false.

   !-------------------------------------------------------------------------
   ! [1.0] calculate grad_v (jo):
   !-------------------------------------------------------------------------
   call da_allocate_y(iv, jo_grad_y)
#if (WRF_CHEM == 1)
   call da_allocate_y_chem(iv, jo_grad_y)
#endif

   cv_dummy = cv

   if (present(re)) then
      call da_calculate_grady(iv, re, jo_grad_y)   
      if ( iter > 0 .and. test_gradient ) jcdf_flag = .true.
      call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
              grid%vp, grid%vv, grid%vp6, grid%vv6, &
              xbx, jo_grad_y, grid, config_flags, jcdf_flag &
#if (WRF_CHEM == 1)
              , grid%vch &
#endif
              )
   else
      !Apply preconditioner:  P * cv
      if (spectral_precon .and. it.gt.1 .and. ( use_lanczos &
#if defined(LAPACK)
            .or. use_randomblock &
#endif
            )) then
         call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 it-1, cv, cv_dummy, 'TL' )
      end if

      call da_transform_vtoy(cv_size, be, grid%ep, cv_dummy, iv, grid%vp, &
              grid%vv, grid%vp6, grid%vv6, xbx, y, grid, config_flags &
#if (WRF_CHEM == 1)
              , grid%vch &
#endif
              )
      call da_calculate_grady(iv, y, jo_grad_y)   
      call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
              grid%vp, grid%vv, grid%vp6, grid%vv6, &
              xbx, jo_grad_y, grid, config_flags, .true. &
#if (WRF_CHEM == 1)
              , grid%vch &
#endif
              )
      grad_jo = - grad_jo    !! Compensate for sign in calculation of grad_v (Jo)
   end if

    
   call da_deallocate_y(jo_grad_y)

#if defined(LAPACK)
   if ( present(re) &
        .or. .not.use_randomblock &
        .or. (spectral_precon .and. it.gt.1) )then
!!        .or. rand_type.eq.3 &
#endif

   !-------------------------------------------------------------------------
   ! [2.0] calculate grad_v (jb):
   !-------------------------------------------------------------------------
   if (cv_size_jb > 0) grad_jb(1:cv_size_jb) = jb_factor * cv_dummy(1:cv_size_jb)

   !-------------------------------------------------------------------------
   ! [3.0] calculate grad_v (je):
   !-------------------------------------------------------------------------
   if (cv_size_je > 0) grad_je(je_start:je_end) = je_factor * cv_dummy(je_start:je_end)
   
   !----------------------------------------------------------------------
   ! [3.1] calculate grad_v (jd):
   !----------------------------------------------------------------------
   if (use_wpec) then

      if (var4d) call da_error(__FILE__,__LINE__,(/'Cannot use 4dvar with dynamic constraint'/))
      if (wpec_factor <= 0) call da_error(__FILE__,__LINE__,(/'"wpec_factor" for dynamic constraint must be greater than zero'/))

      grid%xa%grad_p_x(:,:,:)=0.0
      grid%xa%grad_p_y(:,:,:)=0.0

      call da_transform_vtod_wpec(cv_size, be, grid%ep, cv_dummy, grid%vp, grid%vv, xbx, grid)

      grid%xa%grad_p_x=(grid%xa%grad_p_x)/wpec_factor
      grid%xa%grad_p_y=(grid%xa%grad_p_y)/wpec_factor

      call da_transform_vtod_wpec_adj(cv_size, be, grid%ep, grad_jd, grid%vp, grid%vv, xbx, grid)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] calculate grad_v (jp):
   !-------------------------------------------------------------------------
#if defined(RTTOV) || defined(CRTM)
   if (use_varbc .and. cv_size_jp > 0) then
      do inst = 1, iv % num_inst   
         do ichan = 1, iv%instid(inst)%nchan
            npred    = iv%instid(inst)%varbc(ichan)%npred
            if (npred <= 0) cycle               !! VarBC channels only	 
            do ipred = 1, npred
               id     = iv%instid(inst)%varbc(ichan)%index(ipred)
	       bgerr  = iv%instid(inst)%varbc(ichan)%bgerr(ipred)
	       if (bgerr > 0.0) &
                  grad_jp(id) = (1/sqrt(bgerr)) * &
                     SUM(cv_dummy(id) * iv%instid(inst)%varbc(ichan)%vtox(ipred,1:npred))
	    end do
         end do
      end do
   end if
#endif
      
   !-------------------------------------------------------------------------
   ! [5.0] calculate grad_v (js):
   !-------------------------------------------------------------------------
   if (ANY(use_satcv)) then
      do inst = 1, iv % num_inst   
         do n = iv%instid(inst)%info%n1, iv%instid(inst)%info%n2 ! loop for pixel
         ! Skin Temperature
         !-----------------
	    if (use_satcv(1)) &
            grad_js(iv%instid(inst)%cv_index(n)%ts) = cv_dummy(iv%instid(inst)%cv_index(n)%ts)
	    
         ! Cloud cover(s)
         !---------------
	    if (use_satcv(2)) then
	    grad_js(iv%instid(inst)%cv_index(n)%cc) = cv_dummy(iv%instid(inst)%cv_index(n)%cc)

	    WHERE (cv_dummy(iv%instid(inst)%cv_index(n)%cc) < 0.0 .or.                                &
	           cv_dummy(iv%instid(inst)%cv_index(n)%cc) > 1.0 )                                   &
	    grad_js(iv%instid(inst)%cv_index(n)%cc) = grad_js(iv%instid(inst)%cv_index(n)%cc) + &
                                                       10.0 * cv_dummy(iv%instid(inst)%cv_index(n)%cc)
            end if
	 end do
      end do	      
   end if

   !-------------------------------------------------------------------------
   ! [6.0] calculate grad_v (jl):
   !-------------------------------------------------------------------------
   if (cv_size_jl > 0) grad_jl(jl_start:jl_end) = cv_dummy(jl_start:jl_end)

#if defined(LAPACK)
   end if !use_randomblock test
#endif

   !--------------------------------------------------------------------------------------------------
   ! [7.0] calculate grad_v (j) = grad_v (jb) + grad_v (jo) + grad_v (je) + grad_v (jd) + grad_v (jp) + grad_v (js) + grad_v (jl)
   !--------------------------------------------------------------------------------------------------   
   grad = grad_jo + grad_jb + grad_je + grad_jd + grad_jp + grad_js + grad_jl

   !Apply preconditioner: 
   ! Full Hessian - Lanczos: P^T * (I_n + L^T H^T R^-1 H L * P * cv) = P^T * grad
   ! Obs. Hessian - RSVD:    P^T * (I_n + L^T H^T R^-1 H L * P * cv) - I_n * cv = P^T * grad - cv
   ! gradient: P^T * (cv + L^T H^T R^-1 d) = P^T * grad
   if (spectral_precon .and. it.gt.1 .and. ( use_lanczos &
#if defined(LAPACK)
         .or. use_randomblock &
#endif
         )) then
      cv_dummy = grad
      call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                              it-1, cv_dummy, grad, 'AD' )

#if defined(LAPACK)
      !Account for the fact that use_randblock multiplies vector only by obs. Hessian
!!      if (.not.present(re) .and. use_randomblock .and. rand_type.ne.3) grad = grad - cv
      if (.not.present(re) .and. use_randomblock) grad = grad - cv
#endif
   end if

   !-------------------------------------------------------------------------
   ! [8.0] write Gradient:
   !-------------------------------------------------------------------------
   if (present(re)) then
      gnorm_j  = sqrt(da_dot_cv(cv_size, grad,    grad,    grid, mz, jp_start, jp_end))
      gnorm_jo = sqrt(da_dot_cv(cv_size, grad_jo, grad_jo, grid, mz))
      gnorm_jb = sqrt(da_dot_cv(cv_size, grad_jb, grad_jb, grid, mz))
      gnorm_je = sqrt(da_dot_cv(cv_size, grad_je, grad_je, grid, mz))
      gnorm_jd = sqrt(da_dot_cv(cv_size, grad_jd, grad_jd, grid, mz))
      gnorm_jp = sqrt(da_dot_cv(cv_size, grad_jp, grad_jp, grid, mz, jp_start, jp_end))
      gnorm_js = sqrt(da_dot_cv(cv_size, grad_js, grad_js, grid, mz))
      gnorm_jl = sqrt(da_dot_cv(cv_size, grad_jl, grad_jl, grid, mz))
 
      if (rootproc) &
         write(grad_unit,fmt='(2x,i2,1x,e10.3,2x,i4,8(1x,f10.3))') & 
               it, eps(it), iter, gnorm_j, gnorm_jb, gnorm_jo, gnorm_je, gnorm_jd, gnorm_jp, gnorm_js, gnorm_jl
   end if
   if (trace_use) call da_trace_exit("da_calculate_gradj")

end subroutine da_calculate_gradj


