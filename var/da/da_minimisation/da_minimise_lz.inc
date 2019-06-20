subroutine da_minimise_lz(grid, config_flags,            &
                           it, cv_size, xbx, be, iv,     &
                           j_grad_norm_target, xhat, qhat, cv, &
                           re, y, j_cost, eignvec, eignval, neign, &
                           j0_cost, j0_total, hessian_info)

   !-------------------------------------------------------------------------
   ! Purpose:         Main Lanczos minimisation routine 
   !
   ! Here 
   !    cv   is updated in outer-loop.
   !    xhat is the control variable in inner-loop.
   !
   ! Called from da_solve
   !
   ! History: 07/30/2008  Creation (Tom Auligne)
   !
   ! Reference: Golub and Van Loan 1996 (p493)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: it                            ! external iteration.
   integer, intent(in)               :: cv_size                       ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                           ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                            ! background error structure.
   type (iv_type), intent(inout)     :: iv                            ! ob. increment vector.
   real, intent(inout)               :: j_grad_norm_target            ! Target norm.
   real, intent(out)                 :: xhat(1:cv_size)               ! control variable (local).
   real, intent(out)                 :: qhat(1:cv_size, 0:ntmax(it))  ! cv copy.
   real, intent(inout)               :: cv(1:cv_size)                 ! control variable (local).
   type (y_type), intent(inout)      :: re                            ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                             ! y = H(x_inc) structure.
   type (j_type), intent(out)        :: j_cost                        ! cost function
   type (j_type), intent(out)        :: j0_cost                       ! cost function
   real*8, intent(out)               :: eignvec(ntmax(it), ntmax(it))
   real*8, intent(out)               :: eignval(ntmax(it))
   integer, intent(out)              :: neign
   real, intent(out)                 :: j0_total
   type(hessian_type), intent(inout), optional :: hessian_info
 
   integer                           :: i, j, k, n, inst, iter
   integer                           :: info, nsstwrk
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   integer                           :: mz(13)
#else
   integer                           :: mz(7)
#endif
#else
#ifdef CLOUD_CV
   integer                           :: mz(15)
#else
   integer                           :: mz(9)
#endif
#endif
   integer                           :: npred, ipred
   integer                           :: jp_start, jp_end       ! Start/end indices of Jp.
   real                              :: fhat(1:cv_size)        ! cv copy.
   real                              :: ghat(1:cv_size)        ! cv copy.
   real                              :: ghat0(1:cv_size)       ! cv copy.
   real                              :: shat(1:cv_size)        ! cv copy.
   real                              :: gdot, rho, mu, bndlm
   real                              :: alpha(ntmax(it)), beta(0:ntmax(it))
   real*8                            :: subdiag(ntmax(it))
   real                              :: c(cv_size), d(ntmax(it))
   real*8                            :: sstwrk(2*ntmax(it)-2)
   real                              :: bnds(ntmax(it))
   real                              :: adj_rhs, adj_sum_rhs ! < cv, cv >
   real                              :: j_total
   real                              :: dfs, dfsmax
   character(len=3) :: dummychar=""
   integer          :: ii
!#if defined(LAPACK)
!   real*8                            :: KMAT(ntmax(it),ntmax(it))
!#endif


#if (WRF_CHEM == 1)
   character(len=2) :: iterc,itc
#endif

   real, allocatable                 :: temp_e(:)
   real                              :: temp_l
   logical                           :: swapped

   if (trace_use) call da_trace_entry("da_minimise_lz")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using Lanczos method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
           be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /)
#endif
#else
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
           be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#endif
#endif
!#if defined(LAPACK)
!   KMAT = 0.D0
!#endif
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
   
   call da_calculate_j(it, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat0, grid, config_flags)
 
   j0_total = j_cost%total
   j0_cost = j_cost
   if (j0_total == 0.0) then
      if (trace_use) call da_trace_exit("da_minimise_lz")
      return
   end if
   ghat      = - ghat0
   beta(0)   = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
   qhat(:,0) = 0.D0
   j_grad_norm_target = beta(0)

   write(unit=stdout,fmt='("Starting outer iteration : ",i3)') it
   write(unit=stdout,fmt=11) j0_total, beta(0),eps(it)
11 format('Starting cost function: ' ,1PD15.8,', Gradient: ',1PD15.8,/,&
          'For this outer iteration, the targeted reduction is:       ',1PD15.8)
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'
   write(unit=stdout,fmt='(A)') 'Iter     Cost Function       Info Content'

   !-------------------------------------------------------------------------
   ! [2.0] Iteratively solve for minimum of cost function:
   !-------------------------------------------------------------------------
   do iter=1, ntmax(it)
      qhat(:,iter) = ghat / beta(iter-1)
      
      call da_calculate_gradj(it,iter,cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                              be%cv%size_jl, xbx,be,iv,qhat(:,iter),y,fhat,grid,config_flags)
       
    ! Apply Lanczos recurrence and orthonormalize new gradient (using modified Gramm-Schmidt)
    !----------------------------------------------------------------------------------------
      alpha(iter) = da_dot_cv(cv_size, qhat(:,iter), fhat, grid, mz, jp_start, jp_end)

      ghat        = fhat - alpha(iter)*qhat(:,iter) - beta(iter-1)*qhat(:,iter-1)

      if (grm_shmt_order .eq. 1) then
         do i = iter, 1, -1
            gdot = da_dot_cv(cv_size, ghat, qhat(:,i), grid, mz, jp_start, jp_end)
            ghat = ghat - gdot * qhat(:,i)
         end do
      else
         do j = 1,2 !Second iteration improves orthogonalization for numerical precision (Bjorck, 1994)
         do i = 1, iter !Reordering reduces numerical error (maybe not always?)
            gdot = da_dot_cv(cv_size, ghat, qhat(:,i), grid, mz, jp_start, jp_end)
            ghat = ghat - gdot * qhat(:,i)
         end do
         end do
      end if

      beta(iter)  = SQRT(da_dot_cv (cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
      
    ! Lanczos iteration  
    !------------------
      if (iter == 1) then
         d(1) = alpha(1)
	 c    = qhat(:,1)
	 rho  = beta(0) / alpha(1)
	 xhat = rho * qhat(:,1)
         dfs  = da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end)
         dfsmax = dfs
      else
         mu      = beta(iter-1) / d(iter-1)
	 d(iter) = alpha(iter) - beta(iter-1)*mu
	 c       = qhat(:,iter) - mu*c
	 rho     = - mu*d(iter-1)*rho / d(iter)
	 xhat    = xhat + rho*c
         dfs     = da_dot_cv(cv_size, rho*c, rho*c, grid, mz, jp_start, jp_end)
         if (dfs > dfsmax) dfsmax = dfs
      end if
      
    ! Determine eigenvalues and eigenvectors of the Lanczos tri-diagonal matrix
    !--------------------------------------------------------------------------

      info              = 0

!!#if defined(LAPACK)
!!      KMAT(iter,iter) = alpha(iter)
!!      if (iter .gt. 1) then
!!         KMAT (iter-1,iter) = beta(iter-1)
!!         KMAT (iter,iter-1) = beta(iter-1)
!!      end if
!!      eignvec(1:iter,1:iter)=KMAT(1:iter,1:iter)
!!      call SYEV(eignvec(1:iter,1:iter), eignval(1:iter), JOBZ = 'V' , INFO = info) !F95
!!
!!#else
      eignval(1:iter)   = alpha(1:iter)
      subdiag(1:iter-1) = beta(1:iter-1)
      nsstwrk           = MAX(2*iter-2,1)
#if defined(LAPACK)
      call RSTEQR(eignval(1:iter),subdiag(1:iter-1),eignvec(:,1:iter),'I', INFO = info) !F95
#else
      call DSTEQR('I',iter,eignval(1:iter),subdiag(1:iter-1),eignvec(:,1:iter),ntmax(it),&
	           sstwrk(1:nsstwrk),info)
#endif
      if (info /=0) write(stdout,*) 'Error in Lanczos minimization: SSTEQR returned ',info 

    ! Count converged eigenpairs (using Arnoldi relation)
    !----------------------------------------------------
      bndlm        = eps(it) * eignval(iter)  
      bnds(1:iter) = abs(beta(iter) * eignvec(iter,1:iter))
      neign        = COUNT(bnds(1:iter) <= bndlm)
               
    ! Print Cost Function and Infomation Content
    !-----------------------------------
      j_total = j0_total + 0.5 * da_dot_cv(cv_size,ghat0,xhat,grid,mz,jp_start,jp_end)
      write(unit=stdout,fmt='(A)') 'Iter     Cost Function       Info Content'
      write(unit=stdout,fmt=12)iter, j_total, dfs      	 
12    format(i3,5x,1PD15.8,5x,1PD15.8)

    ! Stop minimization based on Information Content
    ! ----------------------------------------------
      if (sqrt(dfs) < eps(it) * sqrt(dfsmax) .and. info_stop) exit

    ! Option to calculate J Cost Function explicitly
    ! ----------------------------------------------
      if (calculate_cg_cost_fn) &
         call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
	                     be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, fhat, grid, config_flags)

    ! Print current eigenvalues
    ! ----------------------------------------------
      write(unit=stdout,fmt='(A,2I5)') 'it, iter = ',it,iter
      write(dummychar,fmt='(I0)') iter
      write(unit=stdout,fmt='(A,'//trim(dummychar)//'E27.16E3)') 'Ritz eigenvalues: ',eignval(1:iter)     
!      do ii = 1,iter
!         write(unit=stdout,fmt='(A,I3,A,'//trim(dummychar)//'E27.16E3)') 'EIGENVECTOR',ii,': ',eignvec(ii,1:iter) 
!      end do
      write(unit=stdout,fmt='(A,I5)') 'Number of converged eigenpairs: ', neign

      if ( var4d_inc_out &
           .and. present(hessian_info) &
           .and. allocated(hessian_info%ghat) ) then
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

          if (iter .eq. 1) hessian_info%ghat = - ghat0

          hessian_info%qhat(:,iter) = qhat(:,iter)

          allocate(hessian_info%eigs(iter)%eignval(1:iter))
          allocate(hessian_info%eigs(iter)%eignvec(1:iter,1:iter))

          hessian_info%eigs(iter)%eignval = eignval(1:iter)
          hessian_info%eigs(iter)%eignvec = eignvec(1:iter,1:iter)
          hessian_info%eigs(iter)%neign = iter

          hessian_info%nmodes_max = hessian_info%nmodes_max + 1
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: increment_diagnostics1')
#endif
      end if

   end do

   !-------------------------------------------------------------------------
   ! End of the minimization of cost function
   !-------------------------------------------------------------------------
   iter = MIN(iter, ntmax(it))
   
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'
   write(unit=stdout,fmt='(A)') " "
   write(unit=stdout, fmt='("Inner iteration stopped after ",i4," iterations")') iter
   write(unit=stdout,fmt='(A)') " "

   call da_calculate_j(it, iter, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                       be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)

   write(unit=stdout,fmt=15) iter, j_cost%total , &
                             SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
15  format('Final: ',I3,' iter, J=',1PD15.8,', g=',1PD15.8)
   write(unit=stdout,fmt='(A)') &
      '----------------------------------------------------------'

   write(stdout,*) 'Ritz eigenvalues: ',eignval(iter:1:-1)     
   write(stdout,*) 'Number of converged eigenpairs: ', neign
   neign = iter

   if (spectral_precon) then
      !Sort eigenvalues and eigenvectors in descending order
      allocate(temp_e(1:neign))
      do j = neign-1, 1, -1
        swapped = .false.
        do i = 1, j
          if (eignval(i) .lt. eignval(i+1)) then
            temp_l = eignval(i)
            eignval(i) = eignval(i+1)
            eignval(i+1) = temp_l

            temp_e = eignvec(1:neign,i)
            eignvec(1:neign,i) = eignvec(1:neign,i+1)
            eignvec(1:neign,i+1) = temp_e

            swapped = .true.
          end if
        end do
        if (.NOT. swapped) exit
      end do
      deallocate(temp_e)

      call da_hessian_io ('w',grid, cv_size,it,neign,neign,&
                         eignvec(1:neign,1:neign),eignval(1:neign),qhat(:,0:neign),mz)
!      ntused(it) = neign

      if (it.gt.1) &
         call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 it-1, xhat, xhat, 'TL' )
   end if

   if (trace_use) call da_trace_exit("da_minimise_lz")

end subroutine da_minimise_lz

