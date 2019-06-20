subroutine da_output_increments(grid,config_flags,mz, &
                                it,cv_size,xbx,be, &
                                cvt,hessian_info)

   !-------------------------------------------------------------------------
   ! Purpose:   Output increment for all inner loop truncations using Hessian
   !            spectra and qhat basis vectors.  This subroutine is for
   !            diagnostic purposes and should not be used when computational
   !            expense need be minimized.
   !
   ! Called from da_solve
   !
   ! History: 10/31/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   type (grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: mz(:)

   integer,         intent(in)   :: it    ! outer-loop index
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   real*8, intent(in)            :: cvt(1:cv_size)
   type(hessian_type), intent(in) :: hessian_info

!   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
!   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   integer :: i, j, k, ic, imode, nmodes, nmodes_REF, neign, itype, iend, dmode, icount
   real*8, allocatable :: LAM(:)
   real*8, allocatable :: eignval(:)
   real*8, allocatable :: eignvec(:,:)
   real*8, allocatable, target :: qhat(:,:)
   real*8, allocatable :: nuhat(:,:)
   real*8, allocatable :: nuhat_REF(:,:)

   real*8 :: ghat(1:cv_size)
   real*8 :: xhat(1:cv_size)
   real*8 :: cvt_REF(1:cv_size)
   real*8 :: dvnorm_(1:cv_size)
   real*8 :: EUCNORM_xPa(2), EUCNORM_v(2), EUCNORM_x(2), EUCNORM_xt(2), EUCNORM_xtPa(2), SUMLAM, SUMLAM3, temp_real
!   real*8 :: minlam(1:1)
   type(hessian_eig_type) :: hess_REF
   real*8, pointer :: qhat_REF(:,:)
   real*8, allocatable :: LAM_REF(:)

   character(len=10) :: PREFIX
   character(len=256) :: timestr
!   character(len=256) :: dummy_name
   character(len=filename_len)     :: filename
   character(len=filename_len)     :: cvtfile
   character*10                    :: cproc
   character*10                    :: cmode
   character*10                    :: cit

   real, allocatable                 :: temp_e(:)
   real                              :: temp_l
   logical                           :: swapped

   if (trace_use) call da_trace_entry("da_output_increments")

!   je_start   = be % cv % size_jb + 1
!   je_end     = be % cv % size_jb + be % cv % size_je
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
!   jl_start   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + 1
!   jl_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + be % cv % size_jl

   nmodes = hessian_info%nmodes_max

   ghat = hessian_info%ghat

   allocate(qhat(cv_size,0:nmodes))
   allocate(LAM(1:nmodes))

!   dummy_name = trim(config_flags%auxhist18_outname)
!   call nl_set_io_form_auxhist18 ( grid%id, 2 )
!   config_flags%io_form_auxhist18 = 2
!   config_flags%frames_per_auxhist18 = 1

   iend = 2
   if (use_lanczos) iend = 1
!   if (use_lanczos .or. rand_type.eq.3) iend = 1

   do imode = 1, nmodes
      qhat(:,imode) = hessian_info%qhat(:,imode)
   end do

   icount = 0

   write(unit=stdout,fmt='(A)') &
      'EUCNORM_BY_ITER_AND_SAMPLE: imode,neign,min(lambda),EUCNORM_xPa[LRA,LRU],EUCNORM_v[LRA,LRU],EUCNORM_x[LRA,LRU],EUCNORM_xt[LRA,LRU],EUCNORM_xtPa[LRA,LRU] = '

   dmode = min(inc_out_interval,nmodes)
   if (use_randomblock .and. rand_type.eq.3) then
      dmode = nmodes
   end if
   do imode = nmodes, 0, -dmode
      allocate(eignval(1:imode))
      allocate(eignvec(1:imode,1:imode))

      if (imode .gt. 0) then 
         if (use_randomblock .and. rand_type.eq.3) then
            neign = hessian_info%eigs(1)%neign
            eignval = hessian_info%eigs(1)%eignval
            eignvec = hessian_info%eigs(1)%eignvec
         else
            neign = hessian_info%eigs(imode)%neign
            eignval = hessian_info%eigs(imode)%eignval
            eignvec = hessian_info%eigs(imode)%eignvec
         end if
      else
         neign = 0
      end if

      EUCNORM_xPa = 0.D0
      EUCNORM_v = 0.D0
      EUCNORM_x = 0.D0
      EUCNORM_xt = 0.D0
      EUCNORM_xtPa = 0.D0

      do itype = 1, iend
         icount = icount+ 1
         xhat = 0.D0


         if (imode .gt. 0) then
            if (itype .eq. 2) then
               do i = 1, neign
                  LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
               end do

               PREFIX = 'LRU'
            else
               do i = 1, neign
                  LAM(i) = 1.D0 / eignval(i)
               end do

               PREFIX = 'LRA'
               if (use_lanczos) PREFIX = 'CGLZ'
            end if
            write(unit=filename,fmt='(2a,i2.2,a,i4.4)') trim(PREFIX),'inc_it',it,'_nmode',imode

            call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                   imode,imode,neign, &
                                   LAM(1:imode),eignvec,qhat(:,0:imode), &
                                   ghat,xhat)

         end if

         if (itype .eq. 2) then !LRU
            xhat = xhat + ghat
         end if

         if (spectral_precon .and. it.gt.1) &
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                    it-1, xhat, xhat, 'TL' )

         !-----------------------------------------
         ! Calculate ||dx - dx_REF||^2_{Pa_REF^-1}
         !-----------------------------------------
         !Define reference case
         if (icount .eq. 1) then
            cvt_REF = cvt + xhat

            write(cproc,fmt='(i4.4)') myproc
            write(cit,fmt='(i2.2)') it
            write(cmode,fmt='(i4.4)') nmodes
            filename = '../hess_REF_out.nmodes'//trim(adjustl(cmode))//'.it'//trim(adjustl(cit))//'.'//trim(adjustl(cproc))
            call da_hessian_io ('w',grid,cv_size,it,nmodes,neign, &
                         eignvec,eignval,qhat,mz, &
                         cvt_REF,hess_file=filename)

            if (read_cvt_REF.eq.1) then
               cvtfile = 'cvt_REF_in'
               call da_cv_io ('r', cvtfile, cv_size, cvt_REF(1:cv_size), grid, mz)
            end if

            if (read_hess_REF.gt.0) then
               !Allows hess_REF and cvt_REF to be defined in a separate inversion scenario
               ! (read in from file)
               if (read_hess_REF.eq.1 .or. read_hess_REF.eq.3) then
                  filename = 'hess_REF_in'
!               else
!                  !Debug setting to compare against read_hess_REF==0
!                  !Unchanged filename
               end if

               !Read nmodes_REF from file
               call da_hessian_io ('r',grid,cv_size,it,nmodes_REF,hess_REF%neign,hess_file=filename)

               if (read_hess_REF.ne.3) &
                  write(unit=stdout,fmt='(A,I8,2A)') &
                           'Bayes Risk Reference: Utilizing ',nmodes_REF,' modes from ',filename

               allocate(qhat_REF(1:cv_size,0:nmodes_REF))
               allocate(hess_REF%eignval(1:nmodes_REF))
               allocate(hess_REF%eignvec(1:nmodes_REF,1:nmodes_REF))

               qhat_REF(:,0) = 0.D0
               if (read_cvt_REF.eq.2) then
                  call da_hessian_io ('r',grid,cv_size,it,nmodes_REF,hess_REF%neign, &
                               hess_REF%eignvec,hess_REF%eignval,qhat_REF,mz, &
                               cvt_REF,hess_file=filename)
               else
                  call da_hessian_io ('r',grid,cv_size,it,nmodes_REF,hess_REF%neign, &
                               hess_REF%eignvec,hess_REF%eignval,qhat_REF,mz, &
                               hess_file=filename)
               end if
               if (read_hess_REF.eq.3) then
                  deallocate(hess_REF%eignval)
                  deallocate(hess_REF%eignvec)
                  deallocate(qhat_REF)
               end if
            end if
            if (read_hess_REF.le.0 .or. read_hess_REF.eq.3) then
               nmodes_REF = nmodes

               write(unit=stdout,fmt='(A,I8,2A)') &
                        'Bayes Risk Reference: Utilizing ',nmodes_REF,' modes from this outer iteration.'

               qhat_REF => qhat
               allocate(hess_REF%eignval(1:nmodes_REF))
               allocate(hess_REF%eignvec(1:nmodes_REF,1:nmodes_REF))

               hess_REF%neign = neign
               hess_REF%eignval = eignval
               hess_REF%eignvec = eignvec
            end if

            allocate(LAM_REF(1:nmodes_REF))
            do i = 1, hess_REF%neign
               LAM_REF(i) = hess_REF%eignval(i) - 1.D0
            end do

            allocate(nuhat_REF(cv_size,0:nmodes_REF))
            nuhat_REF = 0.D0
            do j = 1, nmodes_REF
               do i = 1, nmodes_REF
                  nuhat_REF(:,j) = nuhat_REF(:,j) + qhat_REF(:,i) * hess_REF%eignvec(i,j)
               end do
            end do

            if (ASSOCIATED(qhat_REF,qhat)) then
               NULLIFY(qhat_REF)
            else
               deallocate(qhat_REF)
            end if
         end if

#if (WRF_CHEM == 1)
         !Only for chem so far
         if (use_chemobs .and. chem_opt(1).gt.0 .and. osse_chem) then
            call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                    xhat(1:be%cv%size_jb),grid%vv,grid%vp &
!#if (WRF_CHEM == 1)
                                    , grid%vch &
!#endif
                                    )

            temp_real = 0.D0
            do ic = PARAM_FIRST_SCALAR, num_scaleant
            do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
               temp_real = temp_real + &
                           (grid%xbch%alpha_ant(i,j,k,ic) + &
                            grid%xach%alpha_ant(i,j,k,ic)) ** 2
            end do
            end do
            end do
            end do

            do ic = PARAM_FIRST_SCALAR, num_scalebb
            do k = 1, num_bb_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + &
                           (grid%xbch%alpha_bb(i,j,k,ic) + &
                            grid%xach%alpha_bb(i,j,k,ic)) ** 2
            end do
            end do
            end do
            end do
#ifdef DM_PARALLEL
            EUCNORM_xt(itype) = wrf_dm_sum_real(temp_real)
#else
            EUCNORM_xt(itype) = temp_real
#endif

!            !Output increment for each imode value
!            call domain_clock_get( grid, current_timestr=timestr )
!            config_flags%auxhist18_outname = trim(filename)
!            call da_transfer_xatowrftl(grid, config_flags, 'DUMM', timestr, .true.)

            !Current increment
!            inc_current = cvt + xhat
!           Convert to ||x-xt||_Pa^-1 somehow
!           Might need cvt_true...right?
!           I know the true value of xbch + xach = 0, but not cvt.  
!           That value is stored in reference directory...
         end if
#endif

         

         !Difference between current increment and REF increment
         xhat = cvt + xhat - cvt_REF

         !Multiply Pav_REF^-1 * (dv - dv_REF) in two steps
         !(1) (L^T H^T R^-1 H L)_REF * diffvec
         dvnorm_ = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                nmodes_REF,nmodes_REF,hess_REF%neign, &
                                LAM_REF,hess_REF%eignvec,nuhat_REF, &
                                xhat,dvnorm_,.true.)

         !(2) (I + L^T H^T R^-1 H L) * diffvec
         dvnorm_ = dvnorm_ + xhat

         !Multiply (dv - dv_REF)^T * Pav_REF^-1 * (dv - dv_REF) = ||dx - dx_REF||^2_{Pa_REF^-1}
         EUCNORM_xPa(itype) = da_dot_cv(cv_size, xhat, dvnorm_, grid, mz, jp_start, jp_end)

         EUCNORM_v(itype) = da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end)

#if (WRF_CHEM == 1)
         !Only for chem so far
         if (use_chemobs .and. chem_opt(1) .gt. 0) then
            call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                    xhat(1:be%cv%size_jb),grid%vv,grid%vp &
!!#if (WRF_CHEM == 1)
                                    , grid%vch &
!!#endif
                                    )

            temp_real = 0.D0
            do ic = PARAM_FIRST_SCALAR, num_scaleant
            do k = 1, num_ant_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + grid%xach%alpha_ant(i,j,k,ic)**2
            end do
            end do
            end do
            end do

            do ic = PARAM_FIRST_SCALAR, num_scalebb
            do k = 1, num_bb_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + grid%xach%alpha_bb(i,j,k,ic)**2
            end do
            end do
            end do
            end do
#ifdef DM_PARALLEL
            EUCNORM_x(itype) = wrf_dm_sum_real(temp_real)
#else
            EUCNORM_x(itype) = temp_real
#endif
         end if
#endif
      end do

      write(unit=stdout,fmt='(A,2I6,11E27.16E3)') &
         'EUCNORM_BY_ITER_AND_SAMPLE: ',&
         imode,neign,minval(eignval),EUCNORM_xPa,EUCNORM_v,EUCNORM_x,EUCNORM_xt,EUCNORM_xtPa

      deallocate(eignval)
      deallocate(eignvec)

   end do


   !Calculate Euclidean norm by modes
   if (nmodes_final .gt. 0) nmodes = nmodes_final
   allocate(eignval(1:nmodes))
   allocate(eignvec(1:nmodes,1:nmodes))

   if (use_randomblock .and. rand_type.eq.3) then
      neign = hessian_info%eigs(1)%neign
      eignval = hessian_info%eigs(1)%eignval
      eignvec = hessian_info%eigs(1)%eignvec
   else
      neign = hessian_info%eigs(nmodes)%neign
      eignval = hessian_info%eigs(nmodes)%eignval
      eignvec = hessian_info%eigs(nmodes)%eignvec
   end if

   !Sort eigenvalues and eigenvectors in descending order
   allocate(temp_e(1:nmodes))
   do j = nmodes-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (eignval(i) .lt. eignval(i+1)) then
         temp_l = eignval(i)
         eignval(i) = eignval(i+1)
         eignval(i+1) = temp_l

         temp_e = eignvec(:,i)
         eignvec(:,i) = eignvec(:,i+1)
         eignvec(:,i+1) = temp_e

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do
   deallocate(temp_e)

   allocate(nuhat(cv_size,0:nmodes))
   nuhat = 0.D0
   do j = 1, nmodes
      do i = 1, nmodes
         nuhat(:,j) = nuhat(:,j) + qhat(:,i) * eignvec(i,j)
      end do
   end do

   write(unit=stdout,fmt='(A)') &
      'EUCNORM_BY_FINAL_MODES: imode,neign,lambda(i),EUCNORM_xPa[LRA,LRU],EUCNORM_v[LRA,LRU],EUCNORM_x[LRA,LRU],EUCNORM_xt[LRA,LRU],EUCNORM_xtPa[LRA,LRU], SUM(LAM), SUM(LAM^3) = '

!   dmode = floor(real(nmodes) / 400.D0) + 1
!   dmode = 1
   dmode = min(inc_out_interval,nmodes)

   do imode = 0, nmodes, dmode
      EUCNORM_xPa = 0.D0
      EUCNORM_v = 0.D0
      EUCNORM_x = 0.D0
      EUCNORM_xt = 0.D0
      EUCNORM_xtPa = 0.D0

      if (imode .gt. 0) then
         SUMLAM = sum(eignval(1:imode) - 1.D0)
         SUMLAM3 = sum((eignval(1:imode) - 1.D0)**3)
      else
         SUMLAM = 0.D0 
         SUMLAM3 = 0.D0 
      end if 
      do itype = 1, iend
         xhat = 0.D0

         if (imode .gt. 0) then
            if (itype .eq. 2) then
               do i = 1, nmodes
                  LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
               end do

               PREFIX = 'LRU'
            else
               do i = 1, nmodes
                  LAM(i) = 1.D0 / eignval(i)
               end do

               PREFIX = 'LRA'
               if (use_lanczos) PREFIX = 'CGLZ'
            end if
            write(unit=filename,fmt='(2a,i2.2,a,i4.4)') trim(PREFIX),'inc_it',it,'_nmode',imode

            call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                   nmodes,nmodes,imode, &
                                   LAM,eignvec,nuhat, &
                                   ghat,xhat,.true.)
         end if

         if (itype .eq. 2) then !LRU
            xhat = xhat + ghat
         end if

         if (spectral_precon .and. it.gt.1) &
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                    it-1, xhat, xhat, 'TL' )

#if (WRF_CHEM == 1)
         !Only for chem so far
         if (use_chemobs .and. chem_opt(1).gt.0 .and. osse_chem) then
            call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                    xhat(1:be%cv%size_jb),grid%vv,grid%vp &
!#if (WRF_CHEM == 1)
                                    , grid%vch &
!#endif
                                    )

            temp_real = 0.D0
            do ic = PARAM_FIRST_SCALAR, num_scaleant
            do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
               temp_real = temp_real + &
                           (grid%xbch%alpha_ant(i,j,k,ic) + &
                            grid%xach%alpha_ant(i,j,k,ic)) ** 2
            end do
            end do
            end do
            end do

            do ic = PARAM_FIRST_SCALAR, num_scalebb
            do k = 1, num_bb_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + &
                           (grid%xbch%alpha_bb(i,j,k,ic) + &
                            grid%xach%alpha_bb(i,j,k,ic)) ** 2
            end do
            end do
            end do
            end do
#ifdef DM_PARALLEL
            EUCNORM_xt(itype) = wrf_dm_sum_real(temp_real)
#else
            EUCNORM_xt(itype) = temp_real
#endif
         end if
#endif

         !Difference between current increment and REF increment
         xhat = cvt + xhat - cvt_REF

         !Multiply Pav_REF^-1 * (dv - dv_REF) in two steps
         !(1) (L^T H^T R^-1 H L)_REF * diffvec
         dvnorm_ = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                nmodes_REF,nmodes_REF,hess_REF%neign, &
                                LAM_REF,hess_REF%eignvec,nuhat_REF, &
                                xhat,dvnorm_,.true.)

         !(2) (I + L^T H^T R^-1 H L) * diffvec
         dvnorm_ = dvnorm_ + xhat

         !Multiply (dv - dv_REF)^T * Pav_REF^-1 * (dv - dv_REF) = ||dx - dx_REF||^2_{Pa_REF^-1}
         EUCNORM_xPa(itype) = da_dot_cv(cv_size, xhat, dvnorm_, grid, mz, jp_start, jp_end)

         EUCNORM_v(itype) = da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end)

#if (WRF_CHEM == 1)
         !Only for chem so far
         if (use_chemobs .and. chem_opt(1).gt.0) then
            call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                    xhat(1:be%cv%size_jb),grid%vv,grid%vp &
!!#if (WRF_CHEM == 1)
                                    , grid%vch &
!!#endif
                                    )

            temp_real = 0.D0
            do ic = PARAM_FIRST_SCALAR, num_scaleant
            do k = 1, num_ant_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + grid%xach%alpha_ant(i,j,k,ic)**2
            end do
            end do
            end do
            end do

            do ic = PARAM_FIRST_SCALAR, num_scalebb
            do k = 1, num_bb_steps
            do j=jts,jte
            do i=its,ite
               temp_real = temp_real + grid%xach%alpha_bb(i,j,k,ic)**2
            end do
            end do
            end do
            end do
#ifdef DM_PARALLEL
            EUCNORM_x(itype) = wrf_dm_sum_real(temp_real)
#else
            EUCNORM_x(itype) = temp_real
#endif
         end if
#endif
      end do

      write(unit=stdout,fmt='(A,2I6,13E27.16E3)') &
         'EUCNORM_BY_FINAL_MODES: ',&
         imode,neign,eignval(max(1,imode)),EUCNORM_xPa,EUCNORM_v,EUCNORM_x,EUCNORM_xt,EUCNORM_xtPa,SUMLAM,SUMLAM3
   end do

   deallocate(nuhat)
   deallocate(eignval)
   deallocate(eignvec)

   deallocate(LAM_REF)
   deallocate(hess_REF%eignval)
   deallocate(hess_REF%eignvec)
   deallocate(nuhat_REF)

   deallocate(LAM)
   deallocate(qhat)


!   config_flags%auxhist18_outname = trim(dummy_name)
!   call nl_set_io_form_auxhist18 ( grid%id, 0 )
!   config_flags%io_form_auxhist18 = 0

   if (trace_use) call da_trace_exit("da_output_increments")

end subroutine da_output_increments