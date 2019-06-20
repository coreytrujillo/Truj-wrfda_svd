subroutine da_setup_be_chem(grid, be)

   !---------------------------------------------------------------------------
   ! Purpose: Define and allocate components of background errors
   !
   ! Updates:
   !
   !       Implementation of chem BE
   !       Jonathan Guerrette,  CU Boulder,  Date: 07/01/2015
   !---------------------------------------------------------------------------

   implicit none
   type (domain),  intent(in)    :: grid
   type (be_type), intent(inout) :: be                    ! Back. errors structure.

   integer              :: ic, k1, k2, ix, jy, mz

   real*8, allocatable  :: ant_rf_lengthscale(:), ant_rf_scale_factor(:)
   real*8, allocatable  :: bb_rf_lengthscale(:), bb_rf_scale_factor(:)

   real*8,  allocatable :: eval(:), evec(:,:), covm(:,:)
   real*8               :: dtd, dth

   ix = ide - ids + 1
   jy = jde - jds + 1

   !! Anthropogenic
   mz = num_ant_steps
   do ic = PARAM_FIRST_SCALAR, num_scaleant
      ! Use all modes
      be % alpha_ant_chem(ic) % mz = mz

      if ( temporal_corr .eq. 0 .or. timescale_h_ant .le. 0 ) then
         ! Eigenvalues
         allocate ( be % alpha_ant_chem(ic) % val( 1:jy, 1 : mz) )

         be % alpha_ant_chem(ic) % val = sigma_ant  !|ln(0.5)| = ln(2) = 0.693 --> 100% uncertainty

         ! Eigenvectors
         allocate ( be % alpha_ant_chem(ic) % evec( 1:jy, 1:mz, 1:mz ) )

         be % alpha_ant_chem(ic) % evec = 0.D0
         do k1 = 1, mz
            be % alpha_ant_chem(ic) % evec(:,k1,k1) = 1.D0
         end do

         ! ant RF lengthscales
         allocate ( be % alpha_ant_chem(ic) % rf_alpha (1:mz))
         be % alpha_ant_chem(ic) % rf_alpha (1:mz) = 1.D0

      
      else if ( temporal_corr .ge. 1 ) then
         !! Temporal alpha correlation

         allocate ( eval(1:num_ant_steps) )
         allocate ( evec(1:num_ant_steps,1:num_ant_steps) )
         allocate ( covm(1:num_ant_steps,1:num_ant_steps) )

         ! Exponential covariance function
         dtd=0.D0
         covm=0.D0
         do k1 = 1, num_ant_steps

            covm(k1,k1) = sigma_ant ** 2

            do k2 = k1+1, num_ant_steps
               dth = abs( REAL((k1 - k2) * scaleant_frq, 8) ) / 60.D0
               if ( multiday_corr ) then
                  dtd = dth / 24.D0
                  dth = mod(dth,24.D0)
                  dth = min(24.D0 - dth, dth)
               end if
!!!               if ( dth / timescale_h_ant .lt. 4.D0 ) then
                  covm(k1,k2) = exp( -dth / timescale_h_ant - dtd / timescale_d_ant ) * sigma_ant ** 2
                  covm(k2,k1) = covm(k1,k2)
!!!               end if
            end do
         end do

         ! Positive-definite eigen decomposition
         ! da_gen_eigen truncates negative eigenvalues
         ! This approach loses accuracy in the U factorization if sigma_ant is too
         ! small relative to the hourly and daily timescales and scaleant_frq, 
         ! due to the introduction of negative eigenvalues. 
         call da_1d_eigendecomposition ( covm, evec, eval)

         if ( minval(eval) .lt. 0.D0 ) call da_gen_eigen (num_ant_steps, covm, evec, eval)

!         ! Negative eigenvalues would require complex be%be_sub%val for U factorization by
!         ! sqrt of diagonal eigenvalue matrix in da_allocate_background_errors_chem
!         call da_1d_eigendecomposition ( covm, evec, eval)

         call da_allocate_background_errors_chem (jy, mz, eval, evec, be % alpha_ant_chem(ic) )

         deallocate ( eval )
         deallocate ( evec )
         deallocate ( covm )

      end if

   end do

   !! Biomass Burning
   mz = num_bb_steps
   do ic = PARAM_FIRST_SCALAR, num_scalebb
      ! Use all modes
      be % alpha_bb_chem(ic) % mz = mz

      if ( temporal_corr .eq. 0 .or. timescale_h_bb .le. 0 ) then
         ! Eigenvalues
         allocate ( be % alpha_bb_chem(ic) % val( 1:jy, 1 : mz) )

         be % alpha_bb_chem(ic) % val = sigma_bb  !|ln(0.5)| = ln(2) = 0.693 --> 100% uncertainty

         ! Eigenvectors
         allocate ( be % alpha_bb_chem(ic) % evec( 1:jy, 1:mz, 1:mz ) )

         be % alpha_bb_chem(ic) % evec = 0.D0
         do k1 = 1, mz
            be % alpha_bb_chem(ic) % evec(:,k1,k1) = 1.D0
         end do

         ! bb RF lengthscales
         allocate ( be % alpha_bb_chem(ic) % rf_alpha (1:mz))
         be % alpha_bb_chem(ic) % rf_alpha (1:mz) = 1.D0

      
      else if ( temporal_corr .ge. 1 ) then
         !! Temporal alpha correlation

         allocate ( eval(1:num_bb_steps) )
         allocate ( evec(1:num_bb_steps,1:num_bb_steps) )
         allocate ( covm(1:num_bb_steps,1:num_bb_steps) )

         ! Exponential covariance function
         dtd=0.D0
         covm=0.D0
         do k1 = 1, num_bb_steps

            covm(k1,k1) = sigma_bb ** 2

            do k2 = k1+1, num_bb_steps
               dth = abs( REAL((k1 - k2) * scalebb_frq, 8) ) / 60.D0
               if ( multiday_corr ) then
                  dtd = dth / 24.D0
                  dth = mod(dth,24.D0)
                  dth = min(24.D0 - dth, dth)
               end if
!!!               if ( dth / timescale_h_bb .lt. 4.D0 ) then
                  covm(k1,k2) = exp( -dth / timescale_h_bb - dtd / timescale_d_bb ) * sigma_bb ** 2
                  covm(k2,k1) = covm(k1,k2)
!!!               end if
            end do
         end do

         ! Positive-definite eigen decomposition 
         ! da_gen_eigen truncates negative eigenvalues
         ! This approach loses accuracy in the U factorization if sigma_bb is too
         ! small relative to the hourly and daily timescales and scalebb_frq, 
         ! due to the introduction of negative eigenvalues. 
         call da_1d_eigendecomposition ( covm, evec, eval)

         if ( minval(eval) .lt. 0.D0 ) call da_gen_eigen (num_bb_steps, covm, evec, eval)

!         ! Negative eigenvalues would require complex be%be_sub%val for U factorization by
!         ! sqrt of diagonal eigenvalue matrix in da_allocate_background_errors_chem
!         call da_1d_eigendecomposition ( covm, evec, eval)

         call da_allocate_background_errors_chem (jy, mz, eval, evec, be % alpha_bb_chem(ic) )

         deallocate ( eval )
         deallocate ( evec )
         deallocate ( covm )

      end if
   end do

   ! As more emission species are defined, will need to create/populate regression table 
   ! describing correlations across species for each source type

   ! Now prepare rf_alpha and val for the recursive filter.
   if ( use_rf ) then
! This code should be equivalent to that below it.
!      allocate ( ant_rf_lengthscale( be % alpha_ant_chem(PARAM_FIRST_SCALAR) % mz) )
!      allocate ( bb_rf_lengthscale( be % alpha_bb_chem(PARAM_FIRST_SCALAR) % mz) )
!
!      ant_rf_lengthscale(:) = rf_lengthscale_ant; !meters
!      bb_rf_lengthscale(:)  = rf_lengthscale_bb; !meters
!
!      do ic = PARAM_FIRST_SCALAR, num_scaleant
!         call da_rescale_background_errors ( var_scaling_ant(1), len_scaling_ant(1),  &
!                                             grid%dx, ant_rf_lengthscale, be % alpha_ant_chem(ic) )
!      end do
!
!      do ic = PARAM_FIRST_SCALAR, num_scalebb
!         call da_rescale_background_errors ( var_scaling_bb(1), len_scaling_bb(1),  &
!                                             grid%dx, bb_rf_lengthscale, be % alpha_bb_chem(ic) )
!      end do
!      deallocate ( ant_rf_lengthscale )
!      deallocate ( bb_rf_lengthscale )


      allocate ( ant_rf_lengthscale( be % alpha_ant_chem(PARAM_FIRST_SCALAR) % mz) )
      allocate ( bb_rf_lengthscale( be % alpha_bb_chem(PARAM_FIRST_SCALAR) % mz) )
      allocate ( ant_rf_scale_factor( be % alpha_ant_chem(PARAM_FIRST_SCALAR) % mz) )
      allocate ( bb_rf_scale_factor( be % alpha_bb_chem(PARAM_FIRST_SCALAR) % mz) )

      ant_rf_lengthscale(:) = rf_lengthscale_ant / grid%dx; !Non-dimensionalized from meters
      bb_rf_lengthscale(:)  = rf_lengthscale_bb / grid%dx;  !Non-dimensionalized from meters

      do ic = PARAM_FIRST_SCALAR, num_scaleant
         call da_calculate_rf_factors ( ant_rf_lengthscale, be % alpha_ant_chem(ic) % rf_alpha (:), &
                                        ant_rf_scale_factor (:) )

         mz = be % alpha_ant_chem(ic) % mz
         do k1 = 1, mz
            be % alpha_ant_chem(ic) % val(:,k1) = ant_rf_scale_factor(k1) *  be % alpha_ant_chem(ic) % val(:,k1)
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         call da_calculate_rf_factors ( bb_rf_lengthscale, be % alpha_bb_chem(ic) % rf_alpha (:), &
                                        bb_rf_scale_factor (:) )
         mz = be % alpha_bb_chem(ic) % mz
         do k1 = 1, mz
            be % alpha_bb_chem(ic) % val(:,k1) = bb_rf_scale_factor(k1) *  be % alpha_bb_chem(ic) % val(:,k1)
         end do
      end do

      deallocate ( ant_rf_lengthscale )
      deallocate ( bb_rf_lengthscale )
      deallocate ( ant_rf_scale_factor )
      deallocate ( bb_rf_scale_factor )

   end if

   if (trace_use) call da_trace_exit("da_setup_be_chem")

end subroutine da_setup_be_chem


