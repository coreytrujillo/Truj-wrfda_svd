#if defined(LAPACK)
subroutine da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega, yhat, qhat, cv_size, &
                     eignval, eignvec, nsamp, print_messages, hess_trace_approx)

   !-------------------------------------------------------------------------
   ! Purpose:         Second stage of RSVD 5.6
   !
   ! Called from da_randomise_svd and da_compare_decomp_methods
   !
   ! Solve nearly-symmetric low-dimensional system A * x = b, &
   !  where KMAT * F1^T = F2^T, F1^T = Q^T * OMEGA, F2^T = Q^T * Y
   !  or    F1 * KMAT^T = F1 * KMAT = F2
   !  D = F1
   !  b = F2(:,iens)
   !  x = KMAT(iens,:)
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   ! References: Halko et al. (2011)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(in) :: grid
   integer, intent(in)      :: jp_start, jp_end
   integer, intent(in)      :: mz(:)

   integer, intent(in)      :: cv_size        ! Total cv size

   integer, intent(in)     :: nsamp                     ! Number of total eigenpairs
   real, intent(in)        :: yhat(1:cv_size, 1:nsamp)  ! yhat = A * omega
   real, intent(in)        :: omega(1:cv_size, 1:nsamp) ! Gaussian draw
   real, intent(in)        :: qhat(1:cv_size, 0:nsamp)  ! orthonormal basis set 

   real, intent(out)       :: eignvec(nsamp, nsamp)
   real, intent(out)       :: eignval(nsamp)
   real, optional, intent(out) :: hess_trace_approx(nsamp)

   logical, optional, intent(in) :: print_messages

   real(wp)                     :: LAM(nsamp)
   integer                      :: iens, jens, i, j, p1
   real(wp)                     :: F1(nsamp, nsamp) ! intermediate matrix
   real(wp)                     :: F2(nsamp, nsamp) ! intermediate matrix
   real(wp)                     :: U(nsamp, nsamp)
   real(wp)                     :: VT(nsamp, nsamp)
   real(wp)                     :: SIG(nsamp)
   real(wp)                     :: UTF2(nsamp,nsamp), VSIGinvT(nsamp,nsamp)
   real(wp)                     :: KMAT(nsamp, nsamp)  ! KMAT = Q^T * A * Q

   integer                      :: piv(nsamp)
   integer                      :: info

   real(wp)                     :: temp_l
   real(wp)                     :: temp_e(nsamp)
   character(len=4) :: cnsamp=""
   logical          :: swapped
  
   if (trace_use) call da_trace_entry("da_rsvd56")

   if(.not.present(print_messages)) write(unit=stdout,fmt='(A,I5,A)') 'Performing Stage 2 of RSVD5.6 for ', nsamp ,' ensembles'
   if(.not.present(print_messages)) write(unit=stdout,fmt=*) ' '

   do iens = 1, nsamp
      do jens = 1, nsamp
         !F1 = OMEGA^T * Q
         F1(iens,jens) = da_dot_cv(cv_size, omega(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)

         !F2 = Y^T * Q
         F2(iens,jens) = da_dot_cv(cv_size, yhat(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
      end do
      if (present(hess_trace_approx)) then
         if (iens .eq. 1) then
            hess_trace_approx(iens) = &
               da_dot_cv(cv_size, omega(:,iens), yhat(:,iens), grid, mz, jp_start, jp_end)
         else
            hess_trace_approx(iens) = ( &
               hess_trace_approx(iens-1) * real(iens-1,8) + &
               da_dot_cv(cv_size, omega(:,iens), yhat(:,iens), grid, mz, jp_start, jp_end) &
               )  / real(iens,8)
         end if
      end if
   end do

   write(cnsamp,fmt='(I0)') nsamp

   info = 0

   !! Solve F1*KMAT = F2

   ! DECOMPOSE F1 = U * SIG * VT 
   call GESVD(F1(1:nsamp,1:nsamp), &
              SIG(1:nsamp), &
              U = U(1:nsamp,1:nsamp), &
              VT = VT(1:nsamp,1:nsamp), &
              INFO = info)
#ifdef DM_PARALLEL
   call mpi_barrier(comm, ierr)
#endif
   if (info .ne. 0) then
      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'F1: ',F1(iens,1:nsamp)
      end do

      write(stdout,*) 'Error in randomised SVD: GESVD returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if

   ! KMAT = VT^T * SIG^-1 * U^T * F2 
   do iens = 1, nsamp
      VSIGinvT(iens,1:nsamp) = VT(iens,1:nsamp) / SIG(iens)
      do jens = 1, nsamp
         UTF2(iens,jens) = da_dot(nsamp,U(1:nsamp,iens),F2(1:nsamp,jens))
      end do
   end do

   do iens = 1, nsamp
      do jens = 1, nsamp
         KMAT(jens,iens) = da_dot(nsamp,VSIGinvT(:,iens),UTF2(:,jens))
      end do
   end do


   ! Force symmetry in KMAT
   KMAT = (KMAT + TRANSPOSE(KMAT)) / 2.D0

 ! Determine eigenvalues and eigenvectors of the symmetric KMAT matrix
 !-----------------------------------------------------------------
   info = 0
   call SYEV(KMAT, LAM, JOBZ = 'V' , INFO = info) !F95
#ifdef DM_PARALLEL
   call mpi_barrier(comm, ierr)
#endif
   if (info .ne. 0) then
      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'F2: ',F2(iens,1:nsamp)
      end do

      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'KMAT: ',KMAT(iens,1:nsamp)
      end do

      write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if

   do iens = 1, nsamp
      do jens = 1, nsamp 
         eignvec(iens,jens) = KMAT(iens,jens)
      end do
   end do

   temp_e = 1.D0
   if (all(LAM .le. 0.D0)) then
      if(.not.present(print_messages)) write(unit=stdout,fmt='(A)') &
         'WARNING: Preconditioned Hessian: all(eigs <= 0), setting them positive.'
      LAM = - LAM
   else if ( sum(sign(temp_e,LAM)) .lt. 0.D0 ) then
      if(.not.present(print_messages)) write(unit=stdout,fmt='(A,F10.1,A)') &
         'WARNING: Preconditioned Hessian: More negative than positive eigenvalues,', sum(sign(temp_e,LAM)), ', setting LAM=-LAM'
      LAM = - LAM
   end if
  
   !Sort eigenvalues and eigenvectors in descending order
   do j = nsamp-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (LAM(i) .lt. LAM(i+1)) then
         temp_l = LAM(i)
         LAM(i) = LAM(i+1)
         LAM(i+1) = temp_l

         temp_e = eignvec(:,i)
         eignvec(:,i) = eignvec(:,i+1)
         eignvec(:,i+1) = temp_e

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do

   !Preconditioned observation Hessian should be positive semi-definite
!   neign = nsamp
   do iens = 1, nsamp
      LAM(iens) = max(LAM(iens),0.D0)
   end do
   eignval = LAM

   if (trace_use) call da_trace_exit("da_rsvd56")

end subroutine da_rsvd56
#endif
