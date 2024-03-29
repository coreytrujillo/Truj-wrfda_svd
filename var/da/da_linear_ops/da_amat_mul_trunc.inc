subroutine da_amat_mul_trunc(grid, mz, jp_start, jp_end, cv_size, &
                             nmodes, neign, ntrunc, &
                             eignval, eignvec, qhat, &
                             shat, xhat, nuflag)

   !-------------------------------------------------------------------------
   ! Purpose:  Multiply a control vector by the Analysis Error Cov Matrix A 
   !
   ! Called from da_output_increments, da_randomise_svd, and others
   !
   ! History: 09/01/2016  Creation (JJ Guerrette)
   !
   ! Modified from da_amat_mul
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: nmodes
   integer,        intent(in) :: neign
   integer,        intent(in) :: ntrunc
   real,         intent(in) :: eignvec(nmodes, nmodes)  
   real,         intent(in) :: eignval(nmodes)
   real,         intent(in) :: qhat(cv_size, 0:nmodes) ! Ritz vectors
   real,         intent(in) :: shat(cv_size)          ! Input vector to multiply by A
   real,      intent(inout) :: xhat(cv_size)          ! Output vector: xhat = A.shat
   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)
   logical, intent(in), optional :: nuflag

   real                     :: nuhat(1:cv_size,1:ntrunc)
   real                     :: gam(ntrunc), gam_(ntrunc), temp_l, temp_e(cv_size)
   integer                    :: i, j, order(ntrunc), temp_i
!   integer                    :: start, fin, step
   real                     :: dot_cv

   logical                    :: swapped
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   if (trace_use) call da_trace_entry("da_amat_mul_trunc")

   write(cproc,fmt='(i4.4)') myproc

   if (present(nuflag)) then
      nuhat = qhat(:,1:ntrunc)
   else
      nuhat = 0.D0
      do j = 1, ntrunc
         ! Is this slow due to non-parallel da_dot across cv_size?
!         do i = 1, cv_size
!            nuhat(i,j) = da_dot(neign, qhat(i,1:neign), eignvec(1:neign,j))

!          ! Is this faster? Same precision?
!          nuhat(i,j) = SUM( qhat(i,1:neign) * eignvec(1:neign,j) )

!         nuhat(:,j) = 0.D0
         ! Is this faster? Same precision?
         do i = 1, neign
            nuhat(:,j) = nuhat(:,j) + qhat(:,i) * eignvec(i,j)
         end do
      end do
   end if

   call mpi_barrier(comm, ierr)

   !Sort gamma values in descending order to improve mat-vec precision
   do i = 1, ntrunc 
      dot_cv = da_dot_cv(cv_size, shat, nuhat(:,i), grid, mz, jp_start, jp_end)
      gam(i) = eignval(i) * dot_cv
      order(i) = i
   end do

   ! This summation reordering significantly improves numerical precision
   gam_ = gam
   do j = ntrunc-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (abs(gam_(i)) .lt. abs(gam_(i+1))) then
         temp_l = gam_(i)
         gam_(i) = gam_(i+1)
         gam_(i+1) = temp_l

         temp_i = order(i)
         order(i) = order(i+1)
         order(i+1) = temp_i

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do

   !!Initial xhat value defined in calling subroutine
   do i = ntrunc, 1, -1 !From smallest to largest gamma
      xhat = xhat + gam(order(i)) * nuhat(:,order(i))
   end do

!   !This summation reordering eliminates some precision issues when lam_k>lam_k+1 or lam_k<lam_k+1. Not as precise as above, but slightly faster.
!   if (eignval(1) .lt. eignval(ntrunc)) then
!      start = 1
!      fin = ntrunc
!      step = 1
!   else
!      start = ntrunc
!      fin = 1
!      step = -1
!   end if
!
!   do i = start, fin, step 
!      dot_cv = da_dot_cv(cv_size, nuhat(:,i), shat, grid, mz, jp_start, jp_end)
!      xhat = xhat + eignval(i) * nuhat(:,i) * dot_cv
!
!   end do

   if (trace_use) call da_trace_exit ("da_amat_mul_trunc")

end subroutine da_amat_mul_trunc
