subroutine da_vv_to_cv(vv, xp, mzs, cv_size, rcv &
#if (WRF_CHEM == 1)
                        , vch &
#endif
                        )

   !---------------------------------------------------------------------------
   ! Purpose: Fill (local) 1D cv array from 3D (local) vv arrays.
   !---------------------------------------------------------------------------

   implicit none

   type (vp_type), intent(in)    :: vv          ! Grdipt/EOF rcv.
#if (WRF_CHEM == 1)
   type (xch_type), optional, intent(inout) :: vch          ! Grdipt/EOF rcv for chem
#endif
   type (xpose_type), intent(in) :: xp          ! Dimensions and xpose buffers.
   integer,        intent(in)    :: mzs(:)      ! Background error structure levels.
   integer,        intent(in)    :: cv_size     ! Length of CV array.
   real,           intent(inout) :: rcv(1:cv_size) ! Control variables v.

   integer   :: is,ie       ! Local grid range in y coordinate.
   integer   :: js,je       ! Local grid range in x coordinate.
   integer   :: ix,jy       ! Local grid horizontal dimensions.
   integer   :: mz          ! Max vertical coordinate for v1 through v5 arrays.
   integer   :: ne          ! Ensemble size.
   integer   :: cv_s,cv_e   ! Starting and ending indices into CV array.
   integer   :: ijm         ! Size of interior of v1 through v5 arrays.
   integer   :: ijmn        ! Size of interior of alpha cv arrays.
#if (WRF_CHEM == 1)
   integer   :: ma, mb      ! Max anthro and bioburn time coordinates
   integer   :: ic          ! chem loop iterator
#endif

#ifdef NORESHAPE
   integer   :: i,j,k,ijk,m,n
#endif

   if (trace_use) call da_trace_entry("da_vv_to_cv")

   if( use_rf )then
      is = xp % its
      ie = xp % ite
      js = xp % jts
      je = xp % jte
   else
      call da_error(__FILE__,__LINE__,(/"This subroutine should not be called for use_rf = .false."/))
   endif
   ix = ie-is+1
   jy = je-js+1
   cv_e = 0
#ifdef NORESHAPE
   ijk=0
#endif

#if (WRF_CHEM == 1)
   if (cv_options .gt. 0) then
#endif
   ! Store v1
   mz = mzs(1)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v1(i,j,k)
            end do
        end do
     end do
#else
     rcv(cv_s:cv_e) = RESHAPE(vv % v1(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v2
   mz = mzs(2)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v2(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v2(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v3
   mz = mzs(3)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1

#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v3(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v3(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v4
   mz = mzs(4)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v4(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v4(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v5
   mz = mzs(5)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v5(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v5(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
#ifdef CLOUD_CV
   ! Store v6
   mz = mzs(6)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v6(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v6(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v7
   mz = mzs(7)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v7(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v7(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v8
   mz = mzs(8)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v8(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v8(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v9
   mz = mzs(9)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v9(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v9(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v10
   mz = mzs(10)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v10(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v10(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
   ! Store v11
   mz = mzs(11)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               rcv(ijk) = vv%v11(i,j,k)
            end do
         end do
      end do
#else
      rcv(cv_s:cv_e) = RESHAPE(vv % v11(is:ie,js:je,1:mz), (/ijm/))
#endif
   end if
#endif
#ifdef CLOUD_CV
   ! Store alpha:
   mz = mzs(12)
   ne = mzs(13)
#else
   ! Store alpha:
   mz = mzs(6)
   ne = mzs(7)
#endif
   if ( ne > 0 ) then
      ix = ite_int - its_int + 1
      jy = jte_int - jts_int + 1
      ijmn = ix * jy * mz * ne
      cv_s = cv_e + 1
      cv_e = cv_s + ijmn - 1
#ifdef NORESHAPE
      do n = 1, ne
         do m = 1, mz
            do j = jts_int, jte_int !js, je
               do i = its_int, ite_int ! is, ie
                  ijk = ijk + 1
                  rcv(ijk) = vv%alpha(i,j,m,n)
               end do
            end do
         end do
      end do
#else
!     rcv(cv_s:cv_e) = RESHAPE(vv % alpha(is:ie,js:je,1:mz,1:ne), (/ijmn/))
      rcv(cv_s:cv_e) = RESHAPE(vv % alpha(its_int:ite_int,jts_int:jte_int,1:mz,1:ne), (/ijmn/))
#endif
   end if
#if (WRF_CHEM == 1)
   end if !cv_options > 0

   if (cv_options_chem.gt.0 .and. present(vch)) then
      ix = ie-is+1
      jy = je-js+1
#ifdef CLOUD_CV
      ma = mzs(14)
      mb = mzs(15)
#else
      ma = mzs(8)
      mb = mzs(9)
#endif
      if ( ma > 0 ) then
         ijm = ix * jy * ma
         do ic = PARAM_FIRST_SCALAR, num_scaleant
            cv_s = cv_e + 1
            cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
            do k=1,ma
               do j=js,je
                  do i=is,ie
                     ijk = ijk + 1
                     rcv(ijk) = vch % alpha_ant(i,j,k,ic)
                  end do
               end do
            end do
#else
            rcv(cv_s:cv_e) = RESHAPE(vch % alpha_ant(is:ie,js:je,1:ma,ic))
#endif
         end do
      end if

      if ( mb > 0 ) then
         ijm = ix * jy * mb
         do ic = PARAM_FIRST_SCALAR, num_scalebb
            cv_s = cv_e + 1
            cv_e = cv_s + ijm - 1
#ifdef NORESHAPE
            do k=1,mb
               do j=js,je
                  do i=is,ie
                     ijk = ijk + 1
                     rcv(ijk) = vch % alpha_bb(i,j,k,ic)
                  end do
               end do
            end do
#else
            rcv(cv_s:cv_e) = RESHAPE(vch % alpha_bb(is:ie,js:je,1:mb,ic))
#endif
         end do
      end if
   end if

#endif

   if (trace_use) call da_trace_exit("da_vv_to_cv")

end subroutine da_vv_to_cv

