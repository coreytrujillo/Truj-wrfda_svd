subroutine da_cv_io (io_config, filename, cv_size, rcv, grid, mz, ascii_bool)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write yhat, omega, qhat, and cv vectors
   !
   ! Called from da_solve and da_randomise_svd
   !
   ! History: 09/01/2016  Creation (JJ Guerrette), modifed from da_lanczos_io.inc
   !
   !-------------------------------------------------------------------------

   implicit none

   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   character(*),     intent(in)    :: filename               ! I/O filename
   integer,          intent(in)    :: cv_size
   real*8,           intent(inout) :: rcv(1:cv_size)
   type(domain), optional, intent(inout) :: grid
   integer,      optional, intent(in)    :: mz(:)
   logical,      optional, intent(in)    :: ascii_bool

   integer                         :: cv_unit, dummy

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

   if (trace_use) call da_trace_entry("da_cv_io")

   if (use_global_cv_io .and. present(grid) .and. present(mz)) then
      if (present(ascii_bool)) then
         call da_cv_io_global (io_config, grid, filename, cv_size, rcv, mz, ascii_bool)
      else
         call da_cv_io_global (io_config, grid, filename, cv_size, rcv, mz)
      end if
   else 
      call da_get_unit (cv_unit)

      if (io_config == 'r') then
         write(*,*) 'Reading cv vector from '//trim(filename)
         if (present(ascii_bool)) then
            open (unit=cv_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
            read(unit=cv_unit,fmt='(I10)') dummy
            read(unit=cv_unit,fmt='(E27.18E3)') rcv(1:cv_size)
         else
            open (unit=cv_unit, file = trim(filename), form = 'unformatted', status = 'old')
            read(unit=cv_unit) dummy
            read(unit=cv_unit) rcv(1:cv_size)
         end if
         close(unit=cv_unit)
      else if (io_config == 'w') then
         write(*,*) 'Writing cv vector to '//trim(filename)
         if (present(ascii_bool)) then
            open (unit=cv_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
            write(unit=cv_unit,fmt='(I10)') cv_size
            write(unit=cv_unit,fmt='(E27.18E3)') rcv(1:cv_size)
         else
            open (unit=cv_unit, file = trim(filename), form = 'unformatted', status = 'replace')
            write(unit=cv_unit) cv_size
            write(unit=cv_unit) rcv(1:cv_size)
         end if
         close(unit=cv_unit)
      else
         write(*,*) 'Unknow configuration for CV I/O routine'
      end if

      call da_free_unit (cv_unit)
   end if

   if (trace_use) call da_trace_exit ("da_cv_io")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: IO_cv_generic')
#endif

end subroutine da_cv_io