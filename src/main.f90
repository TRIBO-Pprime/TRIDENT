!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: december, 31 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Trident application**
!<  </span>
program main
!$ use omp_lib
use script,          only : read_job
use data_arch,       only : I4
use miscellaneous,   only : get_unit
use files,           only : make_path, dir_separator, mkdir

implicit none

integer(kind=I4)   :: ter, spy, istat, ind

character(len=512) :: wkd
character(len=  1) :: sep

   sep = dir_separator()

   ! find where command is run
   call get_command( wkd )
   ind = max( index(wkd, "/main"), index(wkd, "\main") )

   ! set the new working dir
   call chdir( trim(wkd(:ind - 1) ) )
   call getcwd( wkd )
   call mkdir(wkd = trim(wkd), directory = "out", sep = sep, exit_status = istat)

   call get_unit( spy )                                           ! file where read keywords of job file
   open(unit = spy, file = "out"//sep//"job_reading.txt")         ! and other stuff are copied (check purposes)
   write(spy, *) 'executable directory: ', trim(wkd)

   call get_unit( ter )
   open(unit = ter, file = "out"//sep//"execution_output.txt")

   call prg_surf( out_unit = ter, spy_unit = spy )

   close( spy )
   close( ter )

   contains

   subroutine prg_surf(out_unit, spy_unit)
   implicit none
   integer(kind=I4), intent(in) :: out_unit
   integer(kind=I4), intent(in) :: spy_unit

      character(len=128) :: arg_prg ! program arguments
      character(len=014) :: d_and_t
      character(len=008) :: chara_d
      character(len=010) :: chara_t
      character(len=512) :: job_file

      integer(kind=I4)   :: var_i, prg_repeat, i_repeat

      call date_and_time(date = chara_d, time = chara_t)
      d_and_t = chara_d//chara_t(1:6)

      ! string initialisation
      arg_prg  = repeat(' ', len(arg_prg))
      job_file = repeat(' ', len(job_file))

      var_i = 1
      !...............................................................................................
      call get_command_argument(var_i, arg_prg) ! argument one: argument string

      if (len_trim(arg_prg) == 0) then          ! if there is no job file, stop

         write(out_unit,*) 'no job file, stop'
         stop

      else

         job_file = trim(arg_prg)

      endif

      var_i = 2
      !...............................................................................................
      call get_command_argument(var_i, arg_prg) ! argument two: number of repetitions

      if (len_trim(arg_prg) == 0) then

         prg_repeat = 1                         ! if no argument, no repetition

      else

         read(arg_prg,*) prg_repeat

      endif

      do i_repeat = 1, prg_repeat
         call read_job( irep     = i_repeat,       &  !
                        time     = d_and_t,        &  !
                        job_file = job_file,       &  !
                        spy_unit = spy_unit )         ! the program executes 'prg_repeat' times
      enddo

      write(out_unit,*) 'Program completed'

   return
   endsubroutine prg_surf

endprogram main
