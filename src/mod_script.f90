!< author: Arthur Francisco
!<  version: 1.0.0
!<  date: december, 31 2024
!<
!<  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<        **Trident batch application**
!<  </span>
module script
!$ use omp_lib
use data_arch,       only : I4, R8
use miscellaneous,   only : get_unit
use surfile,         only : scale_surf, read_surf, write_surf, unit2IUf
use filter,          only : median_filter, fft_filter
use intpl,           only : tborne
use abbott,          only : abbott_param
use morpho,          only : topology, calcul_normales, surf_area
use grad_curv,       only : peaks_and_pits_curvatures
use stat_mom,        only : moment_stat, calc_moments, calc_median
use asfc,            only : calcul_asfc_hermite, indice_fractal
use anisotropy,      only : correlation_parameters, multiple_anisotropy, PAD_FFT_ANI, APO_FFT_ANI
use fftw3,           only : init_fftw3, fftw_plan_with_nthreads, tab_init_fftw3, end_fftw3, tab_end_fftw3, NB_THREADS_FFT, PAD_FFT, FFTW_MEASURE
use files,           only : make_path, path2vec, vec2path, filename, dir_separator, mkdir, dirname
use tchebychev,      only : least_squares_tcheby

implicit none

private

integer(kind=I4) :: JOB, SPY, STA
integer(kind=I4) :: LINE_READ

integer(kind=I4) :: NB_THREADS

integer(kind=I4) :: NB_ITER_DEB, NB_ITER_FIN, SAVE_LINE_READ

logical(kind=I4) :: GLOBAL_OMP
logical(kind=I4) :: WITH_SAMPLING

logical(kind=I4) :: OUT_CMD

character(len=512) :: NOM_SUR, SAV_NOM_SUR

character(len=  1) :: SEP

public :: read_job

contains

   subroutine read_job(irep, time, job_file, spy_unit)
   !================================================================================================
   !< @note
   !<
   !< The principle of the present function is to read a user-made batch file. The batch file, or
   !< job file, contains a sequence of keywords.
   !<
   !< Each keyword is an instruction - a macro - for a specific action:
   !<
   !< + 'READ_SUR', followed by the path a surface file, makes the program store the surface in
   !<               an array
   !< + 'LSSQ_IMG', followed by two integers, makes the program subtract a 2D polynomial, which
   !<               degrees are the provides integers.
   !< + 'SMOOTH__' makes the program smooth the surface.
   !< + etc.
   !<
   !< Some instruction keywords are present in the subroutine but not active; they will be implemented
   !< in the future.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   integer  (kind=I4), intent(in) :: irep       !! *repetition number when the program is ran multiple times*
   integer  (kind=I4), intent(in) :: spy_unit   !! *spy file unit*
   character(len=014), intent(in) :: time       !! *execution time*
   character(len=512), intent(in) :: job_file   !! *job file with macros to execute*

      real(kind=R8), allocatable, dimension(:,:) :: surf             ! surface heights array
      type(tborne) , allocatable, dimension(:)   :: sample_bounds    ! array bounds

      type(scale_surf) :: surf_prop                                  ! surface properties -- Digital Surf format
      type(scale_surf) :: samp_prop                                  ! sub surface properties

      integer(kind=I4) :: vide                                       ! file read status

      integer(kind=I4) :: nb_samples                                 ! number of samples

      character(len=512) :: mot_clef                                 ! keyword read in job file

      SEP = dir_separator()
      SPY = spy_unit

      call get_unit( JOB )                                           ! job file
      open(unit = JOB, file = trim(job_file), status = 'old')

      PAD_FFT_ANI = PAD_FFT

      APO_FFT_ANI = "no_apo"

      LINE_READ = 0
o:    do

         mot_clef = repeat( ' ', len(mot_clef) )

         read(JOB, *, iostat = vide) mot_clef ; LINE_READ = LINE_READ + 1

         write(SPY,*) LINE_READ, trim(mot_clef)

         selectcase( mot_clef(1:8) )

            case('STRT_JOB')

               call strt_job()

            case('ANALYSES')

               call analyses( tab        = surf,            &  ! IN
                              scal       = surf_prop,       &  ! IN
                              scal_samp  = samp_prop,       &  ! IN
                              tab_bounds = sample_bounds )     ! IN

            case('CORR_ECH')


            case('CROP_IMG')

               call crop_img( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! INOUT

            case('FILTER__')


            case('FIND_DEG')


            case('FT_GAUSS')

               call ft_gauss( tab  = surf,         &  ! OUT
                              scal = surf_prop )      ! OUT

            case('GET_NAME')

               NOM_SUR = SAV_NOM_SUR
               call chdir("../")

            case('HISTORY_')


            case('LOW_PASS')


            case('LSSQ_IMG')

               call lssq_img( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! INOUT

            case('NB_PROCS')

               call nb_procs()

            case('READ_SUR')

               call read_sur( tab  = surf,         &  ! OUT
                              scal = surf_prop )      ! OUT

            case('RESTRICT')


            case('SAMPLING')

               call sampling( scal       = surf_prop,          &  ! IN
                              scal_samp  = samp_prop,          &  ! OUT
                              nb_samp    = nb_samples,         &  ! OUT
                              tab_bounds = sample_bounds )        ! OUT

            case('SAV_NAME')

               SAV_NOM_SUR = NOM_SUR

            case('SAVE_SUR')

               call save_sur( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! INOUT

            case('SMOOTH__')

               call smooth__( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! IN

            case('STA_LOOP')

               call sta_loop()

            case('VERBOSES')

               OUT_CMD = .TRUE.

            case('END_LOOP')

               call end_loop( tab_bounds = sample_bounds )        ! INOUT

            case('END__JOB')

               close(JOB)

               if (OUT_CMD) write(*,*) 'Program completed'

               exit o

         endselect

      enddo o

   return
   endsubroutine read_job


   subroutine strt_job()
   !================================================================================================
   !! Some initializations: 'GLOBAL_OMP' (tasks parallelization), 'OUT_CMD' (verbose mode), etc.
   !------------------------------------------------------------------------------------------------
   implicit none

      GLOBAL_OMP = .false.

      OUT_CMD = .false.

      call random_init(.true., .true.)

   return
   endsubroutine strt_job


   subroutine nb_procs()
   !================================================================================================
   !! Define the number of concurrent threads (1 or all available threads)
   !------------------------------------------------------------------------------------------------
   implicit none

      integer(kind=I4) :: nb_th

      read(JOB,*) nb_th ; LINE_READ = LINE_READ + 1 ; write(SPY, *) LINE_READ, 'nb_procs', nb_th

      select case( nb_th )

         case( -1) ! no multihreading
            GLOBAL_OMP     = .false.
            NB_THREADS     = 1
            NB_THREADS_FFT = 1

         case(  0) ! determined by system
            GLOBAL_OMP     = .true.
            NB_THREADS     = omp_get_num_procs()
            NB_THREADS_FFT = NB_THREADS

         case default
            stop 'Bad choice "nb_procs" in "mod_script"'

      endselect

   return
   endsubroutine nb_procs


   subroutine read_sur(tab, scal)
   implicit none
   !================================================================================================
   !! Read a surface file
   !------------------------------------------------------------------------------------------------
   real(kind=R8),    intent(out), allocatable, dimension(:,:) :: tab    !! *array for the surface file*
   type(scale_surf), intent(out)                              :: scal   !! *[[scale_surf]] object*

      NOM_SUR = repeat(' ', len(NOM_SUR) )

      read(JOB,*) NOM_SUR ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, 'nom_sur ', trim(NOM_SUR)

      write(SPY, *) 'New working directory, after surface read: ', dirname(NOM_SUR)

      call read_surf( nom_fic = trim(NOM_SUR),  &  !  in; Digital Surf format
                        tab_s = tab,            &  ! out; array containing the surface
                         scal = scal )             ! out; surface type containing some informations like length, width, etc.

      call chdir( dirname(NOM_SUR) )

   return
   endsubroutine read_sur


   subroutine smooth__(tab, scal)
   !================================================================================================
   !! Smooth out a surface according [[median_filter]] procedure
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(in   )                 :: scal  !! *[[scale_surf]] object*
   real(kind=R8),    intent(inout), dimension(:,:) :: tab   !! *array containing the surface*

      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      call median_filter( tab    = tab(1:long, 1:larg),  &  ! INOUT
                          long   = long,                 &  ! IN
                          larg   = larg,                 &  ! IN
                          snb    = 10,                   &  ! IN
                          kernel = 5,                    &  ! IN
                          sig    = 3._R8,                &  ! IN
                          omp    = GLOBAL_OMP )             ! IN

   return
   endsubroutine smooth__


   subroutine crop_img(tab, scal)
   !================================================================================================
   !! Crop the surface contained in *tab* to the dimensions specified in the job file and update
   !! the object *scal*
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(inout)                              :: scal  !! *[[scale_surf]] object*
   real(kind=R8),    intent(inout), allocatable, dimension(:,:) :: tab   !! *array containing the surface*

      integer(kind=I4) :: long, larg, new_x, new_y

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      read(JOB,*) new_x, new_y ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, new_x, new_y

      long = scal%xres
      larg = scal%yres

      if (long < new_x .or. larg < new_y) stop 'Error crop surface in crop_img'

      allocate( tab_tmp(1:new_x, 1:new_y) )

      tab_tmp(1:new_x, 1:new_y) = tab(1:new_x, 1:new_y)

      deallocate( tab ) ; allocate( tab(1:new_x, 1:new_y) )

      tab(1:new_x, 1:new_y) = tab_tmp(1:new_x, 1:new_y)

      deallocate( tab_tmp )

      scal%xres = new_x
      scal%yres = new_y

      scal%nofpoints = new_x * new_y

      call write_surf(nom_fic = trim(NOM_SUR), tab_s = tab(1:new_x, 1:new_y), scal = scal)

   return
   endsubroutine crop_img


   subroutine save_sur(tab, scal)
   !================================================================================================
   !! Save the surface in the file specified in the job file
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(inout)                 :: scal  !! *[[scale_surf]] object*
   real(kind=R8),    intent(inout), dimension(:,:) :: tab   !! *array containing the surface*

      character(len=512) :: f_name

      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      f_name = repeat(' ', len(f_name) )

      read(JOB,*) f_name ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, 'NOM_SUR ', trim(f_name)

      call write_surf(nom_fic = trim(f_name), tab_s = tab(1:long, 1:larg), scal = scal)

   return
   endsubroutine save_sur


   subroutine sampling(scal, scal_samp, nb_samp, tab_bounds)
   !================================================================================================
   !< @note
   !<
   !< The function *sampling* allow for the user to sample the surface.
   !<
   !< The user is asked for the number of samples as well as the sample dimensions. It should be
   !< noted that:
   !<
   !< + the number of samples *nb_samp* must be a square,
   !< + the samples can overlap.
   !<
   !< The output is just an array that contains the samples bound.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(in )                            :: scal        !! *[[scale_surf]] object of the input surface*
   type(scale_surf), intent(out)                            :: scal_samp   !! *[[scale_surf]] object of a sample*
   integer(kind=I4), intent(out)                            :: nb_samp     !! *number of samples*
   type(tborne),     intent(out), allocatable, dimension(:) :: tab_bounds  !! *array of the sample bounds*

      integer(kind=I4) :: nn_samp, pp_samp
      integer(kind=I4) :: snb, lb1, ub1, lb2, ub2
      integer(kind=I4) :: i, j, k
      integer(kind=I4) :: long, larg

      real(kind=R8) :: wd, ht

      read(JOB,*) nn_samp, pp_samp, nb_samp ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, 'nn_samp ', nn_samp,  &  !
                                                                                                  'pp_samp ', pp_samp,  &  !
                                                                                                  'nb_samp ', nb_samp
      scal_samp = scal

      long = scal%xres
      larg = scal%yres

      if ( nb_samp <= 1 ) then

         WITH_SAMPLING = .false.

      else

         snb = nint( sqrt( 1._R8 * nb_samp ) )
         if ( nb_samp /= snb**2 ) stop 'Bad sampling number in "tab_boot"'

         WITH_SAMPLING = .true.

         scal_samp%xres = nn_samp
         scal_samp%yres = pp_samp

         scal_samp%lx = ( scal_samp%lx / long ) * nn_samp
         scal_samp%ly = ( scal_samp%ly / larg ) * pp_samp

         allocate( tab_bounds(1:nb_samp) )

         wd = real( (long - nn_samp), kind=R8 ) / (snb - 1)
         ht = real( (larg - pp_samp), kind=R8 ) / (snb - 1)

         k = 0
         do i = 1, snb
         do j = 1, snb

            k = k + 1
            lb1 = int( (i - 1) * wd + 1. )
            lb2 = int( (j - 1) * ht + 1. )

            ub1 = lb1 + nn_samp - 1
            ub2 = lb2 + pp_samp - 1

            if ( ub1 > long .or. ub2 > larg ) stop 'Error in sample bounds, in "tab_boot"'

            tab_bounds(k) = tborne( lb1, ub1, lb2, ub2 )

         enddo
         enddo

      endif

   return
   endsubroutine sampling


   subroutine sta_loop()
   !================================================================================================
   !! Start sample loop
   !------------------------------------------------------------------------------------------------
   implicit none

      read(JOB,*) NB_ITER_DEB, NB_ITER_FIN ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, 'nb_iter_deb ', NB_ITER_DEB, &  !
                                                                                                 'nb_iter_fin ', NB_ITER_FIN     !
      SAVE_LINE_READ = LINE_READ - 1 ! save the line where STA_LOOP is found

   return
   endsubroutine sta_loop


   subroutine analyses(tab, scal, scal_samp, tab_bounds)
   !================================================================================================
   !< @note
   !<
   !< The function *analyses* reads, in the Job file, the analysis to perform. The analysis is then
   !< performed.
   !<
   !< The analysis is always performed on the whole surface, and on the samples if it is asked. The
   !< results are written in the file which unit is *STA*.
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(in)                 :: scal           !! *[[scale_surf]] object of the input surface*
   type(scale_surf), intent(in)                 :: scal_samp      !! *[[scale_surf]] object of a sample*
   real(kind=R8),    intent(in), dimension(:,:) :: tab            !! *array of the surface*
   type(tborne),     intent(in), dimension(:)   :: tab_bounds     !! *array of the sample bounds*

      integer(kind=I4) :: lb1, ub1, lb2, ub2
      integer(kind=I4) :: k
      integer(kind=I4) :: lf
      integer(kind=I4) :: nn, pp
      integer(kind=I4) :: ns, ps
      integer(kind=I4) :: ibatch
      integer(kind=I4) :: exit_status

      character(len=512) :: cwd
      character(len=512) :: ana_file, sf, surf_filename, csv_filename, header
      character(len=128) :: ana_type
      character(len=006) :: degxy

      type(tborne) :: bound

      real(kind=R8), allocatable, dimension(:,:) :: tab_samp

      read(JOB,*) ana_type ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, 'ana_type ', trim(ana_type)

      nn = scal%xres
      pp = scal%yres

      call get_unit(STA) ! STA is the file unit for the analyses output

      call getcwd( cwd ) ! what is the working directory?

      ! surface filename without the folders, with .csv instead of .sur
      ! ex: "surface_02_01.sur" -> csv_filename = "surface_02_01.csv"
      surf_filename = filename( NOM_SUR )
      lf = len_trim(surf_filename)
      csv_filename = surf_filename(1:lf - 4)//".csv"

      ! surface filename without the folders and the polynomial degrees and degrees
      ! ex: "surface_02_01.sur" -> degxy = "02_01"
      !                               sf = "surface.sur"
      sf = filename( NOM_SUR )
      degxy = sf(lf - 9:lf - 3)

      sf = sf(1:lf - 10)//".sur"

      ! create result folders/file
      ! ex: resu_topology/surface_02_01.csv
      call make_ana_file(ana_f = ana_file, anal = ana_type(1:8))

      ! "resu_topology/surface_02_01.csv" -> "resu_topoglob/surface_02_01.csv"
      ! contains whole surface results
      ana_file(10:13) = 'glob'

      call fftw_plan_with_nthreads( nthreads = NB_THREADS_FFT )

      call init_fftw3( long = 2 * ( nint(PAD_FFT * nn)/2 ),    &  !
                       larg = 2 * ( nint(PAD_FFT * pp)/2 ) )      ! because of 0 padding

      ! make csv header
      call make_header(head = header, anal = ana_type(1:8))

      ! absolute path of the result file
      call make_path(wkd = trim(cwd), file_path = trim(ana_file), exit_status = exit_status)

      !=======================================================
      ! always perform analyses on the whole surface

      open( unit = STA, file = trim(ana_file) )

         write(STA,'(a)') trim(header)

         call analyses_stat( tab      = tab(1:nn, 1:pp),   &  ! IN
                             sub_samp = .false.,           &  ! IN
                             scal     = scal,              &  ! IN
                             anal     = ana_type,          &  ! IN
                             omp      = GLOBAL_OMP )          ! IN

         if (OUT_CMD) write(*,*) trim(sf), ' ', trim(degxy(2:)), ' ', ana_type(1:8), ' ', 'glob'

      close( STA )

      call end_fftw3()

      !=======================================================
      ! if SAMPLING is asked, perform analyses

      if ( WITH_SAMPLING .and. ana_type(1:8) /= 'abbott__') then

         bound = tab_bounds( 1 )                   ! retrieve 1st surface sample bounds

         lb1 = bound%lb1 ; lb2 = bound%lb2         ! lower bounds
         ub1 = bound%ub1 ; ub2 = bound%ub2         ! upper bounds

         ns = ub1 - lb1 + 1 ; ps = ub2 - lb2 + 1   ! sample size (the same for all samples)

         allocate( tab_samp(1:ns, 1:ps) )          ! surface sample array (used for all samples)

         ! batch size for parallel computing
         ibatch = max( ( NB_ITER_FIN - NB_ITER_DEB + 1 ) / NB_THREADS, 1 )

         ! file unit for analysis results
         call get_unit(STA)

         ! working directory
         call getcwd( cwd )

         ! result file
         call make_ana_file(ana_f = ana_file, anal = ana_type(1:8))

         ! one thread per fft calculus
         call fftw_plan_with_nthreads( nthreads = 1 )

         ! to prevent fft folding, use 0 padding
         call tab_init_fftw3(      long = 2 * ( nint(PAD_FFT * ns)/2 ),    &  !
                                   larg = 2 * ( nint(PAD_FFT * ps)/2 ),    &  ! because of 0 padding
                              plan_flag = FFTW_MEASURE )                      !

         ! absolute path to result file
         call make_path(wkd = trim(cwd), file_path = trim(ana_file), exit_status = exit_status)

         open( unit = STA, file = trim(ana_file), share = 'DENYRW' )

         write(STA,'(a)') trim(header)

         !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(NB_THREADS) IF (GLOBAL_OMP)
         !$OMP DO SCHEDULE (STATIC, ibatch) PRIVATE(k, bound, tab_samp, lb1, lb2, ub1, ub2)

         do k = NB_ITER_DEB, NB_ITER_FIN

            bound = tab_bounds(k)

            lb1 = bound%lb1 ; lb2 = bound%lb2
            ub1 = bound%ub1 ; ub2 = bound%ub2

            tab_samp( 1:ns, 1:ps ) = tab( lb1:ub1, lb2:ub2 )

            call analyses_stat( tab      = tab_samp(1:ns, 1:ps),    &  ! IN
                                sub_samp = .true.,                  &  ! IN
                                scal     = scal_samp,               &  ! IN
                                anal     = ana_type(1:8),           &  ! IN
                                omp      = .false. )                   ! IN

         enddo
         !$OMP END DO
         !$OMP END PARALLEL

         close( STA )

         call tab_end_fftw3()

         deallocate( tab_samp )

         if (OUT_CMD) write(*,*) trim(sf), ' ', trim(degxy(2:)), ' ', ana_type(1:8), ' ', 'samp'

      endif

   return

   contains

      subroutine make_ana_file(ana_f, anal)
      !-------------------------------------------------------------------------
      !! create the path of the analysis result file
      !-------------------------------------------------------------------------

      implicit none
      character(len=*), intent(out) :: ana_f    !! *path of the analysis result file*
      character(len=8), intent(in ) :: anal     !! *analysis asked to be performed*

         ana_f = repeat( ' ', len(ana_f) )

         select case( anal )

            case('abbott__')
               ana_f = "resu_abbott__"//SEP//trim( csv_filename )

            case('complexi')
               ana_f = "resu_complexi"//SEP//trim( csv_filename )

            case('elli_acv')
               ana_f = "resu_elli_acv"//SEP//trim( csv_filename )

            case('facettes')
               ana_f = "resu_facettes"//SEP//trim( csv_filename )

            case('ind_frac')
               ana_f = "resu_ind_frac"//SEP//trim( csv_filename )

            case('statisti')
               ana_f = "resu_statisti"//SEP//trim( csv_filename )

            case('topology')
               ana_f = "resu_topology"//SEP//trim( csv_filename )

            case default
               stop 'Bad choice of analysis in subroutine analyses'

         endselect

      return
      endsubroutine make_ana_file

      subroutine make_header(head, anal)
      !-------------------------------------------------------------------------
      !! create the header for the analysis result csv file
      !-------------------------------------------------------------------------
      implicit none
      character(len=*), intent(out) :: head  !! *header as a string*
      character(len=8), intent(in ) :: anal  !! *analysis asked to be performed*

         head = repeat( ' ', len(head) )

         select case( anal )

            case('abbott__')
               head = 'surf'//',smrk1'//degxy//       &  ! (*) ISO
                              ',smrk2'//degxy//       &  ! (*) ISO
                              ',spk__'//degxy//       &  ! (*) ISO
                              ',svk__'//degxy//       &  ! (*) ISO
                              ',Sk___'//degxy//       &  ! (*) ISO
                              ',pente'//degxy//       &  !
                              ',residus'//degxy//     &  !
                              ',coeffa_tan'//degxy//  &  !
                              ',coeffb_tan'//degxy       !

            case('complexi')
               head = 'surf'//',Sasfc'//degxy//',R2adj'//degxy ! (*) ISO

            case('elli_acv')
               head = 'surf'//',Rmax_'//degxy//         &  ! ellipsis big axis
                              ',Sal__'//degxy//         &  ! ellipsis small axis           (*) ISO
                              ',Stri_'//degxy//         &  ! another anisotropy factor     (*) ISO almost 1/Str
                              ',Std__'//degxy//         &  ! main texture orientation      (*) ISO
                              ',d.sl_'//degxy//         &  ! radius of greatest slope
                              ',b.sl_'//degxy//         &  ! greatest slope
                              ',r.sl_'//degxy//         &  ! slope anisotropy factor
                              ',r.cv_'//degxy//         &  ! curvature anisotropy factor
                              ',bmp__'//degxy//         &  ! maximum over [0,179°] of the peaks mean width'
                              ',smp__'//degxy//         &  ! minimum over [0,179°] of the peaks mean width'
                              ',rmp__'//degxy//         &  ! ratio bmp/smp'
                              ',bml__'//degxy//         &  ! maximum over [0,179°] of the path length'
                              ',sml__'//degxy//         &  ! minimum over [0,179°] of the path length'
                              ',rml__'//degxy//         &  ! ratio bml/sml'
                              ',bms__'//degxy//         &  ! maximum over [0,179°] of the standard deviation of slope'
                              ',sms__'//degxy//         &  ! minimum over [0,179°] of the standard deviation of slope'
                              ',rms__'//degxy              ! ratio bms/sms'

            case('facettes')
               head = 'surf'//',Sh___'//degxy//',Sdr__'//degxy ! %faces in a 5° cone, Sdr: (*) ISO

            case('ind_frac')
               head = 'surf'//',Smbd_'//degxy//',ord_orig'//degxy//',R2adj'//degxy ! Smbd: Minkowski–Bouligand dimension

            case('statisti')
               head = 'surf'//',Sv___'//degxy//         &  ! (*) ISO
                              ',Sp___'//degxy//         &  ! (*) ISO
                              ',Smd__'//degxy//         &  ! median
                              ',Sa___'//degxy//         &  ! (*) ISO
                              ',Sm___'//degxy//         &  ! mean
                              ',Sq___'//degxy//         &  ! (*) ISO
                              ',Ssk__'//degxy//         &  ! (*) ISO
                              ',Sku__'//degxy//         &  ! (*) ISO
                              ',Sks__'//degxy              ! sku/(ssk**2 + 1)

            case('topology')
               head = 'surf'//',Snb1_'//degxy//         &  ! nb_cells_1
                              ',Smc1_'//degxy//         &  ! median_size_1
                              ',Sk1__'//degxy//         &  ! fraction_of_surface_1
                              ',Snb2_'//degxy//         &  ! nb_cells_2
                              ',Smc2_'//degxy//         &  ! median_size_2
                              ',Sk2__'//degxy//         &  ! fraction_of_surface_2
                              ',Sdq__'//degxy//         &  ! gradient quadratic mean
                              ',Scq__'//degxy//         &  ! curvature quadratic mean
                              ',Sh3z_'//degxy//         &  ! 3 highest curvature mean
                              ',Sv3z_'//degxy              ! 3 deepest curvature mean

            case default
               stop 'Bad choice of analysis in subroutine analyses'

         endselect

      return
      endsubroutine make_header

      subroutine analyses_stat(tab, sub_samp, scal, anal, omp)
      !-------------------------------------------------------------------------
      !! Perform analyses on the surface *tab*. The results are written in the file
      !! of unit *STA*.
      !-------------------------------------------------------------------------
      implicit none
      type(scale_surf), intent(in )                 :: scal       !! *[[scale_surf]] object*
      character(len=8), intent(in )                 :: anal       !! *analysis asked to be performed*
      logical(kind=I4), intent(in )                 :: omp        !! *parallel computing?*
      logical(kind=I4), intent(in )                 :: sub_samp   !! *surface sampling?*
      real(kind=R8),    intent(in ), dimension(:,:) :: tab        !! *surface height array*

         integer (kind=I4) :: nx, ny

         real(kind=R8)     :: ra_t, md
         real(kind=R8)     :: dx, dy, si, fft_cutoff

         type(MOMENT_STAT) :: mx

         real(kind=R8), dimension(:,:), allocatable :: tab_tmp

         real(kind=R8), dimension(:), allocatable :: vec_heights
         real(kind=R8), dimension(1:20)           :: ana_res

         ana_res = 0

         nx = scal%xres
         ny = scal%yres

         dx = scal%dx * unit2IUf( scal%dx_unit )   !
         dy = scal%dy * unit2IUf( scal%dy_unit )   !
         si = 1                                    ! scal%si

         select case( anal )

            case( 'abbott__' )

               allocate( vec_heights(1:nx*ny) )

               vec_heights(1:nx*ny) = reshape( tab(1:nx, 1:ny), [ nx*ny ] )

               call abbott_param( tab     = vec_heights(1:nx * ny),           &  !
                                  lg      = nx * ny,                          &  !
                                  nom     = "out"//SEP//"abbott_res",         &  !
                                  curves  = [.false., .false., .false.],      &  !
                                  results = ana_res(1:11),                    &  !
                                  omp     = omp )                                !

               !$omp critical
               write(STA,'(9(a,E18.6))')     trim(sf)//',', &  !
                                               ana_res( 1), &  ! smrk1, iso 25178
                                          ',', ana_res( 2), &  ! smrk2, iso 25178
                                          ',', ana_res( 3), &  ! spk  , iso 25178
                                          ',', ana_res( 4), &  ! svk  , iso 25178
                                                               ! 5 et 6 pour off1 et off2
                                          ',', ana_res( 7), &  ! sk   , iso 25178
                                          ',', ana_res( 8), &  ! core slope
                                          ',', ana_res( 9), &  ! adjustment factor (tangent fit)
                                          ',', ana_res(10), &  ! coeffa_tan        (tangent fit)
                                          ',', ana_res(11)     ! coeffb_tan        (tangent fit)
               !$omp end critical

               deallocate( vec_heights )

            case( 'topology' )

               call topology( tab  = tab(1:nx, 1:ny),    &  !
                              long = nx,                 &  !
                              larg = ny,                 &  !
                              res  = ana_res(1:6) )         !

               fft_cutoff = dx / 5.e-6 ! 5.e-6 = 5 µm

               allocate( tab_tmp(1:nx, 1:ny) )

               call fft_filter(tab       = tab(1:nx, 1:ny),      & ! in
                               long      = nx,                   & ! in
                               larg      = ny,                   & ! in
                               cutoff    = fft_cutoff,           & ! in
                               bf_tab    = tab_tmp(1:nx, 1:ny),  & ! out
                               multi_fft = sub_samp)               ! in

               call peaks_and_pits_curvatures( heights      = tab_tmp(1:nx, 1:ny),  &  !
                                               nx           = nx,                   &  !
                                               ny           = ny,                   &  !
                                               dx           = dx,                   &  !
                                               dy           = dy,                   &  !
                                               S_param_grad = ana_res(07),          &  !
                                               S_param_curv = ana_res(08),          &  !
                                               peak_curv    = ana_res(09),          &  !
                                               pits_curv    = ana_res(10) )            !

               deallocate( tab_tmp )

               !$omp critical
               write(STA,'(10(a,E18.6))')  trim(sf)//',',      &  !
                                              ana_res(01),     &  ! number of cells 0.15/0.85
                                         ',', ana_res(02),     &  ! median size (percentage of surface)
                                         ',', ana_res(03),     &  ! percentage of surface above 0.15/0.85
                                         ',', ana_res(04),     &  ! number of cells 0.15/0.95
                                         ',', ana_res(05),     &  ! median size (percentage of surface)
                                         ',', ana_res(06),     &  ! percentage of surface above 0.15/0.95
                                         ',', ana_res(07),     &  ! Sdq__ gradient  quadratic mean
                                         ',', ana_res(08),     &  ! Scq__ curvature quadratic mean
                                         ',', ana_res(09),     &  ! Sh3z_ 3 highest curvature mean
                                         ',', ana_res(10)         ! Sv3z_ 3 deepest curvature mean
               !$omp end critical

            case( 'statisti' )

               call calc_moments( tab    = reshape( tab(1:nx, 1:ny), [nx * ny] ),   &  !
                                  mx     = mx,                                      &  !
                                  nb_mom = 4 )                                         !

               call calc_median( tab = reshape( tab(1:nx, 1:ny), [nx * ny] ),       &  !
                                 md  = md )                                            !

               ra_t = sum( abs(tab(1:nx, 1:ny) - mx%mu ) / (nx * ny) )

               ana_res(1:8) = [ minval( tab(1:nx, 1:ny) ) - mx%mu,   &  !
                                maxval( tab(1:nx, 1:ny) ) - mx%mu,   &  !
                                                       md - mx%mu,   &  !
                              ra_t, mx%mu, mx%si, mx%sk, mx%ku ]        !

               !$omp critical
               write(STA,'(9(a,E18.6))')  trim(sf)//',',                         &  !
                                               ana_res(1),                       &  ! max depth
                                          ',', ana_res(2),                       &  ! max height
                                          ',', ana_res(3),                       &  ! median
                                          ',', ana_res(4),                       &  ! absolute roughness
                                          ',', ana_res(5),                       &  ! mean height
                                          ',', ana_res(6),                       &  ! standard deviation
                                          ',', ana_res(7),                       &  ! skewness
                                          ',', ana_res(8),                       &  ! kurtosis
                                          ',', ana_res(8)/( ana_res(7)**2 + 1 )     ! kind of kurtosis excess

               !$omp end critical

            case( 'ind_frac' )

               call indice_fractal( tab_in = tab(1:nx, 1:ny),  &  !
                                    long   = nx,               &  !
                                    larg   = ny,               &  !
                                    indf   = ana_res(1:3) )       !

               !$omp critical
               write(STA, '(3(a,E18.6))')  trim(sf)//',',     &  !
                                                ana_res(1),   &  ! fractal dimension (slope)
                                           ',', ana_res(2),   &  ! ordinate at origin
                                           ',', ana_res(3)       ! adjustment factor
               !$omp end critical

            case( 'facettes' )

               call calcul_normales( tab_in     = tab(1:nx, 1:ny),   &  !
                                     long       = nx,                &  !
                                     larg       = ny,                &  !
                                     scale_xyz  = [ dx, dy, si ],    &  !
                                     cone_angle = 5._R8,             &  !
                                     hori       = ana_res(1) )          !

               call surf_area( tab_in     = tab(1:nx, 1:ny),         &  !
                               long       = nx,                      &  !
                               larg       = ny,                      &  !
                               scale_xyz  = [ dx, dy, si ],          &  !
                               aire       = ana_res(2) )                !

               !$omp critical
               write(STA, '(2(a,E18.6))') trim(sf)//',',     &  !
                                               ana_res(1),   &  ! horizontality
                                          ',', ana_res(2)       ! relative area
               !$omp end critical

            case( 'complexi' )

               call calcul_asfc_hermite( tab_in   = tab(1:nx, 1:ny),    &  !
                                         scal     = scal,               &  !
                                         asfc_res = ana_res(1:2),       &  !
                                         omp      = omp )                  !

               !$omp critical
               write(STA, '(2(a,E18.6))') trim(sf)//',',     &  !
                                               ana_res(1),   &  ! asfc
                                          ',', ana_res(2)       ! adjustment factor
               !$omp end critical

            case( 'elli_acv' )

               call correlation_parameters( tab       = tab(1:nx, 1:ny),   &  ! IN
                                            long      = nx,                &  ! IN
                                            larg      = ny,                &  ! IN
                                            res       = ana_res(1:8),      &  ! OUT
                                            cut       = 0.2_R8,            &  ! IN 0.5_R8,            &  ! IN
                                            sub_plane = .true.,            &  ! IN, remove plane?
                                            sub_sampl = sub_samp,          &  ! IN, subsampling?
                                            scale_xy  = [ dx, dy ],        &  ! IN
                                            omp       = omp )                 ! IN

               ana_res(9:17) = 0
               call multiple_anisotropy( tabin     = tab(1:nx, 1:ny),       &  ! IN
                                         long      = nx,                    &  ! IN
                                         larg      = ny,                    &  ! IN
                                         scale_xy  = [ dx, dy ],            &  ! IN
                                         multi_fft = sub_samp,              &  ! IN
                                         vec_ani   = ana_res(9:17) )           ! OUT

               !$omp critical
               write(STA,'(17(a,E18.6))') trim(sf)//',',       &  !
                                              ana_res(01),     &  ! big axis
                                         ',', ana_res(02),     &  ! small axis
                                         ',', ana_res(03),     &  ! big/small -> anisotropy 1
                                         ',', ana_res(04),     &  ! anisotropy angle
                                         ',', ana_res(05),     &  ! radius of greatest slope
                                         ',', ana_res(06),     &  ! greatest slope
                                         ',', ana_res(07),     &  ! slope anisotropy factor
                                         ',', ana_res(08),     &  ! curvature anisotropy factor
                                         ',', ana_res(09),     &  ! maximum over [0,179°] of the peaks mean width'
                                         ',', ana_res(10),     &  ! minimum over [0,179°] of the peaks mean width'
                                         ',', ana_res(11),     &  ! ratio bmp/smp'
                                         ',', ana_res(12),     &  ! maximum over [0,179°] of the path length'
                                         ',', ana_res(13),     &  ! minimum over [0,179°] of the path length'
                                         ',', ana_res(14),     &  ! ratio bml/sml'
                                         ',', ana_res(15),     &  ! maximum over [0,179°] of the standard deviation of slope'
                                         ',', ana_res(16),     &  ! minimum over [0,179°] of the standard deviation of slope'
                                         ',', ana_res(17)         ! ratio bms/sms'
               !$omp end critical

            case default

               stop 'Bad choice of analysis in subroutine analyses_stat'

         endselect

      return
      endsubroutine analyses_stat

   endsubroutine analyses


   subroutine end_loop(tab_bounds)
   !================================================================================================
   !! End of the Job file - some finalizations.
   !------------------------------------------------------------------------------------------------
   implicit none
   type(tborne),  intent(inout), allocatable, dimension(:)   :: tab_bounds !! *array to deallocate*

      if ( allocated( tab_bounds ) )  deallocate( tab_bounds )

      WITH_SAMPLING = .FALSE.

   return
   endsubroutine end_loop


   subroutine ft_gauss(tab, scal)
   !================================================================================================
   !< @note
   !<
   !< The function *ft_gauss* performs a Gaussian filter on the surface *tab*. It:
   !<
   !< + reads the *fft_cutoff* in the Job file
   !< + applies the corresponding Gaussian filter
   !< + creates a new folder
   !< + stores the new surface in the new folder.
   !<
   !< Example:
   !<
   !< if cutoff = 30.5 µm and the working directory `sur`, the resulting surface is
   !< `sur/FLT_030.5/file_FLT_030.5.sur`
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(inout)                 :: scal  !! *[[scale_surf]] object of the input surface*
   real(kind=R8),    intent(inout), dimension(:,:) :: tab   !! *array of the surface*

      integer(kind=I4) :: nx, ny
      integer(kind=I4) :: degx, degy
      integer(kind=I4) :: i, istat

      real(kind=R8) :: dx, fft_cutoff, cutoff

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      character(len=512) :: surf_filename, new_surf_filename
      character(len=512) :: wkd
      character(len=128) :: str

      read(JOB,*) fft_cutoff ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, fft_cutoff

      nx = scal%xres
      ny = scal%yres

      dx = scal%dx * unit2IUf( scal%dx_unit )

      allocate( tab_tmp(1:nx, 1:ny) )

      cutoff = dx / fft_cutoff

      call fftw_plan_with_nthreads( nthreads = NB_THREADS_FFT )

      call init_fftw3( long = 2 * ( nint(PAD_FFT * nx)/2 ),    &  !
                       larg = 2 * ( nint(PAD_FFT * ny)/2 ) )      ! because of 0 padding

      call fft_filter(tab       = tab(1:nx, 1:ny),      & ! in
                      long      = nx,                   & ! in
                      larg      = ny,                   & ! in
                      cutoff    = cutoff,               & ! in
                      bf_tab    = tab_tmp(1:nx, 1:ny),  & ! out
                      multi_fft = .false.)                ! in

      call end_fftw3()

      ! surface filename without the folders
      surf_filename = filename( NOM_SUR )

      ! folder name where the new surface is stored. Ex: FLT_030.5
      write( str,'(i3.3,f0.1)') int(fft_cutoff*1e6), fft_cutoff*1e6 - int(fft_cutoff*1e6) ! writes 030, then .5
      str = 'FLT_'//adjustl( trim(str) )

      ! the folder must be created under the folder where the surface is. Recall that the currect directory
      ! is set in subroutine read_sur
      call getcwd( wkd )
      call mkdir(wkd = trim(wkd), directory = trim(str), sep = SEP, exit_status = istat)

      ! set the new working directory
      call chdir( trim(wkd)//SEP//trim(str) )
      write(SPY, *) 'New working directory, after FT_GAUSS: ', trim(wkd)//SEP//trim(str)

      ! the name of the new surface
      NOM_SUR = repeat( ' ', len(NOM_SUR) )
      NOM_SUR = surf_filename(1:len_trim(surf_filename) - 4)//'_'//trim(str)//'.sur'

      ! the new surface is stored in, as an example: sur/FLT_030.5/file_FLT_030.5.sur
      new_surf_filename = repeat( ' ', len(new_surf_filename) )
      new_surf_filename = trim(str)//SEP//trim(NOM_SUR)

      call write_surf(nom_fic = trim(NOM_SUR), tab_s = tab_tmp(1:nx, 1:ny), scal = scal)

      deallocate( tab_tmp )

   return
   endsubroutine ft_gauss


   subroutine lssq_img(tab, scal)
   !================================================================================================
   !< @note
   !<
   !< The function *lssq_img* subtracts a least square 2D polynomial from *tab*. It:
   !<
   !< + reads the degrees in the Job file
   !< + calculates and subtracts the least square 2D polynomial
   !< + creates a new folder
   !< + stores the new surface in the new folder.
   !<
   !< Example:
   !<
   !< if degrees are (5,8) and the working directory `sur`, the resulting surface is
   !< `sur/DEG_05_08/file_X1_005_Y1_008.sur`
   !<
   !< @endnote
   !------------------------------------------------------------------------------------------------
   implicit none
   type(scale_surf), intent(inout)                 :: scal  !! *[[scale_surf]] object of the input surface*
   real(kind=R8),    intent(inout), dimension(:,:) :: tab   !! *array of the surface*

      integer(kind=I4) :: long, larg
      integer(kind=I4) :: degx, degy
      integer(kind=I4) :: istat

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      character(len=512) :: surf_filename, new_surf_filename
      character(len=512) :: wkd
      character(len=128) :: str

      read(JOB,*) degx, degy ; LINE_READ = LINE_READ + 1 ; write(SPY,*) LINE_READ, degx, degy

      long = scal%xres
      larg = scal%yres

      allocate( tab_tmp(1:long, 1:larg) )

      call least_squares_tcheby( tab_in       = tab(1:long, 1:larg),       &  !
                                 tab_out      = tab_tmp(1:long, 1:larg),   &  !
                                 long1        = long,                      &  !
                                 long2        = larg,                      &  !
                                 nvarx        = degx,                      &  !
                                 nvary        = degy,                      &  !
                                 verif        = .false.,                   &  !
                                 multi_thread = GLOBAL_OMP )                  !

      ! the least square polynomial tab_tmp is ubtracted from tab
      tab(1:long, 1:larg) = tab(1:long, 1:larg) - tab_tmp(1:long, 1:larg)

      ! surface filename without the folders
      surf_filename = filename( NOM_SUR )

      ! folder name where the new surface is stored. Ex: DEG_08_06
      write( str, '(i2.2,a1,i2.2)' ) degx, '_', degy
      str = 'DEG_'//adjustl( trim(str) )

      ! the folder must be created under the folder where the surface is. Recall that the currect directory
      ! is set in subroutine read_sur
      call getcwd( wkd )
      call mkdir(wkd = trim(wkd), directory = trim(str), sep = SEP, exit_status = istat)

      ! set the new working directory
      call chdir( trim(wkd)//SEP//trim(str) )
      write(SPY, *) 'New working directory, after LSSQ: ', trim(wkd)//SEP//trim(str)

      ! the name of the new surface
      NOM_SUR = repeat( ' ', len(NOM_SUR) )
      write( NOM_SUR, '(a,2(a,i2.2),a)' ) surf_filename(1:len_trim(surf_filename) - 4), '_', degx, '_', degy, '.sur'

      ! the new surface is stored in, as an example: sur/DEG_05_08/file_X1_005_Y1_008.sur
      new_surf_filename = repeat( ' ', len(new_surf_filename) )
      new_surf_filename = trim(str)//SEP//trim(NOM_SUR)

      call write_surf(nom_fic = trim(NOM_SUR), tab_s = tab(1:long, 1:larg), scal = scal)

      deallocate( tab_tmp )

   return
   endsubroutine lssq_img

endmodule script

