module script
!$ use omp_lib
use data_arch,       only : I4, R8
use miscellaneous,   only : get_unit
use surfile,         only : scale_surf, read_surf, write_surf, unit2IUf
use filter,          only : median_filter
use intpl,           only : tborne
use abbott,          only : abbott_param
use morpho,          only : topology, calcul_normales, surf_area
use stat_mom,        only : moment_stat, calc_moments, calc_median
use asfc,            only : calcul_asfc_hermite, indice_fractal
use anisotropy,      only : correlation_parameters
use fftw3,           only : init_fftw3, fftw_plan_with_nthreads, tab_init_fftw3, end_fftw3, tab_end_fftw3, NB_THREADS_FFT
use files,           only : make_path, path2vec, vec2path, filename, dir_separator, mkdir, dirname
use tchebychev,      only : least_squares_tcheby

implicit none

private

integer(kind=I4) :: JOB, SPY, STA
integer(kind=I4) :: LIGNE_LUE

integer(kind=I4) :: NB_THREADS

integer(kind=I4) :: NB_ITER_DEB, NB_ITER_FIN, SAVE_LIGNE_LUE

logical(kind=I4) :: GLOBAL_OMP
logical(kind=I4) :: WITH_SAMPLING

character(len=512) :: NOM_SUR

character(len=  1) :: SEP

public :: read_job

contains

   subroutine read_job(irep, time, job_file, spy_unit)
   implicit none
   integer  (kind=I4), intent(in) :: irep       !! repetition number when the program is run multiple times
   integer  (kind=I4), intent(in) :: spy_unit   !! spy file unit
   character(len=014), intent(in) :: time       !! execution time
   character(len=512), intent(in) :: job_file   !! job file with macros to execute

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

      LIGNE_LUE = 0
o:    do

         mot_clef = repeat( ' ', len(mot_clef) )

         read(JOB, *, iostat = vide) mot_clef ; LIGNE_LUE = LIGNE_LUE + 1

         write(SPY,*) LIGNE_LUE, trim(mot_clef)

         selectcase( mot_clef(1:8) )

            case('DEBUT___')

               call debut___()

            case('CORR_ECH')



            case('ANALYSES')

               call analyses( tab        = surf,            &  ! IN
                              scal       = surf_prop,       &  ! IN
                              scal_samp  = samp_prop,       &  ! IN
                              tab_bounds = sample_bounds )     ! IN

            case('HISTORY_')



            case('FILTER__')



            case('FIND_DEG')



            case('LECT_BAS')

               call lect_bas( tab  = surf,         &  ! OUT
                              scal = surf_prop )      ! OUT

            case('LSSQ_IMG')

               call lssq_img( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! INOUT

            case('NB_PROCS')

               call nb_procs()

            case('RESTRICT')



            case('SAVE_SUR')

               call save_sur( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! INOUT

            case('LOW_PASS')



            case('SMOOTH__')

               call smooth__( tab  = surf,         &  ! INOUT
                              scal = surf_prop )      ! IN

            case('STA_LOOP')

               call sta_loop()

            case('TAB_BOOT')

               call tab_boot( scal       = surf_prop,          &  ! IN
                              scal_samp  = samp_prop,          &  ! OUT
                              nb_samp    = nb_samples,         &  ! OUT
                              tab_bounds = sample_bounds )        ! OUT

            case('VERBOSES')



            case('END_LOOP')

               call end_loop( tab        = surf,               &  ! INOUT
                              tab_bounds = sample_bounds )        ! INOUT

            case('FIN_____')

               close(JOB)

               exit o

         endselect

      enddo o

   return
   endsubroutine read_job

   subroutine debut___()
   implicit none

      GLOBAL_OMP = .false.

      call random_init(.true., .true.)

   return
   endsubroutine debut___


   subroutine nb_procs()
   implicit none

      integer(kind=I4) :: nb_th

      read(JOB,*) nb_th ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY, *) LIGNE_LUE, 'nb_procs', nb_th

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


   subroutine lect_bas(tab, scal)
   implicit none
   real(kind=R8),    intent(out), allocatable, dimension(:,:) :: tab
   type(scale_surf), intent(out)                              :: scal

!~       character(len=512), allocatable, dimension(:) :: vpath

!~       character(len = :), allocatable               :: path_surf_dir

      NOM_SUR = repeat(' ', len(NOM_SUR) )

      read(JOB,*) NOM_SUR ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, 'nom_sur ', trim(NOM_SUR)

!~       call path2vec( file_path = trim(NOM_SUR),    &  !
!~                      vec_path  = vpath )              !

!~       call vec2path( file_path = path_surf_dir,                   &  !
!~                      vec_path  = vpath( 1:size(vpath) - 1 ) )        !

      call chdir( dirname(NOM_SUR) )
      write(SPY, *) 'New working directory, after surface read: ', dirname(NOM_SUR)

      call read_surf( nom_fic = trim(NOM_SUR),  &  !  in; Digital Surf format
                        tab_s = tab,            &  ! out; array containing the surface
                         scal = scal )             ! out; surface type containing some informations like length, width, etc.

   return
   endsubroutine lect_bas


   subroutine smooth__(tab, scal)
   implicit none
   type(scale_surf), intent(in   )                 :: scal
   real(kind=R8),    intent(inout), dimension(:,:) :: tab

      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      call median_filter( tab    = tab(1:long, 1:larg),  &  !
                          long   = long,                 &  !
                          larg   = larg,                 &  !
                          snb    = 10,                   &  !
                          kernel = 5,                    &  !
                          sig    = 3._R8,                &  !
                          omp    = GLOBAL_OMP )             !

   return
   endsubroutine smooth__


   subroutine save_sur(tab, scal)
   implicit none
   type(scale_surf), intent(inout)                 :: scal
   real(kind=R8),    intent(inout), dimension(:,:) :: tab

      character(len=512) :: f_name

      integer(kind=I4) :: long, larg

      long = scal%xres
      larg = scal%yres

      f_name = repeat(' ', len(f_name) )

      read(JOB,*) f_name ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, 'NOM_SUR ', trim(f_name)

      call write_surf(nom_fic = trim(f_name), tab_s = tab(1:long, 1:larg), scal = scal)

   return
   endsubroutine save_sur


   subroutine tab_boot(scal, scal_samp, nb_samp, tab_bounds)
   implicit none
   type(scale_surf), intent(in )                            :: scal
   type(scale_surf), intent(out)                            :: scal_samp
   integer(kind=I4), intent(out)                            :: nb_samp
   type(tborne),     intent(out), allocatable, dimension(:) :: tab_bounds

      integer(kind=I4) :: nn_samp, pp_samp
      integer(kind=I4) :: snb, lb1, ub1, lb2, ub2
      integer(kind=I4) :: i, j, k
      integer(kind=I4) :: long, larg

      real(kind=R8) :: wd, ht

      read(JOB,*) nn_samp, pp_samp, nb_samp ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, 'nn_samp ', nn_samp,  &  !
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
   endsubroutine tab_boot


   subroutine sta_loop()
   implicit none

      read(JOB,*) NB_ITER_DEB, NB_ITER_FIN ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, 'nb_iter_deb ', NB_ITER_DEB, &  !
                                                                                                 'nb_iter_fin ', NB_ITER_FIN     !
      SAVE_LIGNE_LUE = LIGNE_LUE - 1

   return
   endsubroutine sta_loop


   subroutine analyses(tab, scal, scal_samp, tab_bounds)
   implicit none
   type(scale_surf), intent(in)                 :: scal
   type(scale_surf), intent(in)                 :: scal_samp
   real(kind=R8),    intent(in), dimension(:,:) :: tab
   type(tborne),     intent(in), dimension(:)   :: tab_bounds

      integer(kind=I4) :: lb1, ub1, lb2, ub2
      integer(kind=I4) :: k
      integer(kind=I4) :: nn, pp
      integer(kind=I4) :: ns, ps
      integer(kind=I4) :: ibatch
      integer(kind=I4) :: exit_status

      character(len=512) :: cwd
      character(len=512) :: ana_file, surf_filename
      character(len=128) :: ana_type

      type(tborne) :: bound

!~       real(kind=R8) :: scal_lx_tmp, scal_ly_tmp

      real(kind=R8), allocatable, dimension(:,:) :: tab_samp

      read(JOB,*) ana_type ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, 'ana_type ', trim(ana_type)

      nn = scal%xres
      pp = scal%yres

      call get_unit(STA)

      call getcwd( cwd )

      ! surface filename without the folders, with .txt instead of .sur
      surf_filename = filename( NOM_SUR )
      surf_filename = surf_filename(1:len_trim(surf_filename) - 4)//".txt"

      ana_file = repeat( ' ', len(ana_file) )
      select case( ana_type(1:8) )

         case('abbott__')
            ana_file = "resu_abboglob"//SEP//trim( surf_filename )

         case('complexi')
            ana_file = "resu_compglob"//SEP//trim( surf_filename )

         case('elli_acv')
            call fftw_plan_with_nthreads( nthreads = NB_THREADS_FFT )
            call init_fftw3( long = nn, larg = pp )

            ana_file = "resu_elliglob"//SEP//trim( surf_filename )

         case('facettes')
            ana_file = "resu_faceglob"//SEP//trim( surf_filename )

         case('ind_frac')
            ana_file = "resu_ind_glob"//SEP//trim( surf_filename )

         case('statisti')
            ana_file = "resu_statglob"//SEP//trim( surf_filename )

         case('topology')
            ana_file = "resu_topoglob"//SEP//trim( surf_filename )

         case default
            stop 'Bad choice of analysis in subroutine analyses'

      endselect

      call make_path(wkd = trim(cwd), file_path = trim(ana_file), exit_status = exit_status)

      open( unit = STA, file = trim(ana_file) )

      call analyses_stat( tab      = tab(1:nn, 1:pp),   &  !
                          sub_samp = .false.,           &  !
                          scal     = scal,              &  !
                          anal     = ana_type,          &  !
                          omp      = GLOBAL_OMP )          !

      close( STA )

      if ( index(ana_type, 'elli_acv') /= 0 ) call end_fftw3()

      if ( WITH_SAMPLING ) then

         bound = tab_bounds( 1 )

         lb1 = bound%lb1 ; lb2 = bound%lb2
         ub1 = bound%ub1 ; ub2 = bound%ub2

         ns = ub1 - lb1 + 1 ; ps = ub2 - lb2 + 1

         allocate( tab_samp(1:ns, 1:ps) )

!~          scal_lx_tmp = scal%lx ; scal%lx = (scal%lx / nn) * ns
!~          scal_ly_tmp = scal%ly ; scal%ly = (scal%ly / pp) * ps

         ibatch = max( ( NB_ITER_FIN - NB_ITER_DEB + 1 ) / NB_THREADS, 1 )

         call get_unit(STA)

         call getcwd( cwd )

         ana_file = repeat( ' ', len(ana_file) )
         select case( ana_type(1:8) )

            case('abbott__')
               ana_file = "resu_abbott__"//SEP//trim( surf_filename )

            case('complexi')
               ana_file = "resu_complexi"//SEP//trim( surf_filename )

            case('elli_acv')
               call fftw_plan_with_nthreads( nthreads = 1 )
               call tab_init_fftw3( long = ns, larg = ns )

               ana_file = "resu_elli_acv"//SEP//trim( surf_filename )

            case('facettes')
               ana_file = "resu_facettes"//SEP//trim( surf_filename )

            case('ind_frac')
               ana_file = "resu_ind_frac"//SEP//trim( surf_filename )

            case('statisti')
               ana_file = "resu_statisti"//SEP//trim( surf_filename )

            case('topology')
               ana_file = "resu_topology"//SEP//trim( surf_filename )

            case default
               stop 'Bad choice of analysis in subroutine analyses'

         endselect

         call make_path(wkd = trim(cwd), file_path = trim(ana_file), exit_status = exit_status)

         open( unit = STA, file = trim(ana_file), share = 'DENYRW' )

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

         if ( ana_type(1:8) == 'elli_acv' ) call tab_end_fftw3()

         deallocate( tab_samp )

!~          scal%lx = scal_lx_tmp
!~          scal%ly = scal_ly_tmp

      endif

   return

   contains

      subroutine analyses_stat(tab, sub_samp, scal, anal, omp)
      implicit none
      type(scale_surf), intent(in )                 :: scal
      character(len=8), intent(in )                 :: anal
      logical(kind=I4), intent(in )                 :: omp
      logical(kind=I4), intent(in )                 :: sub_samp
      real(kind=R8),    intent(in ), dimension(:,:) :: tab

         integer (kind=I4) :: nx, ny

         real(kind=R8)     :: ra_t, md
         real(kind=R8)     :: dx, dy, si

         type(moment_stat) :: mx

         real(kind=R8), dimension(:),   allocatable :: vec_heights
         real(kind=R8), dimension(1:20)             :: ana_res

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
               write(STA,'(9E18.6,T172,a)') ana_res( 1), & ! smr1, iso 25178
                                            ana_res( 2), & ! smr2, iso 25178
                                            ana_res( 3), & ! spk , iso 25178
                                            ana_res( 4), & ! svk , iso 25178
                                            ! 5 et 6 pour off1 et off2
                                            ana_res( 7), & ! sk  , iso 25178
                                            ana_res( 8), & ! core slope
                                            ana_res( 9), & ! adjustment factor (tangent fit)
                                            ana_res(10), & ! coeffa_tan        (tangent fit)
                                            ana_res(11), & ! coeffb_tan        (tangent fit)
                                            '  smr1   smr2   spk   svk   sk   pente   residus'// & !
                                            '  coeffa_tan   coeffb_tan'                            !
               !$omp end critical

               deallocate( vec_heights )

            case( 'topology' )

               call topology( tab  = tab(1:nx, 1:ny),    &  !
                              long = nx,                 &  !
                              larg = ny,                 &  !
                              res  = ana_res(1:6) )         !

               !$omp critical
               write(STA,'(6E18.6,T124,a)')  ana_res(1), & ! number of cells 0.15/0.85
                                             ana_res(2), & ! median size (percentage of surface)
                                             ana_res(3), & ! percentage of surface above 0.15/0.85
                                             ana_res(4), & ! number of cells 0.15/0.95
                                             ana_res(5), & ! median size (percentage of surface)
                                             ana_res(6), & ! percentage of surface above 0.15/0.95
                                             '  nb_cells_1   median_size_1   fraction_of_surface_1'// & !
                                             '  nb_cells_2   median_size_2   fraction_of_surface_2'     !
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
               write(STA,'(8E18.6,T154,a)') ana_res(1),  &  ! max depth
                                            ana_res(2),  &  ! max height
                                            ana_res(3),  &  ! median
                                            ana_res(4),  &  ! absolute roughness
                                            ana_res(5),  &  ! mean height
                                            ana_res(6),  &  ! standard deviation
                                            ana_res(7),  &  ! skewness
                                            ana_res(8),  &  ! kurtosis
                                            '  Sv(signed)   Sp   Md   Sa   Mu   Sq   Ssk   Sku'
               !$omp end critical

            case( 'ind_frac' )

               call indice_fractal( tab_in = tab(1:nx, 1:ny),  &  !
                                    long   = nx,               &  !
                                    larg   = ny,               &  !
                                    indf   = ana_res(1:3) )       !

               !$omp critical
               write(STA, '(3E18.6, T64, a)')   ana_res(1),    &  ! fractal dimension (slope)
                                                ana_res(2),    &  ! ordinate at origin
                                                ana_res(3),    &  ! adjustment factor
                                                '   slope   ord_orig   R2adj'
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
               write(STA, '(2E18.6,T46,a)')  ana_res(1),    &  ! horizontality
                                             ana_res(2),    &  ! relative area
                                             '   %age_fac_5%   rel_area'
               !$omp end critical

            case( 'complexi' )

               call calcul_asfc_hermite( tab_in   = tab(1:nx, 1:ny),    &  !
                                         scal     = scal,               &  !
                                         asfc_res = ana_res(1:2),       &  !
                                         omp      = omp )                  !

               !$omp critical
               write(STA, '(2E18.6,T46,a)')  ana_res(1),    &  ! asfc
                                             ana_res(2),    &  ! adjustment factor
                                             '   asfc2  R2adj'
               !$omp end critical

            case( 'elli_acv' )

               call correlation_parameters( tab       = tab(1:nx, 1:ny),   &  ! IN
                                            long      = nx,                &  ! IN
                                            larg      = ny,                &  ! IN
                                            res       = ana_res(1:8),      &  ! OUT
                                            cut       = 0.5_R8,            &  ! IN
                                            sub_plane = sub_samp,          &  ! IN, remove plane in samplings
                                            scale_xy  = [ dx, dy ],        &  ! IN
                                            multi_fft = sub_samp,          &  ! IN
                                            omp       = omp )                 ! IN

               !$omp critical
               write(STA,'(8E18.6, T154, a)')   ana_res(1), &  ! big axis
                                                ana_res(2), &  ! small axis
                                                ana_res(3), &  ! big/small -> anisotropy 1
                                                ana_res(4), &  ! anisotropy angle
                                                ana_res(5), &  ! radius of greatest slope
                                                ana_res(6), &  ! greatest slope
                                                ana_res(7), &  ! slope anisotropy factor
                                                ana_res(8), &  ! curvature anisotropy factor
                                                '   x_acv   y_acv   x_acv/y_acv   angle(°)'       // & !
                                                '   radius_max_slope   max_slope  max/min_slope'  // & !
                                                '   max/min_curvature'                                 !
               !$omp end critical

            case default

               stop 'Bad choice of analysis in subroutine analyses_stat'

         endselect

      return
      endsubroutine analyses_stat

   endsubroutine analyses


   subroutine end_loop(tab, tab_bounds)
   implicit none
   real(kind=R8), intent(inout), allocatable, dimension(:,:) :: tab
   type(tborne),  intent(inout), allocatable, dimension(:)   :: tab_bounds

      if ( allocated( tab ) )         deallocate( tab )
      if ( allocated( tab_bounds ) )  deallocate( tab_bounds )

   return
   endsubroutine end_loop


   subroutine lssq_img(tab, scal)
   implicit none
   type(scale_surf), intent(inout)                 :: scal
   real(kind=R8),    intent(inout), dimension(:,:) :: tab

      integer(kind=I4) :: long, larg
      integer(kind=I4) :: degx, degy
      integer(kind=I4) :: istat

      real(kind=R8), allocatable, dimension(:,:) :: tab_tmp

      character(len=512) :: surf_filename, new_surf_filename
      character(len=512) :: wkd
      character(len=  6) :: str

      read(JOB,*) degx, degy ; LIGNE_LUE = LIGNE_LUE + 1 ; write(SPY,*) LIGNE_LUE, degx, degy

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

      ! folder name where the new surface is stored. Ex: DEG_8
      write( str, '(i2)' ) degx
      str = 'DEG_'//adjustl( trim(str) )

      ! the folder must be created under the folder where the surface is. Recall that the currect directory
      ! is set in subroutine lect_bas
      call getcwd( wkd )
      call mkdir(wkd = trim(wkd), directory = trim(str), sep = SEP, exit_status = istat)

      ! set the new working directory
      call chdir( trim(wkd)//SEP//trim(str) )
      write(SPY, *) 'New working directory, after LSSQ: ', trim(wkd)//SEP//trim(str)

      ! the name of the new surface
      NOM_SUR = repeat( ' ', len(NOM_SUR) )
      write( NOM_SUR, '(a,2(a,i3.3),a)' ) surf_filename(1:len_trim(surf_filename) - 4), '_X1_', degx, '_Y1_', degy, '.sur'

      ! the new surface is stored in, as an example: sur/DEG8/file_X1_008_Y1_008.sur
      new_surf_filename = repeat( ' ', len(new_surf_filename) )
      new_surf_filename = trim(str)//SEP//trim(NOM_SUR)

      call write_surf(nom_fic = trim(NOM_SUR), tab_s = tab(1:long, 1:larg), scal = scal)

      deallocate( tab_tmp )

   return
   endsubroutine lssq_img

endmodule script

