module sed_run_options
   ! One command-line SYSTEM for every program in sed/ -- one vocabulary, one
   ! parser, one refusal, one filename-tag rule -- so that a word names the same
   ! physics whichever program reads it, and a program that does not read a word
   ! says so in those terms instead of ignoring it or calling it unknown.
   !
   ! The settings are grouped into independent AXES, and each program DECLARES
   ! which axes it has a referent for before it reads a single argument
   ! (declare_run_options).  Three things follow from that declaration:
   !
   !   * a word of a declared axis is read;
   !   * a word of an axis the program did NOT declare is refused with the
   !     reason -- "that axis has no referent in this program" -- rather than
   !     silently dropped, which is the failure this system exists to prevent;
   !   * anything else is handed back to the caller as its own positional
   !     argument (a submodel name, a wavelength floor, a file path).
   !
   ! The axes, and the words that select their values:
   !
   !   subject        the first argument: which model, or which product, the
   !                  program is asked to work on.  The list is the program's
   !                  own (declare_run_options(subjects=...)).
   !   solver         heuristic | draine | equil | qm | qm_dbcon | qm_stati
   !   grid           euv
   !   field          mathis_orig, logU=X, hardfield
   !   emission       induced, photcut
   !   qm size        nstate=N, nisrf=N
   !   graphite       gra_d03_sphere | gra_d16_sphere | gra_d16_spheroid
   !   stage1         c2
   !   pah vintage    dl07 | ld01
   !   zubko sizedist formula | table
   !   zubko optics   zda | mie_d03
   !
   ! Every combination ACROSS axes is a valid run.  Only two things are refused:
   ! two values of one axis, and a setting the chosen solver does not read
   ! (check_run_options).
   !
   ! Parsing and tagging are kept apart on purpose.  Reading an argument only
   ! records WHAT was asked for; the filename tag is assembled afterwards in ONE
   ! fixed order, so that settings given in any order produce the same filename
   ! -- which is what lets a non-default run coexist with the production one
   ! instead of overwriting it.  The order runs from the coarsest distinction to
   ! the finest: grid, field extent, solver, graphite, optics set, enthalpy,
   ! field normalization, emission term, numerical sizes.  A default run carries
   ! no tag at all, so the production filenames carry none either.

   ! The SYSTEM depends on nothing but the working precision.  Which module
   ! state a setting writes is a separate question, answered by
   ! sed_apply_options, so that a program with no radiation field and no solver
   ! -- calc_enthalpy.x, calc_qtable.x, calc_polext.x -- still gets the same
   ! vocabulary, the same refusals and the same filename tag without linking a
   ! solver it does not use.
   use constants, only: wp
   implicit none
   private
   public :: run_options_t
   public :: declare_run_options, widen_run_options
   public :: read_run_subject, read_run_option
   public :: check_run_options
   public :: run_options_tag, report_run_options, write_run_option_usage
   public :: LAM_LYMAN_UM, T_HARD_EUV, W_HARD_EUV

   ! Artificial hard component (the 'hardfield' setting): a diluted 1e5 K
   ! blackbody occupying the band below the Lyman limit, where the Mathis field
   ! is identically zero. 1e5 K puts the Wien cut-off of the added component at
   ! a few hundred eV, well above 13.6 eV and well inside the widest grid any
   ! shipped model carries, so the illuminated band ends strictly between the
   ! Lyman limit and the short end of the grid -- which is what separates a
   ! field-based single-photon bound from a grid-based one and from a fixed
   ! 13.6 eV one in a single run. The dilution is chosen so the EUV band
   ! deposits roughly as much power as the whole Mathis field: enough for the
   ! hard photons to matter, little enough that the small grains still heat
   ! stochastically rather than settling into equilibrium.
   real(wp), parameter :: T_HARD_EUV   = 1.0e5_wp
   real(wp), parameter :: W_HARD_EUV   = 2.0e-19_wp
   real(wp), parameter :: LAM_LYMAN_UM = 0.0912_wp

   integer, parameter :: MAXSUBJ = 12

   type :: run_options_t
      ! ---- what this program understands (declare_run_options) ---------
      character(len=32) :: program = ''
      integer           :: n_subject_word = 0
      character(len=16) :: subject_word(MAXSUBJ) = ''
      logical :: ax_solver         = .false.
      logical :: ax_grid           = .false.
      logical :: ax_field          = .false.
      logical :: ax_emission       = .false.
      logical :: ax_qm_size        = .false.
      logical :: ax_graphite       = .false.
      logical :: ax_stage1         = .false.
      logical :: ax_pah_xsec       = .false.
      logical :: ax_zubko_sizedist = .false.
      logical :: ax_zubko_optics   = .false.

      ! ---- what was asked for ------------------------------------------
      character(len=32) :: subject = ''

      logical  :: euv = .false.

      logical  :: hard_euv_field = .false.
      logical  :: mathis_orig    = .false.
      logical  :: set_logU       = .false.
      real(wp) :: logU           = 0.0_wp

      character(len=16) :: solver     = 'heuristic'
      character(len=5)  :: qm_cooling = 'dbdis'
      integer           :: n_solver   = 0

      logical :: set_nstate = .false.
      integer :: nstate     = 0
      logical :: set_nisrf  = .false.
      integer :: nisrf      = 0

      logical :: induced       = .false.
      logical :: photon_cutoff = .false.

      ! Blank leaves the model's own choice, whatever the program set.
      character(len=16) :: graphite   = ''
      integer           :: n_graphite = 0

      logical :: stage1_density_corrected = .false.

      ! Which published carbonaceous absorption the PAH <-> graphite blend
      ! takes: the DL07 (2007) vintage the model is defined with, or the
      ! earlier LD01 (2001) one.
      character(len=8) :: pah_xsec   = 'dl07'
      integer          :: n_pah_xsec = 0

      ! Zubko: which size distribution, and which of the two stored optics sets.
      logical           :: zubko_formula = .true.
      integer           :: n_zubko_sizedist = 0
      character(len=16) :: zubko_optics = 'zda'
      integer           :: n_zubko_optics = 0
   end type run_options_t

contains

   subroutine declare_run_options(o, program, subjects, solver, grid, field, &
                                  emission, qm_size, graphite, stage1, &
                                  pah_xsec, zubko_sizedist, zubko_optics)
      ! One call, before any argument is read, saying what this program has a
      ! referent for.  Everything not declared is refused by name.
      type(run_options_t), intent(inout) :: o
      character(len=*),    intent(in)    :: program
      ! The words the first argument may take -- which model, or which product.
      ! Omitted, the program takes no subject and its first argument is a
      ! setting like any other.
      character(len=*), optional, intent(in) :: subjects(:)
      logical, optional, intent(in) :: solver, grid, field, emission, qm_size
      logical, optional, intent(in) :: graphite, stage1, pah_xsec
      logical, optional, intent(in) :: zubko_sizedist, zubko_optics
      integer :: i

      o%program = program
      if (present(subjects)) then
         if (size(subjects) > MAXSUBJ) then
            write(*,'(a)') ' sed_run_options: subject list longer than MAXSUBJ'
            stop 1
         end if
         o%n_subject_word = size(subjects)
         do i = 1, size(subjects)
            o%subject_word(i) = subjects(i)
         end do
      end if
      call widen_run_options(o, solver, grid, field, emission, qm_size, &
                             graphite, stage1, pah_xsec, zubko_sizedist, &
                             zubko_optics)
   end subroutine declare_run_options


   subroutine widen_run_options(o, solver, grid, field, emission, qm_size, &
                                graphite, stage1, pah_xsec, zubko_sizedist, &
                                zubko_optics)
      ! Add axes to the whitelist after the subject has been read, for a
      ! program whose axes depend on which model was named -- the graphite of
      ! the xi blend exists for astrodust and DL07 and not for Zubko, and the
      ! astrodust Stage-1 prefactor for astrodust alone.
      type(run_options_t), intent(inout) :: o
      logical, optional, intent(in) :: solver, grid, field, emission, qm_size
      logical, optional, intent(in) :: graphite, stage1, pah_xsec
      logical, optional, intent(in) :: zubko_sizedist, zubko_optics
      if (present(solver))         o%ax_solver         = o%ax_solver         .or. solver
      if (present(grid))           o%ax_grid           = o%ax_grid           .or. grid
      if (present(field))          o%ax_field          = o%ax_field          .or. field
      if (present(emission))       o%ax_emission       = o%ax_emission       .or. emission
      if (present(qm_size))        o%ax_qm_size        = o%ax_qm_size        .or. qm_size
      if (present(graphite))       o%ax_graphite       = o%ax_graphite       .or. graphite
      if (present(stage1))         o%ax_stage1         = o%ax_stage1         .or. stage1
      if (present(pah_xsec))       o%ax_pah_xsec       = o%ax_pah_xsec       .or. pah_xsec
      if (present(zubko_sizedist)) o%ax_zubko_sizedist = o%ax_zubko_sizedist .or. zubko_sizedist
      if (present(zubko_optics))   o%ax_zubko_optics   = o%ax_zubko_optics   .or. zubko_optics
   end subroutine widen_run_options


   subroutine read_run_subject(o, arg, ok)
      ! Match the first argument against this program's subject list.  ok is
      ! .false. for a word that is not on it, which the caller answers with its
      ! usage -- the subject is required wherever one is declared.
      type(run_options_t), intent(inout) :: o
      character(len=*),    intent(in)    :: arg
      logical,             intent(out)   :: ok
      integer :: i
      ok = .false.
      do i = 1, o%n_subject_word
         if (trim(arg) == trim(o%subject_word(i))) then
            o%subject = trim(arg)
            ok = .true.
            return
         end if
      end do
   end subroutine read_run_subject


   subroutine read_run_option(arg, o, taken)
      ! Read ONE command-line word.  taken = .false. means the word is in no
      ! axis at all and belongs to the caller (a submodel name, a wavelength
      ! floor, a file path).  A word that names an axis this program did not
      ! declare stops here, with the axis named, rather than being ignored.
      character(len=*),    intent(in)    :: arg
      type(run_options_t), intent(inout) :: o
      logical,             intent(out)   :: taken
      integer :: ieq

      taken = .true.
      select case (trim(arg))

      ! ---- wavelength grid --------------------------------------------
      case ('euv')
         if (.not. axis_here(o, o%ax_grid, arg, 'wavelength-grid')) return
         o%euv = .true.

      ! ---- radiation field --------------------------------------------
      case ('hardfield')
         if (.not. axis_here(o, o%ax_field, arg, 'radiation-field')) return
         o%hard_euv_field = .true.
         o%euv            = .true.   ! the photons need a grid to sit on
      case ('mathis_orig')
         if (.not. axis_here(o, o%ax_field, arg, 'radiation-field')) return
         o%mathis_orig = .true.

      ! ---- stochastic solver ------------------------------------------
      case ('heuristic')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'heuristic', 'dbdis')
      case ('draine')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'draine', 'dbdis')
      case ('equil')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'equil', 'dbdis')
      case ('qm')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'qm', 'dbdis')
      case ('qm_dbcon')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'qm', 'dbcon')
      case ('qm_stati')
         if (.not. axis_here(o, o%ax_solver, arg, 'stochastic-solver')) return
         call name_solver(o, 'qm', 'stati')

      ! ---- emission term ----------------------------------------------
      case ('induced')
         if (.not. axis_here(o, o%ax_emission, arg, 'emission-term')) return
         o%induced = .true.
      case ('photcut')
         if (.not. axis_here(o, o%ax_emission, arg, 'emission-term')) return
         o%photon_cutoff = .true.

      ! ---- astrodust Stage-1 enthalpy ---------------------------------
      case ('c2')
         if (.not. axis_here(o, o%ax_stage1, arg, 'Stage-1-enthalpy')) return
         o%stage1_density_corrected = .true.

      ! ---- graphite of the PAH xi blend -------------------------------
      case ('gra_d03_sphere')
         if (.not. axis_here(o, o%ax_graphite, arg, 'graphite')) return
         call name_graphite(o, 'd03_sphere')
      case ('gra_d16_sphere')
         if (.not. axis_here(o, o%ax_graphite, arg, 'graphite')) return
         call name_graphite(o, 'd16_sphere')
      case ('gra_d16_spheroid')
         if (.not. axis_here(o, o%ax_graphite, arg, 'graphite')) return
         call name_graphite(o, 'd16_spheroid')

      ! ---- carbonaceous cross-section vintage -------------------------
      case ('dl07')
         if (.not. axis_here(o, o%ax_pah_xsec, arg, 'carbonaceous-vintage')) return
         o%pah_xsec = 'dl07';  o%n_pah_xsec = o%n_pah_xsec + 1
      case ('ld01')
         if (.not. axis_here(o, o%ax_pah_xsec, arg, 'carbonaceous-vintage')) return
         o%pah_xsec = 'ld01';  o%n_pah_xsec = o%n_pah_xsec + 1

      ! ---- Zubko size distribution and optics set ---------------------
      case ('formula')
         if (.not. axis_here(o, o%ax_zubko_sizedist, arg, 'Zubko-size-distribution')) return
         o%zubko_formula = .true.;   o%n_zubko_sizedist = o%n_zubko_sizedist + 1
      case ('table')
         if (.not. axis_here(o, o%ax_zubko_sizedist, arg, 'Zubko-size-distribution')) return
         o%zubko_formula = .false.;  o%n_zubko_sizedist = o%n_zubko_sizedist + 1
      case ('zda')
         if (.not. axis_here(o, o%ax_zubko_optics, arg, 'Zubko-optics')) return
         o%zubko_optics = 'zda';      o%n_zubko_optics = o%n_zubko_optics + 1
      case ('mie_d03')
         if (.not. axis_here(o, o%ax_zubko_optics, arg, 'Zubko-optics')) return
         o%zubko_optics = 'mie_d03';  o%n_zubko_optics = o%n_zubko_optics + 1

      case default
         ieq = index(arg, '=')
         if (ieq > 1) then
            select case (arg(1:ieq))
            case ('logU=')
               if (.not. axis_here(o, o%ax_field, arg, 'radiation-field')) return
               read(arg(ieq+1:), *) o%logU
               o%set_logU = .true.
            case ('nstate=')
               if (.not. axis_here(o, o%ax_qm_size, arg, 'transition-matrix-size')) return
               read(arg(ieq+1:), *) o%nstate
               o%set_nstate = .true.
            case ('nisrf=')
               if (.not. axis_here(o, o%ax_qm_size, arg, 'transition-matrix-size')) return
               read(arg(ieq+1:), *) o%nisrf
               o%set_nisrf = .true.
            case default
               taken = .false.
            end select
         else
            taken = .false.
         end if
      end select
   end subroutine read_run_option


   logical function axis_here(o, declared, arg, axis) result(yes)
      ! Is this word's axis one this program has a referent for?  When it is
      ! not, the run stops HERE, naming the axis: a setting silently dropped
      ! would leave the run doing something other than what was asked, and a
      ! bare "unknown argument" would suggest the word means nothing anywhere.
      type(run_options_t), intent(in) :: o
      logical,             intent(in) :: declared
      character(len=*),    intent(in) :: arg, axis
      yes = declared
      if (yes) return
      ! Name the model when there is one: an axis can exist in the program and
      ! not in the model it was asked to work on -- the astrodust Stage-1
      ! enthalpy is a setting of calc_sed, but not of its Zubko model.
      if (len_trim(o%subject) > 0) then
         write(*,'(a)') ' '//trim(o%program)//': '''//trim(arg)//''' selects the '// &
            axis//' axis, which '//trim(o%program)//'''s '//trim(o%subject)// &
            ' has no referent for.'
      else
         write(*,'(a)') ' '//trim(o%program)//': '''//trim(arg)//''' selects the '// &
            axis//' axis, which '//trim(o%program)//' has no referent for.'
      end if
      call write_run_option_usage(o)
      stop 1
   end function axis_here


   subroutine name_solver(o, solver, cooling)
      type(run_options_t), intent(inout) :: o
      character(len=*),    intent(in)    :: solver, cooling
      o%solver     = solver
      o%qm_cooling = cooling
      o%n_solver   = o%n_solver + 1
   end subroutine name_solver


   subroutine name_graphite(o, source)
      type(run_options_t), intent(inout) :: o
      character(len=*),    intent(in)    :: source
      o%graphite   = source
      o%n_graphite = o%n_graphite + 1
   end subroutine name_graphite


   subroutine check_run_options(o)
      ! Refuse the two kinds of combination that cannot be honored, and only
      ! those: two values of one axis, and a setting the chosen solver does not
      ! read.  Everything else is a valid run.  Keeping the last of two values
      ! of one axis would write a file whose name claims something the run did
      ! not do, and accepting a setting no solver reads would do the same.
      type(run_options_t), intent(in) :: o

      if (o%n_solver > 1) then
         write(*,'(a)') ' '//trim(o%program)//': give at most one solver of'// &
            ' heuristic / draine / equil / qm / qm_dbcon / qm_stati'
         stop 1
      end if

      if (o%n_graphite > 1) then
         write(*,'(a)') ' '//trim(o%program)//': give at most one graphite of'// &
            ' gra_d03_sphere / gra_d16_sphere / gra_d16_spheroid'
         stop 1
      end if

      if (o%n_pah_xsec > 1) then
         write(*,'(a)') ' '//trim(o%program)//': give at most one of dl07 / ld01'
         stop 1
      end if

      if (o%n_zubko_sizedist > 1) then
         write(*,'(a)') ' '//trim(o%program)//': give at most one of formula / table'
         stop 1
      end if

      if (o%n_zubko_optics > 1) then
         write(*,'(a)') ' '//trim(o%program)//': give at most one of zda / mie_d03'
         stop 1
      end if

      ! nstate / nisrf size the energy-space transition matrix and its
      ! downsampled radiation field; no other solver builds one.
      if ((o%set_nstate .or. o%set_nisrf) .and. trim(o%solver) /= 'qm') then
         write(*,'(a)') ' '//trim(o%program)// &
            ': nstate= and nisrf= size the energy-space transition matrix,'
         write(*,'(a)') '   which only qm / qm_dbcon / qm_stati build.'// &
                        '  Solver named: '//trim(o%solver)
         stop 1
      end if

      ! photcut zeroes the emission of an enthalpy bin at photon energies above
      ! that bin's own enthalpy.  It acts inside the temperature-window
      ! narrowing solver's bin sum, which is the heuristic solver; the qm
      ! solvers carry the same bound in their own emission kernel and the
      ! equilibrium branch has no bins at all.
      if (o%photon_cutoff .and. trim(o%solver) /= 'heuristic') then
         write(*,'(a)') ' '//trim(o%program)// &
            ': photcut bounds a bin''s emission by its own enthalpy inside the'
         write(*,'(a)') '   heuristic narrowing solver, the only solver that'// &
                        ' reads it.  Solver named: '//trim(o%solver)
         stop 1
      end if
   end subroutine check_run_options



   function run_options_tag(o) result(tag)
      ! The filename tag: the tokens of every non-default setting, joined by
      ! '_' in one fixed order.  Empty for a default run, so the production
      ! filenames carry no tag at all.
      type(run_options_t), intent(in) :: o
      character(len=160) :: tag
      character(len=24)  :: piece

      tag = ''
      if (o%euv)            call append_token(tag, 'euv')
      if (o%hard_euv_field) call append_token(tag, 'hardfield')

      select case (trim(o%solver))
      case ('draine');  call append_token(tag, 'draine')
      case ('equil');   call append_token(tag, 'equil')
      case ('qm')
         select case (trim(o%qm_cooling))
         case ('dbcon');  call append_token(tag, 'qm_dbcon')
         case ('stati');  call append_token(tag, 'qm_stati')
         case default;    call append_token(tag, 'qm')
         end select
      end select

      select case (trim(o%graphite))
      case ('d03_sphere');   call append_token(tag, 'd03gra')
      case ('d16_sphere');   call append_token(tag, 'd16gra')
      case ('d16_spheroid'); call append_token(tag, 'sphdgra')
      end select

      if (o%ax_pah_xsec .and. trim(o%pah_xsec) /= 'dl07') &
         call append_token(tag, trim(o%pah_xsec))
      if (o%ax_zubko_optics .and. trim(o%zubko_optics) /= 'zda') &
         call append_token(tag, trim(o%zubko_optics))
      if (o%ax_zubko_sizedist .and. .not. o%zubko_formula) &
         call append_token(tag, 'table')

      if (o%stage1_density_corrected) call append_token(tag, 'c2')
      if (o%mathis_orig)              call append_token(tag, 'morig')
      if (o%set_logU) then
         write(piece,'(a,f0.2)') 'logU', o%logU
         call append_token(tag, trim(piece))
      end if
      if (o%induced) call append_token(tag, 'induced')
      if (o%set_nstate) then
         write(piece,'(a,i0)') 'ns', o%nstate
         call append_token(tag, trim(piece))
      end if
      if (o%set_nisrf) then
         write(piece,'(a,i0)') 'nisrf', o%nisrf
         call append_token(tag, trim(piece))
      end if
      if (o%photon_cutoff) call append_token(tag, 'photcut')
   end function run_options_tag


   subroutine append_token(tag, token)
      character(len=*), intent(inout) :: tag
      character(len=*), intent(in)    :: token
      if (len_trim(tag) == 0) then
         tag = token
      else
         tag = trim(tag)//'_'//token
      end if
   end subroutine append_token


   subroutine report_run_options(o)
      ! The settings this run is about to work with, printed before the work so
      ! that a mistyped word is caught now rather than after it.  Only the axes
      ! this program declared are printed; the rest do not exist here.
      type(run_options_t), intent(in) :: o

      if (o%ax_grid) then
         if (o%euv) then
            write(*,'(a)') ' grid          : carried into the ionizing band'
         else
            write(*,'(a)') ' grid          : cut at the Lyman limit'
         end if
      end if
      if (o%ax_field) then
         if (o%mathis_orig) then
            write(*,'(a)') ' field         : Mathis 1983 as published'// &
                           ' (4000 K dilution 1e-13, CMB 2.9 K)'
         else
            write(*,'(a)') ' field         : Mathis 1983, Draine-corrected'// &
                           ' (4000 K dilution 1.65e-13, CMB 2.725 K)'
         end if
         if (o%hard_euv_field) &
            write(*,'(a,es9.2,a,es9.2)') ' hard EUV band : diluted blackbody, T = ', &
               T_HARD_EUV, ' K, W = ', W_HARD_EUV
      end if
      if (o%ax_solver) then
         write(*,'(a,a)') ' solver        : ', trim(o%solver)
         if (trim(o%solver) == 'qm') then
            write(*,'(a,a)') ' qm cooling    : ', trim(o%qm_cooling)
            if (o%set_nstate) then
               write(*,'(a,i0)') ' qm bins       : ', o%nstate
            else
               write(*,'(a)')    ' qm bins       : (solver default)'
            end if
            if (o%set_nisrf) then
               write(*,'(a,i0)') ' qm field pts  : ', o%nisrf
            else
               write(*,'(a)')    ' qm field pts  : (solver default)'
            end if
         end if
      end if
      if (o%ax_emission) then
         write(*,'(a,l1)') ' induced (1+J/B): ', o%induced
         write(*,'(a,l1)') ' photon cutoff : ', o%photon_cutoff
      end if
   end subroutine report_run_options


   subroutine write_run_option_usage(o)
      ! What THIS program takes: its subject words, then the words of the axes
      ! it declared.  Printed by each program's usage and by a refusal.
      type(run_options_t), intent(in) :: o
      character(len=256) :: line
      integer :: i
      if (o%n_subject_word > 0) then
         line = '   subject   : '//trim(o%subject_word(1))
         do i = 2, o%n_subject_word
            line = trim(line)//' | '//trim(o%subject_word(i))
         end do
         write(*,'(a)') trim(line)
      end if
      if (o%ax_solver) &
         write(*,'(a)') '   solver    : heuristic (default) | draine | equil |'// &
                        ' qm | qm_dbcon | qm_stati'
      if (o%ax_grid) &
         write(*,'(a)') '   grid      : euv         solve the ionizing band as well'
      if (o%ax_field) then
         write(*,'(a)') '   field     : mathis_orig literal Mathis 1983'// &
                        ' (4000 K dilution 1e-13, CMB 2.9 K)'
         write(*,'(a)') '               logU=X      scale the Mathis field to U = 10^X'
         write(*,'(a)') '               hardfield   fill the band below the Lyman'// &
                        ' limit with a diluted'
         write(*,'(a)') '                           1e5 K blackbody (implies euv)'
      end if
      if (o%ax_emission) then
         write(*,'(a)') '   emission  : induced     multiply by the'// &
                        ' stimulated-emission factor (1 + J/B)'
         write(*,'(a)') '               photcut     bound a bin''s emission by its'// &
                        ' own enthalpy (heuristic only)'
      end if
      if (o%ax_qm_size) &
         write(*,'(a)') '   qm sizes  : nstate=N    enthalpy bins'// &
                        '        nisrf=N  field wavelengths'
      if (o%ax_graphite) &
         write(*,'(a)') '   graphite  : gra_d03_sphere | gra_d16_sphere |'// &
                        ' gra_d16_spheroid  (PAH xi blend)'
      if (o%ax_stage1) &
         write(*,'(a)') '   enthalpy  : c2          astrodust Stage-1'// &
                        ' density-corrected prefactor'
      if (o%ax_pah_xsec) &
         write(*,'(a)') '   PAH xsec  : dl07 (default) | ld01   the published'// &
                        ' carbonaceous absorption'
      if (o%ax_zubko_sizedist) &
         write(*,'(a)') '   size dist : formula (default) | table'
      if (o%ax_zubko_optics) &
         write(*,'(a)') '   optics    : zda (default) | mie_d03'
   end subroutine write_run_option_usage



end module sed_run_options
