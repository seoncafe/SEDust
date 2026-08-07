module sed_paths
   ! Where the library's own input data lives.
   !
   ! Every file the library opens that the CALLER did not name is resolved
   ! through here: the dielectric functions, the default extinction curves,
   ! the stored cross-section tables.  They used to be compile-time constants
   ! beginning '../data/', i.e. relative to the working directory, so a host
   ! that passed its own data_dir got its model directories from that argument
   ! and its dielectric functions from wherever it happened to be standing.
   ! Two RT hosts reported the same consequence: the model and the optical
   ! constants its optics were computed on could be separated, and the failure
   ! arrived as a runtime abort rather than a build status.
   !
   ! The root defaults to '../data', which is what a driver run from sed/
   ! sees, so nothing that worked before needs changing.  build_dust sets it
   ! from its data_dir argument for the length of the build and restores it
   ! afterwards, which is what makes that entry point relocatable.
   !
   ! Paths handed in by a caller are used as given and never composed with the
   ! root -- an absolute path from a host must stay absolute.
   implicit none
   private
   public :: sed_set_data_root, sed_get_data_root, sed_data_path
   public :: SED_PATHLEN

   integer, parameter :: SED_PATHLEN = 512
   character(len=*), parameter :: DEFAULT_ROOT = '../data'

   character(len=SED_PATHLEN), save :: data_root = DEFAULT_ROOT

contains

   subroutine sed_set_data_root(dir)
      ! Set the root.  A blank argument restores the default rather than
      ! composing every path against an empty string.  One trailing slash is
      ! removed so that sed_data_path has a single composition rule.
      character(len=*), intent(in) :: dir
      integer :: n

      if (len_trim(dir) == 0) then
         data_root = DEFAULT_ROOT
         return
      end if
      data_root = dir
      n = len_trim(data_root)
      if (n > 1) then
         if (data_root(n:n) == '/') data_root(n:n) = ' '
      end if
   end subroutine sed_set_data_root


   function sed_get_data_root() result(d)
      character(len=SED_PATHLEN) :: d
      d = data_root
   end function sed_get_data_root


   function sed_data_path(rel) result(p)
      ! Compose a data-root-relative name.  rel is written without a leading
      ! slash ('dielectric/index_silD03'), so that the two halves join here
      ! and nowhere else.
      character(len=*), intent(in) :: rel
      character(len=SED_PATHLEN)   :: p
      p = trim(data_root)//'/'//trim(rel)
   end function sed_data_path

end module sed_paths
