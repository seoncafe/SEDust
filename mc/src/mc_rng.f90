module mc_rng
   ! Thread-local pseudo-random number generator for OpenMP parallel MC.
   !
   ! Uses xorshift64*, a single-stream PRNG with 2^64 - 1 period.  Each
   ! thread owns its own rng_t state, which removes the need for any
   ! locking in the inner MC loop.  Quality is more than adequate for
   ! Monte Carlo absorption/cooling simulations (it passes BigCrush on
   ! all but a few of the most demanding tests, which are not relevant
   ! here).
   !
   ! Usage:
   !   type(rng_t) :: r
   !   call rng_init(r, seed=42)
   !   u = rng_uniform(r)              ! uniform in (0,1]
   !   tau = -log(u) / lambda          ! exponential deviate
   !
   ! To seed independent streams in an OpenMP parallel region:
   !   !$omp parallel private(r)
   !   call rng_init(r, seed = base_seed + 1009 * omp_get_thread_num())
   !   ...
   !   !$omp end parallel
   !
   ! The state is updated by ieor() and ishft() with constants chosen
   ! per Marsaglia (2003) and refined by Vigna (xorshift64*).

   use constants, only: wp
   implicit none
   private
   public :: rng_t, rng_init, rng_uniform, rng_exp

   type :: rng_t
      integer(kind=8) :: state = int(z'9E3779B97F4A7C15', kind=8)
   end type

   ! xorshift64* multiplier (Vigna 2014)
   integer(kind=8), parameter :: MULT = int(z'2545F4914F6CDD1D', kind=8)
   ! Use only the high 53 bits to map cleanly to double precision in (0,1).
   integer(kind=8), parameter :: MASK53 = int(z'001FFFFFFFFFFFFF', kind=8)
   real(wp),        parameter :: SCALE53 = 1.0_wp / 9007199254740992.0_wp  ! 2^53

contains

   subroutine rng_init(r, seed)
      ! Initialize the RNG state.  seed = 0 is treated as 1 to avoid
      ! the zero fixed point of xorshift.  Negative seeds use the wall
      ! clock (NOT thread-safe; only use from a single thread).
      type(rng_t), intent(out) :: r
      integer,     intent(in)  :: seed
      integer(kind=8) :: s
      integer :: clock, i
      if (seed < 0) then
         call system_clock(count=clock)
         s = int(clock, kind=8) + 1_8
      else if (seed == 0) then
         s = 1_8
      else
         s = int(seed, kind=8)
      end if
      r%state = s
      ! Warm up by discarding a few iterations.
      do i = 1, 8
         call advance(r)
      end do
   end subroutine rng_init


   subroutine advance(r)
      type(rng_t), intent(inout) :: r
      r%state = ieor(r%state, ishft(r%state, 12))
      r%state = ieor(r%state, ishft(r%state, -25))
      r%state = ieor(r%state, ishft(r%state, 27))
   end subroutine advance


   function rng_uniform(r) result(u)
      ! Uniform deviate in (0, 1].  Never returns exactly 0.
      type(rng_t), intent(inout) :: r
      real(wp) :: u
      integer(kind=8) :: x
      call advance(r)
      x = r%state * MULT
      u = real(iand(ishft(x, -11), MASK53), wp) * SCALE53
      if (u <= 0.0_wp) u = SCALE53      ! never 0
   end function rng_uniform


   function rng_exp(r, rate) result(t)
      ! Exponentially distributed deviate with mean 1/rate.
      type(rng_t), intent(inout) :: r
      real(wp),    intent(in)    :: rate
      real(wp) :: t
      t = -log(rng_uniform(r)) / rate
   end function rng_exp

end module mc_rng
