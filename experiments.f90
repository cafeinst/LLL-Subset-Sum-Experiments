!===============================================================================
! experiments.f90
!
! LLL lattice-reduction experiments for randomly generated subset-sum
! problems.
!
! For each experiment, the program:
!
!   1. Generates positive integers s(1), ..., s(n).
!   2. Generates a hidden binary vector y(1), ..., y(n).
!   3. Computes the target
!
!          target = sum(s(i) * y(i)).
!
!   4. Constructs a lattice basis for the subset-sum instance.
!   5. Applies LLL basis reduction.
!   6. Searches the reduced basis columns for an encoded solution.
!   7. Verifies any recovered solution directly.
!
! Implementation notes:
!
!   * Basis vectors are stored as columns of B.
!   * Lattice calculations use double-precision floating-point arithmetic.
!   * Subset-sum integers use 64-bit integer arithmetic.
!   * Gram-Schmidt data are recomputed after every basis modification.
!     This is slower than an optimized LLL implementation, but it is simpler
!     and more reliable for small experimental examples.
!
! Compile:
!
!   gfortran -Wall -Wextra -O2 experiments.f90 -o experiments
!
! Run on Windows:
!
!   experiments.exe
!
! Run on Linux or macOS:
!
!   ./experiments
!===============================================================================


module lll
  use iso_fortran_env, only: dp => real64, int64
  implicit none

contains

  !-----------------------------------------------------------------------------
  ! Compute the Gram-Schmidt orthogonalization of the basis columns.
  !
  ! Input:
  !
  !   B(:,i)  = ith lattice basis vector
  !
  ! Output:
  !
  !   Q(:,i)  = ith Gram-Schmidt orthogonalized vector
  !   mu(i,j) = Gram-Schmidt coefficient of B(:,i) along Q(:,j)
  !   A(i)    = squared Euclidean norm of Q(:,i)
  !-----------------------------------------------------------------------------
  subroutine compute_gs(B, dimension, mu, Q, A)
    implicit none

    integer(int64), intent(in) :: dimension

    real(dp), intent(in)  :: B(dimension,dimension)
    real(dp), intent(out) :: mu(dimension,dimension)
    real(dp), intent(out) :: Q(dimension,dimension)
    real(dp), intent(out) :: A(dimension)

    integer(int64) :: i, j

    real(dp), parameter :: tiny_value = 1.0e-20_dp

    mu = 0.0_dp
    Q  = 0.0_dp
    A  = 0.0_dp

    do i = 1, dimension
       Q(:,i) = B(:,i)

       do j = 1, i - 1
          if (A(j) > tiny_value) then
             mu(i,j) = dot_product(B(:,i), Q(:,j)) / A(j)
          else
             mu(i,j) = 0.0_dp
          end if

          Q(:,i) = Q(:,i) - mu(i,j) * Q(:,j)
       end do

       A(i) = dot_product(Q(:,i), Q(:,i))
    end do

  end subroutine compute_gs


  !-----------------------------------------------------------------------------
  ! Size-reduce basis column k using basis column l.
  !
  ! The nearest integer to mu(k,l) is subtracted as a multiple of B(:,l).
  !-----------------------------------------------------------------------------
  subroutine size_reduce(k, l, mu, B, dimension)
    implicit none

    integer(int64), intent(in) :: k
    integer(int64), intent(in) :: l
    integer(int64), intent(in) :: dimension

    real(dp), intent(in)    :: mu(dimension,dimension)
    real(dp), intent(inout) :: B(dimension,dimension)

    integer(int64) :: reduction_multiple

    if (abs(mu(k,l)) > 0.5_dp) then
       reduction_multiple = nint(mu(k,l), kind=int64)

       B(:,k) = B(:,k) - real(reduction_multiple,dp) * B(:,l)
    end if

  end subroutine size_reduce


  !-----------------------------------------------------------------------------
  ! LLL-reduce the columns of B.
  !
  ! The reduction parameter delta must satisfy
  !
  !                  0.25 < delta < 1.
  !
  ! Values closer to 1 usually produce stronger reduction, but may require
  ! more computation.
  !-----------------------------------------------------------------------------
  subroutine lll_reduction(B, dimension, delta)
    implicit none

    integer(int64), intent(in) :: dimension

    real(dp), intent(in)    :: delta
    real(dp), intent(inout) :: B(dimension,dimension)

    real(dp) :: mu(dimension,dimension)
    real(dp) :: Q(dimension,dimension)
    real(dp) :: A(dimension)
    real(dp) :: temporary_column(dimension)

    integer(int64) :: k, l

    if (delta <= 0.25_dp .or. delta >= 1.0_dp) then
       error stop "LLL error: delta must satisfy 0.25 < delta < 1."
    end if

    call compute_gs(B, dimension, mu, Q, A)

    k = 2

    do while (k <= dimension)

       ! Reduce column k using the immediately preceding basis column.
       call size_reduce(k, k-1, mu, B, dimension)
       call compute_gs(B, dimension, mu, Q, A)

       ! Check the Lovasz condition.
       if (A(k) < (delta - mu(k,k-1)**2) * A(k-1)) then

          ! Swap columns k and k-1.
          temporary_column = B(:,k)
          B(:,k)           = B(:,k-1)
          B(:,k-1)         = temporary_column

          call compute_gs(B, dimension, mu, Q, A)

          k = max(2_int64, k-1)

       else

          ! Finish size-reducing column k using all earlier columns.
          do l = k-2, 1, -1
             call size_reduce(k, l, mu, B, dimension)
             call compute_gs(B, dimension, mu, Q, A)
          end do

          k = k + 1
       end if

    end do

  end subroutine lll_reduction


  !-----------------------------------------------------------------------------
  ! Search the reduced basis columns for a vector encoding a subset solution.
  !
  ! For the embedding used here, a candidate solution vector should have:
  !
  !   * final coordinate approximately zero;
  !   * each of the first n coordinates approximately +1/2 or -1/2.
  !
  ! Both signs are tested because every lattice vector occurs together with
  ! its negative.
  !
  ! Output:
  !
  !   x                recovered binary subset
  !   found            true if a verified solution was found
  !   solution_column  reduced basis column containing the candidate
  !-----------------------------------------------------------------------------
  subroutine find_subset_solution(B, n, s, target, x, found, solution_column)
    implicit none

    integer(int64), intent(in) :: n
    integer(int64), intent(in) :: s(n)
    integer(int64), intent(in) :: target

    real(dp), intent(in) :: B(n+1,n+1)

    integer(int64), intent(out) :: x(n)
    integer(int64), intent(out) :: solution_column

    logical, intent(out) :: found

    integer(int64) :: i, j
    integer(int64) :: candidate_sum

    logical :: candidate_shape

    real(dp), parameter :: zero_tolerance = 1.0e-6_dp
    real(dp), parameter :: half_tolerance = 1.0e-3_dp

    x               = 0
    found           = .false.
    solution_column = 0

    do j = 1, n + 1

       ! A candidate solution must have an approximately zero final coordinate.
       if (abs(B(n+1,j)) >= zero_tolerance) cycle

       candidate_shape = .true.

       ! Check that all first n coordinates are approximately +/- 1/2.
       do i = 1, n
          if (abs(abs(B(i,j)) - 0.5_dp) > half_tolerance) then
             candidate_shape = .false.
             exit
          end if
       end do

       if (.not. candidate_shape) cycle

       ! Interpret the basis column directly.
       do i = 1, n
          x(i) = nint(B(i,j) + 0.5_dp, kind=int64)
       end do

       candidate_sum = sum(s*x)

       if (candidate_sum == target) then
          found           = .true.
          solution_column = j
          return
       end if

       ! Interpret the negative of the basis column.
       do i = 1, n
          x(i) = nint(-B(i,j) + 0.5_dp, kind=int64)
       end do

       candidate_sum = sum(s*x)

       if (candidate_sum == target) then
          found           = .true.
          solution_column = j
          return
       end if

    end do

    x = 0

  end subroutine find_subset_solution


  !-----------------------------------------------------------------------------
  ! Print the Euclidean norm of each reduced basis column.
  !
  ! This compact output is usually easier to inspect than the entire matrix.
  !-----------------------------------------------------------------------------
  subroutine print_column_norms(B, dimension)
    implicit none

    integer(int64), intent(in) :: dimension
    real(dp), intent(in) :: B(dimension,dimension)

    integer(int64) :: j
    real(dp) :: column_norm

    print *
    print *, "Reduced basis column norms:"
    print *, "--------------------------------"
    print *, " Column             Norm"
    print *, "--------------------------------"

    do j = 1, dimension
       column_norm = sqrt(dot_product(B(:,j), B(:,j)))
       write(*,'(I7,3X,ES18.8)') j, column_norm
    end do

    print *, "--------------------------------"

  end subroutine print_column_norms


  !-----------------------------------------------------------------------------
  ! Print a binary subset together with its selected indices and values.
  !-----------------------------------------------------------------------------
  subroutine print_subset(label, subset, s, n)
    implicit none

    character(len=*), intent(in) :: label

    integer(int64), intent(in) :: n
    integer(int64), intent(in) :: subset(n)
    integer(int64), intent(in) :: s(n)

    integer(int64) :: i
    logical :: any_selected

    print *
    print *, trim(label)
    print *, "Binary vector:"
    write(*,'(*(I2,1X))') subset

    print *, "Selected entries:"
    any_selected = .false.

    do i = 1, n
       if (subset(i) == 1) then
          write(*,'("  index ",I3,": ",I14)') i, s(i)
          any_selected = .true.
       end if
    end do

    if (.not. any_selected) then
       print *, "  [empty subset]"
    end if

  end subroutine print_subset

end module lll



program main
  use iso_fortran_env, only: dp => real64, int64
  use lll
  implicit none

  ! Number of integers in each subset-sum problem.
  integer(int64), parameter :: n = 20

  ! Approximate bit length of each randomly generated subset-sum integer.
  !
  ! With n = 20 and bit_length = 40, the expected density is approximately
  !
  !                         20 / 40 = 0.5.
  integer(int64), parameter :: bit_length = 40

  ! Number of random subset-sum instances to test.
  integer(int64), parameter :: number_of_problems = 10

  ! LLL reduction parameter.
  real(dp), parameter :: delta = 0.85_dp

  ! Change to .true. to print every coordinate of every reduced basis vector.
  logical, parameter :: show_full_reduced_basis = .false.

  real(dp) :: B(n+1,n+1)
  real(dp) :: u
  real(dp) :: scaling_factor
  real(dp) :: density
  real(dp) :: column_norm
  real(dp) :: integer_bound

  integer(int64) :: i, j
  integer(int64) :: problem
  integer(int64) :: target
  integer(int64) :: solution_column
  integer(int64) :: successful_problems

  integer(int64) :: s(n)
  integer(int64) :: hidden_subset(n)
  integer(int64) :: recovered_subset(n)

  logical :: found
  logical :: exact_match

  call random_seed()

  successful_problems = 0
  integer_bound       = 2.0_dp**real(bit_length,dp)

  print *
  print *, "============================================================"
  print *, "              LLL SUBSET-SUM EXPERIMENT"
  print *, "============================================================"
  write(*,'("Subset size:             ",I0)') n
  write(*,'("Integer bit length:      ",I0)') bit_length
  write(*,'("Number of experiments:   ",I0)') number_of_problems
  write(*,'("LLL delta:               ",F6.3)') delta
  print *, "============================================================"

  do problem = 1, number_of_problems

     print *
     print *, "============================================================"
     write(*,'(" EXPERIMENT ",I0," OF ",I0)') problem, number_of_problems
     print *, "============================================================"

     !-------------------------------------------------------------------------
     ! Generate random positive subset-sum integers.
     !
     ! Each value lies between 1 and approximately 2^bit_length.
     !-------------------------------------------------------------------------
     do i = 1, n
        call random_number(u)

        s(i) = int(u * integer_bound, kind=int64)

        if (s(i) == 0) s(i) = 1
     end do

     !-------------------------------------------------------------------------
     ! Generate a random hidden binary subset and its target sum.
     !-------------------------------------------------------------------------
     hidden_subset = 0
     target        = 0

     do i = 1, n
        call random_number(u)

        if (u > 0.5_dp) then
           hidden_subset(i) = 1
        else
           hidden_subset(i) = 0
        end if

        target = target + s(i)*hidden_subset(i)
     end do

     !-------------------------------------------------------------------------
     ! Estimate the subset-sum density:
     !
     !                    density = n / log2(max(s)).
     !
     ! Lower-density instances are generally more favorable to lattice-based
     ! subset-sum attacks.
     !-------------------------------------------------------------------------
     density = real(n,dp) /                                           &
               (log(real(maxval(s),dp)) / log(2.0_dp))

     print *
     print *, "Subset-sum entries:"
     write(*,'(*(I14,1X))') s

     call print_subset("Hidden subset:", hidden_subset, s, n)

     print *
     write(*,'("Target sum:            ",I0)') target
     write(*,'("Approximate density:   ",F8.4)') density

     !-------------------------------------------------------------------------
     ! Construct the lattice basis.
     !
     ! The first n columns contain:
     !
     !   * the n-dimensional identity matrix in the first n rows;
     !   * scaled subset-sum entries in the final row.
     !
     ! The final column contains:
     !
     !   * 1/2 in each of the first n rows;
     !   * the scaled target in the final row.
     !-------------------------------------------------------------------------
     scaling_factor = real(n,dp) *                                    &
                      2.0_dp**(real(n,dp)/2.0_dp)

     B = 0.0_dp

     do i = 1, n
        B(i,i) = 1.0_dp
     end do

     do i = 1, n
        B(n+1,i) = scaling_factor * real(s(i),dp)
     end do

     do i = 1, n
        B(i,n+1) = 0.5_dp
     end do

     B(n+1,n+1) = scaling_factor * real(target,dp)

     ! Apply LLL reduction to the basis columns.
     call lll_reduction(B, n+1, delta)

     ! Print only column norms unless full-basis output is requested.
     call print_column_norms(B, n+1)

     if (show_full_reduced_basis) then
        print *
        print *, "Complete reduced basis:"

        do j = 1, n+1
           column_norm = sqrt(dot_product(B(:,j), B(:,j)))

           write(*,'(/,"Column ",I0,"; norm = ",ES16.7)')              &
                j, column_norm

           write(*,'(*(F14.4,1X))') B(:,j)
        end do
     end if

     ! Search the reduced basis columns for an encoded subset solution.
     call find_subset_solution(                                       &
          B, n, s, target, recovered_subset, found, solution_column   &
     )

     print *
     print *, "-------------------- RESULT --------------------"

     if (found) then
        successful_problems = successful_problems + 1

        call print_subset(                                            &
             "Recovered subset:", recovered_subset, s, n              &
        )

        exact_match = all(recovered_subset == hidden_subset)

        print *
        write(*,'("Recovered from column: ",I0)') solution_column
        write(*,'("Recovered sum:         ",I0)') sum(s*recovered_subset)
        write(*,'("Target sum:            ",I0)') target
        write(*,'("Verified solution:     YES")')

        if (exact_match) then
           write(*,'("Same as hidden set:   YES")')
        else
           write(*,'("Same as hidden set:   NO")')
           print *, "A different subset produced the same target."
        end if

     else
        print *, "Verified solution:     NO"
        print *, "No encoded solution appeared as a reduced basis column."
     end if

     print *, "------------------------------------------------"

  end do

  print *
  print *, "============================================================"
  print *, "                    FINAL SUMMARY"
  print *, "============================================================"
  write(*,'("Problems tested:       ",I0)') number_of_problems
  write(*,'("Solutions recovered:  ",I0)') successful_problems
  write(*,'("Success rate:         ",F6.2,"%")')                       &
       100.0_dp * real(successful_problems,dp) /                     &
       real(number_of_problems,dp)
  print *, "============================================================"

end program main
