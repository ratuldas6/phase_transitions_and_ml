program ising

  ! Global declaration block
  implicit none
  integer, parameter :: nsnapshots = 250, L = 100
  integer :: snapshot_counter, end_count, temp_index, j, k
  real(kind=16) :: temp, temp_interval
  integer, dimension(L, L) :: spin_in, spin_out
  character(len=20) :: filename

  ! Set initial temperature
  temp = 0.05d0

  ! Set temperature interval and counter end
  temp_interval = 0.05d0
  end_count = 10

  ! Simulation loop
  do temp_index = 1, end_count
    ! Perform Metropolis algorithm for a given temp
    snapshot_counter = 0
    do k = 1, nsnapshots
      call ising_sim(temp, spin_in, spin_out)
      snapshot_counter = snapshot_counter + 1
      write(filename, '(a,i3.3,a,i3.3,a)') 'dat_', int(temp*100), '_', snapshot_counter, '.dat'
      open(unit=10, file=filename, status='replace')
      do j = 1, L
        write(10, '(100i3)') spin_out(j,:)
      end do
      close(10)
    end do
    temp = temp + temp_interval
  end do

contains
  subroutine ising_sim(T, spin_matrix_in, spin_matrix_out) ! Generate snapshots for a chosen T

    ! Local declaration block
    implicit none

    real :: start, finish  ! Timer

    integer :: i, j, step, x, y
    integer, parameter :: L = 100 ! Lattice size
    ! new codefuckery begins here with rand1 and rand2
    real(kind=16) :: energy, delta_E, rand1, rand2
    real(kind=16), dimension(L, L) :: spin_real_matrix
    integer, dimension(L, L) :: spin_matrix_in
    integer, dimension(L, L), intent(out) :: spin_matrix_out
    
    integer, parameter :: N = L*L ! Total number of spins
    integer, parameter :: steps = 1000 ! Number of simulation steps

    real(kind=16), parameter :: coup = 1.0d0 ! Coupling constant
    real(kind=16), parameter :: kB = 1.0d0 ! Boltzmann constant
    real(kind=16), intent(in) :: T ! Temperature

    call cpu_time(start)

    ! Initialize the spin matrix
    call random_number(spin_real_matrix)
    spin_matrix_out = 2*nint(spin_real_matrix) - 1

    ! Run the simulation
    do step = 1, steps
      do i = 1, N

        ! sequential sitepicking (not used right now)
        !x = modulo(i-1, L) + 1
        !y = (i - x)/L + 1

        ! random sitepicking
        call random_number(rand1)
        call random_number(rand2)

        x = ceiling(L*rand1)
        y = ceiling(L*rand2)

        !delta_E = 2.0d0*coup*spin_matrix(x,y)*(spin_matrix(mod(x,L)+1,y) + &
            !spin_matrix(x,mod(y,L)+1) + spin_matrix(mod(x-2,L)+1,y) + &
            !spin_matrix(x,mod(y-2,L)+1))

        delta_E = 2.0d0*coup*spin_matrix_out(x,y)*(spin_matrix_out(mod(x,L)+1,y) + spin_matrix_out(x,mod(y,L)+1))

        if (delta_E <= 0.0d0) then
          spin_matrix_out(x,y) = -spin_matrix_out(x,y)
        else
          if (exp(-delta_E/(kB*T)) > rand()) then
            spin_matrix_out(x,y) = -spin_matrix_out(x,y)
          end if
        end if
      end do
    end do

    ! Time taken to run a snapshot
    call cpu_time(finish)

    print '("Snapshot generated at T = ",f6.3,"J")',T
    print '("Time = ",f6.3," seconds.")',finish-start
  
  end subroutine ising_sim

end program ising