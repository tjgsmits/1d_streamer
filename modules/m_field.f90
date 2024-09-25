module m_field
  use m_cells

  implicit none
  private

  public :: calc_initial_potential
  public :: calc_potential
  public :: calc_electric_field_fc
  public :: calc_electric_field_cc

        contains
        subroutine calc_initial_potential(Ngrid, V_1, V_2, phi)
                integer, intent(in)    :: Ngrid
                real, intent(in)       :: V_1, V_2
                real, intent(out)      :: phi(Ngrid)  ! Use assumed-size array
                integer                :: i
            
                ! Ensure boundary conditions are applied correctly
                if (Ngrid < 2) then
                    print *, "*ERROR* Ngrid must be greater than 1"
                    return
                end if
            
                ! Set potential values
                do i = 1, Ngrid
                    phi(i) = V_1 + (V_2 - V_1) / (Ngrid - 1) * (i - 1)
                end do
            end subroutine calc_initial_potential

        subroutine calc_potential(delta_x, Ngrid, V_1, V_2, nion, n, phi_old, solve_type, phi_calc)  
            real, intent(in)                :: delta_x
            integer, intent(in)             :: Ngrid
            real, intent(in)                :: V_1, V_2
            real, intent(in)                :: nion(Ngrid), n(Ngrid)
            real, intent(in)                :: phi_old(Ngrid)
            character(len=6), intent(in)    :: solve_type
            real, intent(out)               :: phi_calc(Ngrid)
            real                            :: phi_ghost(Ngrid+4)
            real                            :: a(Ngrid), b(Ngrid), c(Ngrid)
            real                            :: d(Ngrid), matrix(Ngrid,Ngrid)
            real                            :: e_charge = 1.6e-19
            real                            :: epsilon_0 = 8.85e-12
            real                            :: frac
            integer                         :: i,j

            frac = (delta_x**2 * e_charge) / epsilon_0
            a(:) = 1.0
            b(:) = -2.0
            c(:) = 1.0

            ! Enforce BC
            call ghost_cell(Ngrid, delta_x, "dirichlet", V_1, V_2, phi_old, phi_ghost)

            ! Set up A and b to solve A * phi = b
            do i = 1, Ngrid
                    do j=1, Ngrid
                       ! At start of the matrix
                       if (i == 1 .and. j == 1) then
                            ! Ensure proper handling of BC
                            matrix(i, :) = 0.0
                            matrix(i, j) = 1.0
                       ! At end of the matrix
                       else if (i == Ngrid .and. j == Ngrid) then
                            ! Ensure proper handling of BC
                            matrix(i, :) = 0.0
                            matrix(i, j) = 1.0

                       ! Filling the rest of the matrix
                       else if (i == j - 1 .and. i /= 1) then
                            matrix(i, j) = a(i)
                       else if (i == j) then 
                            matrix(i, j) = b(i)
                       else if (i == j+1) then
                            matrix(i, j) = c(i)
                       else
                            matrix(i,j) = 0.0
                       end if 
                    end do
                ! Right-hand side vector
                if (i == 1) then
                        ! Apply left boundary conditions
                        d(i) = V_1
                else if (i == Ngrid) then
                        ! Apply right boundary conditions
                        d(i) = V_2
                else 
                        d(i) = frac * (nion(i) - n(i))
                        
                end if
            end do

            call solve_tridiag(Ngrid, phi_old, matrix, d, solve_type, phi_calc)
        end subroutine calc_potential
        
        subroutine solve_tridiag(Ngrid, n_old, matrix, d, solve_type, n)
            integer, intent(in)             :: Ngrid
            real, intent(in)                :: n_old(Ngrid)
            real,intent(in)                 :: matrix(Ngrid, Ngrid), d(Ngrid)
            character(len=6), intent(in)    :: solve_type
            real, intent(out)               :: n(Ngrid)
            integer                         :: k,i, max_iter
            real                            :: tol, error
            real                            :: a(Ngrid), b(Ngrid), c(Ngrid)
            real                            :: n_old_calc(Ngrid)
            real                            :: cprime(Ngrid), dprime(Ngrid)

            if (solve_type == "Jacobi") then
                ! Jacobi method
                tol = 1.e-5
                 max_iter = 1e6
                n_old_calc = 0.0  ! Initialize old guess for the Jacobi iteration

                do k = 1, max_iter
                    ! Jacobi iteration
                    do i = 1, Ngrid
                        n(i) = (1 / matrix(i, i)) * (d(i) - sum(matrix(i, :) * n_old_calc(:)))
                    end do

                    error = maxval(abs(n - n_old_calc))
                    if (error < tol) then
                        write(*,*) "Convergence is reached after ", k, " iterations"
                        exit
                    end if

                    n_old_calc = n  ! Update the old guess with the current iteration
                    if (k == max_iter) then
                        write(*,*) "Jacobi method did not converge"
                        error stop
                    end if
                end do

            else if (solve_type == "Thomas") then
                ! Thomas algorithm, this part assumes that A is tridiagonal
                ! Extract a, b, and c from the full matrix A for Thomas
                 do i = 2, Ngrid
                    a(i) = matrix(i, i-1)  ! Lower diagonal
                    b(i) = matrix(i, i)    ! Main diagonal
                    c(i) = matrix(i-1, i)  ! Upper diagonal
                end do
                a(1) = 0.0
                b(1) = matrix(1, 1)   ! First main diagonal element
                c(1) = matrix(1, 2)   ! First upper diagonal element
                cprime(1) = c(1) / b(1)
                dprime(1) = d(1) / b(1)

                do i = 2, Ngrid
                    cprime(i) = c(i) / (b(i) - a(i) * cprime(i-1))
                    dprime(i) = (d(i) - a(i) * dprime(i-1)) / (b(i) - a(i) * cprime(i-1))
                end do
            
                ! Backward substitution
                n(Ngrid) = dprime(Ngrid)
                do i = Ngrid-1, 1, -1
                    n(i) = dprime(i) - cprime(i) * n(i+1)
                end do

                else
                    write(*,*) "Choose solve type electric potential"
                    error stop
                end if
        end subroutine solve_tridiag

        subroutine calc_electric_field_fc(Ngrid, delta_x, V_1, V_2, phi, E)
                integer, intent(in)    :: Ngrid
                real, intent(in)       :: delta_x
                real, intent(in)       :: V_1, V_2
                real, intent(in)       :: phi(Ngrid)
                real, intent(out)      :: E(Ngrid+1)
                real                   :: phi_ghost(Ngrid+4)
                integer                :: i
            
                ! Ensure boundary conditions are applied correctly
                if (Ngrid < 2) then
                    print *, "*ERROR* Ngrid must be greater than 1"
                    return
                end if

                ! Construct ghost layer
                call ghost_cell(Ngrid, delta_x, "dirichlet", V_1, V_2, phi, phi_ghost)

                ! Calculate electric field based on potential
                do i = 1, Ngrid+1
                    if (i == 1) then
                        ! At boundary
                        E(i) = -(phi_ghost(2) - phi(2)) / delta_x
                    else if (i == Ngrid .or. i == Ngrid + 1) then
                        ! At boundary
                        E(i) = -(phi(Ngrid-1) - phi_ghost(Ngrid+3)) / delta_x
                    else
                        E(i) = -(phi(i) - phi(i+1)) / delta_x
                    end if
                end do
                
        end subroutine calc_electric_field_fc

        subroutine calc_electric_field_cc(Ngrid, E, E_cc)
            integer, intent(in) :: Ngrid
            real, intent(in)    :: E(Ngrid)
            real, intent(out)   :: E_cc(Ngrid)
            integer             :: i
        
            do i=1, Ngrid
                if (i == 1) then
                    ! Forward difference for the left boundary
                    E_cc(i) = (E(i+1) - E(i)) / 2.0
                else if (i == Ngrid) then
                    ! Backward difference for the right boundary
                    E_cc(i) = (E(i) - E(i-1)) / 2.0
                else
                    E_cc(i) = (E(i-1) + E(i+1)) / 2.0
                end if
            end do
        end subroutine calc_electric_field_cc
end module m_field