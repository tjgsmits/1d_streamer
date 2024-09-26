module m_fluid
  use m_cells

  implicit none
  private

  public :: solve_diffusion
  public :: solve_drift
  public :: source_term
  public :: initial_dens
  public :: gas_density
  public :: CFL_check
  public :: crank_nicolson_matrix

        contains
        ! Solves the diffusion term
        subroutine solve_diffusion(diff, delta_x, Ngrid, BC_type, bc_value1, bc_value2, n, n_diff)
                real, intent(in)    :: diff
                real, intent(in)    :: delta_x
                integer, intent(in) :: Ngrid
                character(len=*), intent(in) :: BC_type
                real, intent(in)             :: bc_value1, bc_value2
                real, intent(in)    :: n(Ngrid)
                real, intent(out)   :: n_diff(Ngrid)
                integer             :: i
                real                :: frac
                real                :: flux(Ngrid+1)
                real                :: n_ghost(Ngrid+4)

                ! Ensure boundary conditions are applied correctly
                if (Ngrid < 2) then
                        print *, "*ERROR* Ngrid must be greater than 1"
                        error stop
                        return
                end if

                frac = diff / delta_x**2
                call ghost_cell(Ngrid, delta_x, BC_type, bc_value1, bc_value2, n, n_ghost)
                do i=1, Ngrid
                        ! At start boundary
                        if (i==1) then
                                ! The flux (f_{i+1/2}) at the cell surface
                                flux(i) =  n(i) - n(i+1)
                        else if (i == Ngrid) then
                                ! The flux (f_{i+1/2}) at the cell surface
                                flux(i) =  n(i) - n_ghost(Ngrid+3)
                        else    
                                ! The flux (f_{i+1/2}) at the cell surface
                                flux(i) =  n(i) - n(i+1)
                        end if
                end do

                do i=2,Ngrid+1
                        ! Calculate the advection flux (f^a = f_{i-1/2} - f_{i+1/2}) in the cell volume            
                        n_diff(i-1) = frac * (flux(i-1) - flux(i+1))
                end do
        end subroutine solve_diffusion 

        subroutine solve_drift(mob, delta_x, Ngrid, BC_type, bc_value1, bc_value2, n, E, n_drift)
                real, intent(in)    :: mob
                real, intent(in)    :: delta_x
                integer, intent(in) :: Ngrid
                character(len=*), intent(in) :: BC_type
                real, intent(in)             :: bc_value1, bc_value2
                real, intent(in)    :: n(Ngrid)
                real, intent(in)    :: E(Ngrid+1)
                real, intent(out)   :: n_drift(Ngrid)
                integer             :: i
                real                :: frac
                real                :: v
                real                :: r
                real                :: flux(Ngrid+1)
                real                :: n_ghost(Ngrid+4)
                
                ! Ensure boundary conditions are applied correctly
                if (Ngrid < 2) then
                        print *, "*ERROR* Ngrid must be greater than 1"
                        error stop
                        return
                end if

                frac = mob / delta_x
                call ghost_cell(Ngrid, delta_x, BC_type, bc_value1, bc_value2, n, n_ghost)

                do i=1, Ngrid
                        ! Calculated the drift velocity
                        v = -mob * E(i)
                        if (i==1) then
                                ! Upwind flux depending on the sign of v
                                if (v >= 0.0) then
                                        r = (n(i+2) - n(i+1)) / (n(i+1) - n(i))
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n(i+1) - koren_limiter(r) * (n(i+1) - n(i)))
                                else if (v < 0.0) then
                                        r = (n(i) - n_ghost(2)) / (n(i+1) - n(i))
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n(i) + koren_limiter(r) * (n(i+1) - n(i)))
                                end if
                        else if (i == Ngrid) then
                                ! Upwind flux depending on the sign of v
                                if (v >= 0.0) then
                                        r = (n_ghost(Ngrid+4) - n_ghost(Ngrid+3)) / (n_ghost(Ngrid+3) - n(i))
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n_ghost(Ngrid+3) - koren_limiter(r) * (n_ghost(Ngrid+3) - n(i)))
                                else if (v < 0.0) then
                                        r = (n(i) - n(i-1)) / (n_ghost(Ngrid+3) - n(i))
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n(i) + koren_limiter(r) * (n_ghost(Ngrid+3) - n(i)))
                                end if     
                        else   
                                ! Upwind flux depending on the sign of v
                                if (v >= 0.0) then
                                        if (i == Ngrid-1) then 
                                                r = (n_ghost(Ngrid+3) - n(i+1)) / (n(i+1) - n(i))
                                        else
                                                r = (n(i+2) - n(i+1)) / (n(i+1) - n(i))
                                        end if
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n(i+1) - koren_limiter(r) * (n(i+1) - n(i)))
                                else if (v < 0.0) then
                                        r = (n(i) - n(i-1)) / (n(i+1) - n(i))
                                        ! The flux (f_{i+1/2}) at the cell surface
                                        flux(i) = v * (n(i) + koren_limiter(r) * (n(i+1) - n(i)))
                                end if     
                        end if
                end do

                do i=2,Ngrid+1
                        ! Calculate the advection flux (f^a = f_{i-1/2} - f_{i+1/2}) in the cell volume            
                        n_drift(i-1) = frac * (flux(i-1) - flux(i+1))
                end do
                
        end subroutine solve_drift

        function koren_limiter(r) result(phi_koren)
                real    :: phi_koren
                real    :: r
                phi_koren = max(0.0, min(1.0, (2.0 + r) / 6.0, r))
        end function koren_limiter
        
        subroutine source_term(Ngrid, n, nion, E, n_gas, source_type, n_source)
                integer, intent(in)          :: Ngrid
                real, intent(in)             :: n(Ngrid), nion(Ngrid)
                real, intent(in)             :: E(Ngrid)
                real , intent(in)            :: n_gas(Ngrid)
                character(len=*), intent(in) :: source_type
                real, intent(out)            :: n_source(Ngrid)
                integer                      :: i
                real                         :: k_att, k_ion

                ! Ensure boundary conditions are applied correctly
                if (Ngrid < 2) then
                        print *, "*ERROR* Ngrid must be greater than 1"
                        error stop
                        return
                end if
                select case(source_type)
                case('analytical')
                        do i=1,Ngrid
                                if (E(i) == 0.0) then
                                        n_source(i) = 0.0
                                else 
                                        n_source(i) = n(i) * exponent(-1.0/E(i)) - nion(i) * exponent(1.0/E(i))
                                end if
                        end do
                case('reaction_rates')
                        k_ion = 1e16
                        k_att = 1e11
                        do i=1,Ngrid
                                if (E(i) == 0.0) then
                                        n_source(i) = 0.0
                                else 
                                        n_source(i) = k_ion * n_gas(i) * n(i) - k_att * n(i) * nion(i)
                                end if
                        end do
                case default
                        write(*,*) "*ERROR* choose a source type"
                        error stop
                end select
        end subroutine source_term

        subroutine initial_dens(Ngrid, N0, sigma, x, dens_type, dens)
                integer, intent(in)             :: Ngrid
                real, intent(in)                :: N0
                real, intent(in)                :: sigma
                real, intent(in)                :: x(Ngrid)
                character(len=*), intent(in)    :: dens_type
                real, intent(out)               :: dens(Ngrid)
                integer                         :: i

                select case(dens_type)
                case("Exponential")
                        do i=1, Ngrid
                                dens(i) = N0 * exponent(-(x(i)**2)/(2 * sigma**2))
                        end do
                case("Constant")
                        dens(:) = N0
                case("Point")
                        dens(:) = 0.0
                        dens(Ngrid/2) = N0
                case("Point_step")
                        if (Ngrid < 5) then
                                write(*,*) "Grid is to small for initial density"
                                error stop
                        end if
                        dens(:) = 0.0
                        dens(Ngrid/2) = N0

                        dens(Ngrid/2+1) = N0/2
                        dens(Ngrid/2-1) = N0/2

                        dens(Ngrid/2+2) = N0/4
                        dens(Ngrid/2-2) = N0/4
                case default
                        write(*,*) "CHOOSE TYPE OF INITIAL DENSITY PROFILE"
                        error stop
                end select
        end subroutine initial_dens

        subroutine gas_density(Ngrid, p, T, n_gas)
                integer, intent(in)     :: Ngrid
                real, intent(in)        :: p
                real, intent(in)        :: T
                real                    :: k_b
                real, intent(out)       :: n_gas(Ngrid)
                ! Assume constant pressure and temperature
                k_b = 1.38e-23
                n_gas(:) = p / (k_b * T)
        end subroutine gas_density

        subroutine CFL_check(Ngrid, velocity, delta_t, delta_x)
                real                :: CFL
                real                :: v_mag
                integer, intent(in) :: Ngrid
                real, intent(in)    :: velocity(Ngrid)
                real, intent(in)    :: delta_t
                real, intent(in)    :: delta_x

                v_mag = SQRT(sum(velocity)**2)
                CFL = (v_mag * delta_t) / delta_x
                if (CFL <= 1.0) then
                        write(*,*) "CFL condition is broken, take smaller time steps"
                        error stop
                end if
        end subroutine CFL_check

        ! IGNORE --> HIER START EEN MATIGE POGING VOOR IMPLICITE SCHEMA'S

        subroutine crank_nicolson_matrix(Ngrid, delta_t, delta_x, diff, mob, E, n_old, n, a, b, c, d, matrix, n_cranck)
                integer, intent(in)     :: Ngrid
                real, intent(in)        :: delta_t
                real, intent(in)        :: delta_x
                real, intent(in)        :: diff
                real, intent(in)        :: mob
                real, intent(in)        :: E(Ngrid)
                real, intent(in)        :: n_old(Ngrid)
                real, intent(in)        :: n(Ngrid)  
                real                    :: nion(Ngrid)  
                real                    :: n_gas(Ngrid), n_source(Ngrid)
                real                    :: pres, Temp
                real,intent(out)        :: a(Ngrid)
                real,intent(out)        :: b(Ngrid)
                real,intent(out)        :: c(Ngrid)
                real,intent(out)        :: d(Ngrid)
                real,intent(out)        :: matrix(Ngrid, Ngrid)
                real,intent(out)        :: n_cranck(Ngrid)
                integer :: i,j
            
                ! Loop through each grid point to set up the matrix
                do i = 1, Ngrid
                        if (i == 1) then 
                                a(i) = 0.0
                        else 
                                a(i) = (delta_t / 4.0) * (diff / delta_x**2 - (mob * E(i-1)) / delta_x)
                        end if
                        b(i) = 1 + delta_t / 2.0 * (diff / delta_x**2)
                        if (i == Ngrid) then 
                                c(i) = 0.0
                        else
                                c(i) = - delta_t / 4.0 * (diff / delta_x**2 + (mob * E(i-1)) / delta_x)
                        end if

                        do j=1, Ngrid
                           ! At start of the matrix
                           if (i == 1 .and. j == 1) then
                                matrix(i, j) = b(i)
                           ! At end of the matrix
                           else if (i == Ngrid .and. j == Ngrid) then
                                matrix(i, j) = b(i)
                           ! Fulling the rest of the matrix
                           else if (i == j - 1) then
                                matrix(i, j) = a(i)
                           else if (i == j) then 
                                matrix(i, j) = b(i)
                           else if (i == j+1) then
                                matrix(i, j) = b(i)
                           else
                                matrix(i,j) = 0.0
                           end if 
                        end do

                        nion(:) = 1.e11
                        pres = 1   !bar
                        Temp = 300 !K
                        call gas_density(Ngrid,pres,Temp,n_gas)
                        call source_term(Ngrid, n, nion, E, n_gas, 'analytical', n_source)
                        ! Right-hand side vector
                        if (i > 1) then 
                                d(i) = n(i) * (1.0 - delta_t / (2 * delta_x**2) * (2 * diff)) &
                                        + (delta_t / (4 * delta_x)) * (diff / (delta_x**2) + (diff * E(i+1) / delta_x)) * n(i+1) &
                                        + (delta_t / (4 * delta_x)) * (diff / (delta_x**2) - (diff * E(i-1) / delta_x)) * n(i-1) &
                                        + delta_t * n_source(i)
                        else 
                                ! Dit is fout
                                d(i) = 0.0
                        end if
                end do
            
                ! Solve the tridiagonal system here (you might want to use a suitable method)
                call solve_tridiag(Ngrid, n_old, a, b, c, d, n_cranck)
            
        end subroutine crank_nicolson_matrix
            
        subroutine solve_tridiag(Ngrid, n_old, a, b, c, d, n)
                integer, intent(in)     :: Ngrid
                real, intent(in)        :: n_old(Ngrid)
                real, intent(in)        :: a(Ngrid)
                real, intent(in)        :: b(Ngrid)
                real, intent(in)        :: c(Ngrid)
                real, intent(in)        :: d(Ngrid)
                real, intent(out)       :: n(Ngrid)

                integer                 :: k,i, max_iter
                real                    :: tol, error
                real                    :: n_old_calc(Ngrid)

                tol = 1.e-3
                max_iter = 1e3

                ! Jacobi Method
                do k=1,max_iter
                        write(*,*) k
                        if (k == 1) then 
                                do i=1, Ngrid
                                        n(i) = (1 / b(i)) * (d(i) - a(i) * n_old(i-1) - c(i) * n_old(i+1))
                                end do
                        else
                                do i=1, Ngrid
                                        n(i) = (1 / b(i)) * (d(i) - a(i) * n_old_calc(i-1) - c(i) * n_old_calc(i+1))
                                end do
                        end if

                        error = maxval(abs(n - n_old))
                        write(*,*) error
                        if (error < tol) then
                                write(*,*) "Convergence is reached"
                                exit
                        end if
                        n_old_calc(:) = n(:)
                        if (k == max_iter) then
                                write(*,*) 'No convergence is reaced'
                                error stop
                        end if
                end do

        end subroutine solve_tridiag
end module m_fluid
