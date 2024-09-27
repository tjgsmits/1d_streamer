module m_fluid
  use m_cells
  use m_types

  implicit none
  private

  public :: get_diffusive_flux
  public :: get_drift_flux
  public :: get_source_term
  public :: initial_dens
  ! public :: gas_density
  ! public :: CFL_check
  ! public :: crank_nicolson_matrix

contains
  ! Solves the diffusion term
  subroutine get_diffusive_flux(mesh, diff, flux_diffusion)
    type(mesh_t), intent(in) :: mesh
    real(dp), intent(in)     :: diff
    real(dp), intent(out)    :: flux_diffusion(mesh%Nx+1)
    integer                  :: i
    real(dp)                 :: frac

    frac = diff / mesh%dx

    do i = 1, mesh%Nx+1
       flux_diffusion(i) = frac * (mesh%u_electron(i-1) - mesh%u_electron(i))
    end do

  end subroutine get_diffusive_flux

  subroutine get_drift_flux(mesh, mob, flux_drift)
    type(mesh_t), intent(in) :: mesh
    real(dp), intent(in)     :: mob
    real(dp), intent(out)    :: flux_drift(mesh%Nx+1)
    integer                  :: i
    real(dp)                 :: v

    do i = 1, mesh%Nx+1
       ! Get velocity at cell face
       v = -mob * mesh%E_face(i)

       if (v > 0) then
          flux_drift(i) = v * mesh%u_electron(i-1)
       else
          flux_drift(i) = v * mesh%u_electron(i)
       end if
    end do

    ! do i=1, Ngrid
    !    ! Calculated the drift velocity
    !    v = -mob * E(i)
    !    if (i==1) then
    !       ! Upwind flux depending on the sign of v
    !       if (v >= 0.0) then
    !          r = (n(i+2) - n(i+1)) / (n(i+1) - n(i))
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n(i+1) - koren_limiter(r) * (n(i+1) - n(i)))
    !       else if (v < 0.0) then
    !          r = (n(i) - n_ghost(2)) / (n(i+1) - n(i))
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n(i) + koren_limiter(r) * (n(i+1) - n(i)))
    !       end if
    !    else if (i == Ngrid) then
    !       ! Upwind flux depending on the sign of v
    !       if (v >= 0.0) then
    !          r = (n_ghost(Ngrid+4) - n_ghost(Ngrid+3)) / (n_ghost(Ngrid+3) - n(i))
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n_ghost(Ngrid+3) - koren_limiter(r) * (n_ghost(Ngrid+3) - n(i)))
    !       else if (v < 0.0) then
    !          r = (n(i) - n(i-1)) / (n_ghost(Ngrid+3) - n(i))
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n(i) + koren_limiter(r) * (n_ghost(Ngrid+3) - n(i)))
    !       end if
    !    else   
    !       ! Upwind flux depending on the sign of v
    !       if (v >= 0.0) then
    !          if (i == Ngrid-1) then 
    !             r = (n_ghost(Ngrid+3) - n(i+1)) / (n(i+1) - n(i))
    !          else
    !             r = (n(i+2) - n(i+1)) / (n(i+1) - n(i))
    !          end if
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n(i+1) - koren_limiter(r) * (n(i+1) - n(i)))
    !       else if (v < 0.0) then
    !          r = (n(i) - n(i-1)) / (n(i+1) - n(i))
    !          ! The flux (f_{i+1/2}) at the cell surface
    !          flux(i) = v * (n(i) + koren_limiter(r) * (n(i+1) - n(i)))
    !       end if
    !    end if
    ! end do

    ! do i=2,Ngrid+1
    !    ! Calculate the advection flux (f^a = f_{i-1/2} - f_{i+1/2}) in the cell volume            
    !    n_drift(i-1) = frac * (flux(i-1) - flux(i+1))
    ! end do

  end subroutine get_drift_flux

  function koren_limiter(r) result(phi_koren)
    real(dp)    :: phi_koren
    real(dp)    :: r
    phi_koren = max(0.0_dp, min(1.0_dp, (2.0_dp + r) / 6.0_dp, r))
  end function koren_limiter

  subroutine get_source_term(mesh, mob, source_type, source_term)
    type(mesh_t), intent(in)     :: mesh
    real(dp), intent(in)         :: mob
    character(len=*), intent(in) :: source_type
    real(dp), intent(out)        :: source_term(mesh%Nx)

    integer  :: i
    real(dp) :: c1, c2, E_abs

    select case(source_type)
    case('analytical')
       c1 = 1340713.2425426769
       c2 = 29558220.53333022

       do i = 1, mesh%Nx
          E_abs = abs(mesh%E_cc(i))

          if (E_abs < 1e-2_dp * c2) then
             source_term(i) = 0.0
          else
             source_term(i) = c1 * exp(-c2/E_abs) * &
                  E_abs * mob * mesh%u_electron(i)
          end if
       end do
    case default
       write(*,*) "*ERROR* choose a source type"
       error stop
    end select
  end subroutine get_source_term

  subroutine initial_dens(mesh, N0, sigma, dens_type)
    type(mesh_t), intent(inout)  :: mesh
    real(dp), intent(in)             :: N0
    real(dp), intent(in)             :: sigma
    character(len=*), intent(in) :: dens_type
    integer                      :: i

    select case(dens_type)
    case("Exponential")
       do i=1, mesh%Nx
          mesh%u_electron(i) = N0 * exp(-(mesh%x(i)**2)/(2 * sigma**2))
          mesh%u_ion(i) = mesh%u_electron(i)
       end do
    case("Constant")
       mesh%u_electron(:) = N0
       mesh%u_ion(:) = mesh%u_electron(:)
    case("Point_step")
       if (mesh%Nx < 5) then
          write(*,*) "Grid is to small for initial density"
          error stop
       end if

       mesh%u_electron(:) = 0.0
       mesh%u_electron(mesh%Nx/2) = N0

       mesh%u_electron(mesh%Nx/2+1) = N0/2
       mesh%u_electron(mesh%Nx/2-1) = N0/2

       mesh%u_electron(mesh%Nx/2+2) = N0/4
       mesh%u_electron(mesh%Nx/2-2) = N0/4

       mesh%u_ion(:) = mesh%u_electron(:)
    case default
       write(*,*) "CHOOSE TYPE OF INITIAL DENSITY PROFILE"
       error stop
    end select
  end subroutine initial_dens

  ! subroutine gas_density(Ngrid, p, T, n_gas)
  !   integer, intent(in)     :: Ngrid
  !   real(dp), intent(in)        :: p
  !   real(dp), intent(in)        :: T
  !   real(dp)                    :: k_b
  !   real(dp), intent(out)       :: n_gas(Ngrid)
  !   ! Assume constant pressure and temperature
  !   k_b = 1.38e-23
  !   n_gas(:) = p / (k_b * T)
  ! end subroutine gas_density

  ! subroutine CFL_check(mesh, mob, diff, delta_t)
  !   type(mesh_t), intent(in) :: mesh
  !   integer, intent(in)      :: Ngrid
  !   real(dp), intent(in)         :: velocity(Ngrid)
  !   real(dp), intent(in)         :: delta_t
  !   real(dp), intent(in)         :: delta_x
  !   real(dp)                     :: CFL
  !   real(dp)                     :: v_mag

  !   v_mag = SQRT(sum(velocity)**2)
  !   CFL = (v_mag * delta_t) / delta_x
  !   if (CFL <= 1.0) then
  !      write(*,*) "CFL condition is broken, take smaller time steps"
  !      error stop
  !   end if
  ! end subroutine CFL_check

  ! IGNORE --> HIER START EEN MATIGE POGING VOOR IMPLICITE SCHEMA'S

  ! subroutine crank_nicolson_matrix(Ngrid, delta_t, delta_x, diff, mob, E, n_old, n, a, b, c, d, matrix, n_cranck)
  !   integer, intent(in)     :: Ngrid
  !   real(dp), intent(in)        :: delta_t
  !   real(dp), intent(in)        :: delta_x
  !   real(dp), intent(in)        :: diff
  !   real(dp), intent(in)        :: mob
  !   real(dp), intent(in)        :: E(Ngrid)
  !   real(dp), intent(in)        :: n_old(Ngrid)
  !   real(dp), intent(in)        :: n(Ngrid)  
  !   real(dp)                    :: nion(Ngrid)  
  !   real(dp)                    :: n_gas(Ngrid), n_source(Ngrid)
  !   real(dp)                    :: pres, Temp
  !   real(dp),intent(out)        :: a(Ngrid)
  !   real(dp),intent(out)        :: b(Ngrid)
  !   real(dp),intent(out)        :: c(Ngrid)
  !   real(dp),intent(out)        :: d(Ngrid)
  !   real(dp),intent(out)        :: matrix(Ngrid, Ngrid)
  !   real(dp),intent(out)        :: n_cranck(Ngrid)
  !   integer :: i,j

  !   ! Loop through each grid point to set up the matrix
  !   do i = 1, Ngrid
  !      if (i == 1) then 
  !         a(i) = 0.0
  !      else 
  !         a(i) = (delta_t / 4.0) * (diff / delta_x**2 - (mob * E(i-1)) / delta_x)
  !      end if
  !      b(i) = 1 + delta_t / 2.0 * (diff / delta_x**2)
  !      if (i == Ngrid) then 
  !         c(i) = 0.0
  !      else
  !         c(i) = - delta_t / 4.0 * (diff / delta_x**2 + (mob * E(i-1)) / delta_x)
  !      end if

  !      do j=1, Ngrid
  !         ! At start of the matrix
  !         if (i == 1 .and. j == 1) then
  !            matrix(i, j) = b(i)
  !            ! At end of the matrix
  !         else if (i == Ngrid .and. j == Ngrid) then
  !            matrix(i, j) = b(i)
  !            ! Fulling the rest of the matrix
  !         else if (i == j - 1) then
  !            matrix(i, j) = a(i)
  !         else if (i == j) then 
  !            matrix(i, j) = b(i)
  !         else if (i == j+1) then
  !            matrix(i, j) = b(i)
  !         else
  !            matrix(i,j) = 0.0
  !         end if
  !      end do

  !      nion(:) = 1.e11
  !      pres = 1   !bar
  !      Temp = 300 !K
  !      call gas_density(Ngrid,pres,Temp,n_gas)
  !      call source_term(Ngrid, n, nion, E, n_gas, 'analytical', n_source)
  !      ! Right-hand side vector
  !      if (i > 1) then 
  !         d(i) = n(i) * (1.0 - delta_t / (2 * delta_x**2) * (2 * diff)) &
  !              + (delta_t / (4 * delta_x)) * (diff / (delta_x**2) + (diff * E(i+1) / delta_x)) * n(i+1) &
  !              + (delta_t / (4 * delta_x)) * (diff / (delta_x**2) - (diff * E(i-1) / delta_x)) * n(i-1) &
  !              + delta_t * n_source(i)
  !      else 
  !         ! Dit is fout
  !         d(i) = 0.0
  !      end if
  !   end do

  !   ! Solve the tridiagonal system here (you might want to use a suitable method)
  !   call solve_tridiag(Ngrid, n_old, a, b, c, d, n_cranck)

  ! end subroutine crank_nicolson_matrix

  ! subroutine solve_tridiag(Ngrid, n_old, a, b, c, d, n)
  !   integer, intent(in)     :: Ngrid
  !   real(dp), intent(in)        :: n_old(Ngrid)
  !   real(dp), intent(in)        :: a(Ngrid)
  !   real(dp), intent(in)        :: b(Ngrid)
  !   real(dp), intent(in)        :: c(Ngrid)
  !   real(dp), intent(in)        :: d(Ngrid)
  !   real(dp), intent(out)       :: n(Ngrid)

  !   integer                 :: k,i, max_iter
  !   real(dp)                    :: tol, error
  !   real(dp)                    :: n_old_calc(Ngrid)

  !   tol = 1.e-3
  !   max_iter = 1e3

  !   ! Jacobi Method
  !   do k=1,max_iter
  !      write(*,*) k
  !      if (k == 1) then 
  !         do i=1, Ngrid
  !            n(i) = (1 / b(i)) * (d(i) - a(i) * n_old(i-1) - c(i) * n_old(i+1))
  !         end do
  !      else
  !         do i=1, Ngrid
  !            n(i) = (1 / b(i)) * (d(i) - a(i) * n_old_calc(i-1) - c(i) * n_old_calc(i+1))
  !         end do
  !      end if

  !      error = maxval(abs(n - n_old))
  !      write(*,*) error
  !      if (error < tol) then
  !         write(*,*) "Convergence is reached"
  !         exit
  !      end if
  !      n_old_calc(:) = n(:)
  !      if (k == max_iter) then
  !         write(*,*) 'No convergence is reaced'
  !         error stop
  !      end if
  !   end do

  ! end subroutine solve_tridiag

end module m_fluid
