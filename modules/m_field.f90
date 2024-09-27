module m_field
  use m_types
  use m_cells

  implicit none
  private

  public :: calc_electric_field

contains

  subroutine calc_electric_field(mesh, V_1, V_2)
    type(mesh_t), intent(inout) :: mesh
    real(dp), intent(in)        :: V_1, V_2
    real(dp)                    :: a(mesh%Nx-1), b(mesh%Nx), c(mesh%Nx-1)
    real(dp)                    :: d(mesh%Nx)
    real(dp)                    :: e_charge  = 1.602176634e-19
    real(dp)                    :: epsilon_0 = 8.8541878176d-12
    real(dp)                    :: frac
    integer                     :: i

    frac = (mesh%dx**2 * e_charge) / epsilon_0

    ! Set stencil for interior points
    a(:) = 1.0_dp
    b(:) = -2.0_dp
    c(:) = 1.0_dp

    ! Compute right-hand side
    d(:) = -frac * (mesh%u_ion(1:mesh%Nx) - mesh%u_electron(1:mesh%Nx))

    ! Modify stencil and right-hand side at boundaries
    b(1) = -3.0_dp
    d(1) = d(1) - 2 * V_1
    b(mesh%Nx) = -3.0_dp
    d(mesh%Nx) = d(mesh%Nx) - 2 * V_2

    call solve_tridiag(mesh%Nx, a, b, c, d, mesh%phi(1:mesh%Nx))

    ! Compute face-centered electric field
    mesh%E_face(1) = (V_1 - mesh%phi(1))/(0.5_dp * mesh%dx)
    mesh%E_face(mesh%Nx+1) = (mesh%phi(mesh%Nx) - V_2)/(0.5_dp * mesh%dx)
    do i = 2, mesh%Nx
       mesh%E_face(i) = (mesh%phi(i-1) - mesh%phi(i))/mesh%dx
    end do

    ! Note that this field still has a sign
    mesh%E_cc = 0.5_dp * (mesh%E_face(1:mesh%Nx) + mesh%E_face(2:mesh%Nx+1))

  end subroutine calc_electric_field

  subroutine solve_tridiag(Nx, a, b, c, d, phi)
    integer, intent(in)          :: Nx
    real(dp), intent(in)         :: a(2:Nx), b(Nx), c(Nx-1), d(Nx)
    real(dp), intent(out)        :: phi(Nx)
    integer                      :: i
    real(dp)                     :: cprime(Nx)

    ! Thomas algorithm, this part assumes that A is tridiagonal
    cprime(1) = c(1) / b(1)
    phi(1) = d(1) / b(1)

    do i = 2, Nx
       if (i < Nx) cprime(i) = c(i) / (b(i) - a(i) * cprime(i-1))
       phi(i) = (d(i) - a(i) * phi(i-1)) / (b(i) - a(i) * cprime(i-1))
    end do

    ! Backward substitution
    do i = Nx-1, 1, -1
       phi(i) = phi(i) - cprime(i) * phi(i+1)
    end do

  end subroutine solve_tridiag

end module m_field
