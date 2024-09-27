module m_cells
  use m_types

  implicit none
  private

  public :: set_ghost_cells_densities
  public :: create_output

contains

  subroutine set_ghost_cells_densities(mesh, BC_type1, BC_type2)
    type(mesh_t), intent(inout) :: mesh
    character(len=*), intent(in) :: BC_type1
    character(len=*), intent(in) :: BC_type2
    integer                      :: i

    select case(BC_type1)
    case('neumann_zero')
       do i = 1, mesh%n_ghostcells
          mesh%u_electron(1-i) = mesh%u_electron(i)
          mesh%u_ion(1-i) = mesh%u_ion(i)
       end do
    case('dirichlet_zero')
       do i = 1, mesh%n_ghostcells
          mesh%u_electron(1-i) = 0.0_dp
          mesh%u_ion(1-i) = 0.0_dp
       end do
    case default
       write(*,*) "*ERROR* Invalid boundary condition type"
    end select

    select case(BC_type2)
    case('neumann_zero')
       do i = 1, mesh%n_ghostcells
          mesh%u_electron(mesh%Nx+i) = mesh%u_electron(mesh%Nx+1-i)
          mesh%u_ion(mesh%Nx+i) = mesh%u_ion(mesh%Nx+1-i)
       end do
    case('dirichlet_zero')
       do i = 1, mesh%n_ghostcells
          mesh%u_electron(mesh%Nx+i) = 0.0_dp
          mesh%u_ion(mesh%Nx+i) = 0.0_dp
       end do
    case default
       write(*,*) "*ERROR* Invalid boundary condition type"
    end select

  end subroutine set_ghost_cells_densities

  subroutine create_output(mesh, i_step)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in)      :: i_step
    character(len=100)       :: filename
    integer                  :: i
    integer                  :: file_unit

    ! Open file
    write(filename, '(A, I6.6, A)') 'output/output_it', i_step, '.txt'

    ! Open file for writing
    open(newunit=file_unit, file=filename, status='replace')

    ! Write header
    write(file_unit, '(A)') '|      x     |       nion     |       phi       |       n      |       E      |'
    write(file_unit, '(A)') '-------------------------------------------------------------'

    ! Write data for each grid point
    do i = 1, mesh%Nx
       write(file_unit, *) &
            mesh%x(i), mesh%u_ion(i), mesh%phi(i), mesh%u_electron(i), mesh%E_cc(i)
    end do

    ! Close the file
    close(file_unit)
  end subroutine create_output
end module m_cells
