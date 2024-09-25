module m_cells
    implicit none
    private
    public :: ghost_cell

    contains
    subroutine boundary_conditions(delta_x, Ngrid, BC_type, bc_value1, bc_value2, u_old, u)
        real, intent(in)             :: delta_x
        integer, intent(in)          :: Ngrid
        character(len=*), intent(in) :: BC_type
        real, intent(in)             :: bc_value1, bc_value2
        real, intent(in)            :: u_old(Ngrid)
        real, intent(out)          :: u(Ngrid+4)

        ! Important to note is that the BC are applied at the ghost cells
        select case(BC_type)
        case('neumann')
                u(2)         = u_old(1) - bc_value1 * delta_x
                u(1)         = u_old(1) - bc_value1 * delta_x
                u(size(u)-1) = u_old(Ngrid) + bc_value2 * delta_x
                u(size(u))   = u_old(Ngrid) + bc_value2 * delta_x
        case('dirichlet')
                u(2)         = bc_value1    
                u(1)         = bc_value1
                u(size(u)-1) = bc_value2
                u(size(u))   = bc_value2
        case default
            write(*,*) "*ERROR* Invalid boundary condition type"
        end select
    end subroutine boundary_conditions

    subroutine ghost_cell(Ngrid, delta_x, BC_type, bc_value1, bc_value2, input_value, ghost_layer)
        integer, intent(in)          :: Ngrid
        real, intent(in)             :: delta_x
        character(len=*), intent(in) :: BC_type
        real, intent(in)             :: bc_value1, bc_value2
        real, intent(in)             :: input_value(Ngrid)
        real, intent(out)            :: ghost_layer(Ngrid+4)
        integer                      :: i

        ! Create 2 ghost cells at the outer boundaries
        do i=3, Ngrid+2  
            if (i >= 3 .and. i <= Ngrid+2) then
                ! Fill real cells inside domain
                ghost_layer(i) = input_value(i-2)
            end if
        end do
        ! Enforce BC on the ghost cells outside the domain
        call boundary_conditions(delta_x, Ngrid, BC_type, bc_value1, bc_value2, input_value, ghost_layer)
    end subroutine ghost_cell

end module m_cells