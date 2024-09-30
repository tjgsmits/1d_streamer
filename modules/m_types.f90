module m_types
        implicit none
        public

        ! Define double precision
        integer, parameter      :: dp = kind(1.0d0)

        type mesh_t
                ! Grid information
                integer                 :: Nx
                real(dp)                :: dx
                real(dp)                :: l_domain
                real(dp), allocatable   :: x(:)

                ! Arrays  that contain the field information
                real(dp), allocatable   :: E_face(:)
                real(dp), allocatable   :: E_cc(:)
                real(dp), allocatable   :: phi(:)
        
                ! Arrays that contains the density information
                real(dp), allocatable   :: u_ion(:)
                real(dp), allocatable   :: u_electron(:)
        
                ! Number of ghostcells
                integer                 :: n_ghostcells
        end type mesh_t
end module m_types
