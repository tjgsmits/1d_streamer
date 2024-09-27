program streamer
  use m_field
  use m_fluid
  use m_cells
  use m_types

  implicit none

  ! Main program
  real(dp)       :: t_sim
  real(dp)       :: dt
  real(dp)       :: V_right
  real(dp)       :: V_left
  integer    :: save_interval, j, n_gc
  type(mesh_t) :: mesh

  t_sim         = 3e-9_dp
  dt            = 1e-11_dp
  V_right       = 0.0_dp
  V_left        = 15e3_dp
  save_interval = 1

  mesh%L_domain = 2e-3_dp
  mesh%Nx = 512
  mesh%dx = mesh%L_domain / mesh%Nx
  mesh%n_ghostcells = 1

  n_gc = mesh%n_ghostcells
  allocate(mesh%u_electron(1-n_gc:mesh%Nx+n_gc))
  allocate(mesh%u_ion(1-n_gc:mesh%Nx+n_gc))
  allocate(mesh%E_face(mesh%Nx+1))
  allocate(mesh%E_cc(mesh%Nx))
  allocate(mesh%phi(mesh%Nx))
  allocate(mesh%x(mesh%Nx))

  do j = 1, mesh%Nx
     mesh%x(j) = (j - 1) * mesh%dx
  end do

  ! Set-up initial densities
  call initial_dens(mesh, 1.0e11_dp, 0.1e-3_dp, "Point_step")

  ! Make the simulation
  call time_integration(mesh, t_sim, dt, V_left, V_right, save_interval)

contains

  subroutine time_integration(mesh, t_end, delta_t, V_1, V_2, save_interval)
    type(mesh_t), intent(inout) :: mesh
    real(dp), intent(in)            :: t_end
    real(dp), intent(in)            :: delta_t
    real(dp), intent(in)            :: V_1, V_2
    integer, intent(in)         :: save_interval

    integer            :: tsteps
    integer            :: i_step, i
    real(dp)           :: diff
    real(dp)           :: mob

    real(dp), allocatable :: flux_diffusion(:), flux_drift(:), source_term(:)

    allocate(flux_diffusion(mesh%Nx+1))
    allocate(flux_drift(mesh%Nx+1))
    allocate(source_term(mesh%Nx))

    ! Determine gridsize and nr of timesteps
    tsteps = INT(t_end / delta_t)

    ! Setup mobility and diffusion as constant
    diff = 0.1_dp
    mob  = 0.03_dp

    ! Determine the initial potential due to the applied voltage
    call calc_electric_field(mesh, V_1, V_2)

    ! Write out initial conditions
    write(*,*) "------------------------------------------------------------------"
    write(*,*) "time_step = ", 0
    write(*,*) "max(u_electron), min(u_electron) = ", maxval(mesh%u_electron), &
         minval(mesh%u_electron), "1/m^3"
    write(*,*) "max(u_ion), min(u_ion) = ", maxval(mesh%u_ion), &
         minval(mesh%u_ion),"1/m^3"
    write(*,*) "max(E), min(E) = ", maxval(mesh%E_cc), minval(mesh%E_cc), "V/m"
    write(*,*) "max(phi), min(phi) = ", maxval(mesh%phi), minval(mesh%phi),"V"
    write(*,*) "------------------------------------------------------------------"

    ! Write out the initial conditions
    call create_output(mesh, 0)

    ! Start time loop
    do i_step = 1, tsteps
       ! call CFL_check(mesh, mob, diff, delta_t)

       call set_ghost_cells_densities(mesh, "dirichlet_zero", "dirichlet_zero")
       call get_diffusive_flux(mesh, diff, flux_diffusion)
       call get_drift_flux(mesh, mob, flux_drift)
       call get_source_term(mesh, mob, 'analytical', source_term)

       ! Forward-Euler update
       do i = 1, mesh%Nx
          mesh%u_electron(i) = mesh%u_electron(i) + delta_t/mesh%dx * &
               (flux_diffusion(i) + flux_drift(i) &
               - flux_diffusion(i+1) - flux_drift(i+1)) + &
               delta_t * source_term(i)
          mesh%u_ion(i) = mesh%u_ion(i) + delta_t * source_term(i)
       end do

       ! Update electric field
       call calc_electric_field(mesh, V_1, V_2)

       ! Save data at specified intervals
       if (MOD(i_step, save_interval) == 0) then
          call create_output(mesh, i_step)
       end if

    end do
    write(*,*) 'Total time =', i_step * delta_t, 's'
  end subroutine time_integration

end program streamer
