program streamer
  use m_field
  use m_fluid
  use m_cells
  use m_types

  implicit none

  ! Main program
  real       :: L
  real       :: t_sim
  real       :: dt
  real       :: V_right
  real       :: V_left
  real       :: dx
  integer    :: save_interval, j, n_gc
  type(mesh_t) :: mesh

  t_sim         = 10e-9
  dt            = 1e-13
  V_right       = 0.0
  V_left        = 500.e1
  save_interval = 1000

  mesh%L_domain = 10e-3
  mesh%Nx = 500
  mesh%dx = mesh%L_domain / mesh%Nx
  mesh%n_ghostcells = 1

  n_gc = mesh%n_ghostcells
  allocate(mesh%n_electron(1-n_gc:Nx+n_gc))
  allocate(mesh%n_ion(1-n_gc:Nx+n_gc))
  allocate(mesh%E_face(Nx+1))
  allocate(mesh%E_cc(Nx))
  allocate(mesh%phi(1-n_gc:Nx+n_gc))
  allocate(mesh%x(Nx))

  dx = L / Nx
  do j = 1, Nx
     mesh%x(j) = (j - 1) * dx
  end do

  ! Set-up initial densities
  call initial_dens(mesh, 1.e11_dp, 0.1e-3_dp, x, "Point_step")

  ! Make the simulation
  call simulation(mesh, t_sim, dt, V_left, V_right, save_interval)

contains

  subroutine time_integration(mesh, t_end, delta_t, V_1, V_2, save_interval)
    type(mesh_t), intent(inout) :: mesh
    real, intent(in)            :: t_end
    real, intent(in)            :: delta_t
    real, intent(in)            :: V_1, V_2
    integer, intent(in)         :: save_interval

    real                   :: phi_old(Ngrid)
    integer                :: t, i
    integer, intent(out)   :: tsteps
    real,intent(out)       :: delta_x
    real                   :: diff
    real                   :: mob
    real                   :: pres, Temp
    real                   :: n_enew(Ngrid), n_inew(Ngrid), n_gnew(Ngrid)
    real                   :: n_diff(Ngrid), n_drift(Ngrid), n_source(Ngrid), n_gas(Ngrid)
    character(len=100)     :: filename, sim_type
    integer                :: file = 1
    integer                :: save_interval
    real                   :: x(Ngrid)
    real                   :: v(Ngrid)
    real                   :: matrix(Ngrid, Ngrid)
    real                   :: a(Ngrid), b(Ngrid), c(Ngrid), d(Ngrid)
    real                   :: n_drift_ghost(Ngrid+4), n_diff_ghost(Ngrid+4)
    real                   :: n_ghost(Ngrid+4), nion_ghost(Ngrid+4), n_gas_ghost(Ngrid+4)

    ! Determine gridsize and nr of timesteps
    tsteps = INT(t_end / delta_t)

    ! Setup mobility and diffusion as constant
    diff = 0.1_dp
    mob  = 0.03_dp

    ! Determine the initial potential due to the applied voltage
    call calc_initial_potential(mesh, V_1, V_2)

    ! Determine electric field based on the initial potential
    call calc_electric_field(mesh)

    ! Set up initial densities based on input values
    n_enew(:) = n(:)
    n_inew(:) = nion(:)

    ! ! Calculate the gas density
    ! pres = 1.0   !bar
    ! Temp = 300.0 !K
    ! call gas_density(Ngrid,pres,Temp,n_gas)

    ! Write out initial conditions
    write(*,*) "------------------------------------------------------------------"
    write(*,*) "time_step = ", 0
    write(*,*) "max(n_electron), min(n_electron) = ", maxval(n), minval(n), "1/m^3"
    write(*,*) "max(n_ion), min(n_ion) = ", maxval(nion), minval(nion),"1/m^3"
    write(*,*) "max(n_gas), min(n_gas) = ", maxval(n_gas), minval(n_gas), "1/m^3"
    write(*,*) "max(E), min(E) = ", maxval(E), minval(E),"V/m"
    write(*,*) "max(phi), min(phi) = ", maxval(phi), minval(phi),"V"
    write(*,*) "------------------------------------------------------------------"

    ! Write out the initial conditions
    call create_output(Ngrid, 0, x, nion, phi, n, E)

    ! Start time loop
    sim_type = 'explicit'
    do t = 1, tsteps
       if (sim_type == "explicit") then
          v = -mob * E

          call CFL_check(Ngrid, v, delta_t, delta_t)

          call set_ghost_cells_densities(mesh, "neumann", 0.0_dp)
          call get_diffusive_flux(mesh, diff, flux_diffusion)
          call get_drift_flux(mesh, mob, flux_drift)
          call get_source_term(mesh, 'analytical', source_term)

          do i = 1, Ngrid
             n_enew(i) = n_ghost(i+2) + delta_t * (n_diff_ghost(i+2) + n_drift_ghost(i+2) + n_source(i))
             n_inew(i) = nion_ghost(i+2) + delta_t * n_source(i)
             n_gnew(i) = n_gas_ghost(i+2) - delta_t * n_source(i)
          end do

          ! Update old density field to new density field
          n(:) = n_enew(:)
          nion(:) = n_inew(:)
          n_gas(:) = n_gnew(:)

          ! Update potential
          phi_old(:) = phi(:)
          call calc_potential(delta_x, Ngrid, V_1, V_2, nion, n, phi_old, "Thomas", phi)
          ! Update electric field
          call calc_electric_field_fc(Ngrid, delta_x, V_1, V_2, phi, E)

          ! Write maxium dens, E, phi per time step to check the simulation
          if (MOD(t, save_interval) == 0 .or. t == 1) then
             write(*,*) "------------------------------------------------------------------"
             write(*,*) "time_step = ", t
             write(*,*) "max(n_electron), min(n_electron) = ", maxval(n), minval(n), "1/m^3"
             write(*,*) "max(n_ion), min(n_ion) = ", maxval(nion), minval(nion),"1/m^3"
             write(*,*) "max(n_gas), min(n_gas) = ", maxval(n_gas), minval(n_gas), "1/m^3"
             write(*,*) "max(E), min(E) = ", maxval(E), minval(E),"V/m"
             write(*,*) "max(phi), min(phi) = ", maxval(phi), minval(phi),"V"
             write(*,*) "------------------------------------------------------------------"
          end if

          ! Save data at specified intervals
          if (MOD(t, save_interval) == 0 .or. t == 1) then
             call create_output(Ngrid, t, x, nion, phi, n, E)
          end if

       ! else if (sim_type == "implicit") then
       !    write(*,*) "Nee, dit hoeft niet. Het werkt nog niet"
       !    error stop
       !    ! Hier start een matige poging tot implicit
       !    call solve_diffusion(diff, delta_x, Ngrid, "neumann", 0.0, 0.0, n, n_diff)
       !    call solve_drift(mob, delta_x, Ngrid, "neumann", 0.0, 0.0, n, E, n_drift)
       !    call source_term(Ngrid, n, nion, E, n_gas, 'analytical', n_source)

       !    ! Set up ghost cells
       !    call ghost_cell(Ngrid, delta_x, "neumann", 0.0, 0.0, n_diff, n_drift_ghost)
       !    call ghost_cell(Ngrid, delta_x, "neumann", 0.0, 0.0, n_drift, n_diff_ghost)
       !    do i = 1, Ngrid
       !       call source_term(Ngrid, n, nion, E, n_gas, 'analytical', n_source)
       !       call crank_nicolson_matrix(Ngrid, delta_t, delta_x, diff, mob, E, n, n, a, b, c, d, matrix, n_enew)
       !       n_inew(i) = n_source(i)
       !       call gas_density(Ngrid,pres,Temp,n_gas)
       !    end do

       !    ! Update old density field to new density field
       !    n(:) = n_enew(:)
       !    nion(:) = n_inew(:)
       !    n_gas(:) = n_gnew(:)

       !    ! Update potential
       !    phi_old(:) = phi(:)
       !    call calc_potential(delta_x, Ngrid, V_1, V_2, nion, n, phi_old, 'Thomas', phi)

       !    ! Update electric field
       !    call calc_electric_field_fc(Ngrid, delta_x, V_1, V_2, phi, E)

       !    ! Write maxium dens, E, phi per time step to check the simulation
       !    if (MOD(t, save_interval) == 0 .or. t == 1) then
       !       ! Write maxium dens, E, phi per time step to check the simulation
       !       write(*,*) "---------------------------------"
       !       write(*,*) "time_step = ", t
       !       write(*,*) "max(n_electron) = ", maxval(n), "1/m^3"
       !       write(*,*) "max(n_ion) = ", maxval(nion), "1/m^3"
       !       write(*,*) "max(n_gas) = ", maxval(n_gas), "1/m^3"
       !       write(*,*) "max(E) = ", maxval(E), "V/m"
       !       write(*,*) "max(phi) = ", maxval(phi), "V"
       !       write(*,*) "---------------------------------"
       !    end if

       !    ! Save data at specified intervals
       !    if (MOD(t, save_interval) == 0 .or. t == 1) then
       !       do i = 1, Ngrid
       !          x(i) = (i - 1) * delta_x
       !       end do
       !       ! Save data at specified intervals
       !       if (MOD(t, save_interval) == 0 .or. t == 1) then
       !          call create_output(Ngrid, t, x, nion, phi, n, E)
       !       end if
       !    end if
       else
          write(*,*) 'Choose time stepping scheme'
          error stop
       end if

    end do
    write(*,*) 'Total time =', t * delta_t, 's'
  end subroutine time_integration

end program streamer
