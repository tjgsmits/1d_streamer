module m_read_in_data

    implicit none
    private
    public :: table_from_file
    public :: lin_interp
  
    contains

    ! Taken from Afivo-streamer
    subroutine table_from_file(file_name, data_name, x_data, y_data)
        character(len=*), intent(in)       :: file_name, data_name
        real, allocatable, intent(out) :: x_data(:), y_data(:)
        integer, parameter :: table_max_rows   = 1500
    
        ! Temporary variables
        integer                   :: ioState, nL
        integer                   :: n_rows
        integer                   :: my_unit
        character(LEN=40)         :: line_fmt
        character(LEN=248)        :: line
        real                      :: temp_table(2, table_max_rows)
        real                      :: factor
    
        nL = 0 ! Set the number of lines to 0
    
        ! Set the line format to read, only depends on string_len currently
        write(line_fmt, FMT = "(I6)") 248
        line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"
    
        ! Open 'file_name' (with error checking)
        open(newunit=my_unit, file = trim(file_name), action = "read", &
             err = 999, iostat = ioState, status="old")
        do
           ! Search for 'data_name' in the file
           do
              read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
              if (line == data_name) exit
           end do
    
           factor = 1.0
    
           ! Now we can check whether there is a comment, while scanning lines until
           ! dashes are found, which indicate the start of the data
           do
              read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
              line = adjustl(line)
              if ( line(1:5) == "-----" ) then
                 exit
              else if (line(1:7) == "FACTOR:") then
                 read(line(8:), *) factor
              else if (line(1:8) == "COMMENT:") then
                 continue
              else
                 print *, "In file ", trim(file_name), " at line", nL
                 print *, trim(line)
                 error stop "Unknown statement in input file"
              end if
           end do
    
           ! Read the data into a temporary array
           n_rows = 0
           do
              read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
              line = adjustl(line)
              if ( line(1:5) == "-----" ) then
                 exit  ! Dashes mark the end of the data
              else if (trim(line) == "" .or. line(1:1) == "#") then
                 cycle ! Ignore whitespace or comments
              else if (n_rows < table_max_rows) then
                 n_rows = n_rows + 1
                 read(line, FMT = *, ERR = 999, end = 777) temp_table(:, n_rows)
              else
                 print *, "CS_read_file error: too many rows in ", &
                      file_name, " at line ", nL
              end if
           end do
    
           ! Store the data in the actual table
           if (allocated(x_data)) deallocate(x_data)
           if (allocated(y_data)) deallocate(y_data)
           allocate(x_data(n_rows))
           allocate(y_data(n_rows))
    
           x_data = temp_table(1, 1:n_rows)
           y_data = factor * temp_table(2, 1:n_rows)
    
           exit                   ! Done
        end do
    
        close(my_unit)
        return
    
    777 continue ! If the end of the file is reached after finding data
        print *, "table_from_file unexpectedly reached end of " // trim(file_name)
        print *, "searching '" // trim(data_name) // "'"
        call print_usage(file_name, data_name)
        error stop
    
    888 continue ! If the end of the file is reached without finding data
        print *, "table_from_file: no data in " // trim(file_name)
        print *, "searching '" // trim(data_name) // "'"
        call print_usage(file_name, data_name)
        error stop
    
    999 continue ! If there was an input error, the routine will end here
        print *, "table_from_file error at line", nL
        print *, "ioState = ", ioState, " in ", trim(file_name)
        print *, "searching '" // trim(data_name) // "'"
        call print_usage(file_name, data_name)
        error stop
    
      contains
    
        subroutine print_usage(file_name, data_name)
          character(len=*), intent(in) :: file_name
          character(len=*), intent(in) :: data_name
    
          print *, ""
          print *, "Expected a file '", trim(file_name), &
               "' with the following structure:"
          print *, trim(data_name)
          print *, "FACTOR: 1.0         [optional: multiply with this factor]"
          print *, "[other lines]"
          print *, "------------------  [at least 5 dashes]"
          print *, "xxx       xxx       [data in two column format]"
          print *, "...       ..."
          print *, "xxx       xxx"
          print *, "------------------"
          print *, ""
        end subroutine print_usage
    
      end subroutine table_from_file

      subroutine lin_interp(x, y, m, x_new, y_new)
        ! Input parameters
        integer, intent(in)  :: m
        real, intent(out)    :: y_new(m)
        real, intent(out)    :: x_new(m)
        real, intent(in)     :: x(:)
        real, intent(in)     :: y(:) 
        real                 :: step
        integer              :: i, j, n

        ! Size of old table
        n = size(x)

        ! Calculate step size
        if (n > 1) then
            step = (x(n) - x(1)) / real(n - 1)
        else
            step = 0.0  ! If n is 1, step will be zero
        end if 
        
        write(*,*) m 
        ! Fill the result array
        do i = 1, m
            write(*,*) "test 2"
            write(*,*) x(1), (i-1), step
            x_new(i) = x(1) + (i - 1) * step
            write(*,*) "test 3"
        end do

        ! Loop over new x values
        j = 1
        do i = 1, m
            ! Find the interval in which x_new(i) lies
            do while (j < n .and. x_new(i) > x(j))
                j = j + 1
            end do

            ! Linear interpolation
            if (j > 1 .and. j <= n) then
                y_new(i) = y(j-1) + (y(j) - y(j-1)) * (x_new(i) - x(j-1)) / (x(j) - x(j-1))
            else if (j == 1) then
                y_new(i) = y(1)  ! Extrapolate
            else
                y_new(i) = y(n)  ! Extrapolate
            end if 
        end do
      end subroutine lin_interp

end module m_read_in_data
