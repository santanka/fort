program main
    implicit none

    integer :: i, j
    CHARACTER(128) :: file_input1, file_input2, file_output
    INTEGER :: particle_number
    DOUBLE PRECISION :: time_simulation, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq

    character(8)  :: date ! yyyymmdd
    character(10) :: time ! hhmmss.fff
    character(5)  :: zone ! shhmm
    integer :: value(8)   ! yyyy mm dd diff hh mm ss fff

   

    file_input1 = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory00.dat'
    file_input2 = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory01.dat'

    OPEN(100, file = file_input1)
    OPEN(200, file = file_input2)


    do i = 1, 102
        call date_and_time(date, time, zone, value)
        print *, i, value

        if (i <= 51) then
            REWIND(100)
            WRITE(file_output, '(A103, I3.3, A4)') '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/&
                &myrank000/particle_trajectory00-', i, '.dat'
            OPEN(50, file = file_output)
        else if (i >= 52) then
            REWIND(200)
            WRITE(file_output, '(A103, I3.3, A4)') '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/&
            &myrank000/particle_trajectory01-', i, '.dat'
            OPEN(50, file = file_output)
        end if

        do
            if (i <= 51) then
                read(100, *, end = 99) particle_number, time_simulation, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq
            else if (i >= 52) then
                read(200, *, end = 99) particle_number, time_simulation, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq
            end if

            if (particle_number == i) then
                WRITE(50, '(I2.2, 7E15.7)') particle_number, &
                    & time_simulation, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq
            end if
        end do

        99 CLOSE(50)
    end do

    CLOSE(100)
    CLOSE(200)

    
end program main
