program main
    implicit none

    integer :: i
    CHARACTER(128) :: file_input
    CHARACTER(128) :: file_output
    INTEGER :: particle_number
    DOUBLE PRECISION :: time, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq

    file_input = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory01.dat'
    file_output = '/home/satanka/Documents/fort/TPS-DAW/TPS-DAW-Artemyev/flux_calculation/myrank000/particle_trajectory01-58.dat'

    OPEN(100, file = file_input)
    OPEN(200, file = file_output)

    i = 0
    do
        i = i + 1
        read(100, *, end = 99) particle_number, time, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq
        if (particle_number == 58) then
            WRITE(200, '(I2.2, 7E15.7)') particle_number, time, z_position, u_para, u_perp, u_phi, energy, pitch_angle_eq
        end if
    end do

    99 CLOSE(100)
    CLOSE(200)

    
end program main
