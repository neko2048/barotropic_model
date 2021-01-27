program finproj
    !IMPLICIT NONE
    integer, parameter :: N = 100 ! # of grids
    real, dimension(N+1) :: x, y ! horizontal domain range
    real, dimension(N, N) :: vor0 ! vorticity fields
    real, dimension(N, N) :: phi0 ! streamline fields
    real, dimension(N, N) :: dvordt ! d (vorticity) / dt
    real :: Jpp, Jxp, Jpx ! Arawaka Jacobian parameter
    real, dimension(8) :: ps, qs ! Arakawa Jacobian parameter
    integer :: i, k, turn, step = 0 ! iteration variables
    integer, parameter :: v_iteration = 1000, timesteps = 9000 ! model parameter
    real, parameter :: max_cor = 1 * 10 ** (3), min_cor = - 1 * 10 ** (3) ! domain range
    real, parameter :: dt = 10. !time step interval
    character(len=100) :: sphi, svor ! output file name

    !##### give coordinates #####
    do i = 1, N+1
        x(i) = min_cor + (max_cor - min_cor) / N * (i - 1) ! km
        y(i) = min_cor + (max_cor - min_cor) / N * (i - 1) ! km
    end do
    d = (x(2) - x(1)) ! m

    !##### give initial vorticity field and streamline field #####
    do i = 1, N
        do k = 1, N
            vor0(k, i) = 10. ** (-3) * exp(- (2 * ((x(k)+x(k+1)) / 500) ** 2 + 0.5 * ((y(i)+y(i+1)) / 500) ** 2))
        end do
    end do

    ! output vor data as .dat
    write(svor, '(a,i4.4,a)') "vors/vors",step,".dat"
    !inquire(iolength=lrec) vor0
    open(10,file=svor, status='replace')
    do i = 1, N
            do k = 1, N
                    write(10, '(E14.4)') real(vor0(k, i))
            end do
    end do 
    !write(*, *) 'output ',svor,' as txt'
    close(10)


    do turn = 1, v_iteration
        ! ##### set initial streamline field #####
        do i = 4, N-4
            do k = 4, N-4
                phi0(k ,i) = (phi0(k+1, i) + phi0(k-1, i) + phi0(k, i+1) + phi0(k, i-1) - (d) ** 2. * vor0(k, i)) / 4.
            end do
        end do
    end do


    ! output phi data as .dat
    write(sphi, '(a,i4.4,a)') "phis/phis",step,".dat"
    !inquire(iolength=lrec) phi0
    open(20,file=sphi,status='replace')
    do i = 1, N
            do k = 1, N
                    write(20, '(E14.4)') real(phi0(k, i))
            end do
    end do
    write(*, *) 'output ',sphi,' as txt'
    close(20)

    do step = 1, timesteps
        write(*, *) step
        ! ##### Arakawa Jacobian to obtain dvor/dt#####
        do i = 2, N - 1
            do k = 2, N - 1
                call p_product(phi0, k, i, ps, N)
                call q_product(vor0, k, i, qs, N)
                Jpp = (ps(1) - ps(3)) * (qs(2) - qs(4)) - (ps(2) - ps(4)) * (qs(1) - qs(3))
                Jxp = (qs(2) * (ps(5) - ps(6)) - qs(4) * (ps(8) - ps(7)) - qs(1) * (ps(5) - ps(8)) + qs(3) * (ps(6) - ps(7)))
                Jpx = (-ps(2) * (qs(5) - qs(6)) + ps(4) * (qs(8) - qs(7)) + ps(1) * (qs(5) - qs(8)) - ps(3) * (qs(6) - qs(7)))
                dvordt(k, i) = - 1 / (12 * d ** 2) * (Jpp + Jxp + Jpx)
                !write(*, *) dvordt(k, i)
            end do
        end do

        ! ##### update
        vor0 = vor0 + dvordt * dt

        ! output vor data as .dat
        write(svor, '(a,i4.4,a)') "vors/vors",step,".dat"
        !inquire(iolength=lrec) vor0
        open(10,file=svor, status='replace')
        do i = 1, N
                do k = 1, N
                        write(10, '(E14.4)') real(vor0(k, i))
                end do
        end do 
        !write(*, *) 'output ',svor,' as txt'
        close(10)

            ! ##### set initial streamline, including boundary condition: u, v = 0 #####
            do turn = 1, v_iteration
                do i = 4, N-4
                    do k = 4, N-4
                        phi0(k, i) = (phi0(k+1, i) + phi0(k-1, i) + phi0(k, i+1) + phi0(k, i-1) - (d) ** 2. * vor0(k, i)) / 4.
                    end do
                end do
            end do

        ! output phi data as .dat
        write(sphi, '(a,i4.4,a)') "phis/phis",step,".dat"
        !inquire(iolength=lrec) phi0
        open(20,file=sphi,status='replace')
        do i = 1, N
                do k = 1, N
                        write(20, '(E14.4)') real(phi0(k, i))
                end do
        end do
        !write(*, *) 'output ',sphi,' as txt'
        close(20)
    end do
close(10)
close(20)
end program finproj

subroutine p_product(phi0, k, i, ps, N)
    ! return p values
    integer, intent(in) :: N
    real, dimension(N, N), intent(in) :: phi0
    real, dimension(8), intent(inout) :: ps

    ps(1) = phi0(k+1, i)
    ps(2) = phi0(k, i+1)
    ps(3) = phi0(k-1, i)
    ps(4) = phi0(k, i-1)
    ps(5) = phi0(k+1, i+1)
    ps(6) = phi0(k-1, i+1)
    ps(7) = phi0(k-1, i-1)
    ps(8) = phi0(k+1, i-1)
    return
end subroutine

subroutine q_product(vor0, k, i, qs, N)
    ! return q values
    integer, intent(in) :: N
    real, dimension(100, 100), intent(in) :: vor0
    real, dimension(8), intent(inout) :: qs
    qs(1) = vor0(k+1, i)
    qs(2) = vor0(k, i+1)
    qs(3) = vor0(k-1, i)
    qs(4) = vor0(k, i-1)
    qs(5) = vor0(k+1, i+1)
    qs(6) = vor0(k-1, i+1)
    qs(7) = vor0(k-1, i-1)
    qs(8) = vor0(k+1, i-1)
    return
end subroutine
