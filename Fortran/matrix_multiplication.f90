! 矩阵相乘。
program matrix_multiplication
    implicit none
    integer,parameter :: l = 2, m = 3, n =4
    integer :: i, j, k
    real :: a(l, m)=(/1, 2, 3, 4, 5, 6/)
    real :: b(m, n) = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    real :: c(l, n)

    do i = 1, l
        do j = 1, n
            c(i, j) = 0.0
            do k = 1, m
                c(i, j) = c(i, j) + a(i, k) * b(k, j)
            end do
        end do
    end do


end program matrix_multiplication