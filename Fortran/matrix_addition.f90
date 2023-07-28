! 矩阵相加。
program matrix_test
    implicit none
    integer :: matrixA(2, 2)
    integer :: matrixB(2, 2)
    integer :: matrixC(2, 2)
    integer i, j
    DATA((matrixA(i, j),i=1, 2),j=1, 2)/1, 0, 0, 1/
    DATA((matrixB(i, j),i=1, 2),j=1, 2)/1, 0, 0, 1/

    matrixC = matrixA + matrixB
    do i = 1, 2
        do j = 1, 2
            write(*, *) matrixC(i, j) 
        end do
    end do

    stop
end program matrix_test
