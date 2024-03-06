module matrix_operation
    implicit none
contains
    subroutine MATINV(A,B,C)
        real(kind=8) :: A(3,3)
        real(kind=8) B(3,3), D(3,3)
        real(kind=8) C
        D(1,1)=A(2,2)*A(3,3)-A(2,3)*A(3,2)
        D(2,2)=A(3,3)*A(1,1)-A(3,1)*A(1,3)
        D(3,3)=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        D(1,2)=A(2,3)*A(3,1)-A(2,1)*A(3,3)
        D(2,3)=A(3,1)*A(1,2)-A(3,2)*A(1,1)
        D(3,1)=A(1,2)*A(2,3)-A(1,3)*A(2,2)
        D(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
        D(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
        D(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
        C=A(1,1)*D(1,1)+A(1,2)*D(1,2)+A(1,3)*D(1,3)
        B(1,1)=D(1,1)/C
        B(2,2)=D(2,2)/C
        B(3,3)=D(3,3)/C
        B(1,2)=D(2,1)/C
        B(2,3)=D(3,2)/C
        B(3,1)=D(1,3)/C
        B(2,1)=D(1,2)/C
        B(3,2)=D(2,3)/C
        B(1,3)=D(3,1)/C
        return
    end subroutine MATINV
    ! SUBROUTINE TO EVALUATE THE INNER PRODUCT OF
    ! TWO VECTORS UNDER ARBITRARY METRIC.
    subroutine INPRODUCT(X, Y, HH, IP)
        integer I,J
        real(kind=8) :: X(3), Y(3), HH(3,3)
        real(kind=8) :: IP
        do I = 1,3
            do J = 1,3
            IP = IP + X(I) * HH(I,J) * Y(J)
            end do
        end do
        return
    end subroutine INPRODUCT
end module matrix_operation