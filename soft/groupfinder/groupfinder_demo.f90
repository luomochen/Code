program groupfinder
    implicit none
    logical UNI, EQU
    integer i, j, k, l, m, n, o, p, NOL, NX, IC, MSEC, Ntrans, Ngroup, nonsym
    integer NT(3), NTYPE(8)
    real(kind = 8) IP1, IP2, Accuracy
    real(kind = 8) :: Basis(3, 3), Metric(3, 3), Direct(3, 8), X(3), RO(3, 3), Utest(3, 3)
    real(kind = 8), allocatable :: SS(:, :), SSS(:, :), XX(:, :, :), ROT(:, :, :), translantion(:, :)
    PARAMETER (MSEC = 3, Accuracy = 1D-4)
    
    !读取布拉菲晶格矢量并构成基矢量矩阵。此处需要补充IO模块用于读取POSCAR或者CONTCAR中的数据。
    !赋初值的时候括号内的先开始循环：(1, 1), (1, 2), (1, 3), (2, 1), ……，晶格矢量以行向量的形式构成矩阵。
    DATA((Basis(i, j), i=1, 3), j=1, 3)/5.5881266594, 0.0000000000, 0.0000000000, 0.0000000000, 5.5881266594, 0.0000000000, 0.0000000000, 0.0000000000, 5.5881266594/
    !则计算的时候，(x, y, z) = (s1, s2, s3)A = s1*i + s2*j + s3*k。
    !计算度量矩阵（metric matrix）. 可用于计算分数坐标下的运算，例如计算模方的时候：
    !+为取复共轭，x+x 变为 x+Mx 并且 M = B+B B是布拉菲晶格矩阵。
    !计算度量矩阵，Mij = (B+B)ij = (B+)ikBkj = BkiBkj
    !但是注意到，我们使用的行向量，实际上相当于取过一次转置（实矩阵取转置相当于取复共轭），要是求模方相当于左边
    !右边才是常用的列向量和(i, j, k)。应当为Mij = (BB+)ij = (B)ik(B+)kj = BikBjk
    do i = 1, 3
        do j = 1, 3
            Metric(i, j) = 0.0
            do k = 1, 3
                Metric(i, j) = Metric(i, j) + Basis(i, k) * Basis(j, k)
            end do
        end do
    end do
    !如果VASP文件使用的是笛卡尔坐标，此处应当进行一次笛卡尔坐标向分数坐标的转换。

    !读取VASP文件中的原子数。Number of atom in Lattice = NOL
    NOL = 8
    !读取VASP文件中的坐标，
    !allocate(Direct(3, NOL))
    !allocate(NTYPE(NOL))
    DATA((Direct(i, j), j=1, 3), i=1, 8)/0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5/
    DATA NTYPE /1, 1, 1, 1, 2, 2, 2, 2/
    !---------------------------------------------------
    !寻找所有可能的转动部分并且判断它们是否可以构成空间群操作
    !---------------------------------------------------

    !步骤一：寻找对任意晶格矢量保长的晶格矢量变换。

    do NX = 1, 3
        !计算三个晶格基矢的模方。
        !NX = 1 (1, 0, 0), 2 (0, 1, 0), 3 (0, 0, 1)。
        do i = 1, 3
            if ( i == NX ) then
                X(i) = 1.0
            else
                X(i) = 0.0
            end if
        end do
        !Inner Product1 = IP1
        IP1 = 0.0
        do i = 1, 3
            do j = 1, 3
                IP1 = IP1 + X(i) * Metric(i,j) * X(j)
            end do
        end do
        !循环指定的壳层范围。Maximum number of shells fully enclose the center = MSEC 
        !取值为-3-3，这会形成一个6*6*6的三维网格网格
        IC = 0
        do i = -MSEC, MSEC
            do j = -MSEC, MSEC
                do k = -MSEC, MSEC
                    X(1) = i
                    X(2) = j
                    X(3) = k
                    !Inner Product2 = IP2，计算网格上点对应矢量的模方。
                    IP2 = 0.0
                    do m = 1, 3
                        do l = 1, 3
                            IP2 = IP2 + X(l) * Metric(l,m) * X(m)
                        end do
                    end do
                    !此处的accuracy用于确定网格上点和晶格基矢之间是否等长的精度。
                    if ( abs(IP1-IP2) < Accuracy ) then
                        IC = IC + 1
                    !NT(3)用于记录分别和基矢a，b，c等长的矢量数目。
                        NT(NX) = IC
                    end if
                end do    
            end do
        end do
    end do
    !注意，由于用于记录网格中和基矢等长矢量的数组是可变数组，这里通过第一次循环确定了该数组的大小。
    !该种方法之后需要做改进，通过一次循环直接录入数据。
    !分配数组大小。
    allocate(XX(3, IC, 3))
    do NX = 1, 3
        do i = 1, 3
            if ( i == NX ) then
                X(i) = 1.0
            else
                X(i) = 0.0
            end if
        end do
        !计算模方。
        IP1 = 0.0
        do i = 1, 3
            do j = 1, 3
                IP1 = IP1 + X(i) * Metric(i,j) * X(j)
            end do
        end do
        !跑遍所有包围中心的壳层。
        IC = 0
        do i = -MSEC, MSEC
            do j = -MSEC, MSEC
                do k = -MSEC, MSEC
                    X(1) = i
                    X(2) = j
                    X(3) = k
                    IP2 = 0.0
                    do l = 1, 3
                        do m = 1, 3
                            IP2 = IP2 + X(l) * Metric(l,m) * X(m)
                        end do
                    end do
                    if ( abs(IP1-IP2) < Accuracy ) then
                        IC = IC + 1
                    !XX(NX, IC, 3)为三维数组，NX指定和哪条基矢（a，b，c）相等，IC为对相等矢量的编号。
                        XX(NX, IC, 1) = real(i)
                        XX(NX, IC, 2) = real(j)
                        XX(NX, IC, 3) = real(k)
                    end if
                end do    
            end do         
        end do
    end do
    write(*, *) NT, XX
    !步骤2 在NT(1)*NT(2)*NT(3)个候选当中选择幺正转换矩阵。
    Ntrans = 0
    !提取出和a等长的矢量。
    do i = 1, NT(1)
        !提取出和b等长的矢量。
        do j = 1, NT(2)
            !提取出和c等长的矢量。
            do k = 1, NT(3)
                do l = 1, 3
                    !注意，此处同样是按照行存储向量。ROtation = RO, 3*3*3 matrix。
                    RO(1, l) = XX(1, i, l)
                    RO(2, l) = XX(2, j, l)
                    RO(3, l) = XX(3, k, l)
                end do
                !此处缺少一个用于判断是否是酉矩阵的函数。
                !判断获得的旋转矩阵是否是酉矩阵。Utest = ROHH+RO+ = HH+, (Utest)ij = (RO)ikMkl(RO+)lj =
                !(RO)ikMkl(RO)jl
                do m = 1, 3
                    do n = 1, 3
                        Utest(m, n) = 0.D0
                        do o = 1, 3
                            do p = 1, 3
                                Utest(m, n) = Utest(m, n) + RO(m, o) * Metric(o, p) * RO(n, p)
                            end do
                        end do
                    end do
                end do
                !判断是否和度规矩阵相等。
                do m = 1, 3
                    do n = 1, 3
                        if ( abs(Utest(m, n) - Metric(m, n)) < Accuracy ) then
                            UNI = .True.
                        end if 
                    end do    
                end do
                !计数有多少酉转换矩阵。
                if ( UNI ) then
                    Ntrans = Ntrans + 1
                end if
            end do
        end do    
    end do
    !write(*, *) Ntrans
    !存储所有的幺正转换矩阵。
    allocate(ROT(Ntrans, 3, 3))
    Ntrans = 0
    do i = 1, NT(1)
        do j = 1, NT(2)
            do k = 1, NT(3)
                do l = 1, 3
                    RO(1, l) = XX(1, i, l)
                    RO(2, l) = XX(2, j, l)
                    RO(3, l) = XX(3, k, l)
                end do
                do m = 1, 3
                    do n = 1, 3
                        Utest(m, n) = 0.D0
                        do o = 1, 3
                            do p = 1, 3
                                Utest(m, n) = Utest(m, n) + RO(m, o) * Metric(o, p) * RO(n, p)
                            end do
                        end do
                    end do
                end do
                do m = 1, 3
                    do n = 1, 3
                        if ( abs(Utest(m, n) - Metric(m, n)) < Accuracy ) then
                            UNI = .True.
                        end if
                    end do
                end do 
                if ( UNI ) then
                    Ntrans = Ntrans + 1
                    do m = 1, 3
                        do n = 1, 3
                            ROT(Ntrans, m, n) = RO(m, n)
                        end do
                    end do
                end if    
            end do
        end do    
    end do
    
    !使用旋转矩阵操作原胞内的原子。
    allocate(SS(NOL, 3), SSS(NOL, 3))
    Ngroup = 0
    nonsym = 0
    do i = 1, Ntrans
        
        do j = 1, NOL
            do k = 1, 3
                SS(j, k) = 0.D0
                do l = 1, 3
                    !对原胞中的原子进行旋转操作。
                    SS(j, k) = SS(j, k) + Direct(j, l) * ROT(i, l, k)
                end do
            end do
        end do
        
        loop1: do j = 1, NOL
            
            do k = 1, NOL
                do l = 1, 3
                    SSS(k, l) = SS(k, l) + Direct(j, l) - SS(1, l)
                end do
            end do

            loop2: do k = 1, NOL
                loop3: do l = 1, NOL
                    if ( NTYPE(k) == NTYPE(l) .and. &
                    (abs(SSS(k, 1) - Direct(l, 1)) - nint(SSS(k, 1) - Direct(l, 1))) < Accuracy .and. &
                    (abs(SSS(k, 2) - Direct(l, 2)) - nint(SSS(k, 2) - Direct(l, 2))) < Accuracy .and. &
                    (abs(SSS(k, 3) - Direct(l, 3)) - nint(SSS(k, 3) - Direct(l, 3))) < Accuracy ) then
                        EQU = .True.
                        cycle loop3
                    else
                        EQU = .False.
                        exit loop3
                        exit loop2
                    end if
                end do loop3
            end do loop2

            if ( EQU ) then
                Ngroup = Ngroup + 1
            end if
        end do loop1
    end do
    write(*, *) Ngroup
end program groupfinder