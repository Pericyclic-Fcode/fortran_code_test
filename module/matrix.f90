module matrix!勉强够用就行，真实用还是得看lapack
        implicit none
private
public::linspace,diag,trace,solve_local,gram_schmidt,gram_schmidt_local,inverse,solve,QR
public::eigen_value_vector,powermethod,householder,hessenberg,cond
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号
        real(kind=8),private,parameter::error=1E-7!误差限制
        real(kind=8),private,parameter::error_0=1E-15!用于判定是否等于0
        integer,private,parameter::max_limit=500!最大迭代次数

interface linspace!linspace(a,b,n)在a,b之间均匀取N个点，包含a,b
        module procedure linspace_dp
        module procedure linspace_sp
end interface

interface diag!diag(array)，返回对角阵
        module procedure diag_dp
        module procedure diag_sp
end interface

interface trace!trace(matrix),返回矩阵迹
        module procedure trace_dp
        module procedure trace_sp
end interface

interface swap!swap(a,b),用来交换元
        module procedure swap_dp
        module procedure swap_sp
end interface

interface solve_local!solve_local(A,b,column(A),column(b)),A是方阵，b可以是矩阵，可以是向量
        module procedure solve_local_dp!直接本地求解矩阵，列主元方法，速度快，尽量用于函数内调用
        module procedure solve_local_sp
end interface

interface gram_schmidt!施密特正交化,gram_schmidt(matrix),返回正交化后的矩阵
        module procedure gram_schmidt_dp
        module procedure gram_schmidt_sp
end interface

interface gram_schmidt_local!施密特正交化，本地操作，用于函数内调用，gram_schmidt_local(matrix)
        module procedure gram_schmidt_local_dp
        module procedure gram_schmidt_local_sp
end interface

interface inverse!逆矩阵,inverse(matrix)
        module procedure inverse_dp
        module procedure inverse_sp
end interface

interface solve!解线性方程组，是函数，非本地操作，速度不如欧拉哥给的代码,solve(matrix,array)
        module procedure solve_dp!不对退化阵负责，返回解
        module procedure solve_sp
end interface

interface QR!QR分解，第一个参数是待分解矩阵，第二个是Q输出，第三个是R输出，QR(matrix,Q,R)
        module procedure QR_dp
        module procedure QR_sp
end interface

interface eigen_value_vector!求特征值特征向量，第一个参数是待求矩阵，第二个是特征值组成的列向量
        module procedure eigen_value_vector_dp!第三个是对应特征列向量组成的矩阵
        module procedure eigen_value_vector_sp!主要是考虑到单纯求特征向量不求特征值的情况不多
end interface!传统QR分解求特征值，收敛性不太好,eigen_value_vector(matrix,array[,matrix])

interface powermethod!幂法求主特征值和对应特征向量，第一个是待求阵，第二个是特征值，第三个是对应特征向量
        module procedure powermethod_dp!powermethod(matrix,scalar,array)
        module procedure powermethod_sp
end interface

interface householder!求householder变换矩阵,householder(array) result(matrix)
        module procedure householder_dp!输入为映射操作平面的法向量，不用归一化
        module procedure householder_sp
end interface

interface hessenberg!本地操作，用householder方法化上海森堡型矩阵
        module procedure hessenberg_dp!解出来结果不唯一
        module procedure hessenberg_sp!效果一般
end interface

interface cond!求矩阵条件数，用的是2-范数，cond(A),A是要求条件数的矩阵
        module procedure cond_dp
        module procedure cond_sp
end interface

contains

pure function linspace_dp(a,b,n) result(ans)
        implicit none
        real(kind=dp),intent(in)::a,b
        integer,intent(in)::n
        real(kind=dp)::ans(n)
        integer::counter
        ans=(/(a+real(counter)*(b-a)/real(n-1),counter=0,n-1,1)/)
        return
end function linspace_dp

pure function linspace_sp(a,b,n) result(ans)
        implicit none
        real(kind=sp),intent(in)::a,b
        integer,intent(in)::n
        real(kind=sp)::ans(n)
        integer::counter
        ans=(/(a+real(counter)*(b-a)/real(n-1),counter=0,n-1,1)/)
        return
end function linspace_sp

pure function diag_dp(x) result(ans)!对角阵
        implicit none
        real(kind=dp),intent(in)::x(:)
        real(kind=dp)::ans(size(x),size(x))!用自动数组代替allocate
        integer::n,counter1,counter2
        n=size(x)
        ans=reshape((/(x(counter2),(0D0,counter1=1,n,1),counter2=1,n-1,1),x(n)/),(/n,n/))
        return
end function diag_dp

pure function diag_sp(x) result(ans)
        implicit none
        real(kind=sp),intent(in)::x(:)
        real(kind=sp)::ans(size(x),size(x))!用自动数组代替allocate
        integer::n,counter1,counter2
        n=size(x)
        ans=reshape((/(x(counter2),(0E0,counter1=1,n,1),counter2=1,n-1,1),x(n)/),(/n,n/))
        return
end function diag_sp

pure function trace_dp(A) result(ans)!矩阵迹
        implicit none
        real(kind=dp),intent(in)::A(:,:)!A是方阵
        real(kind=dp)::ans
        integer::n,counter

        n=size(A,dim=1)
        ans=sum((/(A(counter,counter),counter=1,n,1)/))
        return
end function trace_dp

pure function trace_sp(A) result(ans)
        implicit none
        real(kind=sp),intent(in)::A(:,:)!A是方阵
        real(kind=sp)::ans
        integer::n,counter

        n=size(A,dim=1)
        ans=sum((/(A(counter,counter),counter=1,n,1)/))
        return
end function trace_sp

pure elemental subroutine swap_dp(a,b)!用于交换元
    real(kind=dp),intent(inout)::a
    real(kind=dp),intent(inout)::b
    real(kind=dp)::tmp
    tmp=a;a=b;b=tmp
end subroutine swap_dp

pure subroutine solve_local_dp(a,b,n,m)!欧拉哥给的代码，直接本地操作，实测效率极高
    ! n*n * n*m = n*n
    integer,intent(in)::n
    integer,intent(in)::m
    real(kind=dp),intent(inout)::a(n,n)
    real(kind=dp),intent(inout)::b(n,m)
    real(kind=dp)::temp,temp2
    integer::i,j,k,major
    ! column major method
    do i=1,n-1
        ! 寻找第i列主元 find major position
        major=maxloc(abs(a(i:n,i)),dim=1)+i-1
        !数据交换,cahce miss
        call swap(a(i,i:n),a(major,i:n))
        call swap(b(i,:),b(major,:))
        if(abs(a(i,i))<epsilon(b(1,1))) return
        temp=1.d0/a(i,i)
        !计算
        do j=i+1,n
            a(i,j)=a(i,j)*temp!归一化
            a(i+1:n,j)=a(i+1:n,j)-a(i+1:n,i)*a(i,j)
        end do
        do j=1,m
            b(i,j)=b(i,j)*temp!归一化
            b(i+1:n,j)=b(i+1:n,j)-a(i+1:n,i)*b(i,j)
        end do
        !第i列
        a(i,i)=1.d0
        a(i+1:n,i)=0.d0
    end do
    temp=1.d0/a(n,n)
    a(n,n)=1.d0
    ! 回代
    do j=1,m
        !归一化
        b(n,j)=b(n,j)*temp
        do k=1,n-1
            temp2=b(n-k+1,j)
            b(1:n-k,j)=b(1:n-k,j)-temp2*a(1:n-k,n-k+1)
        end do
    end do
end subroutine solve_local_dp

pure elemental subroutine swap_sp(a,b)!用于交换元
    real(kind=sp),intent(inout)::a
    real(kind=sp),intent(inout)::b
    real(kind=sp)::tmp
    tmp=a;a=b;b=tmp
end subroutine swap_sp

pure subroutine solve_local_sp(a,b,n,m)!欧拉哥给的代码，直接本地操作，实测效率极高
    ! n*n * n*m = n*n
    integer,intent(in)::n
    integer,intent(in)::m
    real(kind=sp),intent(inout)::a(n,n)
    real(kind=sp),intent(inout)::b(n,m)
    real(kind=sp)::temp,temp2
    integer::i,j,k,major
    ! column major method
    do i=1,n-1
        ! 寻找第i列主元 find major position
        major=maxloc(abs(a(i:n,i)),dim=1)+i-1
        !数据交换,cahce miss
        call swap(a(i,i:n),a(major,i:n))
        call swap(b(i,:),b(major,:))
        if(abs(a(i,i))<epsilon(b(1,1))) return
        temp=1.d0/a(i,i)
        !计算
        do j=i+1,n
            a(i,j)=a(i,j)*temp!归一化
            a(i+1:n,j)=a(i+1:n,j)-a(i+1:n,i)*a(i,j)
        end do
        do j=1,m
            b(i,j)=b(i,j)*temp!归一化
            b(i+1:n,j)=b(i+1:n,j)-a(i+1:n,i)*b(i,j)
        end do
        !第i列
        a(i,i)=1.d0
        a(i+1:n,i)=0.d0
    end do
    temp=1.d0/a(n,n)
    a(n,n)=1.d0
    ! 回代
    do j=1,m
        !归一化
        b(n,j)=b(n,j)*temp
        do k=1,n-1
            temp2=b(n-k+1,j)
            b(1:n-k,j)=b(1:n-k,j)-temp2*a(1:n-k,n-k+1)
        end do
    end do
end subroutine solve_local_sp

pure subroutine gram_schmidt_local_dp(A) !对列满秩阵施密特正交化，不对非满秩阵的结果负责
        implicit none
        real(kind=dp),intent(inout)::A(:,:)
        integer::column,counter
        column=size(A,dim=2)

        do counter=2,column,1
                A(:,counter)=A(:,counter)-matmul(A(:,:counter-1),matmul(A(:,counter) &
                ,A(:,:counter-1))/sum(A(:,:counter-1)*A(:,:counter-1),dim=1))
        end do 
        forall(counter=1:column) A(:,counter)=A(:,counter)/sqrt(dot_product(A(:,counter),A(:,counter)))
        return
end subroutine gram_schmidt_local_dp

pure subroutine gram_schmidt_local_sp(A) !对列满秩阵施密特正交化，不对非满秩阵的结果负责
        implicit none
        real(kind=sp),intent(inout)::A(:,:)
        integer::column,counter
        column=size(A,dim=2)

        do counter=2,column,1
                A(:,counter)=A(:,counter)-matmul(A(:,:counter-1),matmul(A(:,counter) &
                ,A(:,:counter-1))/sum(A(:,:counter-1)*A(:,:counter-1),dim=1))
        end do 
        forall(counter=1:column) A(:,counter)=A(:,counter)/sqrt(dot_product(A(:,counter),A(:,counter)))
        return
end subroutine gram_schmidt_local_sp

pure function gram_schmidt_dp(A) result(gs)!对列满秩阵施密特正交化，不对非满秩阵的结果负责
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp)::gs(size(A,dim=1),size(A,dim=2))
        integer::column,row,counter
        column=size(A,dim=2)
        row=size(A,dim=1)

        gs(:,1)=A(:,1)
        do counter=2,column,1
                gs(:,counter)=A(:,counter)-matmul(gs(:,:counter-1),matmul(A(:,counter) &
                ,gs(:,:counter-1))/sum(gs(:,:counter-1)*gs(:,:counter-1),dim=1))
        end do 
        forall(counter=1:column) gs(:,counter)=gs(:,counter)/sqrt(dot_product(gs(:,counter),gs(:,counter)))
        return
end function gram_schmidt_dp

pure function gram_schmidt_sp(A) result(gs)!对列满秩阵施密特正交化，不对非满秩阵的结果负责
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp)::gs(size(A,dim=1),size(A,dim=2))
        integer::column,row,counter
        column=size(A,dim=2)
        row=size(A,dim=1)

        gs(:,1)=A(:,1)
        do counter=2,column,1
                gs(:,counter)=A(:,counter)-matmul(gs(:,:counter-1),matmul(A(:,counter) &
                ,gs(:,:counter-1))/sum(gs(:,:counter-1)*gs(:,:counter-1),dim=1))
        end do 
        forall(counter=1:column) gs(:,counter)=gs(:,counter)/sqrt(dot_product(gs(:,counter),gs(:,counter)))
        return
end function gram_schmidt_sp

pure function inverse_dp(A) result(ans)!求逆阵，不对退化的矩阵负责
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp)::ans(size(A,dim=2),size(A,dim=1)),media(size(A,dim=1),2*size(A,dim=2))
        integer::n,counter,i
        
        n=size(A,dim=1)
        
        media=0!构建增广阵
        media(:,:n)=A
        forall (counter=1:n) media(counter,counter+n)=1!右半截是单位阵
        
        do counter=1,n-1,1    !初等行变换化上三角阵
                do i=1,n-counter,1
                        if (abs(media(counter,counter)) .GE. error_0) exit         !cshift两侧开口可循环，eoshift两端闭口
                        media(counter:counter+i,:)=cshift(media(counter:counter+i,:),shift=1,dim=1)
                end do
                if (i .EQ. n-counter+1) cycle
                media(counter,:)=media(counter,:)/media(counter,counter)
                media(counter+1:,:)=media(counter+1:,:)-matmul(media(counter+1:,counter:counter),media(counter:counter,:))
        end do
        
        if (abs(media(n,n)) .GE. error_0) media(n,:)=media(n,:)/media(n,n)!最后一行归一化
        
        do counter=n,2,-1     !化行最简
                if (abs(media(counter,counter)) .LE. error_0) cycle
                media(:counter-1,:)=media(:counter-1,:)-matmul(media(:counter-1,counter:counter),media(counter:counter,:))
        end do
        
        ans=media(:,n+1:)
        return
end function inverse_dp

pure function inverse_sp(A) result(ans)
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp)::ans(size(A,dim=2),size(A,dim=1)),media(size(A,dim=1),2*size(A,dim=2))
        integer::n,counter,i
        
        n=size(A,dim=1)
        
        media=0!构建增广阵
        media(:,:n)=A
        forall (counter=1:n) media(counter,counter+n)=1!右半截是单位阵
        
        do counter=1,n-1,1    !初等行变换化上三角阵
                do i=1,n-counter,1
                        if (abs(media(counter,counter)) .GE. error_0) exit         !cshift两侧开口可循环，eoshift两端闭口
                        media(counter:counter+i,:)=cshift(media(counter:counter+i,:),shift=1,dim=1)
                end do
                if (i .EQ. n-counter+1) cycle
                media(counter,:)=media(counter,:)/media(counter,counter)
                media(counter+1:,:)=media(counter+1:,:)-matmul(media(counter+1:,counter:counter),media(counter:counter,:))
        end do
        
        if (abs(media(n,n)) .GE. error_0) media(n,:)=media(n,:)/media(n,n)!最后一行归一化
        
        do counter=n,2,-1     !化行最简
                if (abs(media(counter,counter)) .LE. error_0) cycle
                media(:counter-1,:)=media(:counter-1,:)-matmul(media(:counter-1,counter:counter),media(counter:counter,:))
        end do
        
        ans=media(:,n+1:)
        return
end function inverse_sp

pure function solve_dp(A,b) result(ans)!用来求解方程组，自行避免矩阵退化，实测列变换快
        implicit none                       !个人非常满意的一个代码，很灵性，很有想法
        real(kind=dp),intent(in)::A(:,:),b(:)
        real(kind=dp)::ans(size(b))
        real(kind=dp)::M(size(A,dim=2)+ 1,size(A,dim=1))
        logical::mask(size(A,dim=1))
        integer::row,column,counter,i,n,order(size(A,dim=1))

        row    = size(A,dim=1)
        column = size(A,dim=2)
        mask   = .TRUE.!掩码数组，用来记录是否被提取过主元
        order  = 0!记录被提取主元的位置
        M(:row,:) = transpose(A)!构造增广阵，耗时，用隐循环反而慢，这里是决速步
        M(row+1,:)= b

        do counter= 1,row,1
                i=maxloc(abs(M(counter,:)),dim=1,mask=mask)!maxloc不加dim返回的是数组
                mask(i)=.FALSE.
                order(counter)=i
                if (abs(M(counter,i)) .LE. epsilon(b(1))) cycle!防止算出一堆NAN，不好看
                M(:,i)=M(:,i)/M(counter,i)!归一化
                forall(n=1:row,mask(n)) M(:,n)=M(:,n)-M(counter,n)*M(:,i)
        end do

        mask= .not.mask
        do counter=count(mask), 1, -1
            mask(order(counter))=.FALSE.    !回代
            forall(n=1:row,mask(n)) M(:,n)=M(:,n)-M(counter,n)*M(:,order(counter))
        end do
        ans=matmul(M(:column,:),M(column+1,:))
        return
end function solve_dp

pure function solve_sp(A,b) result(ans)!用来求解方程组，自行避免矩阵退化，实测列变换快
        implicit none                       !个人非常满意的一个代码，很灵性，很有想法
        real(kind=sp),intent(in)::A(:,:),b(:)
        real(kind=sp)::ans(size(b))
        real(kind=sp)::M(size(A,dim=2)+ 1,size(A,dim=1))
        logical::mask(size(A,dim=1))
        integer::row,column,counter,i,n,order(size(A,dim=1))

        row    = size(A,dim=1)
        column = size(A,dim=2)
        mask   = .TRUE.!掩码数组，用来记录是否被提取过主元
        order  = 0!记录被提取主元的位置
        M(:row,:) = transpose(A)!构造增广阵，耗时，用隐循环反而慢，这里是决速步
        M(row+1,:)= b

        do counter= 1,row,1
                i=maxloc(abs(M(counter,:)),dim=1,mask=mask)!maxloc不加dim返回的是数组
                mask(i)=.FALSE.
                order(counter)=i
                if (abs(M(counter,i)) .LE. epsilon(b(1))) cycle!防止算出一堆NAN，不好看
                M(:,i)=M(:,i)/M(counter,i)!归一化
                forall(n=1:row,mask(n)) M(:,n)=M(:,n)-M(counter,n)*M(:,i)
        end do

        mask= .not.mask
        do counter=count(mask), 1, -1
            mask(order(counter))=.FALSE.    !回代
            forall(n=1:row,mask(n)) M(:,n)=M(:,n)-M(counter,n)*M(:,order(counter))
        end do
        ans=matmul(M(:column,:),M(column+1,:))
        return
end function solve_sp

pure subroutine QR_dp(A,Q,R)
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp),intent(out)::Q(:,:),R(:,:)
        Q=gram_schmidt(A)
        R=matmul(transpose(Q),A)
        return
end subroutine QR_dp

pure subroutine QR_sp(A,Q,R)
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp),intent(out)::Q(:,:),R(:,:)
        Q=gram_schmidt(A)!格拉姆施密特正交化求Q
        R=matmul(transpose(Q),A)!回代求R
        return
end subroutine QR_sp

pure subroutine eigen_value_vector_dp(A,eigen_value,eigen_vector)!求特征值，传统QR分解求特征值，收敛性不太好
        implicit none
        real(kind=dp),intent(in)::A(:,:)!A是满秩方阵
        real(kind=dp),intent(out)::eigen_value(:)
        real(kind=dp),optional,intent(out)::eigen_vector(:,:)!不输入这项就不求特征向量
        real(kind=dp)::B(size(A,dim=1),size(A,dim=2)),Q(size(A,dim=1),size(A,dim=2)),R(size(A,dim=1),size(A,dim=2))
        integer::counter
        B=A
        if(present(eigen_vector)) then
                eigen_vector=diag((/(1E0,counter=1,size(A,dim=1),1)/))
                do counter=1,max_limit,1!设置一个最大迭代次数，否则用if加goto也行
                        call QR(B,Q,R)
                        B=matmul(R,Q)
                        eigen_vector=matmul(eigen_vector,Q)
                        if(abs(sum(B)-trace(B)) .LE. error) exit
                end do
        else
                do counter=1,max_limit,1!设置一个最大迭代次数，否则用if加goto也行
                        call QR(B,Q,R)
                        B=matmul(R,Q)
                        if(abs(sum(B)-trace(B)) .LE. error) exit
                end do
        end if
        eigen_value=(/(B(counter,counter),counter=1,size(A,dim=1),1)/)
        return
end subroutine eigen_value_vector_dp

pure subroutine eigen_value_vector_sp(A,eigen_value,eigen_vector)!求特征值，传统QR分解求特征值，收敛性不太好
        implicit none
        real(kind=sp),intent(in)::A(:,:)!A是满秩方阵
        real(kind=sp),intent(out)::eigen_value(:)
        real(kind=sp),optional,intent(out)::eigen_vector(:,:)!不输入这项就不求特征向量
        real(kind=sp)::B(size(A,dim=1),size(A,dim=2)),Q(size(A,dim=1),size(A,dim=2)),R(size(A,dim=1),size(A,dim=2))
        integer::counter
        B=A
        if(present(eigen_vector)) then
                eigen_vector=diag((/(1E0,counter=1,size(A,dim=1),1)/))
                do counter=1,max_limit,1!设置一个最大迭代次数，否则用if加goto也行
                        call QR(B,Q,R)
                        B=matmul(R,Q)
                        eigen_vector=matmul(eigen_vector,Q)
                        if(abs(sum(B)-trace(B)) .LE. error) exit
                end do
        else
                do counter=1,max_limit,1!设置一个最大迭代次数，否则用if加goto也行
                        call QR(B,Q,R)
                        B=matmul(R,Q)
                        if(abs(sum(B)-trace(B)) .LE. error) exit
                end do
        end if
        eigen_value=(/(B(counter,counter),counter=1,size(A,dim=1),1)/)
        return
end subroutine eigen_value_vector_sp

pure subroutine powermethod_dp(A,lambda,x)!幂法求主特征值和对应特征向量
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp),intent(out)::lambda
        real(kind=dp),intent(inout)::x(:)
        integer::counter

        if(all(abs(x) .LT. error_0)) x(1)=1E0!没初始化给个初值
        lambda=0E0
        do counter=1,max_limit,1
                x=matmul(A,x)
                if(abs(maxval(x)-lambda) .LT. error) exit!收敛了，跳出
                lambda=maxval(x)
                x=x/lambda
        end do
        x=x/sqrt(dot_product(x,x))!归一化
        return
end subroutine powermethod_dp

pure subroutine powermethod_sp(A,lambda,x)!幂法求主特征值和对应特征向量
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp),intent(out)::lambda
        real(kind=sp),intent(inout)::x(:)
        integer::counter

        if(all(abs(x) .LT. error_0)) x(1)=1E0!没初始化给个初值
        lambda=0E0
        do counter=1,max_limit,1
                x=matmul(A,x)
                if(abs(maxval(x)-lambda) .LT. error) exit!收敛了，跳出
                lambda=maxval(x)
                x=x/lambda
        end do
        x=x/sqrt(dot_product(x,x))
        return
end subroutine powermethod_sp

pure function det_dp(A) result(ans)!求行列式!正负号没法确定
        implicit none
        real(kind=dp)::ans
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp)::Q(size(A,dim=1),size(A,dim=2)),R(size(A,dim=1),size(A,dim=2))
        integer::counter

        call QR(A,Q,R)
        ans=product((/(R(counter,counter),counter=1,size(A,dim=1),1)/))
        return
end function det_dp

pure function det_sp(A) result(ans)!求行列式!正负号没法确定
        implicit none
        real(kind=sp)::ans
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp)::Q(size(A,dim=1),size(A,dim=2)),R(size(A,dim=1),size(A,dim=2))
        integer::counter

        call QR(A,Q,R)
        ans=product((/(R(counter,counter),counter=1,size(A,dim=1),1)/))
        return
end function det_sp

pure function householder_dp(x) result(ans)!求householder变换矩阵
        implicit none
        real(kind=dp),intent(in)::x(:)
        real(kind=dp)::ans(size(x),size(x))
        integer::counter
        ans=reshape((/(x(counter)*x,counter=1,size(x),1)/),(/size(x),size(x)/))
        ans=diag((/(1E0,counter=1,size(x),1)/))-2*ans/dot_product(x,x)
        return
end function householder_dp
        
pure function householder_sp(x) result(ans)
        implicit none
        real(kind=sp),intent(in)::x(:)
        real(kind=sp)::ans(size(x),size(x))
        integer::counter
        ans=reshape((/(x(counter)*x,counter=1,size(x),1)/),(/size(x),size(x)/))
        ans=diag((/(1E0,counter=1,size(x),1)/))-2*ans/dot_product(x,x)
        return
end function householder_sp
        
pure subroutine hessenberg_dp(A)
        implicit none
        real(kind=dp),intent(inout)::A(:,:)
        integer::n,counter
        real(kind=dp)::tmp(size(A,dim=1))

        n=size(A,dim=2)
        do counter=1,n-1,1
                tmp=A(:,counter)
                tmp(counter+1)=tmp(counter+1)+sign(1E0,A(counter+1,counter))*sqrt(sum(A(counter+1:,counter)**2))
                tmp(:counter)=0E0
                A=matmul(householder(tmp),A)
        end do
        return
end subroutine hessenberg_dp

pure subroutine hessenberg_sp(A)
        implicit none
        real(kind=sp),intent(inout)::A(:,:)
        integer::n,counter
        real(kind=sp)::tmp(size(A,dim=1))

        n=size(A,dim=2)
        do counter=1,n-1,1
                tmp=A(:,counter)
                tmp(counter+1)=tmp(counter+1)+sign(1E0,A(counter+1,counter))*sqrt(sum(A(counter+1:,counter)**2))
                tmp(:counter)=0E0
                A=matmul(householder(tmp),A)
        end do
        return
end subroutine hessenberg_sp

pure function cond_dp(A) result(ans)
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        real(kind=dp)::lambda(size(A,dim=2)),ans
        call eigen_value_vector(matmul(transpose(A),A),lambda)
        ans=sqrt(maxval(abs(lambda))/minval(abs(lambda)))
        return
end function cond_dp

pure function cond_sp(A) result(ans)
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        real(kind=sp)::lambda(size(A,dim=2)),ans
        call eigen_value_vector(matmul(transpose(A),A),lambda)
        ans=sqrt(maxval(abs(lambda))/minval(abs(lambda)))
        return
end function cond_sp

end module
