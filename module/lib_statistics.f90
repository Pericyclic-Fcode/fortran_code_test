include "matrix.f90"

module lib_statistics
        use matrix
        implicit none
private
public::cov,varience,correlation,regression
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号

interface cov
        module procedure cov_matrix_dp!求协方差矩阵
        module procedure cov_matrix_sp
        module procedure cov_scalar_dp!求协方差值
        module procedure cov_scalar_sp
end interface

interface varience!方差
        module procedure varience_dp
        module procedure varience_sp
end interface

interface correlation!相关系数矩阵
        module procedure correlation_dp
        module procedure correlation_sp
end interface

interface regression!线性回归
        module procedure regression_dp
        module procedure regression_sp
end interface

contains

pure function cov_matrix_dp(a) result(cov)
        !求协方差矩阵，输入矩阵为(X_1，X_2.....X_column),row行column列，X是一个变量多次观测值组成的列向量
        implicit none
        integer::column
        real(kind=dp),intent(in)::a(:,:)
        real(kind=dp)::cov(size(a,dim=2),size(a,dim=2))
        integer::counter1,counter2
        real(kind=dp)::n
        
        column=size(a,dim=2)
        n=real(size(a,dim=1))!numpy库里这个函数是除的n-1,不影响做PCA
        
        !forall语句和隐式循环功能完全一样，用哪个看编译器，支持forall并行的用forall
        forall(counter1=1:column,counter2=1:column) cov(counter1,counter2)=dot_product&
        (a(:,counter1),a(:,counter2))/n-(sum(a(:,counter1))/n)*(sum(a(:,counter2))/n)
        
        !cov=reshape((/((dot_product(a(:,counter1),a(:,counter2))/n-(sum(a(:,counter1))/n)&
        !*(sum(a(:,counter2))/n),counter1=1,column,1),counter2=1,column,1)/),(/column,column/))
        !多重隐式循环其实生成的也是一维向量，只能在给多维数组赋初值时直接用，因为这时还没有确定数组形状，
        !还只是一块连续内存空间，可以按列优先的规则赋值进去，如real(kind=dp)::x(5,2)=(/1,2,3,4,5,4,5,3,2,4/)，
        !但是程序运行起来之后数组形状固定，不能再直接用，会报错，要用reshape函数
        
        return
end function

pure function cov_scalar_dp(x,y) result(cov)!求协方差值
        implicit none
        real(kind=dp),intent(in)::x(:),y(:)
        real(kind=dp)::cov
        real(kind=dp)::n
        n=real(size(x))
        cov=dot_product(x,y)/n-(sum(x)/n)*(sum(y)/n)
        return
end function

pure function varience_dp(x) result(D)!算方差
        implicit none
        real(kind=dp),intent(in)::x(:)
        real(kind=dp)::D
        real(kind=dp)::n
        n=real(size(x))
        D=dot_product(x,x)/n-(sum(x)/n)**2
        return
end function

pure function correlation_dp(a) result(R)!相关系数矩阵
        implicit none
        integer::column
        real(kind=dp),intent(in)::a(:,:)
        real(kind=dp)::R(size(a,dim=2),size(a,dim=2))
        integer::counter1,counter2
        
        column=size(a,dim=2)
        
        forall(counter1=1:column,counter2=1:column) R(counter1,counter2)=&
        sqrt(varience(a(:,counter1))*varience(a(:,counter2)))
        R=cov(a)/R
        return
end function

pure function regression_dp(x,y) result(theta)!x矩阵的列数是n，y是一个n元函数
        !线性回归,y=x_T*a, x=(x_1 x_2 x_3)_T
        implicit none
        real(kind=dp),intent(in)::x(:,:),y(:)!输入x矩阵的行数，是进行试验的次数
        real(kind=dp)::theta(size(x,dim=2)),A(size(x,dim=2),size(x,dim=2))
        
        A=matmul(transpose(x),x)
        A=inverse(A)!A的逆
        theta=matmul(matmul(A,transpose(x)),y)
        
        return
end function

pure function cov_matrix_sp(a) result(cov)
        !求协方差矩阵，输入矩阵为(X_1，X_2.....X_column),row行column列，X是一个变量多次观测值组成的列向量
        implicit none
        integer::column
        real(kind=sp),intent(in)::a(:,:)
        real(kind=sp)::cov(size(a,dim=2),size(a,dim=2))
        integer::counter1,counter2
        real(kind=sp)::n
        
        column=size(a,dim=2)
        n=real(size(a,dim=1))
        
        forall(counter1=1:column,counter2=1:column) cov(counter1,counter2)=dot_product&
        (a(:,counter1),a(:,counter2))/n-(sum(a(:,counter1))/n)*(sum(a(:,counter2))/n)
        return
end function

pure function cov_scalar_sp(x,y) result(cov)!求协方差值
        implicit none
        real(kind=sp),intent(in)::x(:),y(:)
        real(kind=sp)::cov
        real(kind=sp)::n
        n=real(size(x))
        cov=dot_product(x,y)/n-(sum(x)/n)*(sum(y)/n)
        return
end function

pure function varience_sp(x) result(D)!算方差
        implicit none
        real(kind=sp),intent(in)::x(:)
        real(kind=sp)::D
        real(kind=sp)::n
        n=real(size(x))
        D=dot_product(x,x)/n-(sum(x)/n)**2
        return
end function

pure function correlation_sp(a) result(R)!相关系数矩阵
        implicit none
        integer::column
        real(kind=sp),intent(in)::a(:,:)
        real(kind=sp)::R(size(a,dim=2),size(a,dim=2))
        integer::counter1,counter2
        
        column=size(a,dim=2)
        
        forall(counter1=1:column,counter2=1:column) R(counter1,counter2)=&
        sqrt(varience(a(:,counter1))*varience(a(:,counter2)))
        R=cov(a)/R
        return
end function

pure function regression_sp(x,y) result(theta)!x矩阵的列数是n，y是一个n元函数
        !线性回归,y=x_T*a, x=(x_1 x_2 x_3)_T
        implicit none
        real(kind=sp),intent(in)::x(:,:),y(:)!输入x矩阵的行数，是进行试验的次数
        real(kind=sp)::theta(size(x,dim=2)),A(size(x,dim=2),size(x,dim=2))
        
        A=matmul(transpose(x),x)
        A=inverse(A)!A的逆
        theta=matmul(matmul(A,transpose(x)),y)
        
        return
end function

end module
