include "lib_sort.f90"

module lib_interpolate
        use lib_sort
        implicit none
private
public::initialize,lagrange,spline_initialize,spline!算拉格朗日插值和样条插值的函数默认逐元
        !integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号，想改精度在这里调
        real(kind=dp),private,allocatable,save::x(:),y(:)
        real(kind=dp),private,allocatable,save::m(:),h(:)!用于三次样条插值
        real(kind=dp),private,save::cycle_T!记录周期，用于周期边界条件的三次样条插值

interface spline_initialize!重载
        module procedure spline_cubic_initialize_natural!自然边界，无输入,不对外插负责
        module procedure spline_cubic_initialize_clamped!夹持边界，输入两个夹持条件
        module procedure spline_cubic_initialize_cycle!周期边界，输入周期，调用之后函数自动对输入变量加减周期
end interface spline_initialize

contains

subroutine initialize(a,b)!初始化，调用插值函数之前先用这个传入数据
        implicit none
        real(kind=dp),intent(in)::a(:),b(:)!a是x方向数据，b是y方向的
        integer::n
        n=size(a)-1
        
        if (allocated(x)) deallocate(x,y)
        allocate(x(0:n),y(0:n))!输入n+1组数据，返回n次多项式
        x=a
        y=b
        
        return
end subroutine initialize

elemental function lagrange(x_0) result(ans)!拉格朗日插值
        implicit none
        real(kind=dp),intent(in)::x_0
        integer::limit,counter
        real(kind=dp)::ans
        real(kind=dp)::l_vector(0:size(x)-1)!使用自动数组，比用allocate更快
        
        limit=size(x)-1
        l_vector=(/(l(x_0,counter),counter=0,limit,1)/)
        ans=dot_product(l_vector,y)
        
        return
        contains!内部过程
                pure function l(x_0,counter)!插值基函数
                implicit none
                real(kind=dp),intent(in)::x_0
                real(kind=dp)::l,media(limit)!limit可以直接用
                integer,intent(in)::counter
                logical::mask(0:limit)
        
                mask=.TRUE.!内部过程定义的变量不能在声明时初始化
                mask(counter)=.FALSE.!掩码数组，保证i/=j
                media=pack(x,mask)!用media(:counter)=x(:counter-1)不好，容易有下标越界的问题
                l=product((x_0-media)/(x(counter)-media))
                return
                end function l
end function lagrange

subroutine spline_cubic_initialize_natural()!自由边界的三次样条插值，用之前先调用这个函数初始化
!算法来源zhuanlan.zhihu.com/p/62860859
        implicit none
        real(kind=dp)::A(size(x),size(x))
        real(kind=dp),parameter::one=1.0,zero=0.0
        integer::limit,counter,counter1,arr(size(x))!block是2008的语法，还不支持
        
        limit=size(x)-1
        if (allocated(m)) deallocate(m,h)
        allocate(h(limit),m(0:limit))
100     h=x(1:limit)-x(0:limit-1)
        
        if (any(h .LE. zero)) then!如果不按顺序，进行排序
!        block
!                integer::arr(size(x))
                arr=get_order(x)
                x=permutation(x,arr)
                y=permutation(y,arr)
                goto 100
!        end block
        end if
                
        !构建系数矩阵
        A(:,:limit+1)=reshape((/one,(zero,counter=1,limit,1),((h(counter1),&
        2*(h(counter1)+h(counter1+1)),h(counter1+1),(zero,counter=1,limit-1,1)),counter1=1,limit-1,1),zero,one/),(/limit+1,limit+1/))
        call swap(A(1,2),A(2,1))
        call swap(A(limit+1,limit),A(limit,limit+1))
        !等号右侧
        m=6*(/zero,((y(counter+2)-y(counter+1))/h(counter+2)-(y(counter+1)-y(counter+0))/h(counter+1),counter=0,limit-2,1),zero/)
        !调用欧拉哥给的函数，那个速度比较快
        call solve(A,m,limit+1,1)
        return
end subroutine spline_cubic_initialize_natural

subroutine spline_cubic_initialize_clamped(bound_A,bound_B)!夹持边界的三次样条插值
!算法来源zhuanlan.zhihu.com/p/62860859
        implicit none
        real(kind=dp),intent(in)::bound_A,bound_B
        real(kind=dp)::A(size(x),size(x))
        real(kind=dp),parameter::one=1.0,zero=0.0
        integer::limit,counter,counter1,arr(size(x))
        
        limit=size(x)-1
        if (allocated(m)) deallocate(m,h)
        allocate(h(limit),m(0:limit))
100     h=x(1:limit)-x(0:limit-1)
        
        if (any(h .LE. zero)) then!如果不按顺序，进行排序
                arr=get_order(x)
                x=permutation(x,arr)
                y=permutation(y,arr)
                goto 100
        end if
                
        !构建系数矩阵
        A(:,:limit+1)=reshape((/2*h(1),h(1),(zero,counter=2,limit,1),((h(counter1),&
        2*(h(counter1)+h(counter1+1)),h(counter1+1),(zero,counter=1,limit-1,1)),counter1=1,limit-1,1),h(limit),2*h(limit)/),(/limit+1,limit+1/))
        !等号右侧
        m=6*(/(y(1)-y(0))/h(1)-bound_A,((y(counter+2)-y(counter+1))/h(counter+2)-(y(counter+1)-y(counter+0))/&
        h(counter+1),counter=0,limit-2,1),bound_B-(y(limit)-y(limit-1))/h(limit)/)
        !调用欧拉哥给的函数，那个速度比较快
        call solve(A,m,limit+1,1)
        return
end subroutine spline_cubic_initialize_clamped

subroutine spline_cubic_initialize_cycle(T)!周期边界的三次样条插值
        implicit none
        real(kind=dp),intent(in)::T!输入周期长度
        real(kind=dp)::A(size(x)+1,size(x)+1)
        real(kind=dp),parameter::one=1.0,zero=0.0
        integer::limit,counter,counter1,arr(size(x))
        
        cycle_T=T
        limit=size(x)-1

        if (allocated(m)) deallocate(m,h)
        allocate(h(limit+1),m(0:limit+1))

100     h(:limit)=x(1:limit)-x(0:limit-1)
        h(limit+1)=T+x(0)-x(limit)
        if (any(h .LE. zero)) then!如果不按顺序，进行排序
                arr=get_order(x)
                x=permutation(x,arr)
                y=permutation(y,arr)
                goto 100
        end if
                
        !构建系数矩阵
        A=reshape((/2*(h(1)+h(limit+1)),h(1),(zero,counter=2,limit+1,1),((h(counter1),&
        2*(h(counter1)+h(counter1+1)),h(counter1+1),(zero,counter=1,limit,1)),counter1=1,limit,1),h(limit+1),2*(h(limit+1)+h(1))/),(/limit+2,limit+2/))
        !等号右侧
        m=6*(/(y(1)-y(0))/h(1)-(y(0)-y(limit))/h(limit+1),((y(counter+2)-y(counter+1))/h(counter+2)-(y(counter+1)-y(counter+0))/&
        h(counter+1),counter=0,limit-2,1),(y(0)-y(limit))/h(limit+1)-(y(limit)-y(limit-1))/h(limit),(y(1)-y(0))/h(1)-(y(0)-y(limit))/h(limit+1)/)
        !调用欧拉哥给的函数，那个速度比较快
        call solve(A,m,limit+2,1)
        return
end subroutine spline_cubic_initialize_cycle

elemental function spline(x_0) result(ans)!三次样条插值，只能内插
        implicit none
        real(kind=dp),intent(in)::x_0
        real(kind=dp)::a,b,c,d,ans,tmp
        integer::counter,limit
        
        limit=size(x)-1
        tmp=x_0
500     do counter=0,limit-1,1
                if ((tmp .GE. x(counter)) .AND. (tmp .LE. x(counter+1))) goto 100
        end do
        
        if (size(m) .EQ. limit+2) then!是在周期边界条件下
                do
                        tmp=tmp+sign(1E0,x(0)-tmp)*cycle_T
                        if (tmp .GT. x(limit) .AND. tmp .LE. cycle_T+x(0)) goto 200
                        if ((tmp .GE. x(0)) .AND. (tmp .LE. x(limit))) goto 500
                end do
        else
                ans=0!非周期边界条件，落在范围之外，随便返回个数
                return
        end if
200     a=y(limit)!一周期最后一点到下周期第一点之间的区域
        b=(y(0)-y(limit))/h(limit+1)-h(limit+1)*m(limit)/2-h(limit+1)*(m(limit+1)-m(limit))/6
        c=m(limit)/2
        d=(m(limit+1)-m(limit))/(6*h(limit+1))
        ans=a+b*(tmp-x(limit))+c*(tmp-x(limit))**2+d*(tmp-x(limit))**3
        return

100     a=y(counter)!一般情况
        b=(y(counter+1)-y(counter))/h(counter+1)-h(counter+1)*m(counter)/2-h(counter+1)*(m(counter+1)-m(counter))/6
        c=m(counter)/2
        d=(m(counter+1)-m(counter))/(6*h(counter+1))
        ans=a+b*(tmp-x(counter))+c*(tmp-x(counter))**2+d*(tmp-x(counter))**3
        return
end function spline

pure elemental subroutine swap(a,b)!用于交换元
    real(kind=dp),intent(inout)::a
    real(kind=dp),intent(inout)::b
    real(kind=dp)::tmp
    tmp=a;a=b;b=tmp
end subroutine swap

pure subroutine solve(a,b,n,m)!欧拉哥给的代码，直接本地操作，实测效率极高
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
end subroutine solve

end module
