module lib_fixedvector!定义了固定向量的类           !Fortran 2003
implicit none

private!声明所有私有项，除了操作符和固定向量本身
public::vector_fixed,assignment(=),operator(+),operator(-),operator(*),operator(/),polyxpoly

type::vector_fixed(n,precisio)!固定向量，用来表示线段
!        sequence!按顺序储存!绑定操作的类不能有sequence属性
        integer,len::n=3!空间维度                           !给数组长度赋值的变量要加len关键字
        integer,kind::precisio=8!精度                       !给精度赋值的变量要加kind
        real(kind=precisio),dimension(n)::location!固着位置
        real(kind=precisio)::direction(n)!方向，是一个经线段长度为模长的向量
contains
!        generic::parallel => parallel_sp,parallel_dp!自定义类型的函数重载!这个没啥用
!                procedure,pass::parallel_sp!用来判段是否平行，返回布尔值，平行为真
!                procedure,pass::parallel_dp!声明为pass，自动获得作为调用参数的变量
        generic::move => move_sp,move_dp!用于移动位置
                procedure,pass::move_sp
                procedure,pass::move_dp
        generic::mode => mode_sp,mode_dp!用于算模长
                procedure,pass::mode_sp
                procedure,pass::mode_dp
end type

interface assignment(=)!等于号的操作符拓展比较特殊，而且必须指向子程序而不能是函数
        module procedure array_to_vector_sp!输入数组第一列是位置，第二列是向量
        module procedure array_to_vector_dp
        module procedure vector_to_vector_sp
        module procedure vector_to_vector_dp
end interface

interface operator(+)!操作符扩展，加减法结果依第一个向量的位置
        module procedure vector_add_sp!单精度
        module procedure vector_add_dp!双精度
        module procedure vector_add_array_sp!单精度，固定向量加一自由向量
        module procedure vector_add_array_dp!双精度
        module procedure array_add_vector_sp!单精度
        module procedure array_add_vector_dp!双精度
end interface

interface operator(-)
        module procedure vector_minus_sp
        module procedure vector_minus_dp
        module procedure vector_minus_array_sp!单精度，固定向量加一自由向量
        module procedure vector_minus_array_dp!双精度
end interface

interface operator(*)!向量数乘
        module procedure vector_times_scalar_sp
        module procedure vector_times_scalar_dp
        module procedure scalar_times_vector_sp
        module procedure scalar_times_vector_dp
end interface
        
interface operator(/)
        module procedure vector_div_scalar_sp
        module procedure vector_div_scalar_dp
end interface

contains

pure subroutine array_to_vector_sp(vector,array)!An entity declared with the
!CLASS keyword shall be a dummy argument or have the ALLOCATABLE or POINTER
!attribut
        implicit none
        type(vector_fixed(*,precisio=4)),intent(out)::vector
        real(kind=4),intent(in)::array(:,:)

!        vector%n=size(array,dim=1)!type中的常量应该在传入子程序前就予以确定
        vector%location=array(:,1)       !!!!!!!!!!!!!error #5533: Feature found on this line is not yet supported in ifx 
        vector%direction=array(:,2)
        return
end subroutine

pure subroutine array_to_vector_dp(vector,array)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(out)::vector
        real(kind=8),intent(in)::array(:,:)

!        vector%n=size(array,dim=1)
        vector%location=array(:,1)
        vector%direction=array(:,2)
        return
end subroutine

pure subroutine vector_to_vector_sp(v1,v2)
        implicit none
        type(vector_fixed(*,precisio=4)),intent(out)::v1
        type(vector_fixed(*,precisio=4)),intent(in)::v2

        v1%location=v2%location
        v1%direction=v2%direction
        return
end subroutine

pure subroutine vector_to_vector_dp(v1,v2)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(out)::v1
        type(vector_fixed(*,precisio=8)),intent(in)::v2

        v1%location=v2%location
        v1%direction=v2%direction
        return
end subroutine

elemental function vector_add_sp(v1,v2) result(ans)!使用allocatable的参数就不能用elemental
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v1,v2
        type(vector_fixed(v1%n,precisio=4))::ans!那就可以这样分配

!        allocate(vector_fixed(v1%n,4)::ans)        !这样分配参数化派生数据类型
        ans%location=v1%location
        ans%direction=v1%direction+v2%direction
        return
end function

elemental function vector_add_dp(v1,v2) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v1,v2
        type(vector_fixed(v1%n,precisio=8))::ans

!        allocate(vector_fixed(v1%n,8)::ans)        !这样分配参数化派生数据类型
        ans%location=v1%location
        ans%direction=v1%direction+v2%direction
        return
end function

pure function vector_add_array_sp(v,array) result(ans)
        implicit none
        type(vector_fixed(:,precisio=4)),allocatable::ans
        type(vector_fixed(*,precisio=4)),intent(in)::v
        real(kind=4),intent(in)::array(:)

        allocate(vector_fixed(v%n,4)::ans)        !这样分配参数化派生数据类型
        ans=v
        ans%direction=ans%direction+array
        return
end function

pure function vector_add_array_dp(v,array) result(ans)
        implicit none
        type(vector_fixed(:,precisio=8)),allocatable::ans
        type(vector_fixed(*,precisio=8)),intent(in)::v
        real(kind=8),intent(in)::array(:)

        allocate(vector_fixed(v%n,8)::ans)        !这样分配参数化派生数据类型
        ans=v
        ans%direction=ans%direction+array
        return
end function

pure function array_add_vector_sp(array,v) result(ans)
        implicit none
        type(vector_fixed(:,precisio=4)),allocatable::ans
        type(vector_fixed(*,precisio=4)),intent(in)::v
        real(kind=4),intent(in)::array(:)

        allocate(vector_fixed(v%n,4)::ans)        !这样分配参数化派生数据类型
        ans=v
        ans%direction=ans%direction+array
        return
end function

pure function array_add_vector_dp(array,v) result(ans)
        implicit none
        type(vector_fixed(:,precisio=8)),allocatable::ans
        type(vector_fixed(*,precisio=8)),intent(in)::v
        real(kind=8),intent(in)::array(:)

        allocate(vector_fixed(v%n,8)::ans)        !这样分配参数化派生数据类型
        ans=v
        ans%direction=ans%direction+array
        return
end function

elemental function vector_minus_sp(v1,v2) result(ans)
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v1,v2
        type(vector_fixed(v1%n,precisio=4))::ans

!        ans%n=v1%n!n作为派生变量的可声明常量，不能这样赋值，要么声明时就给定，要么用allocate
        ans%location=v1%location
        ans%direction=v1%direction-v2%direction
        return
end function

elemental function vector_minus_dp(v1,v2) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v1,v2
        type(vector_fixed(v1%n,precisio=8))::ans

!        ans%n=v1%n
        ans%location=v1%location
        ans%direction=v1%direction-v2%direction
        return
end function

pure function vector_minus_array_sp(v,array) result(ans)
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v
        type(vector_fixed(v%n,precisio=4))::ans
        real(kind=4),intent(in)::array(:)

        ans=v
        ans%direction=ans%direction-array
        return
end function

pure function vector_minus_array_dp(v,array) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v
        type(vector_fixed(v%n,precisio=8))::ans
        real(kind=8),intent(in)::array(:)

        ans=v
        ans%direction=ans%direction-array
        return
end function

elemental function vector_times_scalar_sp(v,s) result(ans)!向量数乘
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=4))::ans
        real(kind=4),intent(in)::s!标量

        ans=v
        ans%direction=ans%direction*s
        return
end function
        
elemental function vector_times_scalar_dp(v,s) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=8))::ans
        real(kind=8),intent(in)::s

        ans=v
        ans%direction=ans%direction*s
        return
end function
        
elemental function scalar_times_vector_sp(s,v) result(ans)
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=4))::ans
        real(kind=4),intent(in)::s

        ans=v
        ans%direction=ans%direction*s
        return
end function
        
elemental function scalar_times_vector_dp(s,v) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=8))::ans
        real(kind=8),intent(in)::s

        ans=v
        ans%direction=ans%direction*s
        return
end function
        
elemental function vector_div_scalar_sp(v,s) result(ans)
        implicit none
        type(vector_fixed(*,precisio=4)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=4))::ans
        real(kind=4),intent(in)::s!标量

        ans=v
        ans%direction=ans%direction/s
        return
end function
        
elemental function vector_div_scalar_dp(v,s) result(ans)
        implicit none
        type(vector_fixed(*,precisio=8)),intent(in)::v!向量
        type(vector_fixed(v%n,precisio=8))::ans
        real(kind=8),intent(in)::s

        ans=v
        ans%direction=ans%direction/s
        return
end function

pure function move_sp(vector,array) result(ans)
        implicit none
        class(vector_fixed(*,4)),intent(in)::vector
        type(vector_fixed(vector%n,precisio=4))::ans
        real(kind=4),intent(in)::array(:)!移动的方向

        ans=vector
        ans%location=ans%location+array
        return
end function

pure function move_dp(vector,array) result(ans)
        implicit none
        class(vector_fixed(*,8)),intent(in)::vector
        type(vector_fixed(vector%n,precisio=8))::ans
        real(kind=8),intent(in)::array(:)!移动的方向

        ans=vector
        ans%location=ans%location+array
        return
end function

elemental function mode_sp(vector) result(ans)
        implicit none
        real(kind=4)::ans
        class(vector_fixed(*,precisio=4)),intent(in)::vector

        ans=sqrt(dot_product(vector%direction,vector%direction))
        return
end function

elemental function mode_dp(vector) result(ans)
        implicit none
        real(kind=8)::ans
        class(vector_fixed(*,precisio=8)),intent(in)::vector

        ans=sqrt(dot_product(vector%direction,vector%direction))
        return
end function

pure function polyxpoly(segment1,segment2)!判断二维平面上线段是否有交点，历史遗留代码
        implicit none
        logical::polyxpoly!相交为真
        class(vector_fixed(2,8)),intent(in)::segment1,segment2
        real(kind=8)::x_a(2),y_a(2),x_b(2),y_b(2)
        !用associate语句简化名称，方便使用
        associate(x_a=>segment1%location,y_a=>segment1%location+segment1%direction,&!a线段两端点坐标
                  x_b=>segment2%location,y_b=>segment2%location+segment2%direction) !b线段两端点坐标
        
        polyxpoly=.FALSE.
        if (dot_product(line(x_a,y_a),(/x_b,1.0_8/))*dot_product(line(x_a,y_a),(/y_b,1.0_8/)) .GT. 0) return
        if (dot_product(line(x_b,y_b),(/x_a,1.0_8/))*dot_product(line(x_b,y_b),(/y_a,1.0_8/)) .GT. 0) return
        polyxpoly=.TRUE.
        end associate
        return
        
        contains
                pure function line(m_1,m_2)!已知两点求直线方程一般式的三个参数
                implicit none
                real(kind=8)::line(3)
                real(kind=8),intent(in)::m_1(2),m_2(2)
                line(1)=m_1(2)-m_2(2)
                line(2)=-(m_1(1)-m_2(1))
                line(3)=-dot_product(line(1:2),m_1)
                return
                end function
                
end function

end module
