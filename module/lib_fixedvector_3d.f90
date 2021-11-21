module lib_fixedvector_3d
implicit none

private!声明所有私有项，除了操作符和固定向量本身
public::vector_fixed_3d,assignment(=),operator(+),operator(-),operator(*),operator(/)

integer,parameter::kind=4       !在这里手动配置精度

type::vector_fixed_3d           !固定向量
!        sequence               !有类型绑定过程不能有sequence声明
        real(kind=kind)::location(3)
        real(kind=kind)::direction(3)
contains
        procedure,pass::move       !移动
        procedure,pass::mode       !求模
        procedure,pass::plane      !平面定向
        procedure,pass::normalize  !归一化
end type

interface assignment(=)!赋值符号
        module procedure array_to_vector!第一列是位置，第二列是方向
        module procedure vector_to_vector
end interface

interface operator(+)
        module procedure vector_add!位置依第一个向量
        module procedure vector_add_array!数组加到方向上去
        module procedure array_add_vector        
end interface

interface operator(-)
        module procedure vector_minus
        module procedure vector_minus_array        
end interface

interface operator(*)
        module procedure vector_times_scalar!数乘
        module procedure scalar_times_vector
end interface
        
interface operator(/)
        module procedure vector_div_scalar
end interface

contains

pure subroutine array_to_vector(vector,array)
        implicit none
        type(vector_fixed_3d),intent(out)::vector
        real(kind=kind),intent(in)::array(3,2)

        vector%location=array(:,1)
        vector%direction=array(:,2)
        return
end subroutine

elemental subroutine vector_to_vector(v2,v1)
        implicit none
        type(vector_fixed_3d),intent(out)::v2
        type(vector_fixed_3d),intent(in)::v1

        v2%location=v1%location
        v2%direction=v1%direction
        return
end subroutine

elemental function vector_add(v1,v2) result(ans)
        implicit none
        type(vector_fixed_3d),intent(in)::v1,v2
        type(vector_fixed_3d)::ans

        ans%location=v1%location
        ans%direction=v1%direction+v2%direction
        return
end function

pure function vector_add_array(v,array) result(ans)
        implicit none
        type(vector_fixed_3d)::ans
        type(vector_fixed_3d),intent(in)::v
        real(kind=kind),intent(in)::array(3)

        ans=v
        ans%direction=ans%direction+array
        return
end function

pure function array_add_vector(array,v) result(ans)
        implicit none
        type(vector_fixed_3d)::ans
        type(vector_fixed_3d),intent(in)::v
        real(kind=kind),intent(in)::array(3)

        ans=v
        ans%direction=ans%direction+array
        return
end function

elemental function vector_minus(v1,v2) result(ans)
        implicit none
        type(vector_fixed_3d),intent(in)::v1,v2
        type(vector_fixed_3d)::ans

        ans%location=v1%location
        ans%direction=v1%direction-v2%direction
        return
end function

pure function vector_minus_array(v,array) result(ans)
        implicit none
        type(vector_fixed_3d)::ans
        type(vector_fixed_3d),intent(in)::v
        real(kind=kind),intent(in)::array(3)

        ans=v
        ans%direction=ans%direction-array
        return
end function
 
elemental function vector_times_scalar(v,s) result(ans)
        implicit none
        type(vector_fixed_3d),intent(in)::v
        type(vector_fixed_3d)::ans
        real(kind=kind),intent(in)::s

        ans=v
        ans%direction=ans%direction*s
        return
end function
        
elemental function scalar_times_vector(s,v) result(ans)
        implicit none
        type(vector_fixed_3d),intent(in)::v
        type(vector_fixed_3d)::ans
        real(kind=kind),intent(in)::s

        ans=v
        ans%direction=ans%direction*s
        return
end function
 
elemental function vector_div_scalar(v,s) result(ans)
        implicit none
        type(vector_fixed_3d),intent(in)::v
        type(vector_fixed_3d)::ans
        real(kind=kind),intent(in)::s

        ans=v
        ans%direction=ans%direction/s
        return
end function

pure function move(vector,array) result(ans)
        implicit none
        class(vector_fixed_3d),intent(in)::vector!类型绑定过程的虚参必须用class声明
        type(vector_fixed_3d)::ans
        real(kind=kind),intent(in)::array(3)!移动的方向

        ans=vector
        ans%location=ans%location+array
        return
end function

elemental function mode(vector) result(ans)
        implicit none
        real(kind=kind)::ans
        class(vector_fixed_3d),intent(in)::vector

        ans=sqrt(dot_product(vector%direction,vector%direction))
        return
end function

pure function plane(vector,dot) result(ans)
        integer::ans
        class(vector_fixed_3d),intent(in)::vector
        real(kind=kind),intent(in)::dot(3)
        real::media

        media=dot_product(vector%direction,dot-vector%location)
        if (media) 100,200,300
100             ans=1           !点和向量终点在向量确定的平面同侧
                return
200             ans=0           !点在向量确定的平面上
                return
300             ans=-1          !异侧
                return
end function

elemental function normalize(vector) result(ans)
        implicit none
        class(vector_fixed_3d),intent(in)::vector
        type(vector_fixed_3d)::ans

        ans%location=vector%location
        ans%direction=vector%direction/sqrt(dot_product(vector%direction,vector%direction))
        return
end function

end module
