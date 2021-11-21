module lib_fraction
        implicit none

private!reshape,transpose等函数可直接用于派生类型，matmul,sum,dot_product不行
public::fractio,assignment(=),operator(+),operator(-),operator(*),operator(/),operator(**)
public::operator(.x.),operator(.i.),operator(.over.)!自定义运算符优先级最低
        
        integer,parameter,private::p=selected_int_kind(10)!精度，此时p=8，为了应付希尔伯特阵设置得特别大
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号

type::fractio!分数
        integer(kind=p)::numerator!分子
        integer(kind=p)::denominator!分母
contains
        procedure,pass::dec!化小数function
        procedure,pass::symplify!化最简subroutine
        procedure,pass::check!检查分母是否为零，为零返回.FALSE.function
        procedure,pass::reciprocal!求倒数function
end type

interface assignment(=)!等号重载
        module procedure int_to_frac!整数变成分母为1的分数
        module procedure arr_to_frac!二员数组第一个变成分子，第二个分母
        module procedure frac_to_frac!分数复制
end interface

interface operator(+)!加法
        module procedure frac_add_frac
        module procedure int_add_frac
        module procedure frac_add_int
end interface

interface operator(-)!减法
        module procedure frac_minus_frac
        module procedure int_minus_frac
        module procedure frac_minus_int
end interface

interface operator(*)!乘法
        module procedure frac_times_frac
        module procedure int_times_frac
        module procedure frac_times_int
end interface

interface operator(/)!除法
        module procedure frac_div_frac!分数除分数
        module procedure int_div_frac!整数除分数
        module procedure frac_div_int!分数除整数
end interface

interface operator(**)!乘非负整数次方
        module procedure frac_int
end interface

interface operator(.i.)!矩阵求逆
        module procedure matrix_inverse
end interface

interface operator(.x.)!矩阵相乘
        module procedure matrix_multiply_frac!两个都是分数阵
        module procedure matrix_multiply_int1!整数阵在后
        module procedure matrix_multiply_int2!整数阵在前
end interface

interface operator(.over.)!a.over.b=a/b分数
        module procedure over
end interface

contains

elemental function gcd(a,b) result(ans)!求最大公约数
        implicit none
        integer(kind=p),intent(in)::a,b!只有这里和type定义处要加(kind=p)
        integer::m(3),ans

        m(1:2)=abs((/a,b/))
        if(any(m(1:2) .EQ. 0)) then
                ans=1
                return
        end if

        if (m(1) .LT. m(2)) m(1:2)=cshift(m(1:2),shift=1)
        do
                m(3)=mod(m(1),m(2))
                if(m(3) .EQ. 0) exit
                m=eoshift(m,shift=1)
        end do

        ans=m(2)
        return
end function gcd

elemental subroutine int_to_frac(frac,inte)
        implicit none
        integer,intent(in)::inte
        type(fractio),intent(out)::frac

        frac%numerator=inte
        frac%denominator=1
        return
end subroutine

pure subroutine arr_to_frac(frac,arr)
        implicit none
        integer,intent(in)::arr(2)
        type(fractio),intent(out)::frac

        frac%numerator=arr(1)
        frac%denominator=arr(2)
        if(frac%numerator .EQ. 0) frac%denominator=1
        return
end subroutine

elemental subroutine frac_to_frac(f2,f1)
        implicit none
        type(fractio),intent(in)::f1
        type(fractio),intent(out)::f2

        f2%numerator=f1%numerator
        f2%denominator=f1%denominator
        if(f2%numerator .EQ. 0) f2%denominator=1
        return
end subroutine

elemental function frac_add_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%denominator+f1%denominator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function int_add_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%denominator+f2%numerator
        ans%denominator=f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_add_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator+f1%denominator*i
        ans%denominator=f1%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_minus_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%denominator-f1%denominator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function int_minus_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%denominator-f2%numerator
        ans%denominator=f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_minus_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator-f1%denominator*i
        ans%denominator=f1%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_times_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_times_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator*i
        ans%denominator=f1%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function int_times_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%numerator
        ans%denominator=f2%denominator
        if(ans%numerator .EQ. 0) ans%denominator=1
        return
end function

elemental function frac_div_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%denominator
        ans%denominator=f1%denominator*f2%numerator
        return
end function

elemental function frac_div_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator
        ans%denominator=f1%denominator*i
        return
end function

elemental function int_div_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%denominator
        ans%denominator=f2%numerator
        return
end function

elemental function frac_int(frac,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::frac
        integer,intent(in)::i

        if(i .LT. 0) return
        ans%numerator=frac%numerator**i
        ans%denominator=frac%denominator**i
        return
end function

elemental function dec(frac)
        implicit none
        class(fractio),intent(in)::frac
        real(kind=dp)::dec

        dec=real(frac%numerator)/real(frac%denominator)
        return
end function dec

elemental subroutine symplify(frac)
        implicit none
        class(fractio),intent(inout)::frac
        integer::tmp

        if(.not.frac%check()) then
                frac%numerator=1
                return
        end if
        if(frac%numerator .EQ. 0) then
                frac%denominator=1
                return
        end if

        tmp=gcd(frac%numerator,frac%denominator) 
        frac%numerator=frac%numerator/tmp
        frac%denominator=frac%denominator/tmp
        
        return
end subroutine symplify

elemental function check(frac) result(ans)
        implicit none
        class(fractio),intent(in)::frac
        logical::ans
        ans=.TRUE.
        if(frac%denominator .EQ. 0) ans=.FALSE.
        return
end function check

elemental function reciprocal(frac) result(ans)
        implicit none
        class(fractio),intent(in)::frac
        type(fractio)::ans
        ans%numerator=frac%denominator
        ans%denominator=frac%numerator
        return
end function reciprocal

pure function matrix_multiply_frac(A,B) result(ans)
        implicit none
        type(fractio),intent(in)::A(:,:),B(:,:)
        type(fractio)::ans(size(A,dim=1),size(B,dim=2))
        integer::i,counter

        ans=0
        do counter=1,size(B,dim=1),1
                ans=ans+reshape((/(B(counter,i)*A(:,counter),i=1,size(B,dim=2),1)/),(/size(A,dim=1),size(B,dim=2)/))
                call ans%symplify()!为了应付极端病态的希尔伯特矩阵每步都要约分
        end do
        return
end function matrix_multiply_frac

pure function matrix_multiply_int1(A,B) result(ans)
        implicit none
        type(fractio),intent(in)::A(:,:)
        integer,intent(in)::B(:,:)
        type(fractio)::ans(size(A,dim=1),size(B,dim=2))
        integer::i,counter

        ans=0
        do counter=1,size(B,dim=1),1
                ans=ans+reshape((/(B(counter,i)*A(:,counter),i=1,size(B,dim=2),1)/),(/size(A,dim=1),size(B,dim=2)/))
                call ans%symplify()!为了应付极端病态的希尔伯特矩阵每步都要约分
        end do
        return
end function matrix_multiply_int1

pure function matrix_multiply_int2(A,B) result(ans)
        implicit none
        integer,intent(in)::A(:,:)
        type(fractio),intent(in)::B(:,:)
        type(fractio)::ans(size(A,dim=1),size(B,dim=2))
        integer::i,counter

        ans=0
        do counter=1,size(B,dim=1),1
                ans=ans+reshape((/(B(counter,i)*A(:,counter),i=1,size(B,dim=2),1)/),(/size(A,dim=1),size(B,dim=2)/))
                call ans%symplify()!为了应付极端病态的希尔伯特矩阵每步都要约分
        end do
        return
end function matrix_multiply_int2

pure function matrix_inverse(A) result(ans)!求逆阵，不对退化的矩阵负责
        implicit none
        type(fractio),intent(in)::A(:,:)
        type(fractio)::ans(size(A,dim=1),size(A,dim=2)),media(size(A,dim=1),2*size(A,dim=2))
        integer::n,counter,i

        n=size(A,dim=1)
        media=0!构建增广阵
        media(:,:n)=A
        forall (counter=1:n) media(counter,counter+n)=1!右半截是单位阵

        do counter=1,n-1,1    !初等行变换化上三角阵
                i=1
                if(all(media(counter:,counter)%numerator .EQ. 0)) cycle
                do while(media(counter,counter)%numerator .EQ. 0)
                        call swap(media(counter,:),media(counter+i,:))                     !防止对角线上是0
                        i=i+1
                end do
                media(counter,:)=media(counter,:)/media(counter,counter)
                media(counter+1:,:)=media(counter+1:,:)-(media(counter+1:,counter:counter) .x. media(counter:counter,:))
                call media%symplify()
        end do
        if (media(n,n)%numerator .NE. 0) media(n,:)=media(n,:)/media(n,n)!最后一行归一化
        do counter=n,2,-1     !化行最简
                if (media(counter,counter)%numerator .EQ. 0) cycle
                media(:counter-1,:)=media(:counter-1,:)-(media(:counter-1,counter:counter) .x. media(counter:counter,:))!自定义运算符优先级最低
                call media%symplify()!为了应付极端病态的希尔伯特矩阵每步都要约分
        end do
        ans=media(:,n+1:)
        return
end function matrix_inverse

elemental function over(a,b) result(ans)
        implicit none
        integer,intent(in)::a,b
        type(fractio)::ans
        ans=(/a,b/)
        return
end function over

elemental subroutine swap(a,b)!用于交换元
    type(fractio),intent(inout)::a,b
    type(fractio)::tmp
    tmp=a;a=b;b=tmp
end subroutine swap

end module
