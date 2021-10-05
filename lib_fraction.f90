module lib_fraction
implicit none

private
public::fractio,assignment(=),operator(+),operator(-),operator(*),operator(/),operator(**)

type::fractio!分数
        integer::numerator!分子
        integer::denominator!分母
contains
        procedure,pass::dec!化小数
        procedure,pass::symplify!化最简
end type

interface assignment(=)
        module procedure int_to_frac
        module procedure arr_to_frac
        module procedure frac_to_frac
end interface

interface operator(+)
        module procedure frac_add_frac
        module procedure int_add_frac
        module procedure frac_add_int
end interface

interface operator(-)
        module procedure frac_minus_frac
        module procedure int_minus_frac
        module procedure frac_minus_int
end interface

interface operator(*)
        module procedure frac_times_frac
        module procedure int_times_frac
        module procedure frac_times_int
end interface

interface operator(/)
        module procedure frac_div_frac
        module procedure int_div_frac
        module procedure frac_div_int
end interface

interface operator(**)
        module procedure frac_int
end interface

contains

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
        return
end subroutine

elemental subroutine frac_to_frac(f2,f1)
        implicit none
        type(fractio),intent(in)::f1
        type(fractio),intent(out)::f2

        f2%numerator=f1%numerator
        f2%denominator=f1%denominator
        return
end subroutine

elemental function frac_add_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%denominator+f1%denominator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        return
end function

elemental function int_add_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%denominator+f2%numerator
        ans%denominator=f2%denominator
        return
end function

elemental function frac_add_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator+f1%denominator*i
        ans%denominator=f1%denominator
        return
end function

elemental function frac_minus_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%denominator-f1%denominator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        return
end function

elemental function int_minus_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%denominator-f2%numerator
        ans%denominator=f2%denominator
        return
end function

elemental function frac_minus_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator-f1%denominator*i
        ans%denominator=f1%denominator
        return
end function

elemental function frac_times_frac(f1,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1,f2

        ans%numerator=f1%numerator*f2%numerator
        ans%denominator=f1%denominator*f2%denominator
        return
end function

elemental function frac_times_int(f1,i) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f1
        integer,intent(in)::i

        ans%numerator=f1%numerator*i
        ans%denominator=f1%denominator
        return
end function

elemental function int_times_frac(i,f2) result(ans)
        implicit none
        type(fractio)::ans
        type(fractio),intent(in)::f2
        integer,intent(in)::i

        ans%numerator=i*f2%numerator
        ans%denominator=f2%denominator
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

        ans%numerator=frac%numerator**i
        ans%denominator=frac%denominator**i
        return
end function

elemental function dec(frac)
        implicit none
        class(fractio),intent(in)::frac
        real::dec

        dec=real(frac%numerator)/real(frac%denominator)
        return
end function

elemental subroutine symplify(frac)
        implicit none
        class(fractio),intent(inout)::frac
        integer::i1,i2

        associate(i1=>frac%numerator,i2=>frac%denominator)!这里不能用数组
                i1=i1/gcd(abs(i1),abs(i2))        
                i2=i2/gcd(abs(i1),abs(i2))        
        end associate
        
        return

        contains
                elemental function gcd(a,b)
                        implicit none
                        integer,intent(in)::a,b
                        integer::m(3),gcd

                        m(1)=a
                        m(2)=b

                        if (m(1) .LT. m(2)) m(1:2)=cshift(m(1:2),shift=1)

                        do
                                m(3)=mod(m(1),m(2))
                                if(m(3) .EQ. 0) exit
                                m=eoshift(m,shift=1)
                        end do

                        gcd=m(2)
                        return
                end function
end subroutine

end module
