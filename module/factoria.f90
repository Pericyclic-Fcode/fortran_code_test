module factoria

contains

!测试fortran自带的sum函数和用dot_product(x,unus)内积求和的速度
!一千个数据相加没区别，一百万个数据相加sum用时比内积少将近一半
!用sum改写代码，还是自带的函数优化得好
elemental function factorial(x)!阶乘
        integer,intent(in)::x
        integer::factorial,counter
        if(x .EQ. 0) then
                factorial=1
                return
        end if
        factorial=product((/(counter,counter=1,x,1)/))
        return
end function

elemental function ln_factorial(x)!阶乘的对数更常用
        !  x!=exp(ln(PIx_i))=exp(SIGMA(ln(x_i))),向量化，求和用内积
        implicit none
        integer,intent(in)::x
        real(kind=8)::ln_factorial,media(x)
        integer::counter
        
        media=(/(counter,counter=1,x,1)/)!不能直接用x充当limit在声明时赋值
        ln_factorial=sum(log(media))
        
        return
end function

elemental function combination(m,n) result(ans)!组合数
        !C(m,n)=m!/n!(m-n)!=exp(SIGMA(i=m-n+1,m)(ln(m_i))-SIGMA(ln(n)))
        implicit none
        integer,intent(in)::m,n
        integer::counter,ans
        real(kind=8)::top(n),bottom(n)

        if(n .EQ. 0) then
                ans=1
                return
        end if
        
        top=(/(counter,counter=m-n+1,m,1)/)
        bottom=(/(counter,counter=1,n,1)/)
        ans=nint(exp(sum(dlog(top))-sum(dlog(bottom))))
        
        return
end function

elemental function arrangement(m,n)!排列数
        !A(m,n)=m!/(m-n)!=exp(SIGMA(i=m-n+1,m)(ln(m_i)))
        implicit none
        integer,intent(in)::m,n
        integer::arrangement
        integer::counter
        real(kind=8)::media(n)
        
        media=(/(counter,counter=m-n+1,m,1)/)
        arrangement=nint(exp(sum(dlog(media))))
        
        return
end function

end module
