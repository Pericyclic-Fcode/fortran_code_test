module decomposition
implicit none
        real(kind=8),private,parameter::pi=3.14159265258979
        complex(kind=16),private,parameter::i=(0,1)
        real(kind=8),private,parameter::dt=0.0001!积分步长
        real(kind=8),private,parameter::error=1.0E-7
        real(kind=8),private,allocatable,save::t(:)
        real(kind=8),private,allocatable,save::f(:)

contains

pure function fourier(f,n) result(ans)!这样写更清晰
        !傅里叶分解，第一列是分解到余弦上的值，第二列是分解到正弦上的值
        implicit none
        real(kind=8),intent(in)::f(:)!待分解的函数
        integer,intent(in)::n!截断的谐波次数
        integer::counter,limit
        real(kind=8)::location(size(f))
        real(kind=8)::ans(0:n,2)
        
        limit=size(f)!输入的f视作一周期2\pi上的函数
        location=(/(counter*2*pi/real(limit),counter=1,limit,1)/)!取样点的位置
        forall(counter=0:n) ans(counter,1)=dot_product(cos(counter*location),f)/(0.5*limit)
        forall(counter=0:n) ans(counter,2)=dot_product(sin(counter*location),f)/(0.5*limit)
        ans(0,1)=0.5*ans(0,1) !第一项归一化
        return
end function

subroutine func_initialize(func,t_max)
        implicit none
        real(kind=8),external::func!函数作为虚参传入时应加external声明
        real(kind=8),optional,intent(in)::t_max
        integer::n,counter

        if(present(t_max)) then    !没给截断值的话自己算一个出来
                n=nint(t_max/dt)
        else
                n=get_n_max(func)
        end if

        if(allocated(t)) deallocate(t)
        if(allocated(f)) deallocate(f)

        allocate(t(-n:n),f(-n:n))
        t=(/(counter*dt,counter=-n,n,1)/)
        f=(/(func(counter*dt),counter=-n,n,1)/)
        !f=func(t)   !如果func是逐元的
        return

        contains
        function get_n_max(func) result(ans)!确定一个截断
        !无法确定func是纯函数
                implicit none
                integer::ans,counter
                real(kind=8),external::func
                
                counter=2

                do
                        if (abs(func((2**counter)*dt)-func((2**(counter-1))*dt)) .LE. error) exit
                        counter=counter+1
                end do
                
                ans=2**counter
                return
        end function
end subroutine

elemental function ft(omega)!傅里叶变换 F(\omega)
        implicit none
        real(kind=8),intent(in)::omega
        complex(kind=16)::ft

        ft=dot_product(f,exp(-i*omega*t))*dt
        return
end function

end module
