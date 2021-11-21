module lib_legendre
!函数默认逐元
        implicit none
private
public::legendre,spherical,mod_legendre!求模函数输出为双精度
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号
        real(kind=dp),parameter::pi=3.14159265358979
        complex(kind=2*dp),parameter::i=(0D0,1D0)

interface legendre!legendre(x,l,m)
        module procedure legendre_sp!勒让德多项式，输入输出均为单精度
        module procedure legendre_dp!双精度
        module procedure associated_legendre_sp!单精度连带勒让德
        module procedure associated_legendre_dp!双精度
end interface
        
interface spherical!球函数spherical(\theta,\phi,l,m)
        module procedure Y_dp
        module procedure Y_sp
end interface

contains

elemental function combination(m,n) result(ans)!组合数
        !C(m,n)=m!/n!(m-n)!=exp(SIGMA(i=m-n+1,m)(ln(m_i))-SIGMA(ln(n)))
        implicit none
        integer,intent(in)::m,n
        integer::counter,ans
        real(kind=dp),allocatable::top(:),bottom(:)!在逐元函数中，虚参不能用于声明变量，除非是查询函数

        if(n .EQ. 0) then
                ans=1
                return
        end if

        allocate(top(n),bottom(n))!只能这样手动配置
        top=(/(counter,counter=m-n+1,m,1)/)
        bottom=(/(counter,counter=1,n,1)/)
        ans=nint(exp(sum(dlog(top))-sum(dlog(bottom))))
        return
end function

pure function mod_legendre(l,m) result(ans)!勒让德多项式的模
        implicit none
        integer,intent(in)::l
        integer,intent(in),optional::m
        integer::counter
        real(kind=dp)::ans

        if (present(m)) then
                ans=sqrt(product((/(counter,counter=l-m+1,l+m,1)/))*2/real(2*l+1))
        else
                ans=sqrt(2/real(2*l+1))
        end if

        return
end function

elemental function legendre_sp(x,n) result(ans)!P_{n}(x)!n阶勒让德函数，输入-1<x<1
        implicit none
        real(kind=sp),intent(in)::x
        integer,intent(in)::n
        integer::limit,counter
        real(kind=sp),allocatable::x_vector(:),l_vector(:)
        integer,allocatable::m(:)
        real(kind=sp)::ans

        limit=floor(real(n)/2)
        allocate(x_vector(0:limit),l_vector(0:limit),m(0:limit))

        m=(/(counter,counter=0,limit,1)/)   !求和指标
        l_vector=real(((-1)**m)*combination(n,m)*combination(2*(n-m),n))/real(2**n)!系数，注意要用real函数把整型变实型
        x_vector=x**(n-2*m)!自变量
        ans=dot_product(l_vector,x_vector)
        return
end function

elemental function associated_legendre_sp(x,n,m) result(ans) !P_{n}^{m}(x)=-1^{m}(1-x^{2})^{m/2}\frac{d^{m}}{dx^{m}}P_{n}(x)
        implicit none
        real(kind=sp)::ans          !n阶连带勒让德
        real(kind=sp),intent(in)::x
        integer,intent(in)::n,m
        integer::limit,counter
        real(kind=sp),allocatable::x_vector(:),l_vector(:)
        integer,allocatable::k(:)

        if (m .EQ. 0) then
                ans=legendre_sp(x,n)
                return
        end if

        limit=floor(real(n)/2)
        allocate(x_vector(0:limit),l_vector(0:limit),k(0:limit))

        k=(/(counter,counter=0,limit,1)/)   !求和指标
        l_vector=real(((-1)**k)*combination(n,k)*combination(2*(n-k),n))/real(2**n)!系数
        x_vector=product(reshape((/(n-2*k-counter,counter=0,m-1,1)/),(/limit+1,m/)),dim=2)*(x**(n-2*k-m))!自变量
        where (n-2*k-m .LT. 0) x_vector=0
        ans=dot_product(l_vector,x_vector)*((1-x**2)**(real(m)/2))*((-1)**m)
        return
end function

elemental function legendre_dp(x,n) result(ans)!P_{n}(x)!n阶勒让德函数，输入-1<x<1
        implicit none
        real(kind=dp),intent(in)::x
        integer,intent(in)::n
        integer::limit,counter
        real(kind=dp),allocatable::x_vector(:),l_vector(:)
        integer,allocatable::m(:)
        real(kind=dp)::ans

        limit=floor(real(n)/2)
        allocate(x_vector(0:limit),l_vector(0:limit),m(0:limit))

        m=(/(counter,counter=0,limit,1)/)   !求和指标
        l_vector=real(((-1)**m)*combination(n,m)*combination(2*(n-m),n))/real(2**n)!系数，注意要用real函数把整型变实型
        x_vector=x**(n-2*m)!自变量
        ans=dot_product(l_vector,x_vector)
        return
end function

elemental function associated_legendre_dp(x,n,m) result(ans) !P_{n}^{m}(x)=-1^{m}(1-x^{2})^{m/2}\frac{d^{m}}{dx^{m}}P_{n}(x)
        implicit none
        real(kind=dp)::ans          !n阶连带勒让德
        real(kind=dp),intent(in)::x
        integer,intent(in)::n,m
        integer::limit,counter
        real(kind=dp),allocatable::x_vector(:),l_vector(:)
        integer,allocatable::k(:)

        if (m .EQ. 0) then
                ans=legendre_dp(x,n)
                return
        end if

        limit=floor(real(n)/2)
        allocate(x_vector(0:limit),l_vector(0:limit),k(0:limit))

        k=(/(counter,counter=0,limit,1)/)   !求和指标
        l_vector=real(((-1)**k)*combination(n,k)*combination(2*(n-k),n))/real(2**n)!系数
        x_vector=product(reshape((/(n-2*k-counter,counter=0,m-1,1)/),(/limit+1,m/)),dim=2)*(x**(n-2*k-m))!自变量
        where (n-2*k-m .LT. 0) x_vector=0
        ans=dot_product(l_vector,x_vector)*((1-x**2)**(real(m)/2))*((-1)**m)
        return
end function

elemental function Y_dp(theta,phi,l,m) result(ans)
        implicit none
        real(kind=dp),intent(in)::theta,phi
        integer,intent(in)::l,m
        complex(kind=2*dp)::ans
        ans=((-1)**m)*legendre(cos(theta),l,m)*exp(-1*i*m*phi)/(sqrt(2*pi)*mod_legendre(l,m))
        return
end function Y_dp

elemental function Y_sp(theta,phi,l,m) result(ans)
        implicit none
        real(kind=sp),intent(in)::theta,phi
        integer,intent(in)::l,m
        complex(kind=2*sp)::ans
        ans=((-1)**m)*legendre(cos(theta),l,m)*exp(-1*i*m*phi)/(sqrt(2*pi)*mod_legendre(l,m))
        return
end function Y_sp

end module
