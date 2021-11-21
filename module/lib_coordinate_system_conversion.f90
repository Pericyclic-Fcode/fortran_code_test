module lib_coordinate_system_conversion
        implicit none
private
public::polar_to_cartesian,cartesian_to_polar
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号
        real(kind=dp),private,parameter::pi=3.14159265358979

interface polar_to_cartesian!\rho的位置变x，\theta的位置变y
        module procedure polar_to_cartesian_sp
        module procedure polar_to_cartesian_dp
end interface

interface cartesian_to_polar!x的位置变\rho，y的位置变\theta,返回[0,2\pi)
        module procedure cartesian_to_polar_sp
        module procedure cartesian_to_polar_dp
end interface

contains

elemental subroutine polar_to_cartesian_sp(rho,theta)
        implicit none
        real(kind=sp),intent(inout)::rho,theta
        real(kind=sp)::tmp
        tmp=rho;rho=tmp*cos(theta);theta=tmp*sin(theta)
        return
end subroutine polar_to_cartesian_sp

elemental subroutine polar_to_cartesian_dp(rho,theta)
        implicit none
        real(kind=dp),intent(inout)::rho,theta
        real(kind=dp)::tmp
        tmp=rho;rho=tmp*cos(theta);theta=tmp*sin(theta)
        return
end subroutine polar_to_cartesian_dp

elemental subroutine cartesian_to_polar_sp(x,y)
        implicit none
        real(kind=sp),intent(inout)::x,y
        real(kind=sp)::tmp
        tmp=x;x=sqrt(x**2+y**2);y=atan2(y,tmp)
        if(y .LT. 0.0) y=y+2*pi
        return
end subroutine cartesian_to_polar_sp

elemental subroutine cartesian_to_polar_dp(x,y)
        implicit none
        real(kind=dp),intent(inout)::x,y
        real(kind=dp)::tmp
        tmp=x;x=sqrt(x**2+y**2);y=atan2(y,tmp)
        if(y .LT. 0.0) y=y+2*pi
        return
end subroutine cartesian_to_polar_dp

end module
