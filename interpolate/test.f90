module lib_io
contains

subroutine matrix_write(A,n)
        implicit none
        real(kind=8),intent(in)::A(:,:)
        integer,optional,intent(in)::n
        integer::counter,row
        row=size(A,dim=1)
        if (present(n)) then
                do counter=1,row,1
                        write(unit=n,fmt=*)A(counter,:)
                end do
        else
                do counter=1,row,1
                        write(unit=6,fmt=*)A(counter,:)
                end do
        end if
        return
end subroutine

end module

include "lib_interpolate.f90"

program main
        use lib_interpolate
        use lib_io
        implicit none
        integer::counter
        real(kind=8)::x(5)=(/(real(counter)+0.5,counter=0,4,1)/),y(5)=(/2,1,3,4,2/)
        real(kind=8)::location(100)=(/(real(counter)/15,counter=-10,89,1)/),ans(100,4)
        call initialize(x,y)
        call spline_initialize()!自由边界
        ans(:,1)=spline(location)
        call spline_initialize(real(5,kind=8))!周期边界
        ans(:,2)=spline(location)
        call spline_initialize(real(0.5,kind=8),real(-0.3,kind=8))!夹持边界
        ans(:,3)=spline(location)
        ans(:,4)=lagrange(location)!拉格朗日
        open(unit=10,file="result.txt")
                call matrix_write(ans,10)
        close(unit=10)
        stop
end program
