module lib_sort
        implicit none
private
public::permutation,get_order
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号

interface permutation
        module procedure permutation_dp
        module procedure permutation_sp
end interface

interface get_order
        module procedure get_order_dp
        module procedure get_order_sp
end interface

contains
pure function permutation_dp(x,arr) result(ans)!按整型数组arr提供的序列对x进行置换
        implicit none
        integer,intent(in)::arr(:)
        real(kind=dp),intent(in)::x(size(arr))
        real(kind=dp)::ans(size(arr))
        integer::counter

        forall (counter=1:size(arr)) ans(arr(counter))=x(counter)
        return
end function

pure function get_order_dp(x) result(ans)!获得由小到大的顺序
        implicit none
        real(kind=dp),intent(in)::x(:)
        integer::ans(size(x)),counter
        logical::mask(size(x))!纯函数不能在这初始化

        mask=.TRUE.
        do counter=1,size(x),1
                ans(minloc(x,mask))=counter
                mask(minloc(x,mask))=.FALSE.
        end do
        return
end function

pure function permutation_sp(x,arr) result(ans)!按整型数组arr提供的序列对x进行置换
        implicit none
        integer,intent(in)::arr(:)
        real(kind=sp),intent(in)::x(size(arr))
        real(kind=sp)::ans(size(arr))
        integer::counter

        forall (counter=1:size(arr)) ans(arr(counter))=x(counter)
        return
end function

pure function get_order_sp(x) result(ans)!获得由小到大的顺序
        implicit none
        real(kind=sp),intent(in)::x(:)
        integer::ans(size(x)),counter
        logical::mask(size(x))!纯函数不能在这初始化

        mask=.TRUE.
        do counter=1,size(x),1
                ans(minloc(x,mask))=counter
                mask(minloc(x,mask))=.FALSE.
        end do
        return
end function

end module
