module lib_sort
        implicit none
private
public::permutation,get_order,qsort
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号

interface permutation!按给的顺序进行置换
        module procedure permutation_dp
        module procedure permutation_sp
end interface

interface get_order!古老低效的获得由小到大顺序的算法
        module procedure get_order_dp
        module procedure get_order_sp
end interface

!代码来自https://fortran-lang.discourse.group/t/modern-fortran-sample-code/2019/2
interface qsort!快速排序算法，大佬写的
        module procedure qsort_dp!qsort(array)
        module procedure qsort_sp
!        module procedure qsort_associatematrix_dp!qsort(array,matrix)
!        module procedure qsort_associatematrix_sp!按array顺序进行对matrix按行同时排序
!        module procedure qsort_associatevector_dp!qsort(array,vector)
!        module procedure qsort_associatevector_sp
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

pure recursive function qsort_dp(data) result(sorted)
        real(kind=dp),intent(in)::data(:)
        real(kind=dp)::sorted(size(data))
        if (size(data) > 1) then
                sorted=(/qsort(pack(data(2:),data(2:)<data(1))), data(1), &
                qsort(pack(data(2:),data(2:)>=data(1)))/)
        else
                sorted = data
        end if
end function qsort_dp

pure recursive function qsort_sp(data) result(sorted)
        real(kind=sp),intent(in)::data(:)
        real(kind=sp)::sorted(size(data))
        if (size(data) > 1) then
                sorted=(/qsort(pack(data(2:),data(2:)<data(1))), data(1), &
                qsort(pack(data(2:),data(2:)>=data(1)))/)
        else
                sorted = data
        end if
end function qsort_sp

end module
