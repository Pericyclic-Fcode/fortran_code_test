module lib_io
        implicit none
private
public::matrix_write,matrix_read,fileline,num_number,tex_matrix_write,loadtxt
        integer,private,parameter::sp=selected_real_kind(p=6,r=37)!获得单精度实型类别号
        integer,private,parameter::dp=selected_real_kind(p=13,r=200)!获得双精度实型类别号
        character(len=2),parameter,private::quote='"'//"'"
        character(len=3),parameter,private::seperator=" "//","//char(9)

interface loadtxt!loadtxt(A,n)A是输入矩阵，要求有allocatable属性,n是文件通道号，文件要求形式规整
        module procedure loadtxt_file_dp!loadtxt(A,"filename")未打开文件时使用
        module procedure loadtxt_file_sp
        module procedure loadtxt_dp
        module procedure loadtxt_sp
end interface

interface matrix_write!(A[,n][,model])A是矩阵，n是通道号，不填默认6，model填C按行优先的顺序输出，F按列优先输出
        module procedure matrix_write_dp
        module procedure matrix_write_sp
        module procedure matrix_write_vector_dp!统一输出成列向量
        module procedure matrix_write_vector_sp
end interface

interface matrix_read!(A,n[,model])A是矩阵，n是通道号，model填C按行优先的顺序输入，F按列优先输入
        module procedure matrix_read_dp
        module procedure matrix_read_sp
        module procedure matrix_read_vector_dp
        module procedure matrix_read_vector_sp
end interface

interface tex_matrix_write!输出成pmatrix的形式，双精度四位小数，单精度三位
        module procedure tex_matrix_write_dp!(A,n[,model]),A是矩阵，n是通道号
        module procedure tex_matrix_write_sp
        module procedure tex_matrix_write_vector_dp!默认输出成列向量
        module procedure tex_matrix_write_vector_sp!model项加row输出成行向量，加column输出成列向量
end interface

contains

subroutine loadtxt_file_dp(A,filename)
        implicit none
        real(kind=dp),intent(inout),allocatable::A(:,:)
        character(len=*),intent(in)::filename
        integer::row,column,n
        character(len=512)::tmp

        open(newunit=n,file=filename)
                read(unit=n,fmt="(A512)")tmp
                row=fileline(n)!矩阵行列数
                column=num_number(tmp)
        
                if(allocated(A)) deallocate(A)
                allocate(A(row,column))
                rewind(unit=n)
                call matrix_read(A,n)
        close(unit=n)
        return
end subroutine loadtxt_file_dp

subroutine loadtxt_file_sp(A,filename)
        implicit none
        real(kind=sp),intent(inout),allocatable::A(:,:)
        character(len=*),intent(in)::filename
        integer::row,column,n
        character(len=512)::tmp

        open(newunit=n,file=filename)
                read(unit=n,fmt="(A512)")tmp
                row=fileline(n)!矩阵行列数
                column=num_number(tmp)
        
                if(allocated(A)) deallocate(A)
                allocate(A(row,column))
                rewind(unit=n)
                call matrix_read(A,n)
        close(unit=n)
        return
end subroutine loadtxt_file_sp

subroutine loadtxt_dp(A,n)!在子程序里可以对输入矩阵重新allocate
        implicit none
        real(kind=dp),intent(inout),allocatable::A(:,:)
        integer,intent(in)::n
        integer::row,column
        character(len=512)::tmp

        read(unit=n,fmt="(A512)")tmp
        row=fileline(n)!矩阵行列数
        column=num_number(tmp)

        if(allocated(A)) deallocate(A)
        allocate(A(row,column))
        rewind(unit=n)
        call matrix_read(A,n)
        rewind(unit=n)
        return
end subroutine loadtxt_dp

subroutine loadtxt_sp(A,n)
        implicit none
        real(kind=sp),intent(inout),allocatable::A(:,:)
        integer,intent(in)::n
        integer::row,column
        character(len=512)::tmp

        read(unit=n,fmt="(A512)")tmp
        row=fileline(n)!矩阵行列数
        column=num_number(tmp)

        if(allocated(A)) deallocate(A)
        allocate(A(row,column))
        rewind(unit=n)
        call matrix_read(A,n)
        rewind(unit=n)
        return
end subroutine loadtxt_sp

function fileline(passageway) result(ans)!获取文件行数
        implicit none
        integer,intent(in)::passageway
        integer::ios,ans
        character(len=1)::cdummy
        ans=0
        rewind(unit=passageway)
        do
                read(unit=passageway,fmt=*,iostat=ios)cdummy
                if(ios .NE. 0) exit
                ans=ans+1
        end do
        rewind(unit=passageway)
        return
end function fileline

pure function num_number(string) result(ans)                    !获取字符串里的数据个数，引号内的也算
        implicit none
        character(len=*),intent(in)::string
        integer::ans,counter
        logical::mask(len_trim(string))       !掩码数组 !两头补上一个空格，空格和引号，逗号变成.FALSE.
        mask=(/((scan(string(counter:counter),seperator//quote) .EQ. 0),counter=1,size(mask),1)/)
        ans=count((/mask,.FALSE./) .NEQV. (/.FALSE.,mask/))/2!异或运算，变号的地方为真，除以二是数据个数
        return
end function

subroutine matrix_write_dp(A,n,model)
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        integer,optional,intent(in)::n
        character(len=1),optional,target,intent(in)::model
        integer::counter,passageway=6

        character(len=1),pointer::mode

        if (present(n)) passageway=n
        if (present(model)) then!控制一下，防止输入奇怪的参数
                mode=>model
        else 
                allocate(mode)
                mode="C"!默认输出成行优先的形式
        end if

        if (mode .EQ."F" .OR. mode .EQ. "f") then
                write(unit=passageway,fmt=100) (A(:,counter),counter=1,size(A,dim=2),1)
        else if(mode .EQ."C" .OR. mode .EQ. "c") then
                write(unit=passageway,fmt=200) (A(counter,:),counter=1,size(A,dim=1),1)
        end if
        return
100     format(<size(A,dim=2)-1>(<size(A,dim=1)>(ES20.12,2X,\),//),<size(A,dim=1)>(ES20.12,2X))
200     format(<size(A,dim=1)-1>(<size(A,dim=2)>(ES20.12,2X,\),//),<size(A,dim=2)>(ES20.12,2X))
end subroutine matrix_write_dp

subroutine matrix_write_sp(A,n,model)
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        integer,optional,intent(in)::n
        character(len=1),optional,target,intent(in)::model
        integer::counter,passageway=6

        character(len=1),pointer::mode

        if (present(n)) passageway=n
        if (present(model)) then!控制一下，防止输入奇怪的参数
                mode=>model
        else 
                allocate(mode)
                mode="C"!默认输出成行优先的形式
        end if

        if (mode .EQ."F" .OR. mode .EQ. "f") then
                write(unit=passageway,fmt=100) (A(:,counter),counter=1,size(A,dim=2),1)
        else if(mode .EQ."C" .OR. mode .EQ. "c") then
                write(unit=passageway,fmt=200) (A(counter,:),counter=1,size(A,dim=1),1)
        end if
        return
100     format(<size(A,dim=2)-1>(<size(A,dim=1)>(ES13.5,2X,\),//),<size(A,dim=1)>(ES13.5,2X))
200     format(<size(A,dim=1)-1>(<size(A,dim=2)>(ES13.5,2X,\),//),<size(A,dim=2)>(ES13.5,2X))
end subroutine matrix_write_sp

subroutine matrix_write_vector_dp(x,n,model)
        implicit none
        real(kind=dp),intent(in)::x(:)
        integer,optional,intent(in)::n
        integer::passageway=6
        character(len=*),optional,intent(in)::model

        if(present(n)) passageway=n
        if(present(model) .AND. model .EQ."row") then
                write(unit=passageway,fmt=100)x!输出成行向量
        else if(present(model) .AND. model .EQ."column") then
                write(unit=passageway,fmt=200)x!输出成列向量
        else if(.NOT.present(model)) then
                write(unit=passageway,fmt=200)x!默认输出成列向量
        end if
        return
100     format(<size(x)-1>(F14.3,\),(F14.3))
200     format(<size(x)-1>(F14.3,//),(F14.3))
end subroutine matrix_write_vector_dp

subroutine matrix_write_vector_sp(x,n,model)
        implicit none
        real(kind=sp),intent(in)::x(:)
        integer,optional,intent(in)::n
        integer::passageway=6
        character(len=*),optional,intent(in)::model

        if(present(n)) passageway=n
        if(present(model) .AND. model .EQ."row") then
                write(unit=passageway,fmt=100)x!输出成行向量
        else if(present(model) .AND. model .EQ."column") then
                write(unit=passageway,fmt=200)x!输出成列向量
        else if(.NOT.present(model)) then
                write(unit=passageway,fmt=200)x!默认输出成列向量
        end if
        return
100     format(<size(x)-1>(F7.3,\),(F7.3))
200     format(<size(x)-1>(F7.3,//),(F7.3))
end subroutine matrix_write_vector_sp

subroutine matrix_read_dp(A,n,model)
        implicit none
        real(kind=dp),intent(out)::A(:,:)
        integer,intent(in)::n
        character(len=1),optional,target,intent(in)::model
        integer::counter

        character(len=1),pointer::mode

        if (present(model)) then!控制一下，防止输入奇怪的参数
                mode=>model
        else 
                allocate(mode)
                mode="C"!默认输入成行优先的形式
        end if

        if (mode .EQ."F" .OR. mode .EQ. "f") then
                read(unit=n,fmt=*) (A(:,counter),counter=1,size(A,dim=2),1)
        else if(mode .EQ."C" .OR. mode .EQ. "c") then
                read(unit=n,fmt=*) (A(counter,:),counter=1,size(A,dim=1),1)
        end if
        return
end subroutine matrix_read_dp

subroutine matrix_read_sp(A,n,model)
        implicit none
        real(kind=sp),intent(out)::A(:,:)
        integer,intent(in)::n
        character(len=1),optional,target,intent(in)::model
        integer::counter

        character(len=1),pointer::mode

        if (present(model)) then!控制一下，防止输入奇怪的参数
                mode=>model
        else 
                allocate(mode)
                mode="C"!默认输入成行优先的形式
        end if

        if (mode .EQ."F" .OR. mode .EQ. "f") then
                read(unit=n,fmt=*) (A(:,counter),counter=1,size(A,dim=2),1)
        else if(mode .EQ."C" .OR. mode .EQ. "c") then
                read(unit=n,fmt=*) (A(counter,:),counter=1,size(A,dim=1),1)
        end if
        return
end subroutine matrix_read_sp

subroutine matrix_read_vector_dp(x,n)
        implicit none
        real(kind=dp),intent(out)::x(:)
        integer,intent(in)::n
        integer::counter

        read(unit=n,fmt=*)(x(counter),counter=1,size(x),1)
        return
end subroutine matrix_read_vector_dp

subroutine matrix_read_vector_sp(x,n)
        implicit none
        real(kind=sp),intent(out)::x(:)
        integer,intent(in)::n
        integer::counter

        read(unit=n,fmt=*)(x(counter),counter=1,size(x),1)
        return
end subroutine matrix_read_vector_sp

subroutine tex_matrix_write_dp(A,n)
        implicit none
        real(kind=dp),intent(in)::A(:,:)
        integer,intent(in)::n
        integer::counter

        write(unit=n,fmt=*) "$$\begin{pmatrix}"
        write(unit=n,fmt=100) (A(counter,:),counter=1,size(A,dim=1),1)
        write(unit=n,fmt=*) "\end{pmatrix} $$"

        return
100     format(<size(A,dim=1)-1>(<size(A,dim=2)>(F14.4,1x,"&",1x,\),"\\",//),<size(A,dim=2)>(F14.4,1X,"&",1x))
end subroutine tex_matrix_write_dp

subroutine tex_matrix_write_sp(A,n)
        implicit none
        real(kind=sp),intent(in)::A(:,:)
        integer,intent(in)::n
        integer::counter

        write(unit=n,fmt=*) "$$\begin{pmatrix}"
        write(unit=n,fmt=100) (A(counter,:),counter=1,size(A,dim=1),1)
        write(unit=n,fmt=*) "\end{pmatrix} $$"

        return
100     format(<size(A,dim=1)-1>(<size(A,dim=2)>(F7.3,1x,"&",1x,\),"\\",//),<size(A,dim=2)>(F7.3,1X,"&",1x))
end subroutine tex_matrix_write_sp

subroutine tex_matrix_write_vector_dp(x,n,model)
        implicit none
        real(kind=dp),intent(in)::x(:)
        integer,intent(in)::n
        character(len=*),optional,intent(in)::model

        write(unit=n,fmt=*) "$$\begin{pmatrix}"

        if(present(model) .AND. model .EQ."row") then
                write(unit=n,fmt=100)x!输出成行向量
        else if(present(model) .AND. model .EQ."column") then
                write(unit=n,fmt=200)x!输出成列向量
        else if(.NOT.present(model)) then
                write(unit=n,fmt=200)x!默认输出成列向量
        end if

        write(unit=n,fmt=*) "\end{pmatrix} $$"
        return
100     format(<size(x)-1>(F14.4,"&"),(F14.4))
200     format(<size(x)-1>(F14.4,"\\"),(F14.4))
end subroutine tex_matrix_write_vector_dp

subroutine tex_matrix_write_vector_sp(x,n,model)
        implicit none
        real(kind=sp),intent(in)::x(:)
        integer,intent(in)::n
        character(len=*),optional,intent(in)::model

        write(unit=n,fmt=*) "$$\begin{pmatrix}"

        if(present(model) .AND. model .EQ."row") then
                write(unit=n,fmt=100)x!输出成行向量
        else if(present(model) .AND. model .EQ."column") then
                write(unit=n,fmt=200)x!输出成列向量
        else if(.NOT.present(model)) then
                write(unit=n,fmt=200)x!默认输出成列向量
        end if

        write(unit=n,fmt=*) "\end{pmatrix} $$"
        return
100     format(<size(x)-1>(F7.3,"&"),(F7.3))
200     format(<size(x)-1>(F7.3,"\\"),(F7.3))
end subroutine tex_matrix_write_vector_sp

end module
