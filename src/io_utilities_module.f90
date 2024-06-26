!Fortran IO utilities module 
!Max Wood - mw16116@bristol.ac.uk 
!Univeristy of Bristol - Department of Aerospace Engineering

!Version 2.2
!Updated 24-03-2024

!Module 
module io_utilities
use flow_data_mod
contains 


!Check if file exists function =========================
function file_exists(filename) result(exist_state)
implicit none 

!Variables - Result
logical :: exist_state

!Variables - Import 
character(*), intent(in) :: filename

!Check for existance
inquire(file=filename,exist=exist_state)
return 
end function file_exists




!Function to get integer command argument =========================
function get_command_argument_n_int(arg_idx) result(argint)
implicit none 

!Variables - Import 
integer(in) :: argint
integer(in32) :: arg_idx

!Variables - Local 
! integer(in32) :: arglen,argstat
character(len=:), allocatable :: argtemp

!Read argument 
argtemp = get_command_argument_n_str(arg_idx)

!Convert to integer
argint = str2int(argtemp)
return 
end function get_command_argument_n_int




!Function to get real command argument =========================
function get_command_argument_n_real(arg_idx) result(argreal)
implicit none 

!Variables - Import 
real(dp) :: argreal
integer(in32) :: arg_idx

!Variables - Local 
! integer(in32) :: arglen,argstat
character(len=:), allocatable :: argtemp

!Read argument 
argtemp = get_command_argument_n_str(arg_idx)

!Convert to integer
argreal = str2real(argtemp)
return 
end function get_command_argument_n_real




!Function to get string command argument =========================
function get_command_argument_n_str(arg_idx) result(argstr)
implicit none 

!Result 
character(len=:), allocatable :: argstr

!Variables - Import 
integer(in32) :: arg_idx

!Variables - Local 
integer(in32) :: arglen,argstat

!Read argument 
call get_command_argument(number=arg_idx, length=arglen)
allocate(character(len=arglen) :: argstr)
call get_command_argument(number=arg_idx, value=argstr, status=argstat)
return 
end function get_command_argument_n_str




!Set integer option subroutine =========================
subroutine set_int_opt(opt_int,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
integer(in) :: opt_int
character(*), intent(in) :: str_tag

!Variables - Local
integer(in) :: opt_int_read

!Scan for option real 
opt_int_read = scan_opt_int(fh,str_tag)

!If item has been found and read then set this value else leave option as the default
if (opt_int_read .NE. -huge(0_in)) then 
    opt_int = opt_int_read
end if 
return 
end subroutine set_int_opt




!Set real option subroutine =========================
subroutine set_real_opt(opt_real,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
real(dp) :: opt_real
character(*), intent(in) :: str_tag

!Variables - Local
real(dp) :: opt_real_read

!Scan for option real 
opt_real_read = scan_opt_real(fh,str_tag)

!If item has been found and read then set this value else leave option as the default
if (.NOT.isnan(opt_real_read)) then 
    opt_real = opt_real_read
end if 
return 
end subroutine set_real_opt




!Set string option subroutine =========================
subroutine set_str_opt(opt_str,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
character(*), intent(in) :: str_tag
character(len=:), allocatable :: opt_str

!Variables - Local
character(len=:), allocatable :: opt_str_read

!Scan for option string 
opt_str_read = scan_opt_str(fh,str_tag)

!If item has been found and read then set this value else leave option as the default
if (opt_str_read .NE. 'na') then 
    if (allocated(opt_str)) then 
        deallocate(opt_str)
    end if 
    allocate(character(len=len_trim(opt_str_read)) :: opt_str)
    opt_str = opt_str_read
end if 
return 
end subroutine set_str_opt




!Set logical option subroutine =========================
subroutine set_log_opt(opt_log,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
logical :: opt_log
character(*), intent(in) :: str_tag

!Variables - Local
character(len=:), allocatable :: opt_str

!Scan for option string 
opt_str = scan_opt_str(fh,str_tag)

!Set logical option 
if (opt_str .NE. 'na') then
    if (opt_str == 'yes') then 
        opt_log = .true.
    else 
        opt_log = .false.
    end if 
end if 
return 
end subroutine set_log_opt




!Set integer array option subroutine =========================
subroutine set_int_opt_arr1(opt_int,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
integer(in), dimension(:), allocatable :: opt_int
character(*), intent(in) :: str_tag

!Variables - Local
integer(in) :: ii,jj
integer(in) :: arrlen,rs,re
character(len=:), allocatable :: opt_str_read

!Scan for option real
opt_str_read = scan_opt_str(fh,str_tag) 

!If item has been found and read then set this value else leave option as the default
if (opt_str_read .NE. 'na') then 

    !Determine length 
    arrlen = 0 
    do ii=1,len_trim(opt_str_read) !count number of white spaces
        if (opt_str_read(ii:ii) == ' ') then 
            arrlen = arrlen + 1
        end if 
    end do 
    arrlen = arrlen + 1 !set to number of items seperated by white spaces

    !Allocate option array 
    if (allocated(opt_int)) then 
        deallocate(opt_int)
    end if 
    allocate(opt_int(arrlen))
    
    !Write items to the array 
    rs = 1
    do ii=1,arrlen

        !Find next white space or end 
        re = 0
        do jj=rs+1,len_trim(opt_str_read)
            if (opt_str_read(jj:jj) == ' ') then 
                re = jj 
                exit 
            elseif (jj == len_trim(opt_str_read)) then 
                re = jj 
            end if 
        end do 
        
        !Write current item into the option array 
        opt_int(ii) = str2int(opt_str_read(rs:re))

        !Update current starting position 
        rs = re 
    end do 
end if 
return 
end subroutine set_int_opt_arr1




!Set real array option subroutine =========================
subroutine set_real_opt_arr1(opt_real,fh,str_tag)
implicit none

!Variables - Import
integer(in32) :: fh
real(dp), dimension(:), allocatable :: opt_real
character(*), intent(in) :: str_tag

!Variables - Local
integer(in) :: ii,jj
integer(in) :: arrlen,rs,re
character(len=:), allocatable :: opt_str_read

!Scan for option real
opt_str_read = scan_opt_str(fh,str_tag) 

!If item has been found and read then set this value else leave option as the default
if (opt_str_read .NE. 'na') then 

    !Determine length 
    arrlen = 0 
    do ii=1,len_trim(opt_str_read) !count number of white spaces
        if (opt_str_read(ii:ii) == ' ') then 
            arrlen = arrlen + 1
        end if 
    end do 
    arrlen = arrlen + 1 !set to number of items seperated by white spaces

    !Allocate option array 
    if (allocated(opt_real)) then 
        deallocate(opt_real)
    end if 
    allocate(opt_real(arrlen))
    
    !Write items to the array 
    rs = 1
    do ii=1,arrlen

        !Find next white space or end 
        re = 0
        do jj=rs+1,len_trim(opt_str_read)
            if (opt_str_read(jj:jj) == ' ') then 
                re = jj 
                exit 
            elseif (jj == len_trim(opt_str_read)) then 
                re = jj 
            end if 
        end do 
        
        !Write current item into the option array 
        opt_real(ii) = str2real(opt_str_read(rs:re))

        !Update current starting position 
        rs = re 
    end do 
end if 
return 
end subroutine set_real_opt_arr1




!Scan for and import integer valued option following tag 'str_tag' function =========================
function scan_opt_int(fh,str_tag) result(opt_int)
implicit none

!Variables - Result
integer(in) :: opt_int

!Variables - Import
integer(in32) :: fh
character(*), intent(in) :: str_tag

!Variables - Local
character(len=:), allocatable :: opt_str

!Read option as string value 
opt_str = scan_opt_str(fh,str_tag)
if (opt_str == 'na') then !item not found 
    opt_int = -huge(0_in)
    return 
end if 

!Write to integer value 
opt_int = str2int(opt_str)
return 
end function scan_opt_int




!Scan for and import real valued option following tag 'str_tag' function =========================
function scan_opt_real(fh,str_tag) result(opt_real)
use ieee_arithmetic, only : ieee_value,IEEE_QUIET_NAN
implicit none

!Variables - Result
real(dp) :: opt_real

!Variables - Import
integer(in32) :: fh
character(*), intent(in) :: str_tag

!Variables - Local
character(len=:), allocatable :: opt_str

!Read option as string value 
opt_str = scan_opt_str(fh,str_tag)
if (opt_str == 'na') then !item not found 
    opt_real = ieee_value(1.0d0,IEEE_QUIET_NAN)
    return 
end if 

!Write to real value 
opt_real = str2real(opt_str)
return 
end function scan_opt_real




!Scan for and import string valued option following tag 'str_tag' function =========================
function scan_opt_str(fh,str_tag) result(opt_str)
implicit none

!Variables - Result
character(len=:), allocatable :: opt_str

!Variables - Import
integer(in32) :: fh
character(*), intent(in) :: str_tag

!Variables - Local
integer(in32) :: ii 
integer(in32) :: iostatus,len_tag,read_init_pos,item_found
character(len=1000) :: rtemp 

!Set tag length 
len_tag = len_trim(str_tag)

!Rewind file 
rewind(fh)

!Scan for target string tag
iostatus = 0 
item_found = 0 
do while (iostatus == 0)
    read(fh,'(A)',iostat=iostatus) rtemp
    if ((rtemp(1:len_tag) == str_tag) .AND. (rtemp(len_tag+1:len_tag+1) == ' ')) then !tag has been located and is the full tag

        !Locate position of the end of '= ' to begin reading the option string value 
        read_init_pos = -1
        do ii=1,len_trim(rtemp)
            if (rtemp(ii:ii+1) == '= ') then 
                read_init_pos = ii+2
                exit 
            end if 
        end do 
        if (read_init_pos == -1) then !input is the wrong format and the '= ' tag has not been found 
            write(*,'(A)') '** could not locate "= ", the option input following tag "'//str_tag//'" is of incorrect format' 
            opt_str = 'na'
            return 
        end if 

        !Set item_found flag
        item_found = 1

        !Read the option string into the return variable and exit search
        opt_str = rtemp(read_init_pos:len_trim(rtemp))
        exit 
    end if 
end do 

!If the target string has not been found then return na
if (item_found == 0) then 
    opt_str = 'na'
end if 
return 
end function scan_opt_str




!Convert array to F0.X format with leading zero subroutine =========================
subroutine array_real2F0_Xstring(array_str,array,X) 
implicit none 

!Variables - Import 
character(*), dimension(:) :: array_str
integer(in) :: X
real(dp), dimension(:,:) :: array

!Variables - Local
integer(in) :: ii,jj
integer(in) :: Nrep
character(len=2*X) :: frmtI,str
character(len=4*X) :: frmt

!Set number of repetitions 
Nrep = size(array,2)

!Construct format descriptor 
write(frmtI,'(I0)') X
frmtI = 'F0.'//trim(frmtI)
frmt = '('//trim(frmtI)//')'

!Convert array
do ii=1,size(array,1)
    do jj=1,Nrep
        write(str,trim(frmt)) array(ii,jj)
        if (str(1:1) == '.') then 
            str = '0'//trim(str)
        elseif (str(1:2) == '-.') then 
            str = '-0.'//str(3:len_trim(str))
        end if
        if (jj == 1) then 
            array_str(ii) = trim(str)
        else
            array_str(ii) = trim(array_str(ii))//' '//trim(str)
        end if 
    end do 
end do 
return 
end subroutine array_real2F0_Xstring




!F0.X format with leading zero function =========================
function real2F0_Xstring(val,X) result(str)
implicit none 

!Result 
character(len=:), allocatable :: str

!Variables - Import 
character(len=10) :: frmtI
character(len=20), allocatable :: frmt
integer(in) :: X,len_str
real(dp) :: val

!Construct format descriptor 
write(frmtI,'(I0)') X
frmt = '(F0.'//trim(frmtI)//')'

!Format
if (isnan(val)) then 
    str = 'nan'
else

    !Find length of result string 
    if (abs(val) .LT. 1.0d0) then 
        len_str = 1
    else
        len_str = floor(log10(abs(val))) + 1
    end if 
    len_str = len_str + X + 2

    !Write to return string 
    allocate(character(len=len_str) :: str)
    write(str,trim(frmt)) val
    str = trim(str)

    !Assign leading zero if required
    if (str(1:1) == '.') then 
        str = '0'//trim(str)
    elseif (str(1:2) == '-.') then 
        len_str = len_trim(str)
        str = '-0.'//str(3:len_str)
    end if 
end if 
return 
end function real2F0_Xstring




!Integer to string function =========================
function int2str(val) result(stri)
implicit none 

!Result 
character(len=:), allocatable :: stri

!Variables - Import 
integer(in) :: val

!Write
allocate(character(len=100) :: stri)
write(stri,'(I0)') val
stri = stri(1:len_trim(stri))
return 
end function int2str




!String to integer function =========================
function str2int(str) result(val)
implicit none 

!Result 
integer(in) :: val

!Variables - Import 
character(*), intent(in) :: str

!Read to integer
read(str,*) val
return 
end function str2int




!String to real function =========================
function str2real(str) result(val)
implicit none 

!Result 
real(dp) :: val

!Variables - Import 
character(*), intent(in) :: str

!Read to integer
read(str,*) val
return 
end function str2real


end module io_utilities