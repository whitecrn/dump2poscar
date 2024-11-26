program lammps2vasp
        implicit none

        character(len=15) :: file_n,filename,write_file
        character(len=6) :: step
        character(len=25) :: string
        integer :: n_atom
        logical :: alive1
        integer :: alive2
        integer :: i,j,k,stat,type_atom,ii,n,nn
        integer :: n_li,n_o
        real :: x(5000),y(5000),z(5000)
        real:: xlo_bound,xhi_bound,ylo_bound,yhi_bound,zlo_bound,zhi_bound,xlo,xhi,ylo,yhi,zlo,zhi
        real :: a,b,c,xy,xz,yz
!-------------------------open file----------------------------
do k = 1,30
        filename = ''
        file_n = ''
        alive1 = .false.
        alive2 = 0
        nn = 1
        write(file_n,'(I4)') k
        !write(*,*) k
        filename = 'dump' // trim(adjustl(file_n)) // '.atom'
        !write(*,*) filename
        inquire(file = filename,exist = alive1)
                
        if (alive1) then   !aiming at ensuring if the file exists.
               open(unit = 100+k,file = filename,status = 'old',action = 'read')
               write(*,*) k,'file exists, enter the loop'
               do while (alive2 == 0) !aiming at ensuring if it's the end of the file.
                !write(*,*) 'star to read the whole file',k
                        n_atom = 0
                        xlo_bound = 0.0
                        xhi_bound = 0.0
                        ylo_bound = 0.0
                        yhi_bound = 0.0
                        zlo_bound = 0.0
                        zhi_bound = 0.0
                        n_li = 0
                        n_o = 0
                        xlo = 0.0
                        xhi = 0.0
                        ylo = 0.0
                        yhi = 0.0
                        zlo = 0.0
                        zhi = 0.0
                        xy = 0.0
                        xz = 0.0
                        yz = 0.0
                        i = 0
                        do i = 1,3
                                read(100+k,*,iostat = alive2) 
                                if (alive2 /= 0) then
                                        goto 55
                                end if

                        end do

                        
                        read(100+k,'(I3)') n_atom

                        !allocate(x(n_atom))
                        !allocate(y(n_atom))
                        !allocate(z(n_atom))
                        x(:) = 0.0
                        y(:) = 0.0
                        z(:) = 0.0

                        read(100+k,"(A25)") string
                        if (string(18:19) == 'xy') then
                                read(100+k,*) xlo_bound,xhi_bound,xy
                                read(100+k,*) ylo_bound,yhi_bound,xz
                                read(100+k,*) zlo_bound,zhi_bound,yz
                        else if (string(18:19) == 'pp') then
                                read(100+k,*) xlo_bound,xhi_bound
                                read(100+k,*) ylo_bound,yhi_bound
                                read(100+k,*) zlo_bound,zhi_bound
                                xy = 0.0
                                xz = 0.0
                                yz = 0.0
                        end if

                        xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
                        xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
                        ylo = ylo_bound - min(0.0,yz)
                        yhi = yhi_bound - max(0.0,yz)
                        zlo = zlo_bound
                        zhi = zhi_bound
                        a = xhi-xlo
                        b = yhi-ylo
                        c = zhi-zlo

                        read(100+k,'(A)')
                                
                        do j = 1,n_atom
                                read(100+k,*,iostat = alive2) n,type_atom,x(n),y(n),z(n)
                               
                                if (type_atom == 1) then
                                        n_li = n_li + 1
                                else if (type_atom == 4) then
                                        n_o = n_o + 1
                                end if
                        end do
                        write(step,'(I5)') nn
                        write_file = 'POSCAR' // '-' // trim(adjustl(step)) // '-' // trim(adjustl(file_n))
                        open(unit = 10000+k,file = write_file,status = 'replace')
                        !write(*,*) 'start to write file:',n
                        write(10000+k,*) "POSCAR file written by Z.F.Wu"
                        write(10000+k,'(F3.1)') 1.0
                        write(10000+k,"(F15.12,1X,F15.12,1X,F15.12)") a,0.0,0.0
                        write(10000+k,"(F15.12,1X,F15.12,1X,F15.12)") xy,b,0.0
                        write(10000+k,"(F15.12,1X,F15.12,1X,F15.12)") xz,yz,c
                        write(10000+k,"(A4,5X,A4)") 'Li','O'
                        write(10000+k,"(I4,5X,I4)") n_li,n_o
                        write(10000+k,"(A9)") 'Cartesian'
                        do ii = 1,n_atom
                                write(10000+k,"(F12.8,1X,F12.8,1X,F12.8)") x(ii)-xlo,y(ii)-ylo,z(ii)-zlo
                        end do
                        close(10000+k)
                        nn = nn+1
                        !deallocate(x(:))
                        !deallocate(y(:))
                        !deallocate(z(:))
                end do
               close(100+k) 
        end if
55      cycle        
end do 
stop
end program
