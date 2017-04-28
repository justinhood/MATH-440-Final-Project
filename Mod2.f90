module Mod2
contains
        subroutine initializeGrid(master_grid, grid_row, grid_col, x_scale, y_scale, x_min, y_min)
                implicit none
                integer, intent(in) :: grid_row, grid_col 
                double precision, intent(in) :: x_scale, y_scale, x_min, y_min
                double precision, dimension(grid_row, grid_col), intent(inout) :: master_grid
                integer :: i, j

                

                do i=1, grid_row
                        do j=1, grid_col
                                if(i .Eq. 3 .and. j .eq. 3) then
                                        master_grid(i,j)=10.0D0
                                else
                                        master_grid(i,j)=0.0D0
                                endif
                        enddo
                enddo
        end subroutine initializeGrid
       
        subroutine doStep(concat_grid, c_grid_row, c_grid_col, t_step, x_scale, y_scale)
                implicit none
                integer :: i,j
                integer, intent(in) :: c_grid_row, c_grid_col
                double precision, intent(in) :: t_step, x_scale, y_scale
                double precision, dimension(c_grid_row, c_grid_col), intent(inout) :: concat_grid
                double precision, dimension(c_grid_row, c_grid_col) :: temp

                temp=concat_grid
                
                do i=2,c_grid_row-1
                        do j=2, c_grid_col-1
                                concat_grid(i,j)=stepFunction(temp(i,j), temp(i-1,j), temp(i+1,j) &
                                        , temp(i,j-1), temp(i,j+1), t_step, x_scale, y_scale)
                        enddo
                enddo

        end subroutine doStep

        function stepFunction(v, vineg, vipos, vjneg, vjpos, t_step, x_scale, y_scale)
                implicit none
                double precision, intent(in):: v, vineg, vipos, vjneg, vjpos, t_step, x_scale, y_scale
                double precision :: stepFunction
                stepFunction=DBLE(v+t_step*(((vipos-2.0D0*v+vineg)/x_scale)+((vjpos-2.0D0*v+vjneg)/y_scale)))
        end function stepFunction


        subroutine doImplicitStep(uold_con, u_con, unew, x_scale, y_scale, tstep, alpha, beta,&
                        c_grid_row, c_grid_col, rank, temp)
                implicit none
                integer :: i,j
                integer, intent(in) :: rank, c_grid_row, c_grid_col
                double precision, intent(in) :: alpha, beta, x_scale, y_scale, tstep
                double precision, dimension(:,:), intent(in) :: uold_con, u_con
		double precision, dimension(c_grid_row,c_grid_col), intent(inout) :: temp
                double precision, dimension(:,:), intent(in) :: unew
                double precision :: uxx, uyy, b

                do i=2, c_grid_row-1
                        do j=2, c_grid_col-1
                                uxx=(uold_con(i+1, j)-2.0D0*uold_con(i,j)+uold_con(i-1,j))/(x_scale**2)
                                uyy=(uold_con(i, j+1)-2.0D0*uold_con(i,j)+uold_con(i,j-1))/(y_scale**2)
                                b=uold_con(i,j)+0.5D0*tstep*(uxx+uyy)
                                temp(i,j)=b-(alpha*u_con(i-1,j)+alpha*u_con(i+1,j)+beta*u_con(i,j+1)&
                                                +beta*u_con(i,j-1))/(1+2*(alpha+beta))
                               
                        enddo
                enddo
        end subroutine doImplicitStep

end module Mod2
