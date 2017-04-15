module initialMod
contains
        subroutine initializeGrid(master_grid, grid_row, grid_col, x_scale, y_scale, x_min, y_min)
                implicit none
                integer, intent(in) :: grid_row, grid_col 
                double precision, intent(in) :: x_scale, y_scale, x_min, y_min
                double precision, dimension(grid_row, grid_col), intent(inout) :: master_grid
                integer :: i, j

                

                do i=0, grid_row-1
                        do j=0, grid_col-1
                                master_grid(i+1,j+1)=(DBLE(j*x_scale)+x_min)+sin((DBLE(i*y_scale)+y_min))
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

        function initialCondition(x,y)
                implicit none
                integer, intent(in) :: x,y
                integer :: initialCondition

                initialCondition=x*y
        end function initialCondition

        function stepFunction(v, vineg, vipos, vjneg, vjpos, t_step, x_scale, y_scale)
                implicit none
                double precision, intent(in):: v, vineg, vipos, vjneg, vjpos, t_step, x_scale, y_scale
                double precision :: stepFunction
                stepFunction=v+t_step*(((vipos-2.0D0*v+vineg)/x_scale)+((vjpos-2.0D0*v+vjneg)/y_scale))
        end function stepFunction


end module initialMod
