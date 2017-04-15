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
                                master_grid(i+1,j+1)=(DBLE(j*x_scale)+x_min)+(DBLE(i*y_scale)+y_min)
                        enddo
                enddo
        end subroutine initializeGrid
       
        subroutine doStep(master_grid, new_grid, grid_row, grid_col)
                implicit none
                integer :: i,j, grid_row, grid_col
                double precision, dimension(grid_row, grid_col), intent(in) :: master_grid
                double precision, dimension(grid_row, grid_col), intent(inout) :: new_grid
                
                do i=2, grid_row-1
                        do j=2, grid_col-1
                                new_grid(i,j)=stepFunction(master_grid(i,j), master_grid(i-1,j), master_grid(i+1,j) &
                                        , master_grid(i,j-1), master_grid(i,j+1))
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
                stepFunction=v+t_step*(((vipos-2*v+vineg)/x_scale)+((vjpos-2*v+vjneg)/y_scale))
        end function stepFunction


end module initialMod
