module initialMod
contains
        subroutine initializeGrid(master_grid)
                implicit none
                double precision, allocatable, dimension(:,:), intent(inout) :: master_grid
                integer :: i, j
                allocate(master_grid(10,10))
                

                do i=1, size(master_grid, 1)
                        do j=1, size(master_grid, 2)
                                if(i .ne. 5 .and. j .ne. 5) then 
                                        master_grid(i,j)=0
                                else
                                        master_grid(i,j)=100
                                endif
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

        function stepFunction(v, vineg, vipos, vjneg, vjpos)
                implicit none
                double precision, intent(in):: v, vineg, vipos, vjneg, vjpos
                integer :: stepFunction
                stepFunction=v+(vipos-2*v+vineg)+(vjpos-2*v+vjneg)
        end function stepFunction


end module initialMod
