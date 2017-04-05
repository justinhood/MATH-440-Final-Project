module initialMod
contains
        subroutine initializeGrid(master_grid)
                implicit none
                double precision, allocatable, dimension(:,:), intent(inout) :: master_grid
                integer :: i, j
                allocate(master_grid(10,10))
                

                do i=1, size(master_grid, 1)
                        do j=1, size(master_grid, 2)
                                master_grid(i,j)=initialCondition(i,j)
                        enddo
                enddo
        end subroutine initializeGrid
        
        function initialCondition(x,y)
                implicit none
                integer, intent(in) :: x,y
                integer :: initialCondition

                initialCondition=x*y
        end function initialCondition


end module initialMod
