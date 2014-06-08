      SUBROUTINE RANK(array,length)
      IMPLICIT NONE
      INTEGER :: array(length)
      INTEGER :: length
      INTEGER :: i,j,k
      j=length-1
      do while ( j >0)
         do i=1,j
            if( array(i) > array(i+1))then
               k=array(i)
               array(i)=array(i+1)
               array(i+1)=k
            end if
         end do
         j=j-1
      end do
      return
      END SUBROUTINE