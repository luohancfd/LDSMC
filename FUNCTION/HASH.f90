      function hash(int1,int2)
   !get hash code, int1 MUST LITTLE THAN in2
      implicit none
      character(len=20) :: str
      integer :: hash,int1,int2
      integer :: i

      hash = 5381
      write(str,"(2i10)") int1,int2
      do i=1,len(str)
         hash = (ishft(hash,5) + hash) + ichar(str(i:i))
      end do
      end function