   MODULE BEDGE_CHAIN
   IMPLICIT NONE  
   
   TYPE :: bedge_link
      INTEGER :: I1      !cell code
      INTEGER :: I2      !edge code
      INTEGER :: I3      !type
      INTEGER :: I4      !boundary code
      TYPE(bedge_link),POINTER :: NEXT
   END TYPE
   
   contains 
   
   SUBROUTINE DEL_BEDGE_CHAIN(head)
   TYPE(bedge_link),POINTER :: head,action1,action2
   action1=>head
   DO WHILE(associated(action1))
      action2=>action1
      action1=>action1%NEXT
      deallocate(action2)
   END DO
   END SUBROUTINE
   
   END MODULE