   MODULE BIN_TREE
   IMPLICIT NONE

   TYPE :: cell_link
      INTEGER :: C  !cell code
      TYPE(cell_link),pointer :: next
   END TYPE

   TYPE :: tree
      INTEGER :: h  !hash code
      INTEGER :: l  !cell_link length
      TYPE(cell_link),POINTER :: head
      TYPE(cell_link),POINTER :: p
      TYPE(tree),POINTER :: left,right
   END TYPE
   
   
   
   contains
   !SUBROUTINE: search the shared edge 
   RECURSIVE SUBROUTINE SCAN_TREE_SEDGE(treenode)
   IMPLICIT NONE
   TYPE(tree),POINTER :: treenode
   !this is a subroutine to scan tree and build ICELL(9:11,*)
   IF( associated(treenode))THEN
      CALL SCAN_TREE_SEDGE(treenode%left)
      CALL SEARCH_SHARED_EDGE(treenode)
      CALL SCAN_TREE_SEDGE(treenode%right)
   END IF
   RETURN
   END SUBROUTINE SCAN_TREE_SEDGE
   
   SUBROUTINE SEARCH_SHARED_EDGE(treenode)
   USE CELLINFO
   IMPLICIT NONE
   TYPE(tree),POINTER :: treenode
   INTEGER :: I,J,K,L1,L2,L3,L4
   INTEGER,ALLOCATABLE :: hashcode(:,:)
   INTEGER,EXTERNAL :: hash

   IF( treenode%l .gt. 1)then
      allocate(hashcode(3,treenode%l))
      treenode%p=>treenode%head
      DO I=1,treenode%l
         DO J=1,3
            K=J+1
            IF(K .EQ. 4)  K=1
            IF(ICELL(J,treenode%p%C) .LT. ICELL(K,treenode%p%C))then
               L1=ICELL(J,treenode%p%C); L2=ICELL(K,treenode%p%C)
            ELSE
               L1=ICELL(K,treenode%p%C); L2=ICELL(J,treenode%p%C)
            END IF
            IF(hash(L1,L2) .EQ. treenode%h)THEN
               hashcode(1,I)=J
               hashcode(2,I)=K
               hashcode(3,I)=treenode%p%C
               EXIT
            END IF
         END DO
         IF(I .NE. treenode%l)THEN
            treenode%p=>treenode%p%next
         END IF
      END DO

      DO I=1,treenode%l
         IF(ICELL(8+hashcode(1,I),hashcode(3,I)) .eq. 0)THEN
            DO J=I+1,treenode%l
               IF(ICELL(8+hashcode(1,J),hashcode(3,J)) .eq. 0)THEN
                  L1=ICELL(hashcode(1,I),hashcode(3,I))
                  L2=ICELL(hashcode(2,I),hashcode(3,I))
                  L3=ICELL(hashcode(2,J),hashcode(3,J))
                  L4=ICELL(hashcode(1,J),hashcode(3,J))
                  IF((L1 .EQ. L3) .AND. (L2 .EQ. L4))THEN
                     ICELL(8+hashcode(1,I),hashcode(3,I))=hashcode(3,J)
                     ICELL(8+hashcode(1,J),hashcode(3,J))=hashcode(3,I)
                     EXIT
                  END IF
               END IF
            END DO
         END IF
      END DO
   END IF
   
      
            

      
   
  
   
   
   
   
   
   RETURN
   END SUBROUTINE SEARCH_SHARED_EDGE
   
   
   !SUBROUTINE: delete the tree
   RECURSIVE SUBROUTINE DEL_TREE(treenode)
   IMPLICIT NONE
   TYPE(tree),POINTER :: treenode
   IF(associated(treenode))THEN
      CALL DEL_TREE(treenode%left)
      CALL DEL_TREE(treenode%right)
      CALL DEL_LEAVE(treenode)
   END IF
   RETURN
   END SUBROUTINE DEL_TREE
   
   SUBROUTINE DEL_LEAVE(treenode)
   TYPE(tree),POINTER :: treenode
   CALL DEL_CELL_LINK(treenode%head)
   deallocate(treenode)
   RETURN
   END SUBROUTINE
   
   SUBROUTINE DEL_CELL_LINK(head)
   IMPLICIT NONE
   TYPE(cell_link),pointer :: head,action1,action2
   INTEGER :: l,i
   action1=>head
   DO WHILE(associated(action1))
      action2=>action1
      action1=>action1%next
      deallocate(action2)
   END DO
   nullify(head,action1,action2)
   END SUBROUTINE DEL_CELL_LINK
   !SUBROUTINE: search the boundary edge
   SUBROUTINE SEARCH_BOUND_EDGE(treenode,K,N1,N2)
   USE CELLINFO
   IMPLICIT NONE
   TYPE(tree),POINTER :: treenode
   INTEGER :: K,N1,N2,bhash,I,II
   INTEGER,EXTERNAL ::hash,ADD_ONE,MINUS_ONE
   !K is the cell number, N1 is the first node, N2 is the second node
   IF(N1 .GT. N2)THEN
      K=N1
      N1=N2
      N2=K
   END IF
   bhash=hash(N1,N2)
   DO WHILE( .true. )
      IF(bhash < treenode%h)THEN
         treenode=>treenode%left
      ELSE IF( bhash > treenode%h)THEN
         treenode=>treenode%right
      ELSE
         IF(treenode%l .eq. 1)THEN
            K=treenode%head%C
            RETURN
         ELSE
            treenode%p=>treenode%head
            DO WHILE( associated(treenode%p))
               II=0
               DO I=1,3
                  IF( ICELL(I,treenode%p%C) .eq. N1)THEN
                     IF( ICELL(ADD_ONE(I),treenode%p%C) .eq. N2)THEN
                        II=1
                        K=treenode%p%C
                     ELSE IF (ICELL(MINUS_ONE(I),treenode%p%C) .eq. N2)THEN
                        II=1
                        K=treenode%p%C
                     END IF
                  END IF
               END DO
               IF( II .EQ. 1)    RETURN
               treenode%p=>treenode%p%next
            END DO
         END IF
      END IF
   END DO

   END SUBROUTINE SEARCH_BOUND_EDGE

   END MODULE BIN_TREE
