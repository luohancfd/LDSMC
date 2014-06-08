      SUBROUTINE PREPROCESSING
!  pre-processing for DSMC
         USE BIN_TREE
         USE BEDGE_CHAIN
         USE NODEINFO
         USE CELLINFO
         USE CALC
         USE GAS
         USE boundinfo
         USE IFPORT  !change to USE DFPORT for CVF compiler
         IMPLICIT NONE
         
         include 'cgnslib_f90.h'
         integer :: index_file,index_base,index_zone,index_section
         integer :: size_zone(3),type_cell
         real*8 :: x(3),y(3),temp              !used to calculate area
         REAL*8,EXTERNAL :: VEC_PRODUCT
         INTEGER,EXTERNAL :: ADD_ONE,MINUS_ONE
! the following user defined data structure is stored in COMPLEXDATA directory
         TYPE(bedge_link),POINTER :: P,HEAD
         TYPE(tree),POINTER :: treeroot,treep,newleave
         TYPE(cell_link),POINTER :: cellp
         INTEGER :: nodecode(3),nodehash(3),ii,jj,numleave,NWALL_BAK
         INTEGER,EXTERNAL :: hash

!   useless,temp use
         integer :: istart,iend,nbndry,iparent_flag,nbocos
         integer :: i,j,k,l,m,n,j1,j2,j3,j4 !j1-j4 count number for boundary edge
         real*8 :: a,b,c
         character*30 :: text
         character*22 :: num2text !this is used for increase number precision
         integer :: ier
         
         
!   create directory to store result
         ier=SYSTEM("md result") !for linux, md->mkdir
         ier=SYSTEM("COPY grid.cgns .\result\flowfield.cgns") !for linux COPY->CP
         
!   open file
         call cg_open_f('grid.cgns',CG_MODE_READ,index_file,ier)
         if (ier .ne. 0) call cg_error_exit_f
         
!   we know there is only one base, once needed, could use cg_nbase_f to check
!   an example is call cg_nbases_f(index_file,base,ier)
         index_base=1
!   we know there is only one zone, once needed, could use cg_nzone_f to check
         index_zone=1
!   we know there is only one zone, once needed, could use cg_nzone_f to check
         index_section=1
!   get the size of grid and provide node storage
         
         
         call cg_zone_read_f(index_file,index_base,index_zone,text,size_zone,ier)
         NCELLS=size_zone(2)
         NNODES=size_zone(1)
         NCCELLS=NCELLS*4
         ALLOCATE(NODES(2,NNODES),CELL(1,NCELLS),IS(4,NCELLS))         
         ALLOCATE (CCELL(5,NCCELLS),ICCELL(3,NCCELLS),ICELL(11,NCELLS))
         ICELL(:,:)=0
!  get node coordinate
         call cg_coord_read_f(index_file,index_base,index_zone,'CoordinateX',RealDouble,1,NNODES,NODES(1,:),ier)
         call cg_coord_read_f(index_file,index_base,index_zone,'CoordinateY',RealDouble,1,NNODES,NODES(2,:),ier)
         do i=1,NNODES
            do j=1,2
               WRITE(num2text,"(E22.15E3)") NODES(j,i)
               READ(num2text,"(E22.15E3)") NODES(j,i)
            end do
         end do
!   get the name of section  "text",element type"itype",start\end of element
         call cg_section_read_f(index_file,index_base,index_zone,index_section,text,type_cell,istart,iend,nbndry,iparent_flag,ier)        
         if (ElementTypeName(type_cell) .ne. 'TRI_3')then
            write(*,*) "ERROR:Mesh is not compatible!!STOP"
            STOP
         end if      
         
!   get the cells composition
         call cg_elements_read_f(index_file,index_base,index_zone,index_section,ICELL(1:3,:),NULL,ier)
         
!   get the area of cell 
         MAXCELL(:)=0.0D00
         MINCELL(:)=10.0D10
         numleave=0
         do i=1,NCELLS
            do j=1,3
               x(j)=NODES(1,ICELL(j,i))
               y(j)=NODES(2,ICELL(j,i))
            end do
            CELL(1,i)=0.5D00*VEC_PRODUCT(x(3)-x(1),y(3)-y(1),x(2)-x(1),y(2)-y(1))
            a=DSQRT((x(1)-x(2))**2.0D00+(y(1)-y(2))**2.0D00)
            b=DSQRT((x(3)-x(2))**2.0D00+(y(3)-y(2))**2.0D00)
            c=DSQRT((x(1)-x(3))**2.0D00+(y(1)-y(3))**2.0D00)
            temp=DMAX1(a,b,c); b=DMIN1(a,b,c)
           !To make the nodes of cell clockwise
            IF (CELL(1,i) < 0 )THEN
               j=ICELL(2,i)
               ICELL(2,i)=ICELL(3,i)
               ICELL(3,i)=j
               CELL(1,i)=-CELL(1,i)
            END IF            
            IF(CELL(1,i) .GT. MAXCELL(1)) MAXCELL(1)=CELL(1,i)
            IF(CELL(1,i) .LT. MINCELL(1)) MINCELL(1)=CELL(1,i)
            IF(temp .GT. MAXCELL(2))      MAXCELL(2)=temp
            IF(b .LT. MINCELL(2))         MINCELL(2)=b
            nodecode(:)=ICELL(1:3,i)
            call rank(nodecode,int(3))
            nodehash(1)=hash(nodecode(1),nodecode(2))
            nodehash(2)=hash(nodecode(1),nodecode(3))
            nodehash(3)=hash(nodecode(2),nodecode(3))
            
            !create the tree
            DO j=1,3
               numleave=numleave+1
               allocate(newleave,stat=ier)
               if( ier .ne. 0)then
                  write(*,*) "The mesh is too large, out of memory"
                  stop
               end if
               newleave%h=nodehash(j)
               nullify(newleave%left,newleave%right)
               if(numleave .eq. 1)then
                  treep=>newleave
                  treeroot=>newleave
                  treep%l=1
                  allocate(treep%head)
                  treep%p=>treep%head
               else
                  jj=1
                  treep=>treeroot
                  do while( jj .eq. 1)
                     if( nodehash(j) > treep%h)then
                        if(associated(treep%right))then
                           treep=>treep%right
                        else
                           treep%right=>newleave
                           treep=>newleave
                           jj=-1   
                        end if
                     else if(nodehash(j) < treep%h)then
                        if(associated(treep%left))then
                           treep=>treep%left
                        else
                           treep%left=>newleave
                           treep=>newleave
                           jj=-1
                        end if
                     else
                        deallocate(newleave)
                        jj=0
                     end if
                  end do
                  if( jj .eq. 0)then
                     treep%l=treep%l+1
                     treep%p=>treep%head
                     do while(associated(treep%p%next))
                        treep%p=>treep%p%next
                     end do
                     allocate(cellp)
                     treep%p%next=>cellp
                     treep%p=>cellp
                  else if(jj .eq. -1)then
                     treep%l=1
                     nullify(treep%right,treep%left)
                     allocate(treep%head)
                     treep%p=>treep%head
                  end if
               end if
               treep%p%C=i
               nullify(treep%p%next)
            END DO  
            ICELL(8,i)=4*(i-1)
            DO j=1,4
               l=ICELL(8,i)+j
               CCELL(1,l)=CELL(1,i)/4.0D00
               CCELL(2,l)=0.0D00
               ICCELL(3,l)=i
            END DO
         end do

!   scan the tree and get shared edge information
         treep=>treeroot
         CALL SCAN_TREE_SEDGE(treep) !subroutine in BIN_TREE module
            
         
!---get the number of boundary conditions-------------------------------------------------------------------------------------------------------
         NWALL_BAK=NWALL
         call cg_nbocos_f(index_file,index_base,index_zone,NBOUND,ier)
         allocate(bound(NBOUND),IWBOUND(NBOUND))
         NBEDGES=0; NWALL=0; IWBOUND=0
         ALLOCATE(HEAD)
         P=>HEAD
!   get the type of BC and BC composition
         do i=1,NBOUND
            call cg_boco_info_f(index_file,index_base,index_zone,i,bound(i)%name,&
                &bound(i)%type_bound,bound(i)%type_points,bound(i)%num_points,&
                &bound(i)%index_normal,bound(i)%size_normal,bound(i)%datatype_normal,bound(i)%ndataset,ier) 
            allocate(bound(i)%points(bound(i)%num_points))
            call cg_boco_read_f(index_file,index_base,index_zone,i,bound(i)%points,bound(i)%ndataset,ier)
            SELECT CASE( BCTypeName(bound(i)%type_bound))
            CASE("BCInflow")
               l=1
            CASE("BCOutflow")
               l=2
               !Actually there is no difference between inflow and out flow
               !whether a boundary edge is inflow or outflow depends on the
               !intersection angel between the velocity of stream and the normal vector of the edge
            CASE("BCSymmetryPlane")
               l=3
            CASE("BCWall")
               l=4
               NWALL=NWALL+1
               IWBOUND(i)=NWALL
            CASE DEFAULT
               WRITE(*,*) "ERROR:Unsupported boundary condition"
               STOP
            END SELECT
            do j=1,bound(i)%num_points-1
               treep=>treeroot
               CALL SEARCH_BOUND_EDGE(treep,k,bound(i)%points(j),bound(i)%points(j+1))
               DO m=1,3
                  IF( bound(i)%points(j) .EQ. ICELL(m,k))THEN
                     IF( bound(i)%points(j+1) .EQ. ICELL(ADD_ONE(m),k))THEN
                        ICELL(3+m,k)=l
                        n=m!n is the edge code of sampling cell
                     ELSE IF( bound(i)%points(j+1) .EQ. ICELL(MINUS_ONE(m),k))THEN
                        ICELL(3+MINUS_ONE(m),k)=l
                        n=MINUS_ONE(m)
                     END IF
                     EXIT
                  END IF
               END DO
               ICELL(7,k)=ICELL(7,k)+1
               NBEDGES(1)=NBEDGES(1)+1
               NBEDGES(l+1)=NBEDGES(l+1)+1
               P%I1=k  !cell code
               P%I2=n  !edge code
               P%I3=l  !type code
               P%I4=i  !BC code
               allocate(P%NEXT)
               P=>P%NEXT
            end do
         end do
         NULLIFY(P%NEXT)
         IF( NWALL .NE. NWALL_BAK)THEN
            WRITE(9,*) "The mesh file and init.txt file don't match"
            STOP
         END IF
         ALLOCATE(IBEDGE1(2,NBEDGES(2)),IBEDGE2(2,NBEDGES(3)),IBEDGE3(2,NBEDGES(4)),IBEDGE4(2,NBEDGES(5)))
         ALLOCATE(SSSEG(9,MSP,NBEDGES(5)),SSEG(5,NBEDGES(5)),ENTRYIN(2,NBEDGES(2)),ENTRYOUT(2,NBEDGES(3)))
         ALLOCATE(IWALL(0:NBEDGES(5),NWALL),SURFACEVEC(2,NBEDGES(5)))
         IWALL=0
         P=>HEAD
         j1=0  !overall code for inflow BC, j1= 1, NBEDGES(2)
         j2=0
         j3=0
         j4=0  !overall code for solid surface, j4= 1, NBEDGES(5)
         DO I=1,NBEDGES(1)
            IF( P%I3 .EQ. 4)THEN
               j4=j4+1
               DO L=1,MSP
                  SSSEG(1:9,L,j4)=SSSEG_INT(1:9,L)
               END DO
               IBEDGE4(1,j4)=P%I1
               IBEDGE4(2,j4)=P%I2
               IS(P%I2,P%I1)=j4
               j=IWBOUND(P%I4)
               SSEG(1:5,j4)=SSEG_INT(1:5,j)
               IWALL(0,j)=IWALL(0,j)+1
               IWALL(IWALL(0,j),j)=j4
               x(1)=NODES(1,ICELL(P%I2,P%I1))
               x(2)=NODES(1,ICELL(ADD_ONE(P%I2),P%I1))
               y(1)=NODES(2,ICELL(P%I2,P%I1))
               y(2)=NODES(2,ICELL(ADD_ONE(P%I2),P%I1))
               temp=DSQRT((x(1)-x(2))**2.0D00+(y(1)-y(2))**2.0D00)
               SURFACEVEC(1,j4)=(y(2)-y(1))/temp
               SURFACEVEC(2,j4)=(x(1)-x(2))/temp
      !----actually there is no difference between inflow boundary and outflow boundary
            ELSE IF(P%I3 .EQ. 1)THEN
               j1=j1+1
               IBEDGE1(1,j1)=P%I1
               IBEDGE1(2,j1)=P%I2
               IS(P%I2,P%I1)=j1
               x(1)=NODES(1,ICELL(P%I2,P%I1))
               x(2)=NODES(1,ICELL(ADD_ONE(P%I2),P%I1))
               y(1)=NODES(2,ICELL(P%I2,P%I1))
               y(2)=NODES(2,ICELL(ADD_ONE(P%I2),P%I1))
               temp=DSQRT((x(1)-x(2))**2.0D00+(y(1)-y(2))**2.0D00)
               ENTRYIN(1,j1)=(y(2)-y(1))/temp
               ENTRYIN(2,j1)=(x(1)-x(2))/temp
            ELSE IF(P%I3 .EQ. 2)THEN
               j2=j2+1
               IBEDGE2(1,j2)=P%I1
               IBEDGE2(2,j2)=P%I2
               IS(P%I2,P%I1)=j2
               x(1)=NODES(1,ICELL(IBEDGE2(2,j2),IBEDGE2(1,j2)))
               x(2)=NODES(1,ICELL(ADD_ONE(IBEDGE2(2,j2)),IBEDGE2(1,j2)))
               y(1)=NODES(2,ICELL(IBEDGE2(2,j2),IBEDGE2(1,j2)))
               y(2)=NODES(2,ICELL(ADD_ONE(IBEDGE2(2,j2)),IBEDGE2(1,j2)))
               temp=DSQRT((x(1)-x(2))**2.0D00+(y(1)-y(2))**2.0D00)
               ENTRYOUT(1,j2)=(y(2)-y(1))/temp
               ENTRYOUT(2,j2)=(x(1)-x(2))/temp
            ELSE IF(P%I3 .EQ. 3)THEN
               j3=j3+1
               IBEDGE3(1,j3)=P%I1
               IBEDGE3(2,j3)=P%I2
               IS(P%I2,P%I1)=j2
            END IF
            P=>P%NEXT
         END DO
         call cg_close_f(index_file,ier)
         
         !del the bin-tree, bedge_link and other temporal variable
         treep=>treeroot
         CALL DEL_TREE(treep)
         
         P=>HEAD
         CALL DEL_BEDGE_CHAIN(P)
         
         deallocate(bound)
         END SUBROUTINE