      SUBROUTINE POST
      USE OUTPUT
      USE CELLINFO
      use ifport
      IMPLICIT NONE
      
      include 'cgnslib_f90.h'
      integer :: index_file,index_base,index_zone,index_section,index_sol,index_field
      integer :: size_zone(3),type_cell,ier

      call cg_open_f('.\result\flowfield.cgns',CG_MODE_MODIFY,index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      index_base=1
      index_zone=1
      call cg_sol_write_f(index_file,index_base,index_zone,"FlowSolution",&
         &CellCenter,index_sol,ier)
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'NumberDensity',VAR(3,:),index_field,ier)
      if (ier .ne. CG_OK) then
         write(*,*) "cgns write error!!"
         call cg_error_exit_f
         stop
      end if
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'Density',VAR(4,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VelocityX',VAR(5,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VelocityY',VAR(6,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VelocityZ',VAR(7,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'TranslationalTemperature',VAR(8,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'RotationalTemperature',VAR(9,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VibrationalTemperature',VAR(10,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'Temperature',VAR(11,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'Mach',VAR(12,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'MoleculePerCell',VAR(13,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'MeanCollisionTime',VAR(14,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'MeanFreePath',VAR(15,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'Mcs/Mfp',VAR(16,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VelocityMagnitude',VAR(17,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'VelocityAngle',VAR(18,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'Pressure',VAR(19,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'TranslationalTemperatureX',VAR(20,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'TranslationalTemperatureY',VAR(21,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'TranslationalTemperatureZ',VAR(22,:),index_field,ier)
      
      call cg_field_write_f(index_file,index_base,index_zone,index_sol,&
         &RealDouble,'NonequilibriumPara',VAR(23,:),index_field,ier)
      
      call cg_close_f(index_file,ier)

      end subroutine


      