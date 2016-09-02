program pre_config   
implicit none
integer, parameter                        :: PARAMETER_LEN =  50
integer, parameter                        :: LINE_LEN      = 100
integer, parameter                        :: INTE_TYPE     =   1
integer, parameter                        :: REAL_TYPE     =   2
type :: input_para
   character ( len = PARAMETER_LEN) :: category
   character ( len = PARAMETER_LEN) :: name
   integer                          :: paratype
   integer                          :: iValue
   real(kind=8)                     :: rValue
end type
character ( len = *), parameter           :: INPUTFILE ='config.cfg'
integer, parameter                        :: NPARAMETER = 11
integer, parameter                        :: CONFIG_VERSION = 1

type(input_para)                          :: input_parameters(NPARAMETER)
integer                                   :: paratype,i
integer                                   :: fu
integer                                   :: status
character ( len = LINE_LEN)               :: line

i = 0

i = i + 1; input_parameters(i)  = input_para("control"     ,"iterations"        ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("output"      ,"solution"          ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("output"      ,"residual"          ,INTE_TYPE,0   ,1.0d0   )

i = i + 1; input_parameters(i)  = input_para("model"       ,"equation"          ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("model"       ,"turbulence"        ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("disc"        ,"space_order"       ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("disc"        ,"riemann_solver"    ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("disc"        ,"timestep_method"   ,INTE_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("disc"        ,"time_order"        ,INTE_TYPE,0   ,1.0d0   )

i = i + 1; input_parameters(i)  = input_para("disc"        ,"CFL"               ,REAL_TYPE,0   ,1.0d0   )
i = i + 1; input_parameters(i)  = input_para("disc"        ,"timestep"          ,REAL_TYPE,0   ,1.0d-8  )

if (i /= NPARAMETER) then
   write(*,*) "Fehler bei der Anzahl der Input-Parameter"
   stop 1
end if

open(newunit=fu, file=trim(INPUTFILE) )

do
   read( fu, '(A)',iostat = status) line
   if (status /= 0 ) then
      write(*,*) "Input Datei Ende erreicht"
      exit
   end if
   write(*,*) line

end do



close(fu)





!
!for section_name in parser.sections():
!    #print ('Section:', section_name)
!#     print ('  Options:', parser.options(section_name))
!    for name, value in parser.items(section_name):
!        #print ('  %s = %s' % (name, value))
!        name = name.lower()
!        value = value.lower()
!        section_name = section_name.lower()
!        for i,var in enumerate(integers):
!            if (var[0].lower() == section_name):
!                if (var[1].lower() == name):
!                    var[2] = value
!                    int_in_file[i] = 1
!
!        for i,var in enumerate(reals):
!            if (var[0].lower() == section_name):
!                if (var[1].lower() == name):
!                    var[2] = value
!                    real_in_file[i] = 1
!
!        for i,var in enumerate(strings):
!            if (var[0].lower() == section_name):
!                if (var[1].lower() == name):
!                    var[2] = value
!                    string_in_file[i] = 1
!    
!for i,var in enumerate(int_in_file):
!    if (var == 0):
!        print (integers[i][1],"not found in file")
!        exit (1)
!#
!#
!#   KOMMANDOZEILEN PARAMETER
!#
!#
!for i in sys.argv[1:]:
!    section_name = i.split("/")[0].strip().lower()
!    name = i.split("/")[1].split("=")[0].strip().lower()
!    value = i.split("/")[1].split("=")[1].strip().lower()
!    for i,var in enumerate(integers):
!        if (var[0].lower() == section_name):
!            if (var[1].lower() == name):
!                var[2] = value
!                int_in_file[i] = 1
!
!    for i,var in enumerate(reals):
!        if (var[0].lower() == section_name):
!            if (var[1].lower() == name):
!                var[2] = value
!                real_in_file[i] = 1
!
!    for i,var in enumerate(strings):
!        if (var[0].lower() == section_name):
!            if (var[1].lower() == name):
!                var[2] = value
!                string_in_file[i] = 1
!# Write binary data to a file
!with open('config.bin', 'wb') as f:
!    f.write(version.to_bytes(4, byteorder='little', signed=True))
!    for var in integers:
!        f.write(int(var[2]).to_bytes(4, byteorder='little', signed=True))
!    for var in reals:
!        f.write(struct.pack('d',float(var[2])))
end program
