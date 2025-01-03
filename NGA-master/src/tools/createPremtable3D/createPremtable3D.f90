program createPremtable3D
  use prem_table3D
  use parser
  use precision
  use string
  implicit none

  integer :: n
  
  ! Filename of input file
  character(len=str_medium) :: name

  ! -------------------------------------------------------

  ! Parse the input file
  call commandline_args(name)
  call parser_init
  call parser_parsefile(name)
  
  ! Read the list of files
  call parser_getsize("List of Flamelets", nfiles)
  allocate(files(nfiles))
  call parser_read("List of Flamelets", files)
  
  ! Initialize the modules
  call prem_flamelet3D_init
  call prem_table3D_init
  
  ! -------------------------------------------------------

  print*,''
  print*,'** Files in the table **'
  do n=1,nfiles
     write(*,'(a)') trim(files(n))
     ! Read the file
     call prem_flamelet3D_readfile(files(n))
     ! Convolute with PDF
     call prem_table3D_convolute(n)
     ! Deallocate the data array
     call prem_flamelet3D_cleanup
  end do
  
  ! Change the variables names
  call prem_table3D_convert_names
  ! Compute the table
  call prem_table3D_setup
  ! Extent the bounds of the chemtable
  !call table_extent
  ! Print some statistics
  call prem_table3D_stats
  ! Write the table
  call prem_table3D_write

end program createPremtable3D



! -----------------------------------------
subroutine commandline_args(input_name)
  implicit none
  integer, external :: iargc
  character(len=*) :: input_name
  external :: getarg
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
