program createMMMChemtable
  use mmm_table
  use parser
  use precision
  use string
  implicit none

  integer :: n
  logical :: triang
  
  ! Filename of input file
  character(len=str_medium) :: name
  character(len=str_medium) :: interp

  ! -------------------------------------------------------

  ! Parse the input file
  call commandline_args(name)
  call parser_init
  call parser_parsefile(name)
  
  ! Read the list of files, get interpolation method (bilinear is default)
  call parser_getsize("List of Flamelets", nfiles)
  allocate(files(nfiles))
  call parser_read("List of Flamelets", files)
  call parser_read("Interpolation method", interp,'bilinear-Z*')
  call parser_read("Triangularize Z* table", triang,.false.)
  
  ! Initialize the modules
  call flamelet_init
  call mmm_table_init

  ! -------------------------------------------------------

  print*,''
  print*,'** Files in the table **'
  do n=1,nfiles
     write(*,'(a)') trim(files(n))
     ! Read the file
     call flamelet_readfile(files(n))
     ! Convolute with PDF
     call mmm_table_convolute(n)
     ! Deallocate the data array
     call flamelet_cleanup
  end do
  
  ! Change the variables names
  call mmm_table_convert_names
  ! Compute the table

  select case(trim(interp))
  case ('delaunay')
     print *, ''
     print *, 'Using Delaunay triangulation'
     call mmm_table_setup_delaunay
  case ('bilinear-Z*')
     print *, ''
     print *, 'Using bilinear interpolation - C and Z* directions'
     call mmm_table_setup_bilinear
  case ('bilinear-F')
     print *, ''
     print *, 'Using bilinear interpolation - C and F directions'
     call mmm_table_setup_bilinear_F
  end select

  if ((trim(combModel) .eq. 'MMFPVA-Z*') .and. triang) then
    print *, ''
    print *, 'Converting to triangular table'
    call mmm_table_triangularize
  end if
  
  ! Extent the bounds of the chemtable
  !call table_extent
  ! Print some statistics
  call mmm_table_stats
  ! Write the table
  call mmm_table_write

end program createMMMChemtable



! -----------------------------------------
subroutine commandline_args(input_name)
  implicit none
  integer, external :: iargc
  character(len=*) :: input_name
  external :: getarg
  
  call getarg(1,input_name)
  
  return
end subroutine commandline_args
