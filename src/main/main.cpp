#include "main.hpp"

int main() {
  hoibc::data_t data;
  
  data.main.frequency = .2;
  data.main.s1 = {0.,1.,.1};
  data.main.s2 = {0.,0.,0.};

  data.material.thickness = {0.05};
  data.material.epsr      = {hoibc::complex(1.,-1.)};
  data.material.mur       = {hoibc::complex(1.,0.)};

  data.hoibc.name         = {"ibc0","ibc0","ibc0"};
  data.hoibc.label        = {"ibc0-1","ibc0-2","ibc0-3"};
  data.hoibc.suc          = {false,false,true};
  data.hoibc.type         = {'P','P','P'};
  data.hoibc.inner_radius = {0.,0.,0.};
  data.hoibc.mode         = {1,1,2};
  data.hoibc.normalised   = {true,false,true};

  // prints the parameters to stdout
  hoibc::disp_data(data);

  // A polymorphic array of all the IBC, that will store their coefficients
  // Unknown IBC are nullptr
  std::vector<hoibc::hoibc_class*> hoibc_list;

  // Compute HOIBC coefficients
  hoibc::main(data,hoibc_list);

  std::cout << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << "HOIBC: All IBC coefficients computed with success" << std::endl;
  std::cout << "=================================================" << std::endl;
  std::cout << std::endl;

  // Write results (impedance, coeff, ...) to screen and csv files
  write_impedance_errors(data, hoibc_list);

  std::cout << "main: ended successfully " << std::endl ;
  return 0;
}

void write_impedance_errors(const hoibc::data_t& data, std::vector<hoibc::hoibc_class*>& hoibc_list){
 /*

    ! Now we will write many files and print results to the screen
    nhoibc = size(hoibc_list)

    allocate(error(5,nhoibc,2))
    
    error(:,:,:) = huge(1.0_wp)

    k0 = 2._wp*pi*data%main%frequency*1.E9_wp/speed_of_light

    ! Set the format character string to write the results in the csv file
    ! and the character format string for IBC coefficient
    call set_backend(data_extended%out%backend)

    ! To print the impédance we will need the value of the Fourier variables
    ! depending on the IBC
    */
  for ( const auto& ibc : hoibc_list ) {
    std::cout << std::string(60,'#') << std::endl;
    std::cout << std::string(60,'#') << std::endl;



      /*

      if (allocated(R_ex)) deallocate(R_ex)
      if (allocated(R_ap)) deallocate(R_ap)
      if (allocated(Z_ex)) deallocate(Z_ex)
      if (allocated(Z_ap)) deallocate(Z_ap)

      ibc = hoibc_list(i)%ibc
*/
    /*
      Set up the scan of incident angle
      For the plane, depends on (kx,ky)
      For the cylinder, depends on kz, incidence is from theta = 0, but expressed as a Fourier serie over n, truncated
      For the sphere, incidence is always from theta, phi = 0, but expressed as a Mie serie, truncated
    */
    std::vector<hoibc::real> f1, f2;
    ibc->set_fourier_variables(data,f1,f2);
    // Display the coefficient to the standard output (stdout)
    ibc->print_coeff();
    // Display the SUC to stdout
    ibc->print_suc();

    /*
      if (data_extended%out%coeff) then
        open(NEWUNIT=newunit,FILE=data_extended%out%basename//"."//ibc%label//'.coeff.txt',ACTION='write')
        ! write the coefficient to the file
        call ibc%print_coeff(UNIT=newunit)
        ! write the SUC to the file
        call ibc%disp_suc(data%optim%acc,UNIT=newunit)
        close(newunit)
      end if

          ! ########################################################################################

          ! Write info depending on the type

          ! Compute reflexion/fourier/mie coefficient
          if (allocated(R_ap)) deallocate(R_ap)
            allocate(R_ap(size(f1),size(f2),2,2))
          R_ap = ibc%get_reflexion(k0,f1,f2)

          ! Compute impedance
          if (allocated(Z_ap)) deallocate(Z_ap)
            allocate(Z_ap(size(f1),size(f2),2,2))
          Z_ap = ibc%get_impedance(k0,f1,f2)

          if (allocated(Z_ex)) deallocate(Z_ex)
            allocate(Z_ex(size(f1),size(f2),2,2))

          if (allocated(R_ex)) deallocate(R_ex)
            allocate(R_ex(size(f1),size(f2),2,2))


          if (data_extended%other%hoppe_imp) then
            ! To write the same reflexion coefficient as HOPPE & RAHMAT SAMII, 1995
          ! we must express cmplx argument in degree between 0 and 360
          ! that will be written in the csv file
          call set_arg(.True.)
          end if

          select case(ibc%type)
          case('P') ! For infinite plane, write reflection coefficients
          ! we're looking for NaN when k^2 - kx^2 = 0, ky = 0
          do i2=1,size(f2)
          do i1=1,size(f1)
          if (all(R_ap(i1,i2,:,:).ne.R_ap(i1,i2,:,:))) then
            R_ap(i1,i2,:,:) = reshape(cmplx([-1.,0.,0.,1.],[0.,0.,0.,0.],wp),[2,2])
          end if
          end do
          end do

          select case (ibc%mode)
          case (1)
          R_ex = reflexion_infinite_plane(f1,f2,&
          k0, &
          data%material%epsr, &
          data%material%mur, &
          data%material%thickness, &
          data%material%initial_impedance, &
          data%material%loss)

          Z_ex = impedance_from_reflexion_plane(f1,f2,k0,R_ex)
          case (2)
          Z_ex = impedance_infinite_plane(f1,f2,&
          k0, &
          data%material%epsr, &
          data%material%mur, &
          data%material%thickness, &
          data%material%initial_impedance, &
          data%material%loss)

          R_ex = reflexion_from_impedance_plane(f1,f2,k0,Z_ex)
          end select
          if (data_extended%out%r_ex) then
            write(str,'(I1)') ibc%mode
          filename = data_extended%out%basename//'.r_ex.MODE_'//trim(str)//"_TYPE_"//ibc%type//".csv"
          ! write(output_unit,'(a,a)') 'Writing exact reflection to ',filename
          if (data_extended%other%reflex_vs_theta) then
            call dump_to_csv(filename,&
          180._wp/pi*asin(s1),&
          180._wp/pi*asin(s2),&
          r_ex,"theta1","theta_2","r_ex",data_extended%out%skip_nan)
          else
            call dump_to_csv(filename,s1,s2,r_ex,"s1","s2","r_ex",data_extended%out%skip_nan)
          end if       
          end if

          if (data_extended%out%r_ibc) then
          filename = data_extended%out%basename//'.r_ibc.'//ibc%label//".csv"
        ! write(output_unit,'(a,a)') 'Writing IBC reflection to ',filename
        if (data_extended%other%reflex_vs_theta) then
          call dump_to_csv(filename,&
            180._wp/pi*asin(s1),&
            180._wp/pi*asin(s2),&
              r_ap,"theta1","theta_2","r_"//ibc%name,data_extended%out%skip_nan)
              else
            call dump_to_csv(filename,s1,s2,r_ap,"s1","s2","r_"//ibc%name,data_extended%out%skip_nan)
              end if       
              end if


              case ('C') ! For the cylinder, write the reflexion matrices that includes Fourier coefficient
              select case (ibc%mode)
              case (1)
              R_ex = reflexion_infinite_cylinder(f1,f2,&
                k0, &
                data%material%epsr, &
                data%material%mur, &
                data%material%thickness, &
                data%material%initial_impedance, &
                data%material%loss, &
                ibc%inner_radius)

              Z_ex = impedance_from_reflexion_cylinder(f1,f2,k0,R_ex,ibc%outer_radius)
              case (2)
              Z_ex = impedance_infinite_cylinder(f1,f2,&
                k0, &
                data%material%epsr, &
                data%material%mur, &
                data%material%thickness, &
                data%material%initial_impedance, &
                data%material%loss, &
                ibc%inner_radius)

              R_ex = reflexion_from_impedance_cylinder(f1,f2,k0,Z_ex,ibc%outer_radius)
              end select

              if (data_extended%out%r_ex) then
                write(str,'(i1)') ibc%mode
          filename = data_extended%out%basename//'.f_ex.MODE_'//trim(str)//"_TYPE_"//ibc%type
          write(str,'(SP,es10.3)') ibc%inner_radius
          filename = filename//'_'//str//"m.csv"
          ! write(output_unit,'(a,a)') 'Writing exact Fourier coefficient to ',filename
          if (data_extended%other%reflex_vs_theta) then
            call dump_to_csv(filename,&
              180._wp/pi*asin(s2),&
              f1,&
              r_ex,"theta","n","f_ex",data_extended%out%skip_nan)
          else
            call dump_to_csv(filename,s2,f1,r_ex,"s","n","f_ex",data_extended%out%skip_nan)
          end if       
          end if

          if (data_extended%out%r_ibc) then
          filename = data_extended%out%basename//'.f_ibc.'//ibc%label//".csv"
        ! write(output_unit,'(a,a)') 'Writing IBC Fourier coefficient to ',filename
        if (data_extended%other%reflex_vs_theta) then
          call dump_to_csv(filename,&
            180._wp/pi*asin(s2),&
            f1,&
              r_ap,"theta","n","f_"//ibc%name,data_extended%out%skip_nan)
              else
            call dump_to_csv(filename,s2,f1,r_ap,"s","n","f_"//ibc%name,data_extended%out%skip_nan)
              end if
              end if

              case ('S')
              select case (ibc%mode)
              case (1)
              R_ex = reflexion_infinite_sphere(f2,&
                k0, &
                data%material%epsr, &
                data%material%mur, &
                data%material%thickness, &
                data%material%initial_impedance, &
                data%material%loss, &
                ibc%inner_radius)

              Z_ex = impedance_from_reflexion_sphere(f2,k0,R_ex,ibc%outer_radius)
              case (2)
              Z_ex = impedance_infinite_sphere(f2,&
                k0, &
                data%material%epsr, &
                data%material%mur, &
                data%material%thickness, &
                data%material%initial_impedance, &
                data%material%loss, &
                ibc%inner_radius)

              R_ex = reflexion_from_impedance_sphere(f2,k0,Z_ex,ibc%outer_radius)
              end select
              ! Write Mie coefficients

              if (data_extended%out%r_ex) then
                write(str,'(I1)') ibc%mode
          filename = data_extended%out%basename//'.m_ex.MODE_'//trim(str)//"_TYPE_"//ibc%type
          write(str,'(SP,es10.3)') ibc%inner_radius
          filename = filename//'_'//str//"m.csv"
          ! write(output_unit,'(a,a)') 'Writing exact Mie coefficient to ',filename
          call dump_to_csv(filename,f2,r_ex(1,:,:,:),"n","m_ex",data_extended%out%skip_nan)
          end if

          if (data_extended%out%r_ibc) then
          filename = data_extended%out%basename//'.m_ibc.'//ibc%label//".csv"
        ! write(output_unit,'(a,a)') 'Writing IBC Mie coefficient to ',filename
            call dump_to_csv(filename,f2,r_ap(1,:,:,:),"n","m_"//ibc%name,data_extended%out%skip_nan)
              end if

              end select

              ! Relative squared error for impedance
              error(1,i,1) = norm(Z_ex(:,:,1,1) - Z_ap(:,:,1,1))**2 / norm(Z_ex(:,:,1,1))**2
              error(2,i,1) = norm(Z_ex(:,:,2,2) - Z_ap(:,:,2,2))**2 / norm(Z_ex(:,:,2,2))**2
              error(3,i,1) = norm(Z_ex(:,:,2,1) - Z_ap(:,:,2,1))**2 / norm(Z_ex(:,:,2,1))**2
              error(4,i,1) = norm(Z_ex(:,:,1,2) - Z_ap(:,:,1,2))**2 / norm(Z_ex(:,:,1,2))**2
              error(5,i,1) = norm(Z_ex - Z_ap)**2 / norm(Z_ex)**2

              ! Relative error for reflexion: same norm as in STUPFEL, IEEE Trans. Ant. v63, n4, 2015
              error(1,i,2) = norm(R_ex(:,:,1,1) - R_ap(:,:,1,1)) / norm(R_ex(:,:,1,1))
              error(2,i,2) = norm(R_ex(:,:,2,2) - R_ap(:,:,2,2)) / norm(R_ex(:,:,2,2))
              error(3,i,2) = norm(R_ex(:,:,2,1) - R_ap(:,:,2,1)) / norm(R_ex(:,:,2,1))
              error(4,i,2) = norm(R_ex(:,:,1,2) - R_ap(:,:,1,2)) / norm(R_ex(:,:,1,2))
              error(5,i,2) = sum(error(1:2,i,2))

              if (data_extended%other%hoppe_imp) then
                ! To get the impedance of HOPPE & RAHMAT SAMII, 1995
              Z_ex = vacuum_impedance*Z_ex
              Z_ex(:,:,1,1) = - Z_ex(:,:,1,1)
              Z_ap = vacuum_impedance*Z_ap
              Z_ap(:,:,1,1) = - Z_ap(:,:,1,1)
              end if

              ! print argument of impedance between -pi and pi
              ! see mod_hoibc_write_csv.myatan function
              call set_arg(.False.)

              if (data_extended%out%z_ex) then
                write(str,'(I1)') ibc%mode
        filename = data_extended%out%basename//'.z_ex.MODE_'//trim(str)//"_TYPE_"//ibc%type
        select case(ibc%type)
        case('C','S')
        write(str,'(SP,es10.3)') ibc%inner_radius
          filename = filename//'_'//str//"m"
          end select
        filename = filename//".csv"
        ! write(output_unit,'(a,a)') 'Writing exact impedance to ',filename
        call dump_to_csv(filename,s1,s2,Z_ex,"s1","s2","z_ex",data_extended%out%skip_nan)
        end if

        if (data_extended%out%z_ibc) then
        filename = data_extended%out%basename//'.z_ibc.'//ibc%label//".csv"
      ! write(output_unit,'(a,a)') 'Writing IBC impedance to ',filename
        call dump_to_csv(filename,s1,s2,Z_ap,"s1","s2","z_"//ibc%name,data_extended%out%skip_nan)
          end if

          if (data_extended%out%z_err) then
        filename = data_extended%out%basename//'.z_err.'//ibc%label//".csv"
      ! write(output_unit,'(a,a)') 'Writing difference of impedance between exact and IBC to ',filename
        call dump_to_csv(filename,s1,s2,Z_ex - Z_ap,"s1","s2","z_"//ibc%name,data_extended%out%skip_nan)
          end if

          deallocate(hoibc_list(i)%ibc)
          hoibc_list(i)%ibc = ibc
          deallocate(ibc)
          end do
*/
  }
      /*
    write(output_unit,*)
    write(output_unit,'(a)') repeat('#',line_width)
    write(output_unit,'(a)') repeat('#',line_width)
    write(output_unit,*)

    write(output_unit,'(a,a)') 'Writing CSV files to ',data_extended%out%basename

    allocate(idx(nhoibc))

    fmt_error_header = '(a40,1x,a10,1x,a4,1x,a3,1x,a4,1x,*(a10,1x))'
    fmt_error_value  = '(a40,1x,a10,1x,a4,1x,l3,1x,i4,1x,*(es10.2,1x))'

    call bubble_sort(error(5,:,1),IDX=idx)
    write(output_unit,*)
    write(output_unit,'(a)') 'For the following tables, NaN values for the antidiagonal terms (21,12) should be expected when theses matrices are diagonal.'
    write(output_unit,'(a)') 'i.e. for the plane when kx or ky = 0, for the cylinder when kz = 0, and always for the sphere.'
    write(output_unit,'(a)') 'Sorted L2 squared relative error of Z_ij and the whole matrix'
    write(output_unit,fmt_error_header) 'LABEL','NAME','TYPE','SUC','MODE','11','22','21','12','Frobenius'
    do i=1,size(idx)
      ibc = hoibc_list(idx(i))%ibc
      write(output_unit,fmt_error_value) ibc%label,ibc%name,ibc%type,ibc%suc,ibc%mode,error(:,idx(i),1)
      deallocate(ibc)
    end do

    call bubble_sort(error(5,:,2),IDX=idx)
    write(output_unit,*)
    write(output_unit,'(a)') 'Sorted L2 relative error of R_ij and Err(R_11) + Err(R_22)'
    write(output_unit,fmt_error_header) 'LABEL','NAME','TYPE','SUC','MODE','11','22','21','12','SUM(11,22)'
    do i=1,size(idx)
      ibc = hoibc_list(idx(i))%ibc
      write(output_unit,fmt_error_value) ibc%label,ibc%name,ibc%type,ibc%suc,ibc%mode,error(:,idx(i),2)
      deallocate(ibc)
    end do

    deallocate(idx)
  end subroutine*/
}
