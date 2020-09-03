#include "write_impedance_errors.hpp"

// https://stackoverflow.com/a/12399290
#include <limits> // std::max
#include <cmath> // std::pow
// When C++20 will be available
// #include <format> // std::format
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::for_each

using error_array = std::array<std::array<hoibc::real,2>,5>;

// Function to sort the errors and returns the sort index to get the corresponding ibc.

std::vector<std::size_t> sort_indexes(const std::vector<error_array> &v, const std::size_t& j, const std::size_t& l) {

    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),[&v,&j,&l](std::size_t i1, std::size_t i2) {return v[i1][j][l] < v[i2][j][l];});

  return idx;
}

#define SEP_WIDTH 82

#include <fstream>
#include "read_json.hpp"
#include "dump_csv.hpp"
#include <valarray> // std::asin

void write_impedance_errors(const data_out_t& data_out, std::vector<hoibc::hoibc_class*>& hoibc_list){

  // Now we will write many files and print error, ibc coeff & suc values to the screen.

  // Array to store the error of the impedance and the reflexion.
  // errors[i=ibc][j=coeff or matrix][l=impedance or reflexion]
  std::vector<error_array> errors;

  const hoibc::data_t& data = data_out.data_t;
  const hoibc::real k0 = hoibc::free_space_wavenumber(data.main.frequency);

    // ! Set the format character string to write the results in the csv file
    // ! and the character format string for IBC coefficient
    // call set_backend(data_extended%out%backend)

    // To print the impedance we will need the value of the Fourier variables
    // depending on the IBC

  for ( const auto& ibc : hoibc_list ){
    std::cout << std::endl;
    std::cout << std::string(SEP_WIDTH,'#') << std::endl;
    std::cout << std::string(SEP_WIDTH,'#') << std::endl;
    std::cout << std::endl;

    // Display the coefficient to the standard output (stdout)
    ibc->print_coeff();
    // Display the SUC to stdout
    ibc->print_suc(data.optim.tol);

    std::cout << std::endl;

    if (data_out.coeff){
      std::ofstream coeff_file;
      coeff_file.open(data_out.basename+"."+ibc->label+".coeff.txt");
      assert(coeff_file.is_open());
      // write the coefficient to the file
      ibc->print_coeff(coeff_file);
      // write the SUC to the file
      ibc->print_suc(data.optim.tol,coeff_file);
      coeff_file.close();
    }
    // ########################################################################################

    // Set up the scan of incident angle
    // For the plane, depends on (kx,ky)
    // For the cylinder, depends on kz, incidence is from theta = 0, but expressed as a Fourier serie over n, truncated
    // For the sphere, incidence is always from theta, phi = 0, but expressed as a Mie serie over n, truncated
    hoibc::array<hoibc::real> f1, f2, s1 ,s2;
    ibc->set_fourier_variables(data,f1,f2,s1,s2);
    // Compute reflexion/fourier/mie coefficient
    const hoibc::big_matrix<hoibc::complex> reflexion_ibc = ibc->get_reflexion(k0,f1,f2);
    const hoibc::big_matrix<hoibc::complex> impedance_ibc = ibc->get_impedance(k0,f1,f2);

    //   if (data_extended%other%hoppe_imp) then
    //   ! To write the same reflexion coefficient as HOPPE & RAHMAT SAMII, 1995
    //   ! we must express cmplx argument in degree between 0 and 360
    //   ! that will be written in the csv file
    //   call set_arg(.True.)
    //   end if

    hoibc::big_matrix<hoibc::complex> reflexion_ex;
    hoibc::big_matrix<hoibc::complex> impedance_ex;

    // Write info depending on the type

    switch (ibc->type){
    case hoibc::type_t::P : // For infinite plane, write reflection coefficients
    // we're looking for NaN when k^2 - kx^2 = 0, ky = 0
    // do i2=1,size(f2)
    // do i1=1,size(f1)
    // if (all(R_ap(i1,i2,:,:).ne.R_ap(i1,i2,:,:))) then
    //   R_ap(i1,i2,:,:) = reshape(cmplx([-1.,0.,0.,1.],[0.,0.,0.,0.],wp),[2,2])
    // end if
    // end do
    // end do

      switch (ibc->mode){
      case hoibc::mode_t::R:
        reflexion_ex = hoibc::plane::reflexion_infinite(f1,f2,k0,data.material);
        impedance_ex = hoibc::plane::impedance_from_reflexion(f1,f2,k0,reflexion_ex);
        break;
      case hoibc::mode_t::Z:
        impedance_ex = hoibc::plane::impedance_infinite(f1,f2,k0,data.material);
        reflexion_ex = hoibc::plane::reflexion_from_impedance(f1,f2,k0,impedance_ex);
        break;
      }

      if (data_out.reflexion_ex) {
        const std::string filename = data_out.basename+".r_ex.MODE_"+std::to_string(mode_to_int(ibc->mode))+"_TYPE_"+type_to_char(ibc->type)+".csv";
        std::cout << "Writing exact reflexion to " << filename << std::endl;
        if (data_out.reflex_vs_theta) {
          dump_to_csv(filename,180./hoibc::pi*std::asin(s1),180./hoibc::pi*std::asin(s2),reflexion_ex,"theta_1","theta_2","r_ex","");
        } else {
          dump_to_csv(filename,s1,s2,reflexion_ex,"s1","s2","r_ex","");
        }
      }

      if (data_out.reflexion_ibc) {
        const std::string filename = data_out.basename+".r_ibc."+ibc->label+".csv";
        std::cout << "Writing IBC reflexion to " << filename << std::endl;
        if (data_out.reflex_vs_theta) {
          dump_to_csv(filename,180./hoibc::pi*std::asin(s1),180./hoibc::pi*std::asin(s2),reflexion_ibc,"theta_1","theta_2","r_"+ibc->name,ibc->label);
        } else {
          dump_to_csv(filename,s1,s2,reflexion_ibc,"s1","s2","r_"+ibc->name,ibc->label);
        }
      }
      break;

    case hoibc::type_t::C : // For the cylinder, write the reflexion matrices that includes Fourier coefficient
      switch (ibc->mode){
      case hoibc::mode_t::R:
//        reflexion_ex = hoibc::cylinder::reflexion_infinite_cylinder(f1,f2,k0,data,ibc->inner_radius);
//        impedance_ex = hoibc::cylinder::impedance_from_reflexion(f1,f2,k0,reflexion_ex,ibc->outer_radius)
        break;
      case hoibc::mode_t::Z:
//        impedance_ex = hoibc::cylinder::impedance_infinite(f1,f2,k0,data,ibc->inner_radius);
//        reflexion_ex = hoibc::cylinder::reflexion_from_impedance(f1,f2,k0,Z_ex,ibc->outer_radius);
        break;
      }

      if (data_out.reflexion_ex){
        const std::string filename = data_out.basename+".f_ex.MODE_"+std::to_string(mode_to_int(ibc->mode))+"_TYPE_"+type_to_char(ibc->type)+"_"+std::to_string(ibc->inner_radius)+"m.csv";
        std::cout << "Writing exact Fourier coefficient to " << filename << std::endl;
        if (data_out.reflex_vs_theta){
          dump_to_csv(filename,180./hoibc::pi*std::asin(s2),f1,reflexion_ex,"theta","n","f_ex","");
        } else {
          dump_to_csv(filename,s2,f1,reflexion_ex,"s","n","f_ex","");
        }
      }

      if (data_out.reflexion_ibc) {
        const std::string filename = data_out.basename+".f_ibc."+ibc->label+".csv";
        std::cout << "Writing IBC Fourier coefficient to " << filename << std::endl;
        if (data_out.reflex_vs_theta) {
          dump_to_csv(filename,180./hoibc::pi*std::asin(s1),180./hoibc::pi*std::asin(s2),reflexion_ibc,"theta_1","theta_2","f_"+ibc->name,ibc->label);
        } else {
          dump_to_csv(filename,s1,s2,reflexion_ibc,"s1","s2","f_"+ibc->name,ibc->label);
        }
      }
      break;

    case hoibc::type_t::S:
      switch(ibc->mode){
      case hoibc::mode_t::R:
//        reflexion_ex = hoibc::math:sphere::reflexion_infinite_sphere(f2,k0,data,ibc->inner_radius);
//        impedance_ex = hoibc::math::sphere::impedance_from_reflexion(f2,k0,reflexion_ex,ibc->outer_radius);
        break;
      case hoibc::mode_t::Z:
//        impedance_ex = hoibc::math:sphere::impedance_infinite(f2,k0,data,ibc->inner_radius);
//        reflexion_ex = hoibc::math:sphere::reflexion_from_impedance(f2,k0,impedance_ex,ibc->outer_radius);
        break;
      }

      // Write Mie coefficients
      if (data_out.reflexion_ex){
        const std::string filename = data_out.basename+".m_ex.MODE_"+std::to_string(mode_to_int(ibc->mode))+"_TYPE_"+type_to_char(ibc->type)+"_"+std::to_string(ibc->inner_radius)+"m.csv";
        std::cout << "Writing exact Mie coefficient to " << filename << std::endl;
        dump_to_csv(filename,f2,reflexion_ex,"n","m_ex","");
      }

      if (data_out.reflexion_ibc) {
        const std::string filename = data_out.basename+".m_ibc."+ibc->label+".csv";
        std::cout << "Writing IBC Mie coefficient to " << filename << std::endl;
        dump_to_csv(filename,f2,reflexion_ibc,"n","m_"+ibc->name,ibc->label);
      }
      break;
    }

    error_array error;
    // Relative squared error for impedance
    using hoibc::operator-;
    {
      const hoibc::big_matrix<hoibc::complex> tmp = impedance_ex - impedance_ibc;
      error[0][0] = std::pow(hoibc::norm(tmp,0,0),2) / std::pow(hoibc::norm(impedance_ex,0,0),2);
      error[1][0] = std::pow(hoibc::norm(tmp,1,1),2) / std::pow(hoibc::norm(impedance_ex,1,1),2);
      error[2][0] = std::pow(hoibc::norm(tmp,1,0),2) / std::pow(hoibc::norm(impedance_ex,1,0),2);
      error[3][0] = std::pow(hoibc::norm(tmp,0,1),2) / std::pow(hoibc::norm(impedance_ex,0,1),2);
      error[4][0] = std::pow(hoibc::norm(tmp),2) / std::pow(hoibc::norm(impedance_ex),2);
    }

    // Relative error for reflexion: same norm as in STUPFEL, IEEE Trans. Ant. v63, n4, 2015
    {
      const hoibc::big_matrix<hoibc::complex> tmp = reflexion_ex - reflexion_ibc;
      error[0][1] = hoibc::norm(tmp,0,0) / hoibc::norm(reflexion_ex,0,0);
      error[1][1] = hoibc::norm(tmp,1,1) / hoibc::norm(reflexion_ex,1,1);
      error[2][1] = hoibc::norm(tmp,1,0) / hoibc::norm(reflexion_ex,1,0);
      error[3][1] = hoibc::norm(tmp,0,1) / hoibc::norm(reflexion_ex,0,1);
      error[4][1] = error[0][1] + error[1][1];
    }
    errors.push_back(error);

  //         if (data_extended%other%hoppe_imp) then
  //           ! To get the impedance of HOPPE & RAHMAT SAMII, 1995
  //         Z_ex = vacuum_impedance*Z_ex
  //         Z_ex(:,:,1,1) = - Z_ex(:,:,1,1)
  //         Z_ap = vacuum_impedance*Z_ap
  //         Z_ap(:,:,1,1) = - Z_ap(:,:,1,1)
  //         end if

  //         ! print argument of impedance between -pi and pi
  //         ! see mod_hoibc_write_csv.myatan function
  //         call set_arg(.False.)

    if (data_out.impedance_ex){
      std::string filename = data_out.basename+".z_ex.MODE_"+std::to_string(mode_to_int(ibc->mode))+"_TYPE_"+type_to_char(ibc->type);
      switch (ibc->type){
        case hoibc::type_t::C:
        case hoibc::type_t::S:
        filename += "_"+std::to_string(ibc->inner_radius)+"m";
        break;
      }
      filename += ".csv";
      std::cout << "Writing exact impedance to " << filename << std::endl;
      dump_to_csv(filename,s1,s2,impedance_ex,"s1","s2","z_ex",ibc->label);
    }

    if (data_out.impedance_ibc){
      std::string filename = data_out.basename+".z_ibc."+ibc->label+".csv";
      std::cout << "Writing IBC impedance to " << filename << std::endl;
      dump_to_csv(filename,s1,s2,impedance_ibc,"s1","s2","z_"+ibc->name,ibc->label);
    }

    if (data_out.impedance_err){
      std::string filename = data_out.basename+".z_err."+ibc->label+".csv";
      std::cout << "Writing difference of impedance between exact and IBC to " << filename << std::endl;
      dump_to_csv(filename,s1,s2,impedance_ex - impedance_ibc,"s1","s2","z_"+ibc->name,"z_ex - z_ibc "+ibc->label);
    }
  }
    
    std::cout << std::endl;
    std::cout << std::string(SEP_WIDTH,'#') << std::endl;
    std::cout << std::string(SEP_WIDTH,'#') << std::endl;
    std::cout << std::endl;

    // write(output_unit,'(a,a)') 'Writing CSV files to ',data_extended%out%basename

    // When C++20 will be available
    // const string fmt_error_header  = "{40s} {10s} {4s} {3s} {4s} {10s} {10s} {10s} {10s} {10s}\n";
    // const string fmt_error_vaue    = "{40s} {10s} {4s} {3s} {4d} {10e} {10e} {10e} {10e} {10e}\n";
    const std::string fmt_error_header  = "%40s %10s %4s %3s %4s %10s %10s %10s %10s %10s\n";
    const std::string fmt_error_value   = "%40s %10s %4c %3c %4d %10.2E %10.2E %10.2E %10.2E %10.2E\n";

    std::cout << "For the following tables, NaN values for the antidiagonal terms (21,12) should be expected when theses matrices are diagonal." << std::endl;
    std::cout << "i.e. for the plane when kx or ky = 0, for the cylinder when kz = 0, and always for the sphere." << std::endl;

    std::cout << std::endl;
    std::cout << "Sorted L2 squared relative error of Z_ij and the whole matrix" << std::endl;
    // When C++20 will be available
    // std::cout << std::format(fmt_error_header,"LABEL","NAME","TYPE","SUC","MODE","11","22","21","12","Frobenius");
    std::cout << string_format(fmt_error_header,"LABEL","NAME","TYPE","SUC","MODE","11","22","21","12","Frobenius");

    for (auto i : sort_indexes(errors,4,0)){
      const hoibc::hoibc_class* ibc = hoibc_list[i];
      // When C++20 will be available
      // std::cout << std::format(fmt_error_value,ibc->label,ibc->name,ibc->type,ibc->suc,ibc->mode,errors[i][0][0],errors[i][1][0],errors[i][2][0],errors[i][3][0],errors[i][4][0]);
      std::cout << string_format(fmt_error_value,ibc->label.c_str(),ibc->name.c_str(),ibc->type,ibc->suc?'T':'F',ibc->mode,errors[i][0][0],errors[i][1][0],errors[i][2][0],errors[i][3][0],errors[i][4][0]);
    }

    std::cout << std::endl;
    std::cout << "Sorted L2 relative error of R_ij and Err(R_11) + Err(R_22)" << std::endl;
    // When C++20 will be available
    // std::cout << std::format(fmt_error_header,"LABEL","NAME","TYPE","SUC","MODE","11","22","21","12","SUM(11,22)");
    std::cout << string_format(fmt_error_header,"LABEL","NAME","TYPE","SUC","MODE","11","22","21","12","SUM(11,22)");

    for (auto i : sort_indexes(errors,4,1)){
      const hoibc::hoibc_class* ibc = hoibc_list[i];
      // When C++20 will be available
      // std::cout << std::format(fmt_error_value,ibc->label,ibc->name,ibc->type,ibc->suc,ibc->mode,errors[i][0][0],errors[i][1][0],errors[i][2][0],errors[i][3][0],errors[i][4][0]);
      std::cout << string_format(fmt_error_value,ibc->label.c_str(),ibc->name.c_str(),ibc->type,ibc->suc?'T':'F',ibc->mode,errors[i][0][1],errors[i][1][1],errors[i][2][1],errors[i][3][1],errors[i][4][1]);
    }
}
