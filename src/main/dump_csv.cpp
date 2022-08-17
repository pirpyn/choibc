#include "dump_csv.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

#ifdef _HOIBC_HAS_CPP17
#include <filesystem>
#else
#warning C++17 is needed to use std::filesystem::create_directories. We will not create output directories.
#endif

static std::string prefix;
static std::string fmt_head,line_head;
static std::string fmt_val,line_val;

void set_backend(const std::string &backend) {
  prefix = "";
  fmt_head = "%s,";
  fmt_val = "%+15.8E,";
  if (
    backend == "csv" ||
    backend == "comma" ||
    backend == "paraview" ||
    backend == "libreoffice"
  ) {
      ;
  }
  else if (
    backend == "matlab" ||
    backend == "octave" ||
    backend == "scilab"
  ) {
    prefix = hoibc::get_cmt();
    fmt_head = "%s ";
    fmt_val = "+%15.8E ";
  }
  else if ( backend == "semicolon" ) {
    prefix = hoibc::get_cmt();
    fmt_head = "%s;";
    fmt_val = "+%15.8E;";
  }
  else {
    std::cerr << "warning: set_backend: backend '" << backend << "' unknown. Standard csv used." << std::endl;
  }
}

#include <iterator>
std::string repeat(const std::string& input, size_t num) {
  std::ostringstream os;
  std::fill_n(std::ostream_iterator<std::string>(os), num, input);
  std::string repeated = os.str();
  repeated.pop_back();
  return repeated;
}

static void set_fmt(const size_t num){
  line_head = repeat(fmt_head,num);
  line_val = repeat(fmt_val,num);
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std::string& header){
  std::ofstream myfile;
  #ifdef _HOIBC_HAS_CPP17
  std::filesystem::path filepath { filename };
  std::filesystem::create_directories(filepath.parent_path());
  #endif
  myfile.open(filename,std::ios_base::out);
  assert(myfile.is_open());

  set_fmt(18);

  myfile << prefix;
  myfile << string_format(line_head,s1.c_str(),s2.c_str(),
    ("Re("+label+".11)").c_str(),("Im("+label+".11)").c_str(),("Abs("+label+".11)").c_str(),("Arg("+label+".11)").c_str(),
    ("Re("+label+".21)").c_str(),("Im("+label+".21)").c_str(),("Abs("+label+".21)").c_str(),("Arg("+label+".21)").c_str(),
    ("Re("+label+".12)").c_str(),("Im("+label+".12)").c_str(),("Abs("+label+".12)").c_str(),("Arg("+label+".12)").c_str(),
    ("Re("+label+".22)").c_str(),("Im("+label+".22)").c_str(),("Abs("+label+".22)").c_str(),("Arg("+label+".22)").c_str()
  );
  myfile << std::endl;

  for (std::size_t j = 0; j < f2.size(); j++){
    for (std::size_t i = 0; i < f1.size(); i++){
      myfile << string_format(line_val,f1[i],f2[j],
        std::real(gex[i][j][0][0]),std::imag(gex[i][j][0][0]),std::abs(gex[i][j][0][0]),std::arg(gex[i][j][0][0]),
        std::real(gex[i][j][1][0]),std::imag(gex[i][j][1][0]),std::abs(gex[i][j][1][0]),std::arg(gex[i][j][1][0]),
        std::real(gex[i][j][0][1]),std::imag(gex[i][j][0][1]),std::abs(gex[i][j][0][1]),std::arg(gex[i][j][0][1]),
        std::real(gex[i][j][1][1]),std::imag(gex[i][j][1][1]),std::abs(gex[i][j][1][1]),std::arg(gex[i][j][1][1])
      );
      myfile << std::endl;
    }
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::real>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std:: string& header){
  std::ofstream myfile;
  #ifdef _HOIBC_HAS_CPP17
  std::filesystem::path filepath { filename };
  std::filesystem::create_directories(filepath.parent_path());
  #endif
  myfile.open(filename,std::ios_base::out);
  assert(myfile.is_open());

  set_fmt(6);

  myfile << prefix;
  myfile << string_format(line_head,s1.c_str(),s2.c_str(),
    (label+".11").c_str(),
    (label+".21").c_str(),
    (label+".12").c_str(),
    (label+".22").c_str()
  );
  myfile << std::endl;

  for (std::size_t j = 0; j < f2.size(); j++){
    for (std::size_t i = 0; i < f1.size(); i++){
      myfile << string_format(line_val,f1[i],f2[j],
        gex[i][j][0][0],
        gex[i][j][1][0],
        gex[i][j][0][1],
        gex[i][j][1][1]
      );
      myfile << std::endl;
    }
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& sx, const std::string& label, const std::string& header){
  std::ofstream myfile;
  #ifdef _HOIBC_HAS_CPP17
  std::filesystem::path filepath { filename };
  std::filesystem::create_directories(filepath.parent_path());
  #endif
  myfile.open(filename,std::ios_base::out);
  assert(myfile.is_open());

  set_fmt(17);

  myfile << prefix;
  myfile << string_format(line_head,sx.c_str(),
    ("Re("+label+".11)").c_str(),("Im("+label+".11)").c_str(),("Abs("+label+".11)").c_str(),("Arg("+label+".11)").c_str(),
    ("Re("+label+".21)").c_str(),("Im("+label+".21)").c_str(),("Abs("+label+".21)").c_str(),("Arg("+label+".21)").c_str(),
    ("Re("+label+".12)").c_str(),("Im("+label+".12)").c_str(),("Abs("+label+".12)").c_str(),("Arg("+label+".12)").c_str(),
    ("Re("+label+".22)").c_str(),("Im("+label+".22)").c_str(),("Abs("+label+".22)").c_str(),("Arg("+label+".22)").c_str()
  );
  myfile << std::endl;

  for (std::size_t j = 0; j < x.size(); j++){
    myfile << string_format(line_val,x[j],
      std::real(gex[0][j][0][0]),std::imag(gex[0][j][0][0]),std::abs(gex[0][j][0][0]),std::arg(gex[0][j][0][0]),
      std::real(gex[0][j][1][0]),std::imag(gex[0][j][1][0]),std::abs(gex[0][j][1][0]),std::arg(gex[0][j][1][0]),
      std::real(gex[0][j][0][1]),std::imag(gex[0][j][0][1]),std::abs(gex[0][j][0][1]),std::arg(gex[0][j][0][1]),
      std::real(gex[0][j][1][1]),std::imag(gex[0][j][1][1]),std::abs(gex[0][j][1][1]),std::arg(gex[0][j][1][1])
    );
    myfile << std::endl;
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::real>& gex, const std::string& sx, const std::string& label, const std:: string& header){
  std::ofstream myfile;
  #ifdef _HOIBC_HAS_CPP17
  std::filesystem::path filepath { filename };
  std::filesystem::create_directories(filepath.parent_path());
  #endif
  myfile.open(filename,std::ios_base::out);
  assert(myfile.is_open());

  set_fmt(5);

  myfile << prefix;
  myfile << string_format(line_head,sx.c_str(),
    (label+".11").c_str(),
    (label+".21").c_str(),
    (label+".12").c_str(),
    (label+".22").c_str()
  );
  myfile << std::endl;

  for (std::size_t j = 0; j < x.size(); j++){
    myfile << string_format(line_val,x[j],
      gex[0][j][0][0],
      gex[0][j][1][0],
      gex[0][j][0][1],
      gex[0][j][1][1]
    );
    myfile << std::endl;
  }
  myfile.close();
}