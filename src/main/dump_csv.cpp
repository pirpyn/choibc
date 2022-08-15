#include "dump_csv.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <filesystem>

void set_backend(const std::string &backend)
{
  hoibc::set_cmt(backend);
  hoibc::set_fmt_cmplx(backend);
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std::string& header){
  std::ofstream myfile;
  std::filesystem::path filepath { filename };
  std::filesystem::create_directories(filepath.parent_path());
  myfile.open(filename,std::ios_base::out);
  assert(myfile.is_open());

  // myfile << "# " << header << std::endl;
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);

  const std::string fmt_head =
    "# %15s;%15s;"
    "%15s;%15s;%15s;%15s;"
    "%15s;%15s;%15s;%15s;"
    "%15s;%15s;%15s;%15s;"
    "%15s;%15s;%15s;%15s";

  myfile << string_format(fmt_head,s1.c_str(),s2.c_str(),
    ("Re("+label+".11)").c_str(), ("Im("+label+".11)").c_str(), ("Abs("+label+".11)").c_str(), ("Arg("+label+".11)").c_str(),
    ("Re("+label+".21)").c_str(), ("Im("+label+".21)").c_str(), ("Abs("+label+".21)").c_str(), ("Arg("+label+".21)").c_str(),
    ("Re("+label+".12)").c_str(), ("Im("+label+".12)").c_str(), ("Abs("+label+".12)").c_str(), ("Arg("+label+".12)").c_str(),
    ("Re("+label+".22)").c_str(), ("Im("+label+".22)").c_str(), ("Abs("+label+".22)").c_str(), ("Arg("+label+".22)").c_str()) << std::endl;

  const std::string fmt_val =
    "  %+15.8E;%+15.8E;"
    "%+15.8E;%+15.8E;%+15.8E;%+15.8E;"
    "%+15.8E;%+15.8E;%+15.8E;%+15.8E;"
    "%+15.8E;%+15.8E;%+15.8E;%+15.8E;"
    "%+15.8E;%+15.8E;%+15.8E;%+15.8E";
  for (std::size_t j = 0; j < f2.size(); j++){
    for (std::size_t i = 0; i < f1.size(); i++){
      myfile << string_format(fmt_val,f1[i],f2[j],
      std::real(gex[i][j][0][0]), std::imag(gex[i][j][0][0]), std::abs(gex[i][j][0][0]), std::arg(gex[i][j][0][0]),
      std::real(gex[i][j][1][0]), std::imag(gex[i][j][1][0]), std::abs(gex[i][j][1][0]), std::arg(gex[i][j][1][0]),
      std::real(gex[i][j][0][1]), std::imag(gex[i][j][0][1]), std::abs(gex[i][j][0][1]), std::arg(gex[i][j][0][1]),
      std::real(gex[i][j][1][1]), std::imag(gex[i][j][1][1]), std::abs(gex[i][j][1][1]), std::arg(gex[i][j][1][1])) << std::endl;
    }
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& f1, const hoibc::array<hoibc::real>& f2, const hoibc::big_matrix<hoibc::real>& gex, const std::string& s1, const std::string& s2, const std::string& label, const std:: string& header){
  std::ofstream myfile;
  myfile.open(filename);
  assert(myfile.is_open());

  // myfile << "# " << header << std::endl;
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);

  const std::string fmt_head =
    "#%15s;%15s;"
    "%15s;%15s;%15s;%15s";
  myfile << string_format(fmt_head,s1.c_str(),s2.c_str(),
    (label+".11").c_str(), (label+".12").c_str(), (label+".21").c_str(), (label+".22").c_str()) << std::endl;

  const std::string fmt_val =
    " %+15.8E;%+15.8E;"
    "%+15.8E;%+15.8E;%+15.8E;%+15.8E";
  for (std::size_t j = 0; j < f2.size(); j++){
    for (std::size_t i = 0; i < f1.size(); i++){
      myfile << string_format(fmt_val,f1[i],f2[j],
      gex[i][j][0][0], gex[i][j][0][1], gex[i][j][1][0], gex[i][j][1][1]) << std::endl;
    }
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& sx, const std::string& label, const std::string& header){
  std::ofstream myfile;
  myfile.open(filename);
  assert(myfile.is_open());

  // myfile << "# " << header << std::endl;
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);

  const std::string fmt_head =
    "# %15s;"
    " %15s; %15s; %15s; %15s;"
    " %15s; %15s; %15s; %15s;"
    " %15s; %15s; %15s; %15s;"
    " %15s; %15s; %15s; %15s;";

  myfile << string_format(fmt_head,sx.c_str(),
    ("Re("+label+".11)").c_str(), ("Im("+label+".11)").c_str(), ("Abs("+label+".11)").c_str(), ("Arg("+label+".11)").c_str(),
    ("Re("+label+".21)").c_str(), ("Im("+label+".21)").c_str(), ("Abs("+label+".21)").c_str(), ("Arg("+label+".21)").c_str(),
    ("Re("+label+".12)").c_str(), ("Im("+label+".12)").c_str(), ("Abs("+label+".12)").c_str(), ("Arg("+label+".12)").c_str(),
    ("Re("+label+".22)").c_str(), ("Im("+label+".22)").c_str(), ("Abs("+label+".22)").c_str(), ("Arg("+label+".22)").c_str()) << std::endl;

    const std::string fmt_val =
    "  %15.8E;"
    " %15.8E; %15.8E; %15.8E; %15.8E;"
    " %15.8E; %15.8E; %15.8E; %15.8E;"
    " %15.8E; %15.8E; %15.8E; %15.8E;"
    " %15.8E; %15.8E; %15.8E; %15.8E;";
  for (std::size_t j = 0; j < x.size(); j++){
    myfile << string_format(fmt_val,x[j],
    std::real(gex[0][j][0][0]), std::imag(gex[0][j][0][0]), std::abs(gex[0][j][0][0]), std::arg(gex[0][j][0][0]),
    std::real(gex[0][j][1][0]), std::imag(gex[0][j][1][0]), std::abs(gex[0][j][1][0]), std::arg(gex[0][j][1][0]),
    std::real(gex[0][j][0][1]), std::imag(gex[0][j][0][1]), std::abs(gex[0][j][0][1]), std::arg(gex[0][j][0][1]),
    std::real(gex[0][j][1][1]), std::imag(gex[0][j][1][1]), std::abs(gex[0][j][1][1]), std::arg(gex[0][j][1][1])) << std::endl;
  }
  myfile.close();
}

void dump_to_csv(const std::string filename, const hoibc::array<hoibc::real>& x, const hoibc::big_matrix<hoibc::real>& gex, const std::string& sx, const std::string& label, const std:: string& header){
  std::ofstream myfile;
  myfile.open(filename);
  assert(myfile.is_open());

  // myfile << "# " << header << std::endl;
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);

  const std::string fmt_head =
    "# %15s;"
    " %15s; %15s; %15s; %15s;";
  myfile << string_format(fmt_head,sx.c_str(),
    (label+".11").c_str(), (label+".12").c_str(), (label+".21").c_str(), (label+".22").c_str()) << std::endl;

  const std::string fmt_val =
    "  %15.8E;"
    " %15.8E; %15.8E; %15.8E; %15.8E;";
  for (std::size_t j = 0; j < x.size(); j++){
    myfile << string_format(fmt_val,x[j],
    gex[0][j][0][0], gex[0][j][0][1], gex[0][j][1][0], gex[0][j][1][1]) << std::endl;
  }
  myfile.close();
}