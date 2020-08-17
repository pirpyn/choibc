#include "dump_csv.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

void dump_to_csv(const std::string filename, const std::vector<hoibc::real>& f1, const std::vector<hoibc::real>& f2, const hoibc::big_matrix<hoibc::complex>& gex, const std::string& s1, const std::string& s2, const std::string& label){
  std::fstream myfile;
  myfile.open(filename);
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);
  const std::string fmt_head = 
  "# %15s; %15s;"
  " %15s; %15s; %15s; %15s"
  " %15s; %15s; %15s; %15s"
  " %15s; %15s; %15s; %15s"
  " %15s; %15s; %15s; %15s";
  myfile << string_format(fmt_head,s1,s2,
  "Re("+label+".11)", "Im("+label+".11)", "Abs("+label+".11)", "Arg("+label+".11)",
  "Re("+label+".21)", "Im("+label+".21)", "Abs("+label+".21)", "Arg("+label+".21)",
  "Re("+label+".12)", "Im("+label+".12)", "Abs("+label+".12)", "Arg("+label+".12)",
  "Re("+label+".22)", "Im("+label+".22)", "Abs("+label+".22)", "Arg("+label+".22)") << std::endl;

  const std::string fmt_val = 
  "# %15.8e; %15.8e;"
  " %15.8e %15.8e; %15.8e; %15.8e"
  " %15.8e %15.8e; %15.8e; %15.8e"
  " %15.8e %15.8e; %15.8e; %15.8e"
  " %15.8e %15.8e; %15.8e; %15.8e";
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

void dump_to_csv(const std::string filename, const std::vector<hoibc::real>& f1, const std::vector<hoibc::real>& f2, const hoibc::big_matrix<hoibc::real>& gex, const std::string& s1, const std::string& s2, const std::string& label){
  std::ofstream myfile;
  myfile.open(filename);
  myfile << std::showpos;
  myfile.flags(std::ios::scientific | std::ios::uppercase | std::ios::right);
  const std::string fmt_head = 
  "# %15s; %15s;"
  " %15s; %15s; %15s; %15s";
  myfile << string_format(fmt_head,s1,s2,
  label+".11", label+".12", label+".21", label+".22)") << std::endl;

  const std::string fmt_val = 
  "# %15.8e; %15.8e;"
  " %15.8e %15.8e; %15.8e; %15.8e";
    for (std::size_t j = 0; j < f2.size(); j++){
        for (std::size_t i = 0; i < f1.size(); i++){
        myfile << string_format(fmt_val,f1[i],f2[j],
        gex[i][j][0][0], gex[i][j][0][1], gex[i][j][1][0], gex[i][j][1][1]) << std::endl;
      }
  }
  myfile.close();
}