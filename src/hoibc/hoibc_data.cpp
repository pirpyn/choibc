#include "hoibc_data.hpp"
#include "hoibc_constants.hpp"
#include <iomanip>
#include <cmath>

using namespace hoibc;

void hoibc::check_data(const data_t& data){
  std::cout << "hoibc::check_data: i do nothing" << std::endl;
}

void hoibc::disp_data(const data_t& data, std::ostream& out) {
  using std::endl;
  out.precision(3);
  std::showpos(out);
  std::scientific(out);
  out << "# Frequency: " << data.main.frequency << " GHz" <<  endl;
  out << "# " + std::string(60,'-') << endl;
  out << "# Layer 0:" << endl;
  out << "#   TM Impedance: (" << std::real(data.material.initial_impedance[0][0]) << "," 
    << std::imag(data.material.initial_impedance[0][0]) << ")" << endl;
  out << "#   TE Impedance: (" << std::real(data.material.initial_impedance[1][1]) << "," 
    << std::imag(data.material.initial_impedance[1][1]) << ")" << endl;
  for (unsigned int i=0;i<data.material.epsr.size();i++) 
  {
    complex epsr = data.material.epsr[i];
    complex mur = data.material.mur[i];
    out << "# " + std::string(60,'-') << endl;
    out << "# Layer " << i+1 << ":" << endl;
    real ratio = std::abs(2*pi/std::sqrt(mur*epsr)/20.96/data.main.frequency)/data.material.thickness[i];
    out << "#   Thickness: " << std::noshowpos << data.material.thickness[i] << " m [ λ/" 
      << std::setprecision(1)  << std::fixed << ratio << ", " << std::setprecision(5)  << std::fixed << 1./ratio << "λ ]" << endl;
    out.precision(3);
    std::showpos(out);
    std::scientific(out);
    out << "#   Epsilon:   (" << std::real(epsr) << "," << std::imag(epsr) << ")" << endl;
    out << "#   Mu:        (" << std::real(mur) << "," << std::imag(mur) << ")" << endl;
    complex nur = std::sqrt(epsr*mur);
    out << "#   Nu:        (" << std::real(nur) << "," << std::imag(nur) 
      << ") [ index = " << std::noshowpos << std::abs(nur) << std::showpos << "]" << endl;
    complex etar = std::sqrt(mur/epsr);
    out << "#   Eta:       (" << std::real(etar) << "," << std::imag(etar) << ")" << endl;
  }
    out << "# " + std::string(60,'-') << endl;

/*  
    write(newunit,'(a,3(es10.3,","))') '# s_1: [start,end,step] = ',data%main%s1
    write(newunit,'(a,3(es10.3,","))') '# s_2: [start,end,step] = ',data%main%s2
*/
}