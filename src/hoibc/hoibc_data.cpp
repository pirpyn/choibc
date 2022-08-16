#include "hoibc_data.hpp"
#include "hoibc_constants.hpp"
#include <iomanip>
#include <cmath>

namespace hoibc {

  void check_data(const data_t& data){
    std::cout << "check_data: i do nothing" << std::endl;
  }

  void disp_data(const data_t& data, std::ostream& out) {
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
      out << "# s_1: [start, end, step] = " << data.main.s1[0] << ", " << data.main.s1[1] << ", " << data.main.s1[2] << std::endl;
      out << "# s_2: [start, end, step] = " << data.main.s2[0] << ", " << data.main.s2[1] << ", " << data.main.s2[2] << std::endl;
  }

  void check_and_set_material(integer layer_index, complex& epsr, complex&mur, complex& etar, complex& nur,const real& loss){
    if ( (std::imag(mur)==0.) && (std::imag(epsr)==0.) ) {
      if (loss >= 0){
        std::cerr << "Layer " << layer_index << ": adding artificial loss of " << loss << std::endl;
      }
      // The minus is very important : in the case of lossless material (ie loss == 0),
      // Then this complex ( 0., -0. ) will garantee numericaly exponential decaying waves at infinity,
      // Due to the way floating complex works.
      mur = complex(std::real(mur),-std::abs(loss));
      epsr = complex(std::real(epsr),-std::abs(loss));
    }

    etar = std::sqrt(mur/epsr);
    nur = std::sqrt(mur*epsr);

    if (std::imag(nur)>0.) {
      std::cerr << "error: impedance_infinite_plane: Im(nur) > 0" << std::endl;
      exit(1);
    }
    if (std::real(etar)<0.) {
      std::cerr << "error: impedance_infinite_plane: Re(etar) < 0" << std::endl;
      exit(1);
    }
  }

}