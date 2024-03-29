#include "hoibc.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric> // accumulate
#include "hoibc_ibc0.hpp"
#include "hoibc_ibc3.hpp"

namespace hoibc {

  enum class hoibc_names {
    ibc0,
    ibc3,
    unknown
  };

  void main(const data_t& data, std::vector<hoibc_class*>& hoibc_list) {
    setup(data,hoibc_list);

    for ( auto&& ibc : hoibc_list ) {
      if (!ibc){ // for unknown ibc, ibc is a nullptr, so we skip
        continue;
      }
      array<real> f1, f2;
      ibc->set_fourier_variables(data,f1,f2);
      ibc->get_coeff(data,f1,f2);
    }
  }

  void setup(const data_t& data, std::vector<hoibc_class*>& hoibc_list) {

    for (std::size_t i = 0; i < data.hoibc.name.size();i++) {
      std::string name = data.hoibc.name[i];
      hoibc_class* ibc = nullptr;
      switch (resolve_names(name)) {
        case hoibc_names::ibc0:
          ibc = new hoibc_ibc0();
          break;
        case hoibc_names::ibc3:
          ibc = new hoibc_ibc3();
          break;
        default :
          std::cerr << "warning: setup: hoibc '" + name + "'' is unknown. Skipping." << std::endl;;
          continue;
          break;
      }
      hoibc_list.push_back(ibc);
      ibc->name         = name;
      ibc->suc          = (data.hoibc.suc.size()>i)?data.hoibc.suc[i]:false;
      ibc->type         = (data.hoibc.type.size()>i)?data.hoibc.type[i]:type_t::P;
      ibc->inner_radius = (data.hoibc.inner_radius.size()>i)?data.hoibc.inner_radius[i]:0.;
      ibc->outer_radius = ibc->inner_radius + std::accumulate(begin(data.material.thickness),end(data.material.thickness),0.);
      ibc->normalised   = (data.hoibc.normalised.size()>i)?data.hoibc.normalised[i]:true;
      ibc->mode         = (data.hoibc.mode.size()>i)?data.hoibc.mode[i]:mode_t::Z;
      ibc->label        = (data.hoibc.label.size()>i)?data.hoibc.label[i]:"";
      if (ibc->label == "") {
        ibc->label += "IBC_"+ibc->name+"_SUC_"+ (ibc->suc ? "T" : "F" )+ "_MODE_" + std::to_string(mode_to_int(ibc->mode)) + "_TYPE_" + type_to_char(ibc->type) ;
        if ( ibc->type != type_t::P ) {
          std::stringstream ss;
          ss << std::setprecision(3) << std::showpos << std::scientific << std::uppercase;
          ss << ibc->inner_radius;
          ibc->label += ss.str();
        }
      }
    }
  }

  hoibc_names resolve_names(const std::string name) {
    if (name == "ibc0") return hoibc_names::ibc0;
    if (name == "ibc3") return hoibc_names::ibc3;
    return hoibc_names::unknown;
  }

  void free_hoibc_list(std::vector<hoibc_class*>& hoibc_list) {
    for (auto&& ibc: hoibc_list) {
      if (ibc)
        delete ibc;
    }
  }

}