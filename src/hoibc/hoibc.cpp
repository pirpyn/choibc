#include "hoibc.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric> // accumulate
#include "hoibc_ibc0.hpp"

using namespace hoibc;

enum class hoibc::hoibc_names {
  ibc0,
  ibc3,
  unknown
};

void hoibc::main(const data_t& data, std::vector<hoibc_class*>& hoibc_list) { 
  setup(data,hoibc_list);

  for ( auto&& ibc : hoibc_list ) {
    if (!ibc){ // for unknown ibc, ibc is a nullptr, so we skip
      continue; 
    }
    std::vector<real> f1, f2;
    ibc->set_fourier_variables(data,f1,f2);
    ibc->get_coeff(data,f1,f2);
  }
}

void hoibc::setup(const data_t& data, std::vector<hoibc_class*>& hoibc_list) {

  for (unsigned int i = 0; i < data.hoibc.name.size();i++) {
    std::string name = data.hoibc.name[i];
    hoibc_class* ibc = nullptr;
    switch (resolve_names(name)) {
      case hoibc_names::ibc0 :
        ibc = new hoibc_ibc0();
        break;
      default :
        std::cerr << "warning: hoibc::setup: hoibc '" + name + "'' is unknown. Skipping." << std::endl;;
        continue;
        break;
    }
    hoibc_list.push_back(ibc);
    ibc->name         = name;
    ibc->label        = data.hoibc.label[i];
    if (ibc->label == "") {
      ibc->label += "IBC_"+name+"_SUC_"+ (ibc->suc ? "T" : "F" )+ "_MODE_" + std::to_string(ibc->mode) + "_TYPE_" + ibc->type ;
      if ( ibc->type != 'P' ) {
        std::stringstream ss;
        ss << std::setprecision(3) << std::showpos << std::scientific << ibc->inner_radius;
        ibc->label += ss.str();
      }
    }
    ibc->suc          = data.hoibc.suc[i];
    ibc->type         = data.hoibc.type[i];
    ibc->inner_radius = data.hoibc.inner_radius[i];
    ibc->outer_radius = ibc->inner_radius + std::accumulate(data.material.thickness.begin(),data.material.thickness.end(),0);
    ibc->normalised   = data.hoibc.normalised[i];
    ibc->mode         = data.hoibc.mode[i];
  }
}

hoibc_names hoibc::resolve_names(const std::string name) {
  if (name == "ibc0") return hoibc_names::ibc0;
  if (name == "ibc3") return hoibc_names::ibc3;
  return hoibc_names::unknown;
}

void hoibc::free_hoibc_list(std::vector<hoibc_class*>& hoibc_list) {
  for (auto&& ibc: hoibc_list) {
    if (ibc)
      delete ibc;
  }
}
