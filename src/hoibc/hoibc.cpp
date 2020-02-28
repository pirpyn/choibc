#include "hoibc.hpp"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include "hoibc_ibc0.hpp"

using namespace hoibc;

enum class hoibc::hoibc_names {
  ibc0,
  ibc3,
  unknown
};

void hoibc::main(const data_t& data, std::vector<hoibc_class*>& hoibc_list) { 
  setup(data,hoibc_list);

  for ( const auto& ibc : hoibc_list ) {
    if (!ibc) continue; // for unknown ibc, ibc is a nullptr, so we skip

    std::vector<real> f1, f2;
    ibc->set_fourier_variables(data,f1,f2);
    ibc->get_coeff(data,f1,f2);
  }

  std::cout << "hoibc::main: ended successfully " << std::endl;
}

void hoibc::setup(const data_t& data, std::vector<hoibc_class*>& hoibc_list) {
  hoibc_list.resize(data.hoibc.name.size());

  for (unsigned int i = 0; i < hoibc_list.size();i++) {
    std::string name = data.hoibc.name[i];
    switch (resolve_names(name)) {
      case hoibc_names::ibc0 :
        hoibc_list[i] = new hoibc_ibc0();
        break;
      default :
        std::cerr << "warning: hoibc::setup: hoibc '" + name + "'' is unknown. Skipping." << std::endl;;
        hoibc_list[i] = nullptr;
        continue;
        break;
    }
    hoibc_class* ibc  = hoibc_list[i];
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
