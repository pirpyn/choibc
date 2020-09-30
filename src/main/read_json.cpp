#include "read_json.hpp"
#include <fstream>
#include "json.hpp"
#include <vector>

using json = nlohmann::json;

// https://github.com/nlohmann/json/issues/2207 Addint support for complex type
namespace std {

    void to_json(json& j, const hoibc::complex& d)
    {
        j = {d.real(), d.imag()};
    }

    void from_json(const json& j, hoibc::complex& d)
    {
        d.real(j.at(0).get<hoibc::real>());
        d.imag(j.at(1).get<hoibc::real>());
    }

    void to_json(json& j, const hoibc::mode_t& d)
    {
        switch (d){
        case hoibc::mode_t::R:
            j = { 1 };
            break;
        case hoibc::mode_t::Z:
            j = { 2 };
            break;
        }
    }

    void from_json(const json& j, hoibc::mode_t& d)
    {
        switch (j.get<int>()){
        case 1:
            d = hoibc::mode_t::R;
            break;
        case 2:
            d = hoibc::mode_t::Z;
            break;
        }
    }

    void to_json(json& j, const hoibc::type_t& d)
    {
        switch (d){
        case hoibc::type_t::P:
            j = { "P" };
            break;
        case hoibc::type_t::C:
            j = { "C" };
            break;
        case hoibc::type_t::S:
            j = { "S" };
            break;
        }
    }

    void from_json(const json& j, hoibc::type_t& d)
    {
        if (j.get<std::string>() == "P"){
            d = hoibc::type_t::P;
        }
        else if (j.get<std::string>() == "C"){
            d = hoibc::type_t::C;
        }
        else if (j.get<std::string>() == "S"){
            d = hoibc::type_t::S;
        }
    }

}

data_out_t read_data_from_json(const std::string& filename){

    // read a JSON file
    std::ifstream json_file(filename);
    assert(json_file.is_open());
    json json_data;
    json_file >> json_data;

    hoibc::data_t data;
    data.main.frequency         = json_data["data"]["main"]["frequency"].get<hoibc::real>();
    data.main.s1                = json_data["data"]["main"]["s1"].get<std::array<hoibc::real,3>>();
    data.main.s2                = json_data["data"]["main"]["s2"].get<std::array<hoibc::real,3>>();

    data.material.thickness     = json_data["data"]["material"]["thickness"].get<hoibc::array<hoibc::real>>();
    data.material.epsr          = json_data["data"]["material"]["epsr"].get<hoibc::array<hoibc::complex>>();
    data.material.mur           = json_data["data"]["material"]["mur"].get<hoibc::array<hoibc::complex>>();

    data.hoibc.name             = json_data["data"]["hoibc"]["name"].get<hoibc::array<std::string>>();
    data.hoibc.suc              = json_data["data"]["hoibc"]["suc"].get<hoibc::array<bool>>();
    data.hoibc.type             = json_data["data"]["hoibc"]["type"].get<hoibc::array<hoibc::type_t>>();
    data.hoibc.inner_radius     = json_data["data"]["hoibc"]["inner_radius"].get<hoibc::array<hoibc::real>>();
    data.hoibc.mode             = json_data["data"]["hoibc"]["mode"].get<hoibc::array<hoibc::mode_t>>();
    data.hoibc.normalised       = json_data["data"]["hoibc"]["normalised"].get<hoibc::array<bool>>();

    data.optim.grad_delta       = json_data["data"]["optim"]["grad_delta"].get<hoibc::real>();
    data.optim.max_iter         = json_data["data"]["optim"]["max_iter"].get<hoibc::integer>();
    data.optim.no_constraints   = json_data["data"]["optim"]["no_constraints"].get<bool>();
    data.optim.show_iter        = json_data["data"]["optim"]["show_iter"].get<bool>();
    data.optim.tol              = json_data["data"]["optim"]["tol"].get<hoibc::real>();
    data.optim.toldx            = json_data["data"]["optim"]["toldx"].get<hoibc::real>();

    data_out_t data_out;
    data_out.data_t = data;

    data_out.basename           = json_data["data"]["out"]["basename"].get<std::string>();
    data_out.impedance_ibc      = json_data["data"]["out"]["impedance_ibc"].get<bool>();
    data_out.impedance_ex       = json_data["data"]["out"]["impedance_ex"].get<bool>();
    data_out.impedance_err      = json_data["data"]["out"]["impedance_err"].get<bool>();
    data_out.coeff              = json_data["data"]["out"]["coeff"].get<bool>();
    data_out.reflexion_ibc      = json_data["data"]["out"]["reflexion_ibc"].get<bool>();
    data_out.reflexion_ex       = json_data["data"]["out"]["reflexion_ex"].get<bool>();
    data_out.reflex_vs_theta    = json_data["data"]["out"]["reflex_vs_theta"].get<bool>();

    return data_out;
}
