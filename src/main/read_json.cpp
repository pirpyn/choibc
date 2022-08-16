#include "read_json.hpp"
#include "json.hpp"
#include <fstream>
#include <vector>
using json = nlohmann::json;

// https://github.com/nlohmann/json/issues/2207 Addint support for complex type
namespace std {

    void from_json(const json& j, hoibc::complex& d)
    {
        d.real(j.at(0).get<hoibc::real>());
        d.imag(j.at(1).get<hoibc::real>());
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

    void from_json(const json& j, hoibc::start_pt_t& d)
    {
        if (j.get<std::string>() == "feasible"){
            d = hoibc::start_pt_t::feasible;
        }
        else if (j.get<std::string>() == "best"){
            d = hoibc::start_pt_t::best;
        }
    }

    #define SINGLE_ARG(...) __VA_ARGS__
    #define _assign_if_contains(j,d,field,type) \
        if (j.contains(#field)) { \
            d.field = j[#field].get<type>(); \
        }

    void from_json(const json& j, hoibc::main_t& d)
    {
        _assign_if_contains(j,d,frequency,hoibc::real);
        _assign_if_contains(j,d,s1,SINGLE_ARG(std::array<hoibc::real,3>));
        _assign_if_contains(j,d,s2,SINGLE_ARG(std::array<hoibc::real,3>));
    }

    void from_json(const json& j, hoibc::material_t& d)
    {
        _assign_if_contains(j,d,thickness,hoibc::array<hoibc::real>);
        _assign_if_contains(j,d,epsr,hoibc::array<hoibc::complex>);
        _assign_if_contains(j,d,mur,hoibc::array<hoibc::complex>);
    }

    void from_json(const json& j, hoibc::hoibc_t& d)
    {
        _assign_if_contains(j,d,name,hoibc::array<std::string>);
        _assign_if_contains(j,d,suc,hoibc::array<bool>);
        _assign_if_contains(j,d,type,hoibc::array<hoibc::type_t>);
        _assign_if_contains(j,d,inner_radius,hoibc::array<hoibc::real>);
        _assign_if_contains(j,d,mode,hoibc::array<hoibc::mode_t>);
        _assign_if_contains(j,d,normalised,hoibc::array<bool>);
    }

    void from_json(const json& j, hoibc::optim_t& d)
    {
        _assign_if_contains(j,d,grad_delta,hoibc::real);
        _assign_if_contains(j,d,max_iter,hoibc::integer);
        _assign_if_contains(j,d,no_constraints,bool);
        _assign_if_contains(j,d,show_iter,bool);
        _assign_if_contains(j,d,tol,hoibc::real);
        _assign_if_contains(j,d,toldx,hoibc::real);
        _assign_if_contains(j,d,starting_point,hoibc::start_pt_t);
    }

    void from_json(const json& j, hoibc::data_t& d)
    {
        _assign_if_contains(j,d,main,hoibc::main_t);
        _assign_if_contains(j,d,material,hoibc::material_t);
        _assign_if_contains(j,d,hoibc,hoibc::hoibc_t);
        _assign_if_contains(j,d,optim,hoibc::optim_t);
    }

    void from_json(const json& j, data_out_t& d)
    {
        d.data_t = j.get<hoibc::data_t>();
        if (j.contains("out")) {
            const json j_out = j["out"];
            _assign_if_contains(j_out,d,basename,std::string);
            _assign_if_contains(j_out,d,impedance_ibc,bool);
            _assign_if_contains(j_out,d,impedance_ex,bool);
            _assign_if_contains(j_out,d,impedance_err,bool);
            _assign_if_contains(j_out,d,coeff,bool);
            _assign_if_contains(j_out,d,reflexion_ibc,bool);
            _assign_if_contains(j_out,d,reflexion_ex,bool);
            _assign_if_contains(j_out,d,reflex_vs_theta,bool);
            _assign_if_contains(j_out,d,backend,std::string);
        }
    }
}

data_out_t read_data_from_json(const std::string& filename){

    // read a JSON file
    std::ifstream json_file(filename);
    assert(json_file.is_open());
    json json_data;
    json_file >> json_data;

    data_out_t data_out = json_data.get<data_out_t>();

    return data_out;
}
