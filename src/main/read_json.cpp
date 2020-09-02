#include "read_json.hpp"
#include <fstream>
#include "json.hpp"
#include <vector>

using json = nlohmann::json;

// https://github.com/nlohmann/json/issues/2207 Addint support for complex type
namespace std {
    template<typename R, typename std::enable_if<std::is_arithmetic<R>::value>::type* = nullptr>
    void to_json(nlohmann::json& j, const complex<R>& d)
    {
        j = {d.real(), d.imag()};
    }

    template<typename R, typename std::enable_if<std::is_arithmetic<R>::value>::type* = nullptr>
    void from_json(const nlohmann::json& j, complex<R>& d)
    {
        d.real(j.at(0).get<R>());
        d.imag(j.at(1).get<R>());
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
    
    data.material.thickness     = json_data["data"]["material"]["thickness"].get<std::vector<hoibc::real>>();
    data.material.epsr          = json_data["data"]["material"]["epsr"].get<std::vector<hoibc::complex>>();
    data.material.mur           = json_data["data"]["material"]["mur"].get<std::vector<hoibc::complex>>();

    data.hoibc.name             = json_data["data"]["hoibc"]["name"].get<std::vector<std::string>>();
    data.hoibc.suc              = json_data["data"]["hoibc"]["suc"].get<std::vector<bool>>();
    data.hoibc.type             = json_data["data"]["hoibc"]["type"].get<std::vector<short>>();
    data.hoibc.inner_radius     = json_data["data"]["hoibc"]["inner_radius"].get<std::vector<hoibc::real>>();
    data.hoibc.mode             = json_data["data"]["hoibc"]["mode"].get<std::vector<char>>();
    data.hoibc.normalised       = json_data["data"]["hoibc"]["normalised"].get<std::vector<bool>>();

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

    return data_out;
}
