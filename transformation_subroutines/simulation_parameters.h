//
// Created by erolsson on 25/09/2019.
//

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include <map>
#include <fstream>
#include <sstream>
#include <string>


class SimulationParameters {
public:
    explicit SimulationParameters(const std::string& parameter_file_name) {
        std::map<std::string, std::string> parameter_map;
        std::ifstream settings_file(parameter_file_name);
        std::string data_string;

        // Read all lines in the file
        auto line_count(1);
        while (getline(settings_file, data_string)) {
            auto del_position = data_string.find('=');

            // If format is not identifier=value, throw exception
            if (del_position > data_string.size()-1) {
                std::stringstream error_ss;
                error_ss << "Settings parameters on line " << line_count << "has to be on the form x=data";
                throw std::invalid_argument(error_ss.str());
            }

            // Split string in identifier, value pair
            std::string key = data_string.substr(0, del_position);
            std::string data = data_string.substr(del_position+1, data_string.size() - 1 - del_position);

            // If the same identifier occurs twice, throw exception
            if (parameter_map.count(key) != 0) {
                std::stringstream error_ss;
                error_ss << "Parameter " << key << " on line " << line_count << " is already defined";
                throw std::invalid_argument(error_ss.str());
            }

            parameter_map.insert(std::pair<std::string, std::string>(key, data));
            ++line_count;
        }
        E = std::stod(parameter_map["E"]);
        v = std::stod(parameter_map["v"]);
    }
    double E;
    double v;
};


#endif // SIMULATION_PARAMETERS_H
