//
// Created by erolsson on 06/11/2019.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>

std::vector<std::tuple<unsigned, unsigned, double> > austenite;
/*
extern "C" void sdvini_(double* statev, const double* coords, const int* nstatev, const int* ncrds, const int* noel,
        const int* npt, const int* layer, const int* kspt) {

}
*/

extern "C" void uexternaldb_(const int* lop, const int* lrestart, const double* time, const double* dtime,
        const int* kstep, const int* kinc) {
    std::cout << "Calling uexternaldb" << std::endl;
    if (*lop == 0) {
        std::cout << "start, lop == 0" << std::endl;
        // Beginning of analysis
        std::vector<std::string> line_data;
        std::string data_line;
        std::ifstream austenite_file("austenite.dat");
        while (getline(austenite_file, data_line)) {
            std::stringstream line(data_line);
            std::string val;
            while (getline(line, val, ',')) {
                line_data.push_back(val);
            }
            austenite.emplace_back(stoi(line_data[0]), stoi(line_data[0]), stod(line_data[0]));
        }
        for (auto line: austenite) {
            auto [elem, gp, val] = line;
            std::cout << elem << ", " << gp << ", " << val << std::endl;
        }
        std::abort();
    }
}