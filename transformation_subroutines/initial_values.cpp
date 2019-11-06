//
// Created by erolsson on 06/11/2019.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>

extern "C" void getoutdir_(char* outdir, int&, int);

struct PairHash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

std::unordered_map<std::pair<unsigned, unsigned>, double, PairHash> austenite;
/*
extern "C" void sdvini_(double* statev, const double* coords, const int* nstatev, const int* ncrds, const int* noel,
        const int* npt, const int* layer, const int* kspt) {

}
*/

std::string get_run_directory() {
    char out_dir_char[256];
    int out_dir_len = 0;
    getoutdir_(out_dir_char, out_dir_len, 256);
    std::string out_dir(out_dir_char, out_dir_char+out_dir_len);
    return out_dir;
}

extern "C" void uexternaldb_(const int* lop, const int* lrestart, const double* time, const double* dtime,
        const int* kstep, const int* kinc) {
    if (*lop == 0) {
        // Beginning of analysis
        std::string data_line;
        std::ifstream austenite_file(get_run_directory() + "/austenite.dat");

        while (getline(austenite_file, data_line)) {
            std::vector<std::string> line_data;
            std::stringstream line(data_line);
            std::string val;
            while (getline(line, val, ',')) {
                line_data.push_back(val);
            }
            austenite.emplace(std::make_pair(stoi(line_data[0]), stoi(line_data[1])), stod(line_data[2]));
        }
        for (auto line: austenite) {
            std::cout << line.first.first << ", " << line.first.second << ", " << line.second << std::endl;
        }
        std::abort();
    }
}