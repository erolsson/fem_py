//
// Created by erolsson on 06/11/2019.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

extern "C" void getoutdir_(char* outdir, int&, int);
extern "C" void getpartinfoc_(char* name, int& name_len, const int& num, const int& jtype, int& user_num, int& error);

struct PairHash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

std::unordered_map<std::pair<unsigned, unsigned>, double, PairHash> austenite;
std::mutex part_info_mutex;

extern "C" void sdvini_(double* statev, const double* coords, const int& nstatev, const int& ncrds, const int& noel,
        const int& npt, const int& layer, const int& kspt) {
    int user_elem_number = 0;
    char part_name_char[80];
    int part_name_len = 0;
    int error = 0;
    const int jtype = 1;
    {
        std::lock_guard<std::mutex> lock(part_info_mutex);
        getpartinfoc_(part_name_char, part_name_len, noel, jtype, user_elem_number, error);
    }
    std::string part_name(part_name_char, part_name_char+part_name_len);
    std::pair<unsigned, unsigned> point_key(user_elem_number, npt);
    double martensite = 1 - austenite[point_key];
    for (unsigned i = 0; i != nstatev; ++i) {
        if (i == 1)  {
            // Martensite state var
            statev[i] = martensite;
        }
        else if (i == 2) {
            statev[i] = 0.48;
        }
        else {
            statev[i] = 0;
        }
    }
}


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
    }
}