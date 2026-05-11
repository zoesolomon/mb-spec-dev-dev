#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <mpi.h>
#include <cmath>
#include "json/json.h"

class System {

public:
    System();
    ~System();

    nlohmann::json GetJsonConfig();

    void SetUpFromJson(std::string func, char *json_file = nullptr);
    void SetUpFromJson(std::string json_text, std::string func);
    void SetUpFromJson(nlohmann::json j, std::string func);

    std::string GetJsonText();

    std::vector<double> GetBoxVector();
    std::vector<double> GetBoxMaxMin();
    std::string GetMb();
    size_t GetImcon();
    double GetDumpFrequency();
    int GetIntMax();
    double GetTemperature();
    int GetTau();
    double GetAlpha();
    bool GetMu();
    double GetCutoff1();
    double GetCutoff2();
    double GetBoxVolume();
    std::string GetIjk();
    double GetZc1();
    double GetZc2();
    bool GetRamanResponse();

    
private:

    std::string mb_;
    size_t imcon_;
    double dump_frequency_;
    int int_max_;
    double temperature_;
    int tau_; // double??
    double alpha_;
    bool mu_;

    double cutoff1_;
    double cutoff2_;
    std::string ijk_;
    double zc1_;
    double zc2_;  

    bool initialized_;
    bool monomer_json_read_;
    bool use_pbc_;
    bool mpi_initialized_;
    std::string raman_response_;

    int mpi_rank_;
    MPI_Comm world_;

    std::vector<double> box_;
    std::vector<double> box_vector_;
    std::vector<double> box_maxmin_;

    nlohmann::json mbx_j_;

    // // ?
    // bool simcell_periodic_;
    // int tag_counter_;

    // ?
    std::vector<double> BoxVectorToMaxMin(const std::vector<double>& box, size_t imcon);
    std::vector<double> BoxMaxMinToVector(const std::vector<double>& box, size_t imcon);
};

#endif
