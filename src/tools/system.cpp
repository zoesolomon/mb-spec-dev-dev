#include "system.h"

#include <iomanip>
#include <sstream>

System::System() {
    initialized_ = false;
    monomer_json_read_ = false;
    mpi_initialized_ = false;
    // simcell_periodic_ = false;
    // tag_counter_ = 1;

    std::cerr << std::setprecision(20);

    cutoff1_ = 5.5;
    cutoff2_ = 7.5;

    mb_ = "nb";
    imcon_ = 2;
    box_ = {0, 0, 0, 0, 0, 0};

    mbx_j_ = {};
}

System::~System() {}

nlohmann::json System::GetJsonConfig() {
    return mbx_j_;
}

void System::SetUpFromJson(nlohmann::json j, std::string func) {

    try {
        imcon_ = j[func]["imcon"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["imcon"] = imcon_;

    try {
        dump_frequency_ = j[func]["dump_frequency"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["dump_frequency"] = dump_frequency_;

    try {
        mb_ = j[func]["mb"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["mb"] = mb_;

    try {
        mu_ = j[func]["mu"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["mu"] = mu_;

    try {
        int_max_ = j[func]["int_max"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["int_max"] = int_max_;

    try {
        temperature_ = j[func]["temperature"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["temperature"] = temperature_;

    try {
        tau_ = j[func]["tau"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["tau"] = tau_;

    try {
        alpha_ = j[func]["alpha"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["alpha"] = alpha_;

    try {
        cutoff1_= j[func]["cutoff_1"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["cutoff_1"] = cutoff1_;

    try {
        cutoff2_= j[func]["cutoff_2"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["cutoff_2"] = cutoff2_;

    try {
        ijk_ = j[func]["ijk"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["ijk"] = ijk_;

    try {
        zc1_ = j[func]["zc1"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["zc1"] = zc1_;

    try {
        zc2_ = j[func]["zc2"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["zc2"] = zc2_;

    try {
        box_ = j[func]["box"].get<std::vector<double>>();
    } catch (...) {
        // optional warning
    }

    // if (box_.size() == 9) {
    //     box_ABCabc_ = BoxVecToBoxABCabc(box_);
    // } 
    // else if (box_.size() == 6) {
    //     box_ABCabc_ = box_;
    //     box_ = BoxABCabcToBoxVec(box_ABCabc_);
    // }


    mbx_j_["MB-SPEC"]["box"] = box_;
}

void System::SetUpFromJson(std::string json_text, std::string func) {
    nlohmann::json j = nlohmann::json::parse(json_text);
    SetUpFromJson(j, func);
}

void System::SetUpFromJson(std::string func, char *json_file) {

    nlohmann::json j_default = {};
    nlohmann::json j;

    if (json_file != nullptr) {
        try {
            std::ifstream ifjson(json_file);
            j = nlohmann::json::parse(ifjson);
        } 
        catch (...) {
            j = j_default;
            std::cerr << "Problem loading json file: "
                      << json_file
                      << " — using defaults\n";
        }
    } 
    else {
        j = j_default;
    }

    SetUpFromJson(j, func);
}

std::string System::GetJsonText() {
    std::ostringstream ss;
    ss << mbx_j_;
    return ss.str();
}

std::vector<double> System::GetBoxInfo() {
    return box_;
    //return box_ABCabc_;
}

std::string System::GetMb() {
    return mb_;
}

size_t System::GetImcon() {
    return imcon_;
}

double System::GetDumpFrequency () {
    return dump_frequency_;
}

int System::GetIntMax() {
    return int_max_;
}

double System::GetTemperature() {
    return temperature_;
}

int System::GetTau() {
    return tau_;
}

double System::GetAlpha() {
    return alpha_;
}

bool System::GetMu() {
    return mu_;
}

double System::GetCutoff1() {
    return cutoff1_;
}

double System::GetCutoff2() {
    return cutoff2_;
}

double System::GetBoxVolume() {
    return (box_[1] - box_[0]) * (box_[3] - box_[2]) * (box_[5] - box_[4]);
}

std::string System::GetIjk() {
        return ijk_;
}

double System::GetZc1() {
    return zc1_;
}

double System::GetZc2() {
    return zc2_;
}

// std::vector<double> System::BoxVecToBoxABCabc(const std::vector<double>& box) {
//     double A, B, C, alpha, beta, gamma;

//     // // Check that first vector (3 first elements of the box) is only in the x axis
//     // if (IsZero(box[0])) {
//     //     std::string text =
//     //         "X component of first vector in box cannot be 0. Please double check your box definition.";  // +
//     //                                                                                                      // std::to_string(nmax)
//     //                                                                                                      // + " is not
//     //                                                                                                      // acceptable.
//     //                                                                                                      // Possible
//     //                                                                                                      // values are 2
//     //                                                                                                      // or 3.";
//     //     throw CUException(__func__, __FILE__, __LINE__, text);
//     // }

//     // if (!IsZero(box[1]) or !IsZero(box[2])) {
//     //     std::string text =
//     //         "Y and Z components of first vector in box must be 0. Please double check your box definition.";
//     //     throw CUException(__func__, __FILE__, __LINE__, text);
//     // }

//     // // Check that y component of second vector is not 0
//     // if (IsZero(box[4])) {
//     //     std::string text = "Y component of second vector in box cannot be 0. Please double check your box definition.";
//     //     throw CUException(__func__, __FILE__, __LINE__, text);
//     // }

//     // // Check that second vector is in XY plane
//     // if (!IsZero(box[5])) {
//     //     std::string text = "Z component of second vector in box must be 0. Please double check your box definition.";
//     //     throw CUException(__func__, __FILE__, __LINE__, text);
//     // }

//     A = box[0];
//     B = sqrt(box[3] * box[3] + box[4] * box[4]);
//     C = sqrt(box[6] * box[6] + box[7] * box[7] + box[8] * box[8]);

//     double AdotB = box[0] * box[3];
//     double AdotC = box[0] * box[6];
//     double BdotC = box[3] * box[6] + box[4] * box[7];

//     gamma = acos(AdotB / A / B);
//     beta = acos(AdotC / A / C);
//     alpha = acos(BdotC / B / C);

//     std::vector<double> box_out = {A, B, C, alpha / M_PI * 180.0, beta / M_PI * 180.0, gamma / M_PI * 180.0};
//     return box_out;
// }

// std::vector<double> System::BoxABCabcToBoxVec(const std::vector<double>& box) {
//     double A, B, C, alpha, beta, gamma;
//     A = box[0];
//     B = box[1];
//     C = box[2];
//     alpha = box[3] / 180.0 * M_PI;
//     beta = box[4] / 180.0 * M_PI;
//     gamma = box[5] / 180.0 * M_PI;

//     std::vector<double> box_out(9, 0.0);
//     box_out[0] = A;
//     box_out[3] = B * cos(gamma);
//     box_out[4] = B * sin(gamma);
//     box_out[6] = C * cos(beta);
//     double tmp = (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma);
//     box_out[7] = C * tmp;
//     box_out[8] = C * sqrt(1.0 - cos(beta) * cos(beta) - tmp * tmp);

//     return box_out;
    
    
// }