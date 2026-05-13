#include "system.h"

#include <iomanip>
#include <sstream>

#define EXIT_FAILURE 1

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
        nskip_ = j[func]["nskip"].get<std::vector<std::tuple<size_t, size_t>>>();
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["nskip"] = nskip_;

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
        raman_response_ = j[func]["raman_response"];
    } catch (...) {

    }

    mbx_j_["MB-SPEC"]["raman_response"] = raman_response_;

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
    //     box_maxmin_ = BoxVectorToMaxMin(box_, imcon_);
    // } 
    // else if (box_.size() == 6) {
    //     box_maxmin_ = box_;
    //     box_ = BoxMaxMinToVector(box_maxmin_, imcon_);
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

std::vector<double> System::GetBoxVector() {
    if (box_.size() == 9) {
        return box_;
    } else if (box_.size() == 6) {
        return BoxMaxMinToVector(box_, imcon_);
    } else {
        throw std::runtime_error("box must be in {xx, xy, xz, yx, yy, yz, zx, zy, zz} or {xmin, xmax, ymin, ymax, zmin, zmax} format");
    }
}

std::vector<double> System::GetBoxMaxMin() {
    if (box_.size() == 6) {
        return box_;
    } else if (box_.size() == 9) {
        return BoxVectorToMaxMin(box_, imcon_);
    } else {
        throw std::runtime_error("box must be in {xx, xy, xz, yx, yy, yz, zx, zy, zz} or {xmin, xmax, ymin, ymax, zmin, zmax} format");
    }
}

std::string System::GetMb() {
    return mb_;
}

size_t System::GetImcon() {
    return imcon_;
}

std::vector<std::tuple<size_t, size_t>> System::GetNskip() {
    return nskip_;
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
// TODO: catch all so 0/1 t/f both work
bool System::GetMu() {
    return mu_;
}

double System::GetCutoff1() {
    return cutoff1_;
}

double System::GetCutoff2() {
    return cutoff2_;
}


// TODO: change box so xx xy xz yx yy yz zx zy zz then can calc for monoclinic
double System::GetBoxVolume() {
    if (box_.size() == 9) {
        if (imcon_ < 3) {
            return box_[0] * box_[4] * box_[8];
        }
        else {
            double cross_x = box_[4]*box_[8] - box_[5]*box_[7];
            double cross_y = -box_[1]*box_[8] + box_[7]*box_[2];
            double cross_z = box_[1]*box_[5] - box_[4]*box_[2];

            return box_[0]*cross_x + box_[3]*cross_y + box_[6]*cross_z;
        }
    }
    if (box_.size() == 6) {
        if (imcon_ < 3) {
            box_vector_ = BoxMaxMinToVector(box_, imcon_);
            double vmm = (box_[1] - box_[0]) * (box_[3] - box_[2]) * (box_[5] - box_[4]);
            double vv = box_vector_[0] * box_vector_[4] * box_vector_[8];
            if (vmm == vv) {
                return vmm;
            } else {
                throw std::runtime_error("Cannot calculate box volume using Min/Max coordinates and imcon = " + imcon_);
            }
        } else {
            throw std::runtime_error("Cannot calculate box volume using Min/Max coordinates and imcon = " + imcon_);
        }
    }
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

// true if depol
// TODO: .tolower() catch all 'depol' 'dep' 'd' t/f 'iso' etc
bool System::GetRamanResponse() {
    if (raman_response_.compare("depolarized") == 0) {
        return true;
    }
    return false;
}

std::vector<double> System::BoxVectorToMaxMin(const std::vector<double>& box, size_t imcon) {
    if (imcon < 3) {
        std::vector<double> box_out = {0, box[0], 0, box[4], 0, box[8]};
        return box_out;
    } else {
        throw std::runtime_error("Cannot calculate box using vector coordinates and imcon = " + imcon);
    }
}

std::vector<double> System::BoxMaxMinToVector(const std::vector<double>& box, size_t imcon) {
    if (imcon < 3) {
        std::vector<double> box_out = {box[1] - box[0], 0, 0, 0, box[3] - box[2], 0, 0, 0, box[5] - box[4]};
        return box_out;
    } else {
        throw std::runtime_error("Cannot calculate box using Min/Max coordinates and imcon = " + imcon);
    }
    
}