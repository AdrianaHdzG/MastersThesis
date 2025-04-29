#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <json.hpp>
using json = nlohmann::json;
extern "C" {
    double* run(int nnode, double* meshX, double* meshY, double* meshZ, double* dispV, double* dispL, double* dispT,
            double Eb, double EoverG, double EsNominal, double nis, double dfoot, double bfoot, double ni_foot, double mu_int, double qz_foot,
            const char* solver,
            const char* output
            );
    void hello_world();
    int main(int argc, char* argv[]) {

        if (argc > 1) {
            std::cout << "Arguments passed: " << std::endl;
            for (int i = 1; i < argc; ++i) { // Start from 1 to skip the program name
                std::cout << "  Argument " << i << ": " << argv[i] << std::endl;
            }
        } else {
            std::cout << "No arguments were passed." << std::endl;
        }
        // 
        std::ifstream file(argv[1]);
        json data = json::parse(file);
        for (auto it = data.begin(); it != data.end(); ++it) {
            std::cout << it.key() << std::endl;
        }
        std::vector<double> meshXArray = data["meshX"].get<std::vector<double>>();
        std::vector<double> meshYArray = data["meshY"].get<std::vector<double>>();
        std::vector<double> meshZArray = data["meshZ"].get<std::vector<double>>();
        std::vector<double> dispVArray = data["dispV"].get<std::vector<double>>();
        std::vector<double> dispLArray = data["dispL"].get<std::vector<double>>();
        std::vector<double> dispTArray = data["dispT"].get<std::vector<double>>();
        int nnode = data["nnode"];
        double Eb = data["Eb"];
        double EoverG = data["EoverG"];
        double EsNominal = data["EsNominal"];
        double nis = data["nis"];
        double dfoot = data["dfoot"];
        double bfoot = data["bfoot"];
        double ni_foot = data["ni_foot"];
        double mu_int = data["mu_int"];
        double qz_foot = data["q_foot"];
        const char* solver = data["solver"].get<std::string>().c_str();
        const char* output = data["output"].get<std::string>().c_str();

        hello_world();

        double* result = run(nnode,
            meshXArray.data(), meshYArray.data(),
            meshZArray.data(), dispVArray.data(),
            dispLArray.data(), dispTArray.data(),
            Eb, EoverG, EsNominal, nis, dfoot,
            bfoot, ni_foot, mu_int, qz_foot,
            solver, output);
        return 0;
    }
}