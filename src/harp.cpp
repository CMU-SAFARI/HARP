/**
 * @file harp.cpp
 *
 * @brief [Adapted from the original BEER project] Routines for evaluating
 * different error profiling schemes
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <array>
#include <limits>
#include <cctype>
#include <sstream>
#include <set>
#include <algorithm>
#include <random>
#include <fstream>
#include <thread>

/* libraries */
#define CXXOPTS_VECTOR_DELIMITER ';' /**< hideous hack to allow backwards compatability */
#include "cxxopts/cxxopts.h"
#include "Eigen/Eigen"
#include "z3++.h"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"

/* project includes */
#include "einsim_interface.h"
#include "beer_utils.h"
#include "gf2_noz3.h"
#include "hamming_code.h"
#include "unk_hamming_code.h"
#include "evaluations.h"

int g_verbosity = 0;

/**
 * @brief returns a usable output filename (i.e, a file that does not already exist)
 * 
 * @param provided a CLI-provided filename to use - if empty, an automatically generated name will be returned
 * @param base_directory directory into which to output the procedurally generated JSON ECC code configuration files 
 * @param k data word length
 * @param rseed random seed used to generate the ECC code (or -1 if unknown)
 * @param uid UID of the ECC code
 * @param nth_run integer to append to filename (e.g., file index)
 * @return std::string usable output filename or empty string if not possible
 */
std::string gen_output_ecc_code_json_cfg_filename(const std::string &provided, const std::string &base_directory, uint64_t k, int64_t rseed, uint64_t uid, uint64_t nth_run)
{
    std::string ret;
    if(!provided.empty())
    {
        ret = provided;
        // if((ret.size() < 5) || (ret.substr(ret.size() - 5, 5) != ".json")))
        //     ret += ".json";
        if(nth_run)
        {
            std::size_t pos = ret.find_last_of('.');
            if(pos == std::string::npos)
                ret += "." + std::to_string(nth_run);
            else
                ret.insert(pos, std::string(".") + std::to_string(nth_run));
        }
    }
    else
        ret = base_directory + "/k" + std::to_string(k) + "_p" + (rseed == -1u ? "UNK" : std::to_string(rseed)) + "_uid" + std::to_string(uid) 
            + (nth_run == 0 ? "" : std::string(".") + std::to_string(nth_run)) + ".json";
    
    // check that the file doesn't already exist
    for(uint64_t i = 0; i < 300; i++)
    {
        const std::string ret_uniq = ret + (i == 0 ? "" : std::string(".") + get_filename_uniqifier_from_time());
        std::ifstream f(ret_uniq.c_str());
        if(f.good())
        {
            std::cout << "[WARN] filename " << ret_uniq << " requested but already exists ... will autogenerate uniquifier" << std::endl; 
            if(i != 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(1)); // wait for time-based filename generator if there is truly a collisoin
        }
        else
            return ret_uniq;
    }
    std::cout << "[ERROR] failed to generate unique filename" << std::endl;
    return ""; // unable to generate filename
}

/**
 * @brief entry point to function that simulates the accuracy of different
 * profiling mechanisms
 *
 * @param einsim_binary path to EINSim executable
 * @param json_output_base_directory directory into which to output the procedurally generated JSON ECC code configuration files 
 * @param num_codes_to_simulate number of randomly-generated ECC codes to simulate
 * @param num_ecc_words_to_simulate number of ECC words to simulate per ECC code
 * @param k data word length
 * @param rseed random seed used to generate the ECC code (or -1 if unknown)
 */
void run_evaluations(const std::string einsim_binary, const std::string &json_output_base_directory, const int num_codes_to_simulate, const int num_ecc_words_to_simulate, const uint64_t k, const uint64_t cseed, const uint64_t wseed)
{
    for(int n = 0; n < num_codes_to_simulate; n++)
    {
        uint64_t effective_cseed = cseed + n;
        std::cout << "[INFO] evaluating ECC code: " << n << " with effective rseed: " << effective_cseed << std::endl;

        // use EINSim to generate an ECC code and write it to a file
        std::string ecc_code_json_data = einsim_generate_ecc_code(einsim_binary, k, effective_cseed);
        if(g_verbosity > 1)
        {
            std::cout << "[DEBUG] generated ECC code:" << std::endl;
            std::cout << "    " << ecc_code_json_data << std::endl;
        }

        // create a temporary file
        rapidjson::Document d;
        d.Parse< 0 >(ecc_code_json_data.c_str());
        uint64_t uid = d["uid"].GetUint64();
        std::string output_fname = gen_output_ecc_code_json_cfg_filename("", json_output_base_directory, k, effective_cseed, uid, 0);
        std::FILE *ecc_cfg_json_file = std::fopen(output_fname.c_str(), "w");
        if(ecc_cfg_json_file == nullptr)
        {
            std::cout << "[ERROR] invalid ECC code configuration filename: " << output_fname << std::endl;
            return;
        }

        std::cout << "[INFO] writing generated ECC code to output file: " << output_fname << std::endl;
        std::fputs(ecc_code_json_data.c_str(), ecc_cfg_json_file);
        std::fclose(ecc_cfg_json_file);

        uint64_t effective_wseed = wseed + n;
        evaluate(output_fname, num_ecc_words_to_simulate, effective_wseed);
    }
}

/**
 * @brief analyzes the post-correction error probabilities for each bit given
 * uniform-random pre-correction errors
 *
 * @param einsim_binary path to EINSim executable
 * @param json_output_base_directory directory into which to output the procedurally generated JSON ECC code configuration files 
 * @param num_codes_to_simulate number of randomly-generated ECC codes to simulate
 * @param num_ecc_words_to_simulate number of ECC words to simulate per ECC code
 * @param k dataword length
 * @param rseed random seed
 */
void analyze_probabilities(const std::string einsim_binary, const std::string &json_output_base_directory, const int num_codes_to_simulate, const int num_ecc_words_to_simulate, const uint64_t k, const int64_t cseed, const int64_t wseed)
{
    std::cout << "[INFO] analyzing probabilities" << std::endl; 

    for(int n = 0; n < num_codes_to_simulate; n++)
    {
        uint64_t effective_cseed = cseed + n;
        std::cout << "[INFO] evaluating ECC code: " << n << " with effective rseed: " << effective_cseed << std::endl;

        // use EINSim to generate an ECC code and write it to a file
        // std::string ecc_code_json_data = einsim_eccgen(k, rseed);
        std::string ecc_code_json_data = einsim_generate_ecc_code(einsim_binary, k, effective_cseed);
        if(g_verbosity > 1)
        {
            std::cout << "[DEBUG] generated ECC code:" << std::endl;
            std::cout << "    " << ecc_code_json_data << std::endl;
        }

        // create a temporary file
        rapidjson::Document d;
        d.Parse< 0 >(ecc_code_json_data.c_str());
        uint64_t uid = d["uid"].GetUint64();
        std::string output_fname = gen_output_ecc_code_json_cfg_filename("", json_output_base_directory, k, effective_cseed, uid, 0);
        std::FILE *ecc_cfg_json_file = std::fopen(output_fname.c_str(), "w");
        if(ecc_cfg_json_file == nullptr)
        {
            std::cout << "[ERROR] invalid ECC code configuration filename: " << output_fname << std::endl;
            return;
        }

        std::cout << "[INFO] writing generated ECC code to output file: " << output_fname << std::endl;
        std::fputs(ecc_code_json_data.c_str(), ecc_cfg_json_file);
        std::fclose(ecc_cfg_json_file);

        uint64_t effective_wseed = wseed + n;
        analyze_ecc_code(output_fname, num_ecc_words_to_simulate, effective_wseed);
    }
}


/**
 * @brief analyzes the expected number of post-correction errors for different numbers of pre-correction errors
 *
 * @param einsim_binary path to EINSim executable
 * @param json_output_base_directory directory into which to output the procedurally generated JSON ECC code configuration files 
 * @param num_codes_to_simulate number of randomly-generated ECC codes to simulate
 * @param num_ecc_words_to_simulate number of ECC words to simulate per ECC code
 * @param k dataword length
 * @param rseed random seed
 */
void analyze_secondary_ecc(const std::string einsim_binary, const std::string &json_output_base_directory, const int num_codes_to_simulate, const int num_ecc_words_to_simulate, const uint64_t k, const int64_t cseed, const int64_t wseed)
{
    std::cout << "[INFO] analyzing secondary ECC" << std::endl; 

    std::map< uint64_t /* n_prec_errors */, std::map< uint64_t /* n_posc_errors */, uint64_t /* sample count */ > > conditional_pmf_samples;
    for(int n = 0; n < num_codes_to_simulate; n++)
    {
        uint64_t effective_cseed = cseed + n;
        std::cout << "[INFO] evaluating ECC code: " << n << " with effective rseed: " << effective_cseed << std::endl;

        // use EINSim to generate an ECC code and write it to a file
        // std::string ecc_code_json_data = einsim_eccgen(k, rseed);
        std::string ecc_code_json_data = einsim_generate_ecc_code(einsim_binary, k, effective_cseed);
        if(g_verbosity > 1)
        {
            std::cout << "[DEBUG] generated ECC code:" << std::endl;
            std::cout << "    " << ecc_code_json_data << std::endl;
        }

        // create a temporary file
        rapidjson::Document d;
        d.Parse< 0 >(ecc_code_json_data.c_str());
        uint64_t uid = d["uid"].GetUint64();
        std::string output_fname = gen_output_ecc_code_json_cfg_filename("", json_output_base_directory, k, effective_cseed, uid, 0);
        std::FILE *ecc_cfg_json_file = std::fopen(output_fname.c_str(), "w");
        if(ecc_cfg_json_file == nullptr)
        {
            std::cout << "[ERROR] invalid ECC code configuration filename: " << output_fname << std::endl;
            return;
        }

        std::cout << "[INFO] writing generated ECC code to output file: " << output_fname << std::endl;
        std::fputs(ecc_code_json_data.c_str(), ecc_cfg_json_file);
        std::fclose(ecc_cfg_json_file);

        uint64_t effective_wseed = wseed + n;
        analyze_secondary_ecc(output_fname, num_ecc_words_to_simulate, effective_wseed, conditional_pmf_samples);
    }

    // print the distribution
    std::cout << "[DATA] {";
    for(const std::pair< uint64_t, std::map< uint64_t, uint64_t > > &kv_pair_prec : conditional_pmf_samples)
    {
        const uint64_t n_prec_errors = kv_pair_prec.first;
        std::cout << n_prec_errors << ":{";
        for(const std::pair< uint64_t, uint64_t > &kv_pair_posc : kv_pair_prec.second)
        {
            const uint64_t n_posc_errors = kv_pair_posc.first;
            const uint64_t sample_count = kv_pair_posc.second;
            std::cout << n_posc_errors << ":" << sample_count << ",";
        }
        std::cout << "},";
    }
    std::cout << "}" << std::endl;
}

/**
 * @brief harp entry point that parses and handles CLI options
 * 
 * @param argc number of command line arguments
 * @param argv array of command line arguments
 * 
 * @return application return value
 */
int main(int argc, char **argv)
{
    // save the comand line as a string
    std::stringstream command_line;
    for(int i = 0; i < argc; i++)
        command_line << ' ' << argv[i];
    std::cout << "[CLI]" << command_line.str() << std::endl;

    // parse the CLI options
    cxxopts::Options option_parser("harp", "Tools for evaluating the HARP profiler");
    option_parser.add_options("Common")
        ("einsim_binary", "executable binary for EINSim", cxxopts::value< std::string >())
        ("analysis", "type of analysis to run: \"probabilities\", \"secondary\", or \"evaluations\"", cxxopts::value< std::string >())
        ("k, data_bits", "dataword length of the code to generate", cxxopts::value< uint64_t >())
        ("j, json_output_dir", "directory into wich to create the procedurally generated JSON ECC code configurtion files", cxxopts::value< std::string >()->default_value(P_tmpdir))
        ("r, codeseed", "random seed for code generation (incremented once per code)", cxxopts::value< int64_t >()->default_value("0"))
        ("s, wordseed", "random seed used for ECC word generation and evaluation (incremented once per nwords)", cxxopts::value< int64_t >()->default_value("0"))
        ("c, ncodes", "number of ECC codes to generate and simulate", cxxopts::value< int64_t >()->default_value("1"))
        ("w, nwords", "number of ECC words to generate and simulate per ECC code", cxxopts::value< int64_t >()->default_value("1"))
        ("v, verbose", "Print non-essential messages")
        ("h, help", "Show help")
        ;
    option_parser.parse_positional({"einsim_binary", "analysis"});
    option_parser.positional_help("<string : einsim_binary path> <string : analysis>");

    bool needs_help = (argc == 1);
    auto options = option_parser.parse(argc, argv);
    if(needs_help || options.count("help"))
    {
        std::cout << option_parser.help({"", "Common", "positional parameters"}) << std::endl;
        return 0;
    }

    // set g_verbosity
    g_verbosity = options.count("verbose");

    // prepare Eigen library to use multithreading
    Eigen::initParallel();

    // get the einsim binary path
    if(options.count("einsim_binary") == 0)
    {
        std::cout << "[ERROR] must provide EINSim binary path" << std::endl; 
        std::cout << option_parser.help({"", "Common", "positional parameters"}) << std::endl;
        return 0;
    }

    std::string einsim_binary = options["einsim_binary"].as< std::string >();
    std::ifstream f(einsim_binary, std::ifstream::in);
    if(f.good())
    {
        std::cout << "[INFO] einsim binary provided as " << einsim_binary << std::endl;
        f.close();
    }
    else
    {
        std::cout << "[ERROR] invalid einsim binary path: " << einsim_binary << std::endl;
        std::cout << option_parser.help({"", "Simulation"}) << std::endl;
        return -1;            
    }

    if(options.count("data_bits") == 0)
    {
        std::cout << "[ERROR] dataword length required in code generation" << std::endl;
        std::cout << option_parser.help({"", "Common"}) << std::endl;
        return -1;
    }
    uint64_t k = options["data_bits"].as< uint64_t >();

    int64_t cseed = options["codeseed"].as< int64_t >();
    int64_t wseed = options["wordseed"].as< int64_t >();
    int num_codes_to_simulate = options["ncodes"].as< int64_t >();
    int num_ecc_words_to_simulate = options["nwords"].as< int64_t >();
    std::cout << "[RSEED] " << cseed << std::endl; 
    std::cout << "[SEED] codeseed: " << cseed << " wordseed: " << wseed << std::endl; 
    
    std::string json_output_base_directory = options["json_output_dir"].as< std::string >();
    std::string analysis = options["analysis"].as< std::string >();
    if(analysis == "probabilities")
    {
        analyze_probabilities(einsim_binary, json_output_base_directory, num_codes_to_simulate, num_ecc_words_to_simulate, k, cseed, wseed);
    }
    else if(analysis == "secondary")
    {
        analyze_secondary_ecc(einsim_binary, json_output_base_directory, num_codes_to_simulate, num_ecc_words_to_simulate, k, cseed, wseed);
    }
    else if(analysis == "evaluations")
    {
        run_evaluations(einsim_binary, json_output_base_directory, num_codes_to_simulate, num_ecc_words_to_simulate, k, cseed, wseed);
    }
    else
    {
        std::cout << "[ERROR] unhandled analysis: " << analysis << std::endl;
        return -1;
    }

    return 0;
}