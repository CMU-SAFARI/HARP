/**
 * @file beer_utils.cpp
 *
 * @brief [Adapted from the original BEER project] General utility functions
 * used by the HARP artifacts. 
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1 /**< used for the beta function calculation */
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <chrono>
#include <iomanip>
#include <algorithm>

/* libraries */
#include "Eigen/Eigen"
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/error/en.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"

/* project includes */
#include "beer_utils.h"

/** forward declaratoins for global variables */
extern int g_verbosity;

/**
 * @brief computes the binomial coefficient of (n, k)
 * 
 * @param n binomial distribution parameter n (size of total set of choices)
 * @param k binomial distribution parameter k (size of chosen set)
 * @return double n choose k
 */
double binom(int n, int k)
{
    return 1 / ((n + 1) * std::beta(n - k + 1, k + 1));
}

/**
 * @brief Routine for extracting the contents of a JSON file as a rapidjson::Document
 * 
 * @param ecc_code_json_cfg_fname filename of the JSON configuration file
 * @return rapidjson::Document rapidjson document representation
 */
rapidjson::Document read_json_cfg_file(const std::string &ecc_code_json_cfg_fname)
{
    std::ifstream ifs(ecc_code_json_cfg_fname);
    rapidjson::IStreamWrapper isw(ifs);

    rapidjson::Document d;
    rapidjson::ParseResult result = d.ParseStream< rapidjson::kParseCommentsFlag >(isw);
    if(!result) 
    {
        std::cout << "[ERROR] JSON parse error while parsing ECC configuration file " <<
            ecc_code_json_cfg_fname << ": " << rapidjson::GetParseError_En(result.Code()) 
            << " (" << result.Offset() << ")" << std::endl;
        exit(-1);
    }
    ifs.close();

    return d;
}

/**
 * @brief Initializes the test pattern with the first permutation in the sequence
 * 
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 */
void test_pattern_generator_k_ones_in_n_bits_initialize(std::vector< bool > &test_pattern, const uint64_t K, const uint64_t k_ones)
{
    test_pattern.resize(K);
    for(uint64_t i = 0; i < K; i++)
        test_pattern[i] = (i >= K - k_ones) ? 1 : 0;
    std::sort(test_pattern.begin(), test_pattern.end());
}

/**
 * @brief Advance the test pattern to the next sequential permutation (using
 * C++'s std::next_permutation) - for use with std::vector< bool >
 *
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 *
 * @return true More permutations exist
 * @return false No more permutations exist 
 */
bool test_pattern_generator_k_ones_in_n_bits_advance(std::vector< bool > &test_pattern, const uint64_t K, const uint64_t k_ones)
{
    (void)K; (void)k_ones; // reserved for future use by a diferent algorithm (e.g., not stdlib)
    return std::next_permutation(test_pattern.begin(), test_pattern.end());
}

/**
 * @brief Simple test routine to sanity-check the test pattern generator
 */
void tpg_test(void)
{
    uint64_t K = 5;

    {
        std::cout << "[INFO] testing std::vector< bool > implementation" << std::endl;
        std::vector< bool > test_pattern;
        for(uint64_t k_ones = 1; k_ones < 4; k_ones++)
        {
            std::cout << "K_ones: " << k_ones << std::endl;
            test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K, k_ones);
            do
            {
                std::cout << "    ";
                for(bool bit : test_pattern) std:: cout << bit; 
                std::cout << std::endl;
            } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K, k_ones));
        }
    }

    {
        std::cout << "[INFO] testing Eigen::Matrix< bool > implementation" << std::endl;
        Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern;
        for(uint64_t k_ones = 1; k_ones < 4; k_ones++)
        {
            std::cout << "K_ones: " << k_ones << std::endl;
            test_pattern_generator_k_ones_in_n_bits_initialize(test_pattern, K, k_ones);
            do
            {
                std::cout << "    " << test_pattern << std::endl;
            } while(test_pattern_generator_k_ones_in_n_bits_advance(test_pattern, K, k_ones));
        }
    }
}

/**
 * @brief prints out a progress bar with custom prefix/suffix text
 * 
 * @param progress progress from 0.0 to 1.0
 * @param prefix prefix text before the bar
 * @param suffix suffix text after the bar
 */
void update_progress(float progress, std::string prefix, std::string suffix)
{
	int barWidth = 40;

    std::cout << prefix << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) 
	{
        if(i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    printf("] %6.2f%% %s\r", progress * 100.0, suffix.c_str());
	if(progress >= 1.0)
		std::cout << '\n';
    std::cout.flush();
}

/**
 * @brief accumulates into a comma-separated string
 * 
 * @tparam  template specialization for accumulating a split string
 * @param s base string
 * @param elem string to append
 * @return std::string comma-separted string
 */
template<>
std::string accumulate_single_element_with_delimiter< std::string >(const std::string& s, const std::string&elem)
{
    const std::string delimiter = ", ";
    return s + (s.empty() ? std::string() : delimiter) + elem;
}

/**
 * @brief generates a unique filename suffix using the current time
 * 
 * @return std::string uniquifier
 */
std::string get_filename_uniqifier_from_time(void)
{
    time_t time = std::time(nullptr);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&time),  "%Y-%m-%d_%H:%M:%S");
    
    // append milliseconds
    auto now = std::chrono::system_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    ss << '-' << std::setfill('0') << std::setw(3) << ms.count();

    std::string s = ss.str();
    std::replace(s.begin(), s.end(), ':', '-');
    return s;
}