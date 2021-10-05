/**
 * @file hamming_code.cpp
 *
 * @brief Implementation of a known hamming code - NOT compatible with Z3
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <array>
#include <chrono>
#include <limits>
#include <cctype>
#include <sstream>
#include <set>
#include <algorithm>
#include <random>
#include <fstream>

/* libraries */
#define CXXOPTS_VECTOR_DELIMITER ';' /**< hideous hack to allow backwards compatability */
#include "cxxopts/cxxopts.h"
#include "Eigen/Eigen"
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

/**
 * @brief computes n choose k using a brute-force calculation
 * 
 * @param n binomial parameter n
 * @param k binomial parameter k
 * @return uint64_t n choose k
 */
uint64_t nChoosek(uint64_t n, uint64_t k)
{
    if (k > n)
        return 0;
    if (k * 2 > n)
        k = n - k;
    if (k == 0)
        return 1;

    uint64_t result = n;
    for (uint64_t i = 2; i <= k; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

/**
 * @brief supporting C++-style printing of the GF2 type
 * 
 * @param os C++ output stream
 * @param m GF2 element to print
 * 
 * @return std::ostream& output stream with GF2 element appended
 */
std::ostream &operator<<(std::ostream &os, gf2_noz3 const &m) 
{ 
    return os << m.e;
}


/**
 * @brief Reads Hamming code parameters from a file
 * 
 * @param json_ecc_code_cfg_file JSON configuration file
 */
void hamming_code::extract_ecc_code_from_cfg_file(const std::string &json_ecc_code_cfg_file)
{
    rapidjson::Document d = read_json_cfg_file(json_ecc_code_cfg_file);
   
    K = d["k"].GetUint64();
    NmK = hamming_code::compute_hamming_code_n_parity_bits(K);
    N = NmK + K;
    uid = d["uid"].GetUint64();
    p = d["p"].GetInt64();

    const rapidjson::Value &G_json_obj = d["G"];
    // std::cout << "Size: " << G_json_obj.Size() << std::endl; 
    if(G_json_obj.Size() == N)
    {
        GT.resize(N, K);
        for(rapidjson::SizeType r = 0; r < G_json_obj.Size(); r++)
            for(rapidjson::SizeType c = 0; c < G_json_obj[r].Size(); c++)
                GT(r, c) = (bool)G_json_obj[r][c].GetUint64(); 
        G = GT.transpose();
    }
    else
    {
        G.resize(K, N);
        for(rapidjson::SizeType r = 0; r < G_json_obj.Size(); r++)
            for(rapidjson::SizeType c = 0; c < G_json_obj[r].Size(); c++)
                G(r, c) = (bool)G_json_obj[r][c].GetUint64(); 
        GT = G.transpose();
    }
    
    R.resize(K, N);
    const rapidjson::Value &R_json_obj = d["R"];
    for(rapidjson::SizeType r = 0; r < R_json_obj.Size(); r++)
        for(rapidjson::SizeType c = 0; c < R_json_obj[r].Size(); c++)
            R(r, c) = (bool)R_json_obj[r][c].GetUint64(); 
    RT = R.transpose();
    
    H.resize(NmK, N);
    const rapidjson::Value &H_json_obj = d["H"];
    for(rapidjson::SizeType r = 0; r < H_json_obj.Size(); r++)
        for(rapidjson::SizeType c = 0; c < H_json_obj[r].Size(); c++)
            H(r, c) = (bool)H_json_obj[r][c].GetUint64(); 
    HT = H.transpose();
        
    std::cout << "[INFO] Read Hamming parity matrices:" << std::endl;
    std::cout << "[INFO]     H: (" << H.rows() << ", "<< H.cols() << ")" << std::endl;
    std::cout << "[INFO]     G: (" << G.rows() << ", "<< G.cols() << ")" << std::endl;
    std::cout << "[INFO]     R: (" << R.rows() << ", "<< R.cols() << ")" << std::endl;

    if(g_verbosity > 0)
    {
        std::cout << "[DEBUG1]     G:" << std::endl << G << std::endl;
        std::cout << "[DEBUG1]     H:" << std::endl << H << std::endl;
        std::cout << "[DEBUG1]     R:" << std::endl << R << std::endl;
    }

    // sanity-check the matrices
    assert((G * HT).cwiseEqual(0).all());
}

void hamming_code::calculate_everything(
      const std::unordered_set< uint64_t > &error_positions
    , const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &corrected_codeword
    , int &corrected_col
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &syndrome
) const
{
    // calculate codeword from data
    encode(dataword, codeword);

    // inject errors
    codeword_p = codeword;
    for(uint64_t epos : error_positions)
        codeword_p(epos) = codeword_p(epos) != (bool)1;

    // calcluate error syndrome and corrected codeword
    decode(codeword_p, corrected_codeword, corrected_col, dataword_p, syndrome);
}

void hamming_code::decode(
      Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &corrected_codeword
    , int &corrected_col
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &syndrome
) const
{
    // calcluate error syndrome and corrected codeword
    syndrome = codeword_p * HT;
    corrected_codeword = codeword_p;
    corrected_col = -1;
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        if(HT.row(bidx) == syndrome)
        {
            corrected_col = (int)bidx; 
            corrected_codeword(corrected_col) = !corrected_codeword(corrected_col);
            break;
        }
    }

    // undo generator matrix
    dataword_p = corrected_codeword * RT;
}

void hamming_code::encode(
      const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword
) const
{
    codeword = dataword * G;
}

void hamming_code::calculate_dataword_p(
      const std::unordered_set< uint64_t > &error_positions
    , const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p) const
{
    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword;
    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword_p;
    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > corrected_codeword;
    int corrected_bit_pos;
    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > syndrome;
    this->calculate_everything(error_positions, dataword, codeword, codeword_p, corrected_codeword, corrected_bit_pos, dataword_p, syndrome);
}

/**
 * @brief computes the exacerbation quotient for the given ECC function
 *
 * @param ecc_code_json_cfg_fname JSON configuration file for the ECC code to simulate
 * @param n_errors the number of errors to calculate the exacerbation quotient for
 * @return map from n_errors to occurrence counts across all possible test patterns with n_errors (i.e., an unnormalized PMF)
 */
std::map< uint64_t, uint64_t > hamming_code::compute_exacerbation_factor(const int n_errors) const
{
    // calculate expected patterns
    uint64_t n_expected_patterns = nChoosek(N, n_errors); // TODO: replace with beta function calculation
    printf("[INFO] computing exacerbation factor for %d errors - testing %" PRIu64 " patterns\n", n_errors, n_expected_patterns);

    // compute all possible error masks of n_erros bits
    std::vector< bool > error_positions;
    Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > all_error_positions;
    test_pattern_generator_k_ones_in_n_bits_initialize(error_positions, N, n_errors);
    do
    {
        Eigen::Matrix< bool, 1, Eigen::Dynamic > test_pattern_mat(N);
        for(uint64_t i = 0; i < N; i++) 
            test_pattern_mat(i) = error_positions.at(i);

        all_error_positions.conservativeResize(all_error_positions.rows() + 1, N);
        all_error_positions.row(all_error_positions.rows() - 1) = test_pattern_mat;
    } while(test_pattern_generator_k_ones_in_n_bits_advance(error_positions, N, n_errors));

    // initialize mask (for the masks) with each possible variant of the pattern
    std::vector< Eigen::Matrix< bool, 1, Eigen::Dynamic > > all_error_masks;
    for(uint64_t i = 0; i < (1u << n_errors); i++)
    {
        all_error_masks.push_back(Eigen::Matrix< bool, 1, Eigen::Dynamic >::Zero(n_errors));
        for(int j = 0; j < n_errors; j++)
        {
            all_error_masks.back()(j) = (i >> j) & 1;
        
        }
    }

    // for every possible pattern of n_errors
    std::map< uint64_t, uint64_t > ef;
    for(Eigen::Index ep_idx = 0; ep_idx < all_error_positions.rows(); ep_idx++)
    {
        std::vector< int > resulting_syndromes;
        std::vector< std::unordered_set< uint64_t > > resulting_precorrection_errors;
        std::vector< std::unordered_set< uint64_t > > resulting_postcorrection_errors;

        // for every possible error pattern caused by these errors
        const Eigen::Matrix< bool, 1, Eigen::Dynamic > error_positions = all_error_positions.row(ep_idx);
        for(size_t em_idx = 0; em_idx < all_error_masks.size(); em_idx++)
        {
            int mask_bit_idx = 0;
            std::unordered_set< uint64_t > injected_errors;
            // Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > error_mask = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Zero(N);
            for(uint64_t j = 0; j < N; j++)
                if(error_positions(j) == true && all_error_masks.at(em_idx)(mask_bit_idx++) == true)
                {
                    injected_errors.insert(j);
                    // error_mask = all_error_masks.at(em_idx)(mask_bit_idx++);
                }
            assert(mask_bit_idx == n_errors && "inconsistency - code bug!");
            if(injected_errors.size() > 1) // n_errs_correctable by the ECC code... no point wasting compute on correctable errors
            {
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Zero(K); // DP DOES NOT MATTER!!!
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword_p;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > corrected_codeword;
                int corrected_bit_pos;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > syndrome;
                calculate_everything(injected_errors, dataword, codeword, codeword_p, corrected_codeword, corrected_bit_pos, dataword_p, syndrome);
                
                // check where errors exist
                std::unordered_set< uint64_t > observed_errors;
                for(uint64_t bit_idx = 0; bit_idx < K; bit_idx++)
                    if(dataword(bit_idx) != dataword_p(bit_idx))
                        observed_errors.insert(bit_idx);
                resulting_syndromes.push_back(corrected_bit_pos);
                resulting_precorrection_errors.push_back(injected_errors);
                resulting_postcorrection_errors.push_back(observed_errors);

                // DEBUG
                if(g_verbosity > 1)
                {
                    std::cout << "[DEBUG2] injecting errors at " << gf2_boolvec_to_str(error_positions) << " X " << gf2_boolvec_to_str(all_error_masks.at(em_idx)) << " {";
                    for(const uint64_t j : injected_errors)
                        std::cout << j << ", ";
                    std::cout << "}:" << std::endl;
                    std::cout << "[DEBUG2]     dataword:           " << gf2_boolvec_to_str(dataword) << std::endl;
                    std::cout << "[DEBUG2]     codeword:           " << gf2_boolvec_to_str(codeword) << std::endl;
                    std::cout << "[DEBUG2]     codeword_p:         " << gf2_boolvec_to_str(codeword_p) << std::endl;
                    std::cout << "[DEBUG2]     corrected_codeword: " << gf2_boolvec_to_str(corrected_codeword) << std::endl;
                    std::cout << "[DEBUG2]     dataword_p:         " << gf2_boolvec_to_str(dataword_p) << std::endl;
                    std::cout << "[DEBUG2]     syndrome:           " << gf2_boolvec_to_str(syndrome) << std::endl;
                    std::cout << "[DEBUG2]     bit pos:            " << corrected_bit_pos << std::endl;
                    std::cout << "[DEBUG2]     observed errors:    " << observed_errors.size() << " : {";
                    for(const uint64_t j : observed_errors)
                        std::cout << j << ", ";
                    std::cout << "}:" << std::endl;
                }
            }
            else
            {
                resulting_syndromes.push_back(-1);
                resulting_precorrection_errors.push_back({});
                resulting_postcorrection_errors.push_back({});
            }
        }
        assert(resulting_syndromes.size() == (size_t)(1 << n_errors) && "inconsistency - code bug!");
        if(((ep_idx % 1000) == 0) || (ep_idx + 1 == all_error_positions.rows()))
            update_progress((ep_idx + 1) / (float)all_error_positions.rows(), "[INFO] "
                , std::string("tested ") + std::to_string(ep_idx + 1) + " / " + std::to_string(all_error_positions.rows()) + " patterns");
    
        // add a histogram entry for the number of post-correction errors that these pre-correction errors can cause
        std::unordered_set< uint64_t > errors;
        for(size_t em_idx = 0; em_idx < all_error_masks.size(); em_idx++)
            errors.insert(resulting_postcorrection_errors.at(em_idx).cbegin(), resulting_postcorrection_errors.at(em_idx).cend());
        if(ef.count(errors.size()) == 0)
            ef[errors.size()] = 0;
        ef[errors.size()]++;
    }

    return ef;
}

// ADAPTED DIRECTLY FROM EINSIM
std::string hamming_code::to_json(void)
{
    rapidjson::Document d;
    d.SetObject();
    rapidjson::Document::AllocatorType &allocator = d.GetAllocator();

    // populate the rapidjson document with necessary fields
    d.AddMember("s", "HSC", allocator);
    d.AddMember("k", (uint64_t)K, allocator);
    d.AddMember("p", (int64_t)p, allocator);
    d.AddMember("uid", (uint64_t)get_uid(), allocator);

    // GHR matrices
    rapidjson::Value G_mat(rapidjson::kArrayType);
    rapidjson::Value H_mat(rapidjson::kArrayType);
    rapidjson::Value R_mat(rapidjson::kArrayType);
    for(int r = 0; r < G.rows(); r++)
    {
        rapidjson::Value G_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < G.cols(); c++)
            G_mat_row.PushBack< int >(G(r, c), allocator);
        G_mat.PushBack(G_mat_row, allocator);
    }
    for(int r = 0; r < H.rows(); r++)
    {
        rapidjson::Value H_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < H.cols(); c++)
            H_mat_row.PushBack< int >(H(r, c), allocator);
        H_mat.PushBack(H_mat_row, allocator);
    }
    for(int r = 0; r < R.rows(); r++)
    {
        rapidjson::Value R_mat_row(rapidjson::kArrayType);
        for(int c = 0; c < R.cols(); c++)
            R_mat_row.PushBack< int >(R(r, c), allocator);
        R_mat.PushBack(R_mat_row, allocator);
    }
    d.AddMember("G", G_mat, allocator);
    d.AddMember("H", H_mat, allocator);
    d.AddMember("R", R_mat, allocator);

    // write out the string
    rapidjson::StringBuffer strbuf;
    rapidjson::PrettyWriter< rapidjson::StringBuffer > writer(strbuf);
    writer.SetFormatOptions(rapidjson::kFormatSingleLineArray);
    d.Accept(writer);
    
    std::string json;
    json = strbuf.GetString();
    
    // this regex_replace DOES NOT WORK on OSX for some reason, so instead we manually iterate
    // json = std::regex_replace(json, std::regex("], \\["), "]\n        , ["); // format arrays nicely
    size_t index = 0;
    while(true) 
    {
        index = json.find("], [", index);
        if(index == std::string::npos) 
            break;

        json.replace(index, 4, "]\n        , [");
        index += 13;
    }
    // std::cout << "JSON: " << json << std::endl;

    return json;
}