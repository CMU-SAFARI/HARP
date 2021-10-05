/**
 * @file beer_utils.h
 *
 * @brief [Adapted from the original BEER project] General utility functions
 * used by the HARP artifacts.
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef BEER_UTILS_H
#define BEER_UTILS_H

#include <vector>
#include <unordered_set>
#include <set>
#include <cmath>
#include <string>
#include <numeric>
#include <random>

#include "Eigen/Eigen"
#include "rapidjson/document.h"
// #include "z3++.h"
#include "crc64/crc64.h"

#include "gf2.h"
#include "gf2_noz3.h"

/**
 * @brief apply an AND-reduction on a vector of types
 * 
 * @tparam T type of each element
 * @param vec vector of elements to AND-reduce
 * 
 * @return T AND-reduction of vec
 */
template< typename T > 
T z3_and_reduce(const std::vector< T > &vec)
{
	assert(vec.size() > 0);

	T ret = vec.at(0);
	for(std::size_t i = 1; i < vec.size(); i++)
		ret = ret && vec.at(i);
	return ret;
}


/**
 * @brief apply an OR-reduction on a vector of types
 * 
 * @tparam T type of each element
 * @param vec vector of elements to OR-reduce
 * 
 * @return T OR-reduction of vec
 */
template< typename T >
T z3_or_reduce(const std::vector< T > &vec)
{
	assert(vec.size() > 0);

	T ret = vec.at(0);
	for(std::size_t i = 1; i < vec.size(); i++)
		ret = ret || vec.at(i);
	return ret;
}

/**
 * @brief hash fucntion used to compute the UID of an ECC code (copied from EINSim)
 * 
 * @tparam T matrix type
 * @param mat_list list of matrices to compute the UID of
 * @return uint64_t computed UID
 */
template< typename T >
uint64_t hash_matrix(std::initializer_list< T > mat_list) 
{
    std::stringstream ss;
    for(const auto &mat : mat_list)
        for(int r = 0; r < mat.rows(); r++)
            for(int c = 0; c < mat.cols(); c++)
                ss << mat(r, c);
    return crc64(-1, ss.str().c_str(), ss.str().length());
};

/**
 * @brief Compute the UID for an ECC code in the same manner that EINSim would
 * 
 * @param G ECC generator matrix
 * @param H ECC parity-check matrix
 * @param R ECC degenerator matrix
 * @return uint64_t computed UID
 */
template< typename T >
uint64_t compute_uid(
      const Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > &G
    , const Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > &H
    , const Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > &R)
{
    Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > GT = G.transpose();
    return hash_matrix({GT, H, R});
}

/**
 * @brief computes the number of 1s in an integer-like type
 * 
 * @tparam T integer-like type
 * @param n value to compute hamming weight over
 * @return T hamming weight
 */
template< typename T >
T hamming_weight(T n)
{
    T count = 0;
    while(n) 
    {
        n &= (n - 1);
        count++;
    }
    return count;
}

/**
 * @brief computes the hamming distance (in GF(2)) between two matrices
 * 
 * @tparam T Eigen matrix rows
 * @tparam R Eigen matrix cols
 * @tparam C Eigen matrix type
 * @param a first matrix
 * @param b second matrix
 * @return uint64_t hamming distance
 */
template< typename T, int R, int C >
uint64_t hamming_distance(const Eigen::Matrix< T, R, C > &a, const Eigen::Matrix< T, R, C > &b)
{
    uint64_t hd = 0;
    for(Eigen::Index r = 0; r < a.rows(); r++)
        for(Eigen::Index c = 0; c < a.cols(); c++)
            hd += a(r, c) != b(r, c);
    return hd;
}

/**
 * @brief pretty-formats an Eigen::Matrix row-vector of boolean values to a string
 * 
 * @tparam T row-vector length
 * @param v vector
 * @return std::string pretty-formatted string
 */
template< typename T >
std::string gf2_boolvec_to_str(const Eigen::Matrix< T, 1, Eigen::Dynamic > &v)
{
    std::string ret;
    ret += "[";
    for(Eigen::Index c = 0; c < v.cols(); c++)
        ret += " " + std::to_string((int)v(c));
    ret += " ]";
    return ret;
}

/**
 * @brief appends elem to string s using a pretty delimiter
 * 
 * @tparam T type of the object to accumulate
 * @param s string to append to
 * @param elem object to append on
 * @return std::string accumulated string
 */
template< typename T >
std::string accumulate_single_element_with_delimiter(const std::string& s, const T&elem)
{
    const std::string delimiter = ", ";
    return s + (s.empty() ? std::string() : delimiter) + std::to_string(elem);
}

template<>
std::string accumulate_single_element_with_delimiter< std::string >(const std::string& s, const std::string&elem);

template< typename T >
std::string stl_unary_type_to_string(const T &v)
{
    return std::accumulate(v.begin(), v.end(), std::string(), accumulate_single_element_with_delimiter< typename T::value_type >);
}

template< typename T >
std::string vec_to_string(const std::vector< T > &v)
{
    return stl_unary_type_to_string(v);
}

template< typename T >
std::string set_to_string(const std::set< T > &v)
{
    return stl_unary_type_to_string(v);
}

template< typename T >
std::string set_to_string(const std::unordered_set< T > &v)
{
    return stl_unary_type_to_string(v);
}

template< typename K, typename V >
std::string map_to_string(const std::map< K, V > &m)
{

    const std::string delimiter = ", ";
    return std::accumulate(m.begin(), m.end(), std::string(),
        [delimiter](const std::string& s, const std::pair< uint64_t, uint64_t >& p) 
        {
            return s + (s.empty() ? std::string() : delimiter) + std::to_string(p.first) + " : " + std::to_string(p.second);
        }
    );
}


/**
 * @brief Initializes the test pattern with the first permutation in the sequence
 * 
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 */
template< typename T >
void test_pattern_generator_k_ones_in_n_bits_initialize(Eigen::Matrix< T, 1, Eigen::Dynamic > &test_pattern, const uint64_t K, const uint64_t k_ones)
{
    test_pattern = Eigen::Matrix< T, 1, Eigen::Dynamic >::Zero(K);
    test_pattern.segment(K - k_ones, k_ones) = Eigen::Matrix< T, 1, Eigen::Dynamic >::Ones(k_ones);
    std::sort(test_pattern.data(), test_pattern.data() + test_pattern.size());
}

/**
 * @brief Advance the test pattern to the next sequential permutation (using
 * C++'s std::next_permutation) - for use with Eigen::Matrix< bool >
 *
 * @param test_pattern Test patern to initialize
 * @param K Desired word length
 * @param k_ones Number of bits that should be set
 *
 * @return true More permutations exist
 * @return false No more permutations exist 
 */
template< typename T >
bool test_pattern_generator_k_ones_in_n_bits_advance(Eigen::Matrix< T, 1, Eigen::Dynamic > &test_pattern, const uint64_t K, const uint64_t k_ones)
{
    (void)K; (void)k_ones; // reserved for future use by a diferent algorithm (e.g., not stdlib)
    return std::next_permutation(test_pattern.data(), test_pattern.data()  + test_pattern.size());
}

/**
 * @brief Return a random word
 *
 * @param word length-bit word to set to random values
 * @param dist C++ random distribution to generate vaues from
 * @param rng C++ random number generator
 * @param length desired word length
 * @param n_samples_remaining number of samples yet to be generated
 *
 * @return true More samples remaining
 * @return false No more samples remaining 
 */
template< typename T, typename D, typename R >
bool test_pattern_generator_random_word(Eigen::Matrix< T, 1, Eigen::Dynamic > &word, D &dist, R &rng, const uint64_t length, const int64_t n_samples_remaining)
{
    word.resize(length);
    for(uint64_t c = 0; c < length; c++)
        word(c) = dist(rng);
    return n_samples_remaining > 0;
}

// read an ECC code configuration out of the provided JSON configuration file
rapidjson::Document read_json_cfg_file(const std::string &ecc_code_json_cfg_fname);

// test pattern generator for creating permutations of a base n-hot test pattern
void test_pattern_generator_k_ones_in_n_bits_initialize(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones);
bool test_pattern_generator_k_ones_in_n_bits_advance(std::vector< bool > &test_pattern, uint64_t K, uint64_t k_ones);
void tpg_test(void);

// misc
std::string get_filename_uniqifier_from_time(void);
void update_progress(float progress, std::string prefix, std::string suffix);
double binom(int n, int k);

#endif /* BEER_UTILS_H */
