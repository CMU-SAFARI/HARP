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
#include <unordered_set>
#include <algorithm>
#include <random>
#include <fstream>
#include <iomanip>
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

/** the SAT problem is algorithmically nontrivial, and certain cases cause the
 * solver extreme difficult. in these cases, we generally discard the sample
 * based on this timeout value - empirically, we find 1 minute to be reasonable */
#define Z3_TIMEOUT_MS 60000

/** the different data patterns used by the profilers */
enum data_pattern
{
      DP_RANDOM    /**< uniform-random values programmed to each bit */
    , DP_ALL_ONES  /**< 0xFF */
    , DP_COLSTRIPE /**< 0xAA/0x55 */
    , DP_BEEP      /**< variable data pattern using the SAT-solver-based methodology proposed in the paper */
    , DP_BEEP2     /**< identical to BEEP, but distinguished because BEEP2 is a separate configuration */
};

/** convenience tool to allow pretty printing of the data patterns */
std::unordered_map< enum data_pattern, std::string > dp_to_str_map({
      {DP_RANDOM, "RANDOM"}
    , {DP_ALL_ONES, "ALL_ONES"}
    , {DP_COLSTRIPE, "COLSTRIPE"}
    , {DP_BEEP, "BEEP"}
    , {DP_BEEP2, "BEEP2"}
});

/** defining the Bernoulli probability of failure for each at-risk pre-correction bit */
enum per_bit_error_model
{
      PBEM_100
    , PBEM_75
    , PBEM_50
    , PBEM_25
};

/** convenience tool to allow pretty printing of the per-bit error model */
std::unordered_map< enum per_bit_error_model, std::string > pbem_to_str_map(
     {{PBEM_100, "PBEM_100"}
    , {PBEM_75, "PBEM_75"}
    , {PBEM_50, "PBEM_50"}
    , {PBEM_25, "PBEM_25"}}
);

/**
 * @brief inject errors into at-risk bits based on the per-bit error probabilities
 * 
 * @param codeword the ECC codeword to inject errors into
 * @param per_bit_probabilities vector of random distributions that implement the per-bit error probabilities
 * @param codeword_p the returned codeword containing errors
 * @param rng pre-seeded random number generator for use with the distributions
 */
void inject(
      const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword
    , std::vector< std::bernoulli_distribution > &per_bit_probabilities
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword_p
    , std::mt19937 &rng)
{
    codeword_p = codeword;
    for(Eigen::Index cw_bit_idx = 0; cw_bit_idx < codeword.cols(); cw_bit_idx++)
        if(codeword(cw_bit_idx) && per_bit_probabilities.at(cw_bit_idx)(rng))
            codeword_p(cw_bit_idx) =  !codeword(cw_bit_idx);
}

/**
 * @brief generates the data pattern to write into an ECC dataword based on the
 * input parameters
 *
 * @param dw the dataword to program
 * @param dp the data pattern to use
 * @param invert_dp whether to invert the data pattern (e.g., in alternating profiling rounds)
 * @param round index of the profiling round (used primiarly for BEEP)
 * @param beep_profile_precorrection the currently known errors (for use with BEEP)
 * @param uhc unknow Hamming code instance for use with BEEP's SAT solving requirements
 * @param rng pre-seeded random number generator
 */
void dpgen(
      Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dw
    , const enum data_pattern dp
    , const bool invert_dp
    , const uint64_t round
    , const std::map< uint64_t /* bit idx */, uint64_t /* round */ > &beep_profile_precorrection
    , const unk_hamming_code &uhc
    , std::mt19937 &rng)
{
    static std::bernoulli_distribution dpgen_bdist_0_5(0.5);
    switch(dp)
    {
        case DP_RANDOM:
            if(invert_dp)
                for(Eigen::Index i = 0; i < dw.cols(); i++)
                    dw(i) = !dw(i);
            else
                for(Eigen::Index i = 0; i < dw.cols(); i++)
                    dw(i) = dpgen_bdist_0_5(rng);
            break;

        case DP_ALL_ONES:
            if(invert_dp)
                dw = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Zero(1, dw.cols());
            else
                dw = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Ones(1, dw.cols());
            break;

        case DP_COLSTRIPE:
            for(Eigen::Index i = 0; i < dw.cols(); i++)
                dw(i) = invert_dp ? !(i % 2) : !!(i % 2);
            break;

        // our BEEP implementation uses a colstripe to begin introducing miscorrections
        case DP_BEEP:
        case DP_BEEP2: // they're technically the exact same from the algorithm prespective, just differentiated for use as different mechanisms
        {
            if(beep_profile_precorrection.size() == 0)
                for(Eigen::Index i = 0; i < dw.cols(); i++)
                    // dw(i) = invert_dp ? !(i % 2) : !!(i % 2);
                    dw(i) = dpgen_bdist_0_5(rng);
            else // miscorrections have occurred!
            {
                const uint64_t N = uhc.get_n_code_bits();
                if(beep_profile_precorrection.size() == N)
                    break; // we are done profiling - everything has failed  -> doesn't matter what we test now

                // don't try to test bits we already know - that would artificially make BEEP look bad
                uint64_t bit_under_test = round % N;
                while(beep_profile_precorrection.count(bit_under_test) > 0)
                    bit_under_test = (bit_under_test + 1) % N;

                std::unordered_set< uint64_t > known_prec_errors;
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_precorrection) 
                    known_prec_errors.insert(kv_pair.first);
                
                bool pattern_found = uhc.find_beep_pattern(known_prec_errors, dw, bit_under_test, true) || uhc.find_beep_pattern(known_prec_errors, dw, bit_under_test, false); // counting on short circuit behavior
                if(!pattern_found)
                {
                    for(Eigen::Index i = 0; i < dw.cols(); i++)
                        if(known_prec_errors.count(i) == 0)
                            dw(i) = dpgen_bdist_0_5(rng);
                    for(const uint64_t error_pos : known_prec_errors)
                        if(error_pos < (uint64_t)dw.cols())
                        {
                            dw(error_pos) = 1; // true cell assumption
                            if(error_pos > 0 && known_prec_errors.count(error_pos - 1) == 0)
                                dw(error_pos - 1) = 1; // true cell assumption
                            if((error_pos < ((uint64_t)dw.cols() - 1)) && known_prec_errors.count(error_pos + 1) == 0)
                                dw(error_pos + 1) = 1; // true cell assumption
                        }
                }
                // else
                //     std::cout << "[INFO] round: " << round << " BEEP pattern found: " << dw << std::endl;
            }
            
            break;
        }

        default:
            assert(0 && "illegal dp");
    }
}

/**
 * @brief creates a vector of random distributions that correspond to the
 * per-bit error model
 *
 * @param error_positions locations of at-risk bits
 * @param pbem per-bit error model under simulation
 * @param N length of an ECC codeword
 * @param per_bit_probabilities returned vector of distributions
 */
void initialize_per_bit_probabilities(
      const std::unordered_set< uint64_t > &error_positions
    , const enum per_bit_error_model pbem
    , const uint64_t N
    , std::vector< std::bernoulli_distribution > &per_bit_probabilities)
{
    for(uint64_t cw_bit_idx = 0; cw_bit_idx < N; cw_bit_idx++)
    {
        if(error_positions.count(cw_bit_idx) > 0)
            switch(pbem)
            {
                case PBEM_100:
                    per_bit_probabilities.emplace_back(1.0);
                    break;

                case PBEM_75:
                    per_bit_probabilities.emplace_back(0.75);
                    break;

                case PBEM_50:
                    per_bit_probabilities.emplace_back(0.5);
                    break;

                case PBEM_25:
                    per_bit_probabilities.emplace_back(0.25);
                    break;

                default:
                    assert(0 && "unsupported pbem");
            }
        else
            per_bit_probabilities.emplace_back(0);
    }
}

// calculate output probabilites as raw pattern counts that then must be normalized
// each (row, col) corresponds to a different truth table as the :
//                                     |     Outcome (effective probability)     | Observation
//         Probability                 | P[000] | P[001] | P[010] | ... | P[111] | -> must convert to k
//     --------------------------------+--------+--------+--------+-----+--------+
//      0: ERRMASK = 000 r^0 * (1-r)^3 |        |        |        |     |        | 
//      1: ERRMASK = 001 r^1 * (1-r)^2 |        |        |        |     |        |
//      2: ERRMASK = 010 r^1 * (1-r)^2 |        |        |        |     |        |
//      3: ERRMASK = 011 r^2 * (1-r)^1 |        |        |        |     |        |
//      4: ERRMASK = 100 r^1 * (1-r)^2 |        |        |        |     |        |
//      5: ERRMASK = 101 r^2 * (1-r)^1 |        |        |        |     |        |
//      6: ERRMASK = 110 r^2 * (1-r)^1 |        |        |        |     |        |
//      7: ERRMASK = 111 r^3 * (1-r)^0 |        |        |        |     |        |
//                           total count across 2D grid = 2^(n_at_risk_bits_prec + K)
// normalization constant should be 2^(n_at_risk_bits_prec + K)
void calculate_probability_of_k_bit_missed_pattern_given_n_raw_errors(
      const hamming_code &hc
    , const std::unordered_set< uint64_t > &precorrection_error_positions
    , const std::unordered_set< uint64_t > &known_postcorrection_errors
    , std::map< uint64_t /* n missed bits postcorrection */, std::map< uint64_t /* n raw errors (outcome) */, std::map< std::pair< int, int > /* binomial probaiblity term */, uint64_t /* wordcount */ > > > &probability_of_k_bit_missed_pattern_given_n_raw_errors
    , std::bernoulli_distribution &bdist
    , std::mt19937 &rng
    , int64_t random_samples = 0) // 0 -> don't take random samples
{
    const uint64_t K = hc.get_n_data_bits();
    // const uint64_t normalization_constant = random_samples > 0 ? (uint64_t)std::log2(random_samples) : K;


    // create ordered list of prec_error_positions
    const std::vector< uint64_t > iterable_precorrection_error_positions(precorrection_error_positions.begin(), precorrection_error_positions.end());

    // Compute lookup table of 'effective error pattern' to 'n postcorrection errors'
    std::size_t n_at_risk_bits_prec = precorrection_error_positions.size();
    std::unordered_map< std::vector< bool >, uint64_t > effective_error_pattern_to_n_missed_posc_errors;
    for(uint64_t n_prec_errors = 0; n_prec_errors <= n_at_risk_bits_prec; n_prec_errors++)
    {
        std::vector< bool > candidate_error_pattern;
        test_pattern_generator_k_ones_in_n_bits_initialize(candidate_error_pattern, n_at_risk_bits_prec, n_prec_errors);
        do
        {
            // calculate post-correction errors with this particular error pattern (irrespective of data pattern -> syndrome is independent of data pattern)
            // assuming systematic encoding, parity-check bits can only affect data bits via the syndrome (which is indep of data pattern :)
            std::unordered_set< uint64_t > precorrection_error_positions_this_ecc_word;
            for(uint64_t c = 0; c < candidate_error_pattern.size(); c++)
                if(candidate_error_pattern.at(c))
                    precorrection_error_positions_this_ecc_word.insert(iterable_precorrection_error_positions.at(c));
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Zero(K);
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p;
            hc.calculate_dataword_p(precorrection_error_positions_this_ecc_word, dataword, dataword_p);
            int n_missed_posc_errors = 0;
            for(uint64_t k = 0; k < K; k++)
                if(dataword_p(k) && known_postcorrection_errors.count(k) == 0) // i.e., count number of unrepaired errors
                    n_missed_posc_errors++; 
            effective_error_pattern_to_n_missed_posc_errors[candidate_error_pattern] = n_missed_posc_errors;                   
        } while(test_pattern_generator_k_ones_in_n_bits_advance(candidate_error_pattern, n_at_risk_bits_prec, n_prec_errors));
    }

    // we assume random data -> all data patterns are uniformly possible
    // for a particular per-bit error probability, the acutal errors themselves are less likely
    //      however, this is accounted for by the split into n_prec_errors -> since all bits are IID, we can say that 
    //      a particular pre-correction error pattern occurs with probability based on that exact nubmer of errors ocurring
    // FORALL(codeword) -> do something
    bool use_random_sampling = random_samples > 0;
    std::map< uint64_t /* n missed posc errors */, std::map< std::vector< bool > /* prob component */, std::map< std::vector< bool > /* outcome */, uint64_t > > > probability_as_raw_counts;
    for(uint64_t n_ones = 0; use_random_sampling ? random_samples > 0 : n_ones <= K; n_ones++)
    {
        // initialize the dataword
        Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword;
        if(use_random_sampling)
            test_pattern_generator_random_word(dataword, bdist, rng, K, random_samples--);
        else
            test_pattern_generator_k_ones_in_n_bits_initialize(dataword, K, n_ones);

        // repeatedly generate datawords until we've exhausted all requested patterns
        do
        {
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword;
            hc.encode(dataword, codeword);
            // std::cout << "[DEBUG] generated dataword: " << dataword << std::endl;
            
            // FORALL(error pattern)
            for(uint64_t n_prec_errors = 0; n_prec_errors <= n_at_risk_bits_prec; n_prec_errors++)
            {
                std::vector< bool > candidate_error_pattern;
                test_pattern_generator_k_ones_in_n_bits_initialize(candidate_error_pattern, n_at_risk_bits_prec, n_prec_errors);
                do
                {
                    // ASSUME only data 1 can fail (i.e., all true cells) to get the 'effective' error pattern component
                    std::vector< bool > effective_error_pattern;
                    for(uint64_t err_idx = 0; err_idx < n_at_risk_bits_prec; err_idx++)
                    {
                        // bit is charged and candidate pattern says this bit failed
                        bool bit_is_charged = codeword(iterable_precorrection_error_positions.at(err_idx));
                        effective_error_pattern.push_back(bit_is_charged && candidate_error_pattern.at(err_idx));
                    }
                    
                    // add to probability table
                    const uint64_t n_missed_posc_errors = effective_error_pattern_to_n_missed_posc_errors.at(effective_error_pattern);
                    if(probability_as_raw_counts.count(n_missed_posc_errors) == 0)
                        probability_as_raw_counts[n_missed_posc_errors] =  std::map< std::vector< bool >, std::map< std::vector< bool >, uint64_t > >();
                    if(probability_as_raw_counts.at(n_missed_posc_errors).count(candidate_error_pattern) == 0)
                        probability_as_raw_counts.at(n_missed_posc_errors)[candidate_error_pattern] =  std::map< std::vector< bool >, uint64_t >();
                    if(probability_as_raw_counts.at(n_missed_posc_errors).at(candidate_error_pattern).count(effective_error_pattern) == 0)
                        probability_as_raw_counts.at(n_missed_posc_errors).at(candidate_error_pattern)[effective_error_pattern] = 0;
                    probability_as_raw_counts.at(n_missed_posc_errors).at(candidate_error_pattern).at(effective_error_pattern)++;
                } while(test_pattern_generator_k_ones_in_n_bits_advance(candidate_error_pattern, n_at_risk_bits_prec, n_prec_errors));
            }
        } while(use_random_sampling
            ? test_pattern_generator_random_word(dataword, bdist, rng, K, random_samples--) 
            : test_pattern_generator_k_ones_in_n_bits_advance(dataword, K, n_ones));
    }

    // prett-print the raw probability counts
    // std::cout << "[DEBUG] P[x-bit missed pattern | " << n_at_risk_bits_prec << " raw errors] -> normalizer should be: 2^(" << normalization_constant << " + " << n_at_risk_bits_prec << ") : {" << std::endl;
    // for(const std::pair< uint64_t, std::map< std::vector< bool >, std::map< std::vector< bool >, uint64_t > > > &kv_pair_nposcerrors : probability_as_raw_counts)
    // {
    //     const uint64_t n_missed_posc_errors = kv_pair_nposcerrors.first;
    //     std::cout << "    " << n_missed_posc_errors << "-bit missed pattern : {" << std::endl;
    //     for(const std::pair< std::vector< bool >, std::map< std::vector< bool >, uint64_t > > &kv_pair_prob : kv_pair_nposcerrors.second)
    //     {
    //         std::cout << "        "  << vec_to_string(kv_pair_prob.first) << " (probability) : " << std::endl;
    //         for(const std::pair< std::vector< bool >, uint64_t > &kv_pair_outcome : kv_pair_prob.second)
    //             std::cout << "            " << vec_to_string(kv_pair_outcome.first) << " (outcome) : " << kv_pair_outcome.second << std::endl;
    //     }
    //     std::cout << "    }" << std::endl;
    // }
    // std::cout << "}" << std::endl;

    // transform the raw probability counts into k-bit error probability distributions
    for(const std::pair< uint64_t, std::map< std::vector< bool >, std::map< std::vector< bool >, uint64_t > > > &kv_pair_nposcerrors : probability_as_raw_counts)
    {
        const uint64_t n_missed_posc_errors = kv_pair_nposcerrors.first;
        if(probability_of_k_bit_missed_pattern_given_n_raw_errors.count(n_missed_posc_errors) == 0)
            probability_of_k_bit_missed_pattern_given_n_raw_errors[n_missed_posc_errors] = std::map< uint64_t, std::map< std::pair< int, int >, uint64_t > >();

        // 'sum' the individual probability terms
        for(const std::pair< std::vector< bool >, std::map< std::vector< bool >, uint64_t > > &kv_pair_prob : kv_pair_nposcerrors.second)
        {
            const std::vector< bool > &probability_component = kv_pair_prob.first;
            const int val1 = std::count(probability_component.begin(), probability_component.end(), true);
            const int val0 = std::count(probability_component.begin(), probability_component.end(), false);
            const std::pair< int, int > this_term(val0, val1);
            for(const std::pair< std::vector< bool >, uint64_t > &kv_pair_outcome : kv_pair_prob.second)
            {
                const std::vector< bool > &outcome = kv_pair_outcome.first;
                const uint64_t n_raw_bit_errors = std::count(outcome.begin(), outcome.end(), true);
                const uint64_t wordcount = kv_pair_outcome.second;

                std::map< uint64_t, std::map< std::pair< int, int >, uint64_t > > &this_probability_table = probability_of_k_bit_missed_pattern_given_n_raw_errors.at(n_missed_posc_errors);
                if(this_probability_table.count(n_raw_bit_errors) == 0)
                    this_probability_table[n_raw_bit_errors] = std::map< std::pair< int, int >, uint64_t >();
                if(this_probability_table.at(n_raw_bit_errors).count(this_term) == 0)
                    this_probability_table.at(n_raw_bit_errors)[this_term] = 0;
                this_probability_table.at(n_raw_bit_errors).at(this_term) += wordcount;
            }
        }
    }

    // pretty-print the final probability table
    // std::cout << "[DEBUG] P[x-bit missed pattern | " << n_at_risk_bits_prec << " raw errors] -> normalizer should be: 2^(" << normalization_constant << " + " << n_at_risk_bits_prec << ") : {" << std::endl;
    // for(const std::pair< uint64_t, std::map< uint64_t, std::map< std::pair< int, int >, uint64_t > > > &kv_pair_nposcerrors : probability_of_k_bit_missed_pattern_given_n_raw_errors)
    // {
    //     const uint64_t n_missed_posc_errors = kv_pair_nposcerrors.first;
    //     std::cout << "    " << n_missed_posc_errors << "-bit missed pattern : {" << std::endl;
        
    //     for(const std::pair< uint64_t, std::map< std::pair< int, int >, uint64_t > > &kv_pair_n_raw_errors : kv_pair_nposcerrors.second)
    //     {
    //         const uint64_t n_raw_bit_errors = kv_pair_n_raw_errors.first;
    //         std::cout << "        " << n_raw_bit_errors << " raw bit errors : [";
    //         for(const std::pair< std::pair< int, int >, uint64_t > &kv_pair_term : kv_pair_n_raw_errors.second)
    //         {
    //             const std::pair< int, int > &binomial_probability_term = kv_pair_term.first;
    //             const uint64_t wordcount = kv_pair_term.second;
    //             std::cout << " z" << binomial_probability_term.first << "o" << binomial_probability_term.second << ":" << wordcount;
    //         }
    //         std::cout << " ]" << std::endl;
    //     }
    // }
    // std::cout << "}" << std::endl;
}

/**
 * @brief run an entire round of evaluations on a single ECC word
 * 
 * @param hc the known hamming code under simulation
 * @param uhc the unknown hamming code used for SAT problem evaluation
 * @param precorrection_error_positions bit-indices of at-risk pre-correction bits
 * @param postcorrection_error_positions bit-indices of at-risk post-correction bits
 * @param rng pre-seeded random number generator 
 */
void evaluate_ecc_word(
      const hamming_code &hc
    , const unk_hamming_code &uhc
    , const std::unordered_set< uint64_t > &precorrection_error_positions
    , const std::unordered_set< uint64_t > &postcorrection_error_positions
    , std::mt19937 &rng)
{
    const uint64_t N = hc.get_n_code_bits();
    const uint64_t K = hc.get_n_data_bits();
    std::bernoulli_distribution bdist(0.5);

    std::vector< std::pair< enum data_pattern, bool /* invertible */ > > dp_list = 
    {
          {DP_RANDOM,       1} // ONLY simulate HARP (not HARP+)
        // , {DP_ALL_ONES,     1}
        // , {DP_COLSTRIPE,    1}
        , {DP_BEEP,         1} 
        , {DP_BEEP2,        1} // ONLY simulate BEEP2, which assumes we *already* have full direct error coverage
    };
    std::vector< enum per_bit_error_model > pbem_list = 
    {
          PBEM_100
        , PBEM_75
        , PBEM_50
        , PBEM_25
    };
    const int n_rounds = 128;
    for(enum per_bit_error_model pbem : pbem_list)
    {
        std::vector< std::bernoulli_distribution > per_bit_probabilities;
        initialize_per_bit_probabilities(precorrection_error_positions, pbem, N, per_bit_probabilities);
            

        for(std::pair< enum data_pattern, bool > dp_info : dp_list)
        {
            enum data_pattern dp = dp_info.first;
            bool invertible = dp_info.second;

            // profiler data
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > beep_profile_prec; // PRECORRECTION
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > beep_profile_posc; // POSTCORRECTION
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > beep2_profile_prec; // PRECORRECTION
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > beep2_profile_posc; // POSTCORRECTION
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > harp_profile_prec; // PRECORRECTION
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > harp_p_profile_posc; // POSTCORRECTION (based on HARP profile)
            std::map< uint64_t /* bit idx */, uint64_t /* round */ > naive_profile_posc; // POSTCORRECTION 

            // BEEP2 *already* knows all *direct errors* - that's the guarantee -> as if we ran perfect HARP+ active profiling
            for(const uint64_t bit_idx : precorrection_error_positions)
                if(bit_idx < hc.get_n_data_bits()) // relying on systematic encoding for correctness of this conditional check!!
                    beep2_profile_prec[bit_idx] = -1; // bit is known before profiling even starts!
            for(const uint64_t bit_idx : postcorrection_error_positions)
                if(precorrection_error_positions.count(bit_idx) != 0)
                    beep2_profile_posc[bit_idx] = -1; // bit is known before profiling even starts!
            { // just a scoping operator for code clarity
                std::unordered_set< uint64_t > prec_errs;
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_prec) 
                    prec_errs.insert(kv_pair.first);
                std::unordered_set< uint64_t > posc_errs;
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_posc) // if we already know some from before, don't recompute them
                    posc_errs.insert(kv_pair.first);

                uhc.solver_reset();
                uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                bool timeout = uhc.find_possible_postcorrection_errors(prec_errs, posc_errs);
                if(timeout)
                    std::cout << "[WARN] timeout computing BEEP2 direct error set" << std::endl;
                
                // add the newly-discovered
                for(const uint64_t k : posc_errs)
                    if(beep2_profile_posc.count(k) == 0)
                        beep2_profile_posc[k] = -1;
            }

            // histogram data
            std::map< uint64_t /* round_idx */, uint64_t /* max_simultaneous_unknown_posc_errors_possible */ > beep_posc_histograms;
            std::map< uint64_t /* round_idx */, uint64_t /* max_simultaneous_unknown_posc_errors_possible */ > beep2_posc_histograms;
            std::map< uint64_t /* round_idx */, uint64_t /* max_simultaneous_unknown_posc_errors_possible */ > harp_prec_histograms;
            std::map< uint64_t /* round_idx */, uint64_t /* max_simultaneous_unknown_posc_errors_possible */ > harp_posc_histograms;
            std::map< uint64_t /* round_idx */, uint64_t /* max_simultaneous_unknown_posc_errors_possible */ > naive_histograms;

            if(g_verbosity >= 1)
                std::cout << "[INFO]";
            else
                std::cout << "[DATA]";
            std::cout << " pbem: " << pbem_to_str_map.at(pbem) << " dp: " << dp_to_str_map.at(dp) << " n_rounds: " << n_rounds;
            if(g_verbosity >= 1)
                std::cout << std::endl;
            
            // dataword must persist across iterations for random/invrandom pattern to work correctly
            bool hist_timeout = false;
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword(1, K);
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword(1, N);
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > codeword_p(1, N);
            for(int r = 0; r < n_rounds; r++)
            {
                // create a dataword with the appropriate test pattern
                bool invert_dp = invertible ? !!(r % 2) : false;
                dpgen(dataword, dp, invert_dp, r, dp == DP_BEEP ? beep_profile_prec : beep2_profile_prec, uhc, rng);

                // encode, inject errors
                hc.encode(dataword, codeword);
                inject(codeword, per_bit_probabilities, codeword_p, rng);
                
                // decode
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > corrected_codeword;
                int corrected_col;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p;
                Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > syndrome;
                hc.decode(codeword_p, corrected_codeword, corrected_col, dataword_p, syndrome);

                // debug
                if(g_verbosity >= 2)
                {
                    std::cout << "[DEBUG2] pbem: " << pbem << " dp: " << dp << " round: " << r << std::endl;
                    std::cout << "[DEBUG2]     dataword:      " << dataword << std::endl;
                    std::cout << "[DEBUG2]     codeword:      " << codeword << std::endl;
                    std::cout << "[DEBUG2]     codeword_p:    " << codeword_p << std::endl;
                    std::cout << "[DEBUG2]     dataword_p:    " << dataword_p << std::endl;
                    std::cout << "[DEBUG2]     syndrome:      " << syndrome << std::endl;
                    std::cout << "[DEBUG2]     corrected_col: " << corrected_col << std::endl;
                }

                // calculate profiling stuffs
                if(dp == DP_BEEP)
                {
                    // updated ONLY ON MISCORRECTION
                    bool miscorrection_exists = false;
                    for(uint64_t k = 0; k < K; k++)
                        if(!dataword(k) && dataword_p(k)) // 0 -> 1 error == miscorrection
                        {
                            miscorrection_exists = true;
                            break;
                        }
                        
                    bool new_beep_bit = false;
                    if(miscorrection_exists)
                    {
                        for(uint64_t n = 0; n < N; n++)
                            if((codeword(n) != codeword_p(n)) && (beep_profile_prec.count(n) == 0))
                            {
                                beep_profile_prec[n] = r;
                                new_beep_bit = true;
                            }

                        // determine BEEP posc - all possible combinations of the new bit added
                        if(new_beep_bit)
                        {
                            std::unordered_set< uint64_t > prec_errs;
                            for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_prec) 
                                prec_errs.insert(kv_pair.first);
                            std::unordered_set< uint64_t > posc_errs;
                            for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_posc) // if we already know some from before, don't recompute them
                                posc_errs.insert(kv_pair.first);

                            uhc.solver_reset();
                            uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                            bool timeout = uhc.find_possible_postcorrection_errors(prec_errs, posc_errs);
                            if(timeout)
                            {
                                if(g_verbosity == 0)
                                    std::cout << std::endl;
                                std::cout << "[WARN] timeout computing HARP+ set" << std::endl;
                                if(g_verbosity == 0)
                                    std::cout << "[DATA]";
                                continue;
                            }
                            
                            // add the BEEP posc bits
                            for(const uint64_t k : posc_errs)
                                if(beep_profile_posc.count(k) == 0)
                                    beep_profile_posc[k] = r;
                        }
                    }

                    // update the BEEP (posc) histogram
                    if(!hist_timeout && (r == 0 || new_beep_bit))
                    {
                        uint64_t hist = -1u;
                        std::unordered_set< uint64_t > cur_profile;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_posc)
                            cur_profile.insert(kv_pair.first); // add the bit index
                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        hist_timeout = uhc.get_histogram_of_n_simultaneous_errors(precorrection_error_positions, postcorrection_error_positions, cur_profile, hist);
                        if(!hist_timeout)
                            beep_posc_histograms[r] = hist;
                    }
                }
                else if(dp == DP_BEEP2)
                {
                    // updated ONLY ON MISCORRECTION
                    bool miscorrection_exists = false;
                    for(uint64_t k = 0; k < K; k++)
                        if(!dataword(k) && dataword_p(k)) // 0 -> 1 error == miscorrection
                        {
                            miscorrection_exists = true;
                            break;
                        }
                        
                    bool new_beep_bit = false;
                    if(miscorrection_exists)
                    {
                        for(uint64_t n = 0; n < N; n++)
                            if((codeword(n) != codeword_p(n)) && (beep2_profile_prec.count(n) == 0))
                            {
                                beep2_profile_prec[n] = r;
                                new_beep_bit = true;
                            }

                        // determine BEEP posc - all possible combinations of the new bit added
                        if(new_beep_bit)
                        {
                            std::unordered_set< uint64_t > prec_errs;
                            for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_prec) 
                                prec_errs.insert(kv_pair.first);
                            std::unordered_set< uint64_t > posc_errs;
                            for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_posc) // if we already know some from before, don't recompute them
                                posc_errs.insert(kv_pair.first);

                            uhc.solver_reset();
                            uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                            bool timeout = uhc.find_possible_postcorrection_errors(prec_errs, posc_errs);
                            if(timeout)
                            {
                                if(g_verbosity == 0)
                                    std::cout << std::endl;
                                std::cout << "[WARN] timeout computing HARP+ set" << std::endl;
                                if(g_verbosity == 0)
                                    std::cout << "[DATA]";
                                continue;
                            }
                            
                            // add the BEEP posc bits
                            for(const uint64_t k : posc_errs)
                                if(beep2_profile_posc.count(k) == 0)
                                    beep2_profile_posc[k] = r;
                        }
                    }

                    // update the BEEP (posc) histogram
                    if(!hist_timeout && (r == 0 || new_beep_bit))
                    {
                        uint64_t hist = -1u;
                        std::unordered_set< uint64_t > cur_profile;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_posc)
                            cur_profile.insert(kv_pair.first); // add the bit index
                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        hist_timeout = uhc.get_histogram_of_n_simultaneous_errors(precorrection_error_positions, postcorrection_error_positions, cur_profile, hist);
                        if(!hist_timeout)
                            beep2_posc_histograms[r] = hist;
                    }
                }
                else
                {
                    // Naive profiling
                    bool new_naive_bit = false;
                    for(uint64_t k = 0; k < K; k++)
                        if(dataword(k) != dataword_p(k) && naive_profile_posc.count(k) == 0)
                        {
                            naive_profile_posc[k] = r;
                            new_naive_bit = true;
                        }

                    // HARP
                    bool new_harp_bit = false;
                    bool new_harp_p_bit = false;
                    for(uint64_t k = 0; k < K; k++)
                        if(dataword(k) != codeword_p(k) && harp_profile_prec.count(k) == 0)
                        {
                            harp_profile_prec[k] = r;
                            new_harp_bit = true;

                            // HARP+ includes all of HARP! Required to avoid miscorrection 2-bit error case where one bit is in data region, another in parity-check
                            if(harp_p_profile_posc.count(k) == 0)
                            {
                                harp_p_profile_posc[k] = r; 
                                new_harp_p_bit = true;
                            }
                        }

                    // determine HARP+ - all possible combinations of the new bit added
                    if(new_harp_bit || new_harp_p_bit)
                    {
                        std::unordered_set< uint64_t > prec_errs;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_profile_prec) 
                            prec_errs.insert(kv_pair.first);
                        std::unordered_set< uint64_t > posc_errs;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_p_profile_posc) // if we already know some from before, don't recompute them
                            posc_errs.insert(kv_pair.first);

                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        bool timeout = uhc.find_possible_postcorrection_errors(prec_errs, posc_errs);
                        if(timeout)
                        {
                            if(g_verbosity == 0)
                                std::cout << std::endl;
                            std::cout << "[WARN] timeout computing HARP+ set" << std::endl;
                            if(g_verbosity == 0)
                                std::cout << "[DATA]";
                            continue;
                        }
                        
                        // add the HARP+ bits
                        for(const uint64_t k : posc_errs)
                            if(harp_p_profile_posc.count(k) == 0)
                            {
                                harp_p_profile_posc[k] = r;
                                new_harp_p_bit = true;
                            }
                    }

                    // update the naive histogram histogram
                    if(!hist_timeout && (r == 0 || new_naive_bit))
                    {
                        uint64_t hist = -1u;
                        std::unordered_set< uint64_t > cur_profile;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : naive_profile_posc)
                            cur_profile.insert(kv_pair.first); // add the bit index
                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        hist_timeout = uhc.get_histogram_of_n_simultaneous_errors(precorrection_error_positions, postcorrection_error_positions, cur_profile, hist);
                        if(!hist_timeout)
                            naive_histograms[r] = hist;
                    }

                    // update the HARP histogram
                    if(!hist_timeout && (r == 0 || new_harp_bit))
                    {
                        uint64_t hist = -1u;
                        std::unordered_set< uint64_t > cur_profile;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_profile_prec)
                            cur_profile.insert(kv_pair.first); // add the bit index
                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        hist_timeout = uhc.get_histogram_of_n_simultaneous_errors(precorrection_error_positions, postcorrection_error_positions, cur_profile, hist);
                        if(!hist_timeout)
                            harp_prec_histograms[r] = hist;
                    }

                    // update the HARP+ histogram
                    if(!hist_timeout && (r == 0 || new_harp_p_bit))
                    {
                        uint64_t hist = -1u;
                        std::unordered_set< uint64_t > cur_profile;
                        for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_p_profile_posc)
                            cur_profile.insert(kv_pair.first); // add the bit index
                        uhc.solver_reset();
                        uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
                        hist_timeout = uhc.get_histogram_of_n_simultaneous_errors(precorrection_error_positions, postcorrection_error_positions, cur_profile, hist);
                        if(!hist_timeout)
                            harp_posc_histograms[r] = hist;
                    }
                }

                // check for early termination if all bits have been found by all profilers
                bool all_found = true;
                if(dp == DP_BEEP)
                {
                    for(const uint64_t bidx : precorrection_error_positions)
                        if(beep_profile_prec.count(bidx) == 0)
                        {
                            all_found = false;
                            break;
                        }
                    if(all_found == true)
                        for(const uint64_t bidx : postcorrection_error_positions)
                            if(beep_profile_posc.count(bidx) == 0)
                            {
                                all_found = false;
                                break;
                            }
                }
                else if(dp == DP_BEEP2)
                {
                    for(const uint64_t bidx : precorrection_error_positions)
                        if(beep2_profile_prec.count(bidx) == 0)
                        {
                            all_found = false;
                            break;
                        }
                    if(all_found == true)
                        for(const uint64_t bidx : postcorrection_error_positions)
                            if(beep2_profile_posc.count(bidx) == 0)
                            {
                                all_found = false;
                                break;
                            }
                }
                else
                {
                    for(const uint64_t bidx : precorrection_error_positions)
                        if(harp_profile_prec.count(bidx) == 0)
                        {
                            all_found = false;
                            break;
                        }
                    if(all_found == true)
                        for(const uint64_t bidx : postcorrection_error_positions)
                            if((harp_p_profile_posc.count(bidx) == 0) || (naive_profile_posc.count(bidx) == 0))
                            {
                                all_found = false;
                                break;
                            }
                }
                if(all_found)
                {
                    // std::cout << "[INFO] early termination" << std::endl;
                    break;
                }
            }

            if(g_verbosity >= 1)
                std::cout << "[DATA]";
            std::cout << " profiles ";
            if(dp == DP_BEEP)
            {
                std::cout << "b: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_prec)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
                std::cout << "} +: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_profile_posc)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
            }
            else if(dp == DP_BEEP2)
            {
                std::cout << "2: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_prec)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
                std::cout << "} +: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_profile_posc)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
            }
            else
            {
                std::cout << "n: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : naive_profile_posc)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
                std::cout << "} h: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_profile_prec)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
                std::cout << "} +: { ";
                for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_p_profile_posc)
                    std::cout << kv_pair.first << "," << kv_pair.second << " ";
            }
            std::cout << "}" << std::endl;

            // output histogram of *how many* simultaneous errors can occur given the profile
            // Q: can we do this without re-runnign the SAT solver for every single step? precomputation necessary?
            // Q: that would be an unreasonable prospect for sth like {5 precorrection, 20 postcorrection} errors
            std::cout << "[DATA] histograms ";
            if(!hist_timeout)
            {
                if(dp == DP_BEEP)
                {
                    std::cout << "+: { ";
                    for(const std::pair< uint64_t, uint64_t > &kv_pair : beep_posc_histograms)
                        std::cout << kv_pair.first << "," << kv_pair.second << " ";
                    std::cout << " }";
                }
                else if(dp == DP_BEEP2)
                {
                    std::cout << "2: { ";
                    for(const std::pair< uint64_t, uint64_t > &kv_pair : beep2_posc_histograms)
                        std::cout << kv_pair.first << "," << kv_pair.second << " ";
                    std::cout << " }";
                }
                else
                {
                    std::cout << "n: { ";
                    for(const std::pair< uint64_t, uint64_t > &kv_pair : naive_histograms)
                        std::cout << kv_pair.first << "," << kv_pair.second << " ";
                    std::cout << "} h: { ";
                    for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_prec_histograms)
                        std::cout << kv_pair.first << "," << kv_pair.second << " ";
                    std::cout << "} h+: { ";
                    for(const std::pair< uint64_t, uint64_t > &kv_pair : harp_posc_histograms)
                        std::cout << kv_pair.first << "," << kv_pair.second << " ";
                    std::cout << " }";
                }
            }
            std::cout << std::endl;

            // P[x-bit error pattern | n raw bit errors] -> used for computing FIT rate
            std::map< std::string, std::map< uint64_t, uint64_t > > profilers_to_print;
            if(dp == DP_BEEP)
                profilers_to_print["b"] = beep_profile_posc;
            else if(dp == DP_BEEP2)
                profilers_to_print["2"] = beep2_profile_posc;
            else
            {
                profilers_to_print["n"] = naive_profile_posc;
                profilers_to_print["h"] = harp_profile_prec;
                profilers_to_print["+"] = harp_p_profile_posc;
            }

            for(const std::pair< std::string, std::map< uint64_t, uint64_t > > &kv_pair_profiler : profilers_to_print)
            {
                const std::string &profiler_name = kv_pair_profiler.first;
                const std::map< uint64_t, uint64_t > &profile_allrounds = kv_pair_profiler.second;

                // for(uint64_t n_rounds : {32, 64, 128})
                std::set< uint64_t > round_indices_to_consider;
                for(const std::pair< uint64_t, uint64_t > &profile_entry : profile_allrounds)
                    round_indices_to_consider.insert(profile_entry.second);
                if(round_indices_to_consider.size() == 0)
                    round_indices_to_consider.insert(127);
                if(round_indices_to_consider.count(0) == 0) // this is necessary to get the round 0 bits for those profilers that don't find the first bit until round_idx>0
                    round_indices_to_consider.insert(0); 
                for(const uint64_t max_round_idx_to_consider : round_indices_to_consider)
                {
                    std::unordered_set< uint64_t > known_postcorrection_at_risk_bits;
                    for(const std::pair< uint64_t, uint64_t > &profile_entry : profile_allrounds)
                    {
                        const uint64_t bit_idx = profile_entry.first;
                        const uint64_t round_idx = profile_entry.second;
                        if(round_idx <= max_round_idx_to_consider || round_idx == (uint64_t)-1) // supports BEEP2 bits known before profiling
                            known_postcorrection_at_risk_bits.insert(bit_idx);
                    }

                    // calculate the *probability* of each *combination* of the raw bit errors
                    std::map< uint64_t /* n missed bits postcorrection */, std::map< uint64_t /* n raw errors (outcome) */, std::map< std::pair< int, int > /* binomial probaiblity term */, uint64_t /* wordcount */ > > > probability_of_k_bit_missed_pattern_given_n_raw_errors;
                    uint64_t random_samples = 131072;
                    uint64_t normalization_constant_log2 = random_samples > 0 ? (uint64_t)std::log2(random_samples) : K;
                    calculate_probability_of_k_bit_missed_pattern_given_n_raw_errors(hc, precorrection_error_positions, known_postcorrection_at_risk_bits, probability_of_k_bit_missed_pattern_given_n_raw_errors, bdist, rng, random_samples);            
                    
                    int n_prec_errors = precorrection_error_positions.size();
                    std::cout << "[DATA] P[x|n] profiler: " << profiler_name << " n_prec_errors: " << n_prec_errors << " within_round_idx: " << max_round_idx_to_consider << " known_bits: {";
                    for(const uint64_t bit_idx : known_postcorrection_at_risk_bits)
                        std::cout << " " << bit_idx;
                    std::cout << " } log2(normalizer): " << normalization_constant_log2 << "+" << n_prec_errors << " data:";
                    for(const std::pair< uint64_t, std::map< uint64_t, std::map< std::pair< int, int >, uint64_t > > > &kv_pair_nposcerrors : probability_of_k_bit_missed_pattern_given_n_raw_errors)
                    {
                        const uint64_t n_missed_posc_errors = kv_pair_nposcerrors.first;
                        std::cout << " x=" << n_missed_posc_errors << "{";
                        
                        for(const std::pair< uint64_t, std::map< std::pair< int, int >, uint64_t > > &kv_pair_n_raw_errors : kv_pair_nposcerrors.second)
                        {
                            const uint64_t n_raw_bit_errors = kv_pair_n_raw_errors.first;
                            std::cout << n_raw_bit_errors << ":[";
                            for(const std::pair< std::pair< int, int >, uint64_t > &kv_pair_term : kv_pair_n_raw_errors.second)
                            {
                                const std::pair< int, int > &binomial_probability_term = kv_pair_term.first;
                                const uint64_t wordcount = kv_pair_term.second;
                                std::cout << "z" << binomial_probability_term.first << "o" << binomial_probability_term.second << ":" << wordcount << ",";
                            }
                            std::cout << "]";
                        }
                        std::cout << "}";
                    }
                    std::cout << std::endl;
                }
            }
        }
    }
}

/**
 * @brief run a full round of evaluations for a single ECC code instance
 * 
 * @param fname name of a JSON configuration file for the ECC code to simulate 
 * @param n_ecc_words number of ECC words to simulate
 * @param rng pre-seeded random number generator
 */
void evaluate(const std::string &fname, const int n_ecc_words, const int64_t random_seed)
{
    std::mt19937 rng(random_seed);
    
    z3::context z3_ctx;
    z3::solver s(z3_ctx); // default SMT tactic
    // z3::solver s = (z3::tactic(z3_ctx, "eq2bv")).mk_solver();
    // z3::solver s = (z3::tactic(z3_ctx, "qe") & z3::tactic(z3_ctx, "smt")).mk_solver();
    gf2::init_with_z3_context(z3_ctx);
    unk_hamming_code uhc(z3_ctx, s, fname);
    
    hamming_code hc(fname);
    const uint64_t N = hc.get_n_code_bits();
    const uint64_t K = hc.get_n_data_bits();
    std::uniform_int_distribution< uint64_t > udist(0, N - 1);

    std::cout << "[INFO] evaluating JSON ECC code configuration file: " << fname << std::endl; 
    std::cout << "[INFO] using random_seed: " << random_seed << std::endl;
    std::chrono::steady_clock::time_point time_start_code = std::chrono::steady_clock::now();

    for(size_t n_errors = 2; n_errors <= 5; n_errors++)
    {
        std::cout << "[INFO] n_errors: " << n_errors << std::endl; 
        std::chrono::steady_clock::time_point time_start_nerrs = std::chrono::steady_clock::now();

        for(int ecc_word = 0; ecc_word < n_ecc_words; ecc_word++)
        {
            std::chrono::steady_clock::time_point time_start_word = std::chrono::steady_clock::now();

            // choose the spatial positoins of the n_errors uniform-randomly
            std::unordered_set< uint64_t > precorrection_error_positions;
            while(precorrection_error_positions.size() != n_errors)
                precorrection_error_positions.insert(udist(rng));
            // for(uint64_t cw_bit_idx = 0; cw_bit_idx < N; cw_bit_idx++)
            //     if(bdist(rng))
            //         precorrection_error_positions.insert(cw_bit_idx);

            // calculate UBER when all these bits fail with a CHARGED test pattern
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Ones(1, K);
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p(1, K);
            hc.calculate_dataword_p(precorrection_error_positions, dataword, dataword_p);
            std::unordered_set< uint64_t > postcorrection_errors_with_charged_dp;
            for(uint64_t k = 0; k < K; k++)
                if(dataword(k) != dataword_p(k))
                    postcorrection_errors_with_charged_dp.insert(k);
            
            if(g_verbosity >= 1)
            {
                std::cout << "[DEBUG1] word: " << ecc_word << std::endl;
                std::cout << "[DEBUG1] precorrection_error_positions:  {";
                for(const uint64_t epos : precorrection_error_positions)
                    std::cout << epos << " ";
                std::cout << "}" << std::endl;
                std::cout << "[DEBUG1] n_postcorrection_errors with CHARGED dataword: {";
                for(const uint64_t epos : postcorrection_errors_with_charged_dp)
                    std::cout << epos << " ";
                std::cout << "}" << std::endl;
            }
            
            // calculate all theoretically possible errors (effective BER, i.e., EBER)
            std::unordered_set< uint64_t > postcorrection_error_positions; // empty - all unknown
            uhc.solver_reset();
            uhc.solver_set_param(":timeout", static_cast< unsigned >(Z3_TIMEOUT_MS)); // in milliseconds
            bool timeout = uhc.find_possible_postcorrection_errors(precorrection_error_positions, postcorrection_error_positions);
            if(timeout)
            {
                std::cout << "[WARN] skipping sample due to assumed timeout" << std::endl;
                continue;
            }

            if(g_verbosity >= 1)
            {
                std::cout << "[DEBUG1] postcorrection_error_positions: {";
                for(const uint64_t epos : postcorrection_error_positions)
                    std::cout << epos << " ";
                std::cout << "}" << std::endl;
            }

            // print out important information
            std::cout << "[DATA] word: " << ecc_word << " prec: { ";
            for(const uint64_t epos : precorrection_error_positions)
                std::cout << epos << " ";
            std::cout << "} posc_charged: { ";
            for(const uint64_t epos : postcorrection_errors_with_charged_dp)
                std::cout << epos << " ";
            std::cout << "} posc_all: { ";
            for(const uint64_t epos : postcorrection_error_positions)
                std::cout << epos << " ";
            std::cout << "}" << std::endl;

            evaluate_ecc_word(hc, uhc, precorrection_error_positions, postcorrection_error_positions, rng);
            
            std::chrono::steady_clock::time_point time_end_word = std::chrono::steady_clock::now();
            std::cout << "[INFO] elapsed time for ecc_word: " << ecc_word << " is " 
                << (std::chrono::duration_cast< std::chrono::microseconds >(time_end_word - time_start_word).count()) / 1000000.0  << std::endl;
        }
        std::chrono::steady_clock::time_point time_end_nerrs = std::chrono::steady_clock::now();
        std::cout << "[INFO] elapsed time for n_errors: " << n_errors << " is " 
            << (std::chrono::duration_cast< std::chrono::microseconds >(time_end_nerrs - time_start_nerrs).count()) / 1000000.0  << std::endl;
    }
    std::chrono::steady_clock::time_point time_end_code = std::chrono::steady_clock::now();
    std::cout << "[INFO] elapsed time for entire code is " 
        << (std::chrono::duration_cast< std::chrono::microseconds >(time_end_code - time_start_code).count()) / 1000000.0  << std::endl;
}

/**
 * @brief prints out the correspondence between pre- and post-correction error
 * probabilities for the given ECC function
 *
 * @param fname JSON configuration file for an ECC code
 * @param rng pre-seeded random number generator
 */
void analyze_ecc_code(const std::string &fname, const int n_ecc_words, const int64_t random_seed)
{
    std::mt19937 rng(random_seed);

    hamming_code hc(fname);
    const uint64_t N = hc.get_n_code_bits();
    const uint64_t K = hc.get_n_data_bits();
    std::cout << "[INFO] evaluating JSON ECC code configuration file: " << fname << std::endl; 
    std::cout << "[INFO] using random_seed: " << random_seed << std::endl;
    
    std::uniform_int_distribution< uint64_t > udist(0, N - 1);
    for(size_t n_weak_bits = 2; n_weak_bits <= 8; n_weak_bits++)
    {
        std::cout << "[INFO] n_weak_bits: " << n_weak_bits << std::endl; 
        std::map< int, int > map_initializer;
        for(uint64_t nw = 2; nw <= n_weak_bits; nw++)
            map_initializer[nw] = 0;
        
        for(int ecc_word = 0; ecc_word < n_ecc_words; ecc_word++)
        {
            std::vector< std::map< int, int > > per_bit_probabilities(K, map_initializer);

            // choose error positions for this word based on RBER
            std::vector< uint64_t > weak_bit_indices;
            while(weak_bit_indices.size() != n_weak_bits)
                weak_bit_indices.push_back(udist(rng));
            std::sort(weak_bit_indices.begin(), weak_bit_indices.end());
            
            // iterate over every combination of errors and accumulate probabilities
            for(size_t simultaneous_errors = 2; simultaneous_errors <= n_weak_bits; simultaneous_errors++)
            {
                // std::cout << "[INFO] simultaneous_errors: " << simultaneous_errors << std::endl;
                std::vector< bool > set_bits(n_weak_bits, false);
                test_pattern_generator_k_ones_in_n_bits_initialize(set_bits, n_weak_bits, simultaneous_errors);
                do
                {
                    // for(size_t i = 0; i < n_weak_bits; i++)
                    //     std::cout << set_bits[i];
                    // std::cout << std::endl;

                    std::unordered_set< uint64_t > precorrection_error_positions;
                    for(size_t i = 0; i < n_weak_bits; i++)
                        if(set_bits.at(i))
                            precorrection_error_positions.insert(weak_bit_indices.at(i));

                    // simulate the word and see which errors occur
                    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Ones(1, K);
                    Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p(1, K);
                    hc.calculate_dataword_p(precorrection_error_positions, dataword, dataword_p);
                    // std::unordered_set< uint64_t > postcorrection_error_positions;
                    for(uint64_t k = 0; k < K; k++)
                        if(dataword(k) != dataword_p(k))
                        {
                            // postcorrection_error_positions.insert(k);
                            per_bit_probabilities.at(k).at(simultaneous_errors) += 1;
                        }

                    // std::cout << "prec: " << set_to_string(precorrection_error_positions) << std::endl;
                    // std::cout << "posc: " << set_to_string(postcorrection_error_positions) << std::endl;
                } while(test_pattern_generator_k_ones_in_n_bits_advance(set_bits, n_weak_bits, simultaneous_errors));
            }

            // for(uint64_t k1 = 0; k1 < K; k1++)
            // {
            //     std::cout << "    " << k1 << " : {";
            //     for(uint64_t simultaneous_errors = 2; simultaneous_errors <= n_weak_bits; simultaneous_errors++)
            //         if(per_bit_probabilities.at(k1).at(simultaneous_errors) != 0)
            //             std::cout << " " << simultaneous_errors << ":" << per_bit_probabilities.at(k1).at(simultaneous_errors);
            //     std::cout << " }" << std::endl;
            // }

            double p = 0.5; // failure probability
            double omp = 1 - p;
            std::vector< double > probability_values;
            for(uint64_t k1 = 0; k1 < K; k1++)
            {
                double probability = 0.0;
                for(uint64_t simultaneous_errors = 2; simultaneous_errors <= n_weak_bits; simultaneous_errors++)
                {
                    int n_patterns = per_bit_probabilities.at(k1).at(simultaneous_errors);
                    probability += n_patterns * std::pow(p, simultaneous_errors) * std::pow(omp, n_weak_bits - simultaneous_errors);
                }
                // std::cout << "    " << k1 << " : " << probability << " {" << map_to_string(per_bit_probabilities.at(k1)) << "}" << std::endl;
                probability_values.push_back(probability);
            }

            // get the distribution of probability values
            std::cout << "[DATA] weak_bits: { " << vec_to_string(weak_bit_indices) << " } post_probabilities: [ " << vec_to_string(probability_values) << " ]" << std::endl;
        }
    }
}


/**
 * @brief prints out the expected post-correction error counts for a given number of pre-correction errors
 *
 * @param fname JSON configuration file for an ECC code
 * @param n_ecc_words number of ECC words to simulate
 * @param random_seed random seed to use for ECC word generation
 * @param conditional_pmf_samples probability distribution to return sample counts into
 */
void analyze_secondary_ecc(const std::string &fname, const int n_ecc_words, const int64_t random_seed
    , std::map< uint64_t /* n_prec_errors */, std::map< uint64_t /* n_posc_errors */, uint64_t /* sample count */ > > &conditional_pmf_samples)
{
    std::mt19937 rng(random_seed);

    hamming_code hc(fname);
    const uint64_t N = hc.get_n_code_bits();
    const uint64_t K = hc.get_n_data_bits();
    std::cout << "[INFO] evaluating JSON ECC code configuration file: " << fname << std::endl; 
    std::cout << "[INFO] using random_seed: " << random_seed << std::endl;
    
    std::uniform_int_distribution< uint64_t > udist(0, N - 1);
    for(size_t n_prec_errors = 0; n_prec_errors <= 10; n_prec_errors++)
    {
        if(conditional_pmf_samples.count(n_prec_errors) == 0)
            conditional_pmf_samples[n_prec_errors] = std::map< uint64_t, uint64_t >();

        for(int ecc_word = 0; ecc_word < n_ecc_words; ecc_word++)
        {
            // choose uniform-random error positions for this word based on seeded RNG
            std::unordered_set< uint64_t > precorrection_error_positions;
            while(precorrection_error_positions.size() != n_prec_errors)
                precorrection_error_positions.insert(udist(rng));
            
            // simulate the word and see which errors occur
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword = Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic >::Ones(1, K);
            Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > dataword_p(1, K);
            hc.calculate_dataword_p(precorrection_error_positions, dataword, dataword_p);
            
            // count errors
            uint64_t n_posc_errors = 0;
            for(uint64_t k = 0; k < K; k++)
                if(dataword(k) != dataword_p(k))
                    n_posc_errors++;
            
            if(conditional_pmf_samples.at(n_prec_errors).count(n_posc_errors) == 0)
                conditional_pmf_samples.at(n_prec_errors)[n_posc_errors] = 0;
            conditional_pmf_samples.at(n_prec_errors)[n_posc_errors]++;
        }
    }
}
