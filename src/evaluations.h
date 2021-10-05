#ifndef EVALUATIONS_H
#define EVALUATIONS_H

#include <string>
#include <random>

void evaluate(const std::string &fname, const int n_ecc_words, const int64_t random_seed);
void analyze_ecc_code(const std::string &fname, const int n_ecc_words, const int64_t random_seed);
void analyze_secondary_ecc(const std::string &fname, const int n_ecc_words, const int64_t random_seed
    , std::map< uint64_t /* n_prec_errors */, std::map< uint64_t /* n_posc_errors */, uint64_t /* sample count */ > > &conditional_pmf_samples);


#endif /* EVALUATIONS_H */