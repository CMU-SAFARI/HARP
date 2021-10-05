/**
 * @file hamming_code.h
 *
 * @brief Representation of an Hamming code with known parameters - NOT
 * compatible with Z3
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */
#ifndef HAMMING_CODE_H
#define HAMMING_CODE_H

/* stdlib includes */
#include <unordered_set>

/* lib includes */
#include "Eigen/Eigen"

/* project includes */
#include "gf2_noz3.h"
#include "einsim_interface.h"

extern int g_verbosity;

uint64_t nChoosek(uint64_t n, uint64_t k);

class hamming_code
{
private:
    uint64_t N;
    uint64_t K;
    uint64_t NmK;
    uint64_t uid;
    int64_t p; // 'permutation' or random seed
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > G;
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > H;
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > R;
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > GT;
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > HT;
    Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > RT;

public:
    /**
     * @brief constructor based on a JSON ECC configuration file
     * 
     * @param json_ecc_code_cfg_file JSON ECC configuration filename
     */
    hamming_code(const std::string &json_ecc_code_cfg_file)
    {
        extract_ecc_code_from_cfg_file(json_ecc_code_cfg_file);
    };

    /**
     * @brief constructor based on G/H/R matrices
     * 
     * @param G_in Hamming code G-matrix
     * @param H_in Hamming code H-matrix
     * @param R_in Hamming code R-matrix
     */
    hamming_code(const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G_in
               , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_in
               , const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R_in)
    {
        G.resize(G_in.rows(), G_in.cols());
        for(Eigen::Index r = 0; r < G_in.rows(); r++)
            for(Eigen::Index c = 0; c < G_in.cols(); c++)
                G(r, c) = G_in(r, c);

        H.resize(H_in.rows(), H_in.cols());
        for(Eigen::Index r = 0; r < H_in.rows(); r++)
            for(Eigen::Index c = 0; c < H_in.cols(); c++)
                H(r, c) = H_in(r, c);

        R.resize(R_in.rows(), R_in.cols());
        for(Eigen::Index r = 0; r < R_in.rows(); r++)
            for(Eigen::Index c = 0; c < R_in.cols(); c++)
                R(r, c) = R_in(r, c);

        GT = G.transpose();
        HT = H.transpose();
        RT = R.transpose();

        N = H.cols();
        K = G.rows();
        NmK = H.rows();
        assert(N - K == NmK);

        uid = compute_uid(G, H, R);
        p = -1; // unknown what random seed can generate this matrix
    };

    /**
     * @brief constructor of a random H-matrix based on locations of parity syndromes
     *
     *  this intitializer enumerates H-matrices as so:
     *  Forall H matrices in systematic form [P | I] s.t. P and I are strictly ordered
     * 
     * Given K, we have N. There are 2^NmK-NmK-1 syndromes -> (2^NmK-NmK-1)
     * 
     * @param k dataword length
     * @param non_parity_syndrome_mask desired locations of non-parity-check bits
     */ 
    hamming_code(const uint64_t k, const std::vector< bool > non_parity_syndrome_mask)
    {
        K = k;
        NmK = compute_hamming_code_n_parity_bits(k);
        N = K + NmK;

        assert(non_parity_syndrome_mask.size() == (1 << NmK) - NmK - 1 && "must specify a mask for each possible data syndrome");

        // compute syndromes to use (sorted)
        // we have {1, 2, ..., 2^NmK-1} \ {1, 2, 4, 8, 16, ...}
        std::vector< int > all_possible_syndromes;
        for(uint64_t i = 1; i < (1u << NmK); i++)
            if((i & (i - 1)) != 0)
                all_possible_syndromes.push_back(i);
        assert(all_possible_syndromes.size() == (1 << NmK) - NmK - 1);

        // trim syndromes corresponding to non_parity_syndromes
        std::vector< int > required_syndromes;
        for(uint64_t i = 0; i < (1 << NmK) - NmK - 1; i++)
            if(non_parity_syndrome_mask.at(i))
                required_syndromes.push_back(all_possible_syndromes.at(i));
        assert(required_syndromes.size() == K); // there should only be K syndromes remaining

        printf("included_s: { ");
        for(const int i: required_syndromes)
            printf("%d ", i);
        printf(" } excluded_s: { ");
        for(uint64_t i = 0; i < (1 << NmK) - NmK - 1; i++)
            if(!non_parity_syndrome_mask.at(i))
                printf("%d ", all_possible_syndromes.at(i));
        printf("}\n");
        
        // generate matrices
        Eigen::Matrix< gf2_noz3, Eigen::Dynamic, Eigen::Dynamic > P(K, NmK);
        for(uint64_t r = 0; r < K; r++)
            for(uint64_t c = 0; c < NmK; c++)
                P(r, c) = (required_syndromes.at(r) >> c) & 1;

        G.resize(K, N);
        G << Eigen::Matrix< gf2_noz3, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), P;
        R.resize(K, N);
        R << Eigen::Matrix< gf2_noz3, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), Eigen::Matrix< gf2_noz3, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(K, NmK);
        H.resize(NmK, N);
        H << P.transpose(), Eigen::Matrix< gf2_noz3, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(NmK, NmK);

        GT = G.transpose();
        RT = R.transpose();
        HT = H.transpose();

        // sanity-check that GH = 0
        assert((G * HT).cwiseEqual(0).all());

        uid = compute_uid(G, H, R);
        p = -1; // unknown what random seed can generate this matrix
        
        std::cout << "[INFO] Created Hamming parity matrices:" << std::endl;
        std::cout << "[INFO]     H: (" << H.rows() << ", "<< H.cols() << ")" << std::endl;
        std::cout << "[INFO]     G: (" << G.rows() << ", "<< G.cols() << ")" << std::endl;
        std::cout << "[INFO]     R: (" << R.rows() << ", "<< R.cols() << ")" << std::endl;
        // std::cout << "[INFO]     H:" << std::endl << H << std::endl;
        // std::cout << "[INFO]     G:" << std::endl << G << std::endl;
        // std::cout << "[INFO]     R:" << std::endl << R << std::endl;
    }
    
    // accessor functions
    uint64_t get_uid(void) const { return uid; }
    uint64_t get_n_code_bits(void) const { return N; }
    uint64_t get_n_data_bits(void) const { return K; }
    uint64_t get_n_parity_bits(void) const { return NmK; }
    uint64_t get_random_seed(void) const { return p; }

    void extract_ecc_code_from_cfg_file(const std::string &json_ecc_code_cfg_file);

    std::map< uint64_t, uint64_t > compute_exacerbation_factor(const int n_errors) const;
    void encode(
      const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword) const;
    void decode(
      Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &corrected_codeword
    , int &corrected_col
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p
    , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &syndrome) const;
    void calculate_everything(
          const std::unordered_set< uint64_t > &error_positions
        , const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &codeword_p
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &corrected_codeword
        , int &corrected_col
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &syndrome) const;
    void calculate_dataword_p(
          const std::unordered_set< uint64_t > &error_positions
        , const Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dataword_p) const;

    std::string to_json(void);

    /**
     * @brief utility function to compute the number of parity-check bits used
     * by a Hamming code with k parity-check bits
     *
     * @param k ECC dataword length
     * @return uint64_t number of ECC parity-check bits
     */
    static uint64_t compute_hamming_code_n_parity_bits(uint64_t k)
    {
        uint64_t np = 0;
        while((1u << np) < (np + k + 1u))
            np += 1;
        assert(np == std::ceil(std::log2(k + np + 1)) && "MISMATCH - INCORRECT LOGIC");
        return np;
    }
};

#endif /* HAMMING_CODE_H */