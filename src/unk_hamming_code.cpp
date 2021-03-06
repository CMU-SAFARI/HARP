/**
 * @file unk_hamming_code.cpp
 *
 * @brief [Adapted from the original BEER project] Embodiment of an unknown
 * hamming code and associated routines as used by Z3
 *
 * @author Minesh Patel (minesh.patelh@gmail.com)
 */

#include <cmath>
#include <chrono>
#include <unordered_set>

/* libraries */
#include "Eigen/Eigen"
#include "z3++.h"

/* project includes */
#include "gf2.h"
#include "beer_utils.h"
#include "unk_hamming_code.h"

/* forward declaratoins for global variables */
extern int g_verbosity;

/** z3 context pointer for the currently active context */
z3::context *gf2::ctx;

/** converting the enum error operator representation to a string */
std::map< enum error_operator, std::string > error_operator_to_string = 
{
      {EO_ALL_T,                   "ALL_T"}
    , {EO_ALL_A,                   "ALL_A"}
    , {EO_ALT_T,                   "ALT_T"}
    , {EO_DATA_ALT_T_PARITY_ALT_A, "DATA_ALT_T_PARITY_ALT_A"}
    , {EO_DATA_ALT_T_PARITY_T,     "DATA_ALT_T_PARITY_T"}
    , {EO_DATA_ALT_T_PARITY_A,     "DATA_ALT_T_PARITY_A"}
    , {EO_ALT_A,                   "ALT_A"}
    , {EO_DATA_ALT_A_PARITY_ALT_T, "DATA_ALT_A_PARITY_ALT_T"}
    , {EO_DATA_ALT_A_PARITY_T,     "DATA_ALT_A_PARITY_T"}
    , {EO_DATA_ALT_A_PARITY_A,     "DATA_ALT_A_PARITY_A"}

    // brute-force
    , {EO_DATA_ALT_T_PARITY_0T,    "DATA_ALT_T_PARITY_0T"}
    , {EO_DATA_ALT_T_PARITY_1T,    "DATA_ALT_T_PARITY_1T"}
    , {EO_DATA_ALT_T_PARITY_2T,    "DATA_ALT_T_PARITY_2T"}
    , {EO_DATA_ALT_T_PARITY_3T,    "DATA_ALT_T_PARITY_3T"}
    , {EO_DATA_ALT_T_PARITY_4T,    "DATA_ALT_T_PARITY_4T"}
    , {EO_DATA_ALT_T_PARITY_5T,    "DATA_ALT_T_PARITY_5T"}
    , {EO_DATA_ALT_T_PARITY_6T,    "DATA_ALT_T_PARITY_6T"}
    , {EO_DATA_ALT_T_PARITY_7T,    "DATA_ALT_T_PARITY_7T"}

    , {EO_DATA_ALT_A_PARITY_0T,    "DATA_ALT_A_PARITY_0T"}
    , {EO_DATA_ALT_A_PARITY_1T,    "DATA_ALT_A_PARITY_1T"}
    , {EO_DATA_ALT_A_PARITY_2T,    "DATA_ALT_A_PARITY_2T"}
    , {EO_DATA_ALT_A_PARITY_3T,    "DATA_ALT_A_PARITY_3T"}
    , {EO_DATA_ALT_A_PARITY_4T,    "DATA_ALT_A_PARITY_4T"}
    , {EO_DATA_ALT_A_PARITY_5T,    "DATA_ALT_A_PARITY_5T"}
    , {EO_DATA_ALT_A_PARITY_6T,    "DATA_ALT_A_PARITY_6T"}
    , {EO_DATA_ALT_A_PARITY_7T,    "DATA_ALT_A_PARITY_7T"}
};

/** converting the string error operator representation to an enum */
std::map< std::string, enum error_operator > string_to_error_operator = 
{
      {"ALL_T",                   EO_ALL_T}
    , {"ALL_A",                   EO_ALL_A}
    , {"ALT_T",                   EO_ALT_T}
    , {"DATA_ALT_T_PARITY_ALT_A", EO_DATA_ALT_T_PARITY_ALT_A}
    , {"DATA_ALT_T_PARITY_T",     EO_DATA_ALT_T_PARITY_T}
    , {"DATA_ALT_T_PARITY_A",     EO_DATA_ALT_T_PARITY_A}
    , {"ALT_A",                   EO_ALT_A}
    , {"DATA_ALT_A_PARITY_ALT_T", EO_DATA_ALT_A_PARITY_ALT_T}
    , {"DATA_ALT_A_PARITY_T",     EO_DATA_ALT_A_PARITY_T}
    , {"DATA_ALT_A_PARITY_A",     EO_DATA_ALT_A_PARITY_A}

    // brute-force
    , {"DATA_ALT_T_PARITY_0T",    EO_DATA_ALT_T_PARITY_0T}
    , {"DATA_ALT_T_PARITY_1T",    EO_DATA_ALT_T_PARITY_1T}
    , {"DATA_ALT_T_PARITY_2T",    EO_DATA_ALT_T_PARITY_2T}
    , {"DATA_ALT_T_PARITY_3T",    EO_DATA_ALT_T_PARITY_3T}
    , {"DATA_ALT_T_PARITY_4T",    EO_DATA_ALT_T_PARITY_4T}
    , {"DATA_ALT_T_PARITY_5T",    EO_DATA_ALT_T_PARITY_5T}
    , {"DATA_ALT_T_PARITY_6T",    EO_DATA_ALT_T_PARITY_6T}
    , {"DATA_ALT_T_PARITY_7T",    EO_DATA_ALT_T_PARITY_7T}

    , {"DATA_ALT_A_PARITY_0T",    EO_DATA_ALT_A_PARITY_0T}
    , {"DATA_ALT_A_PARITY_1T",    EO_DATA_ALT_A_PARITY_1T}
    , {"DATA_ALT_A_PARITY_2T",    EO_DATA_ALT_A_PARITY_2T}
    , {"DATA_ALT_A_PARITY_3T",    EO_DATA_ALT_A_PARITY_3T}
    , {"DATA_ALT_A_PARITY_4T",    EO_DATA_ALT_A_PARITY_4T}
    , {"DATA_ALT_A_PARITY_5T",    EO_DATA_ALT_A_PARITY_5T}
    , {"DATA_ALT_A_PARITY_6T",    EO_DATA_ALT_A_PARITY_6T}
    , {"DATA_ALT_A_PARITY_7T",    EO_DATA_ALT_A_PARITY_7T}
};

/**
 * @brief supporting C++-style printing of the GF2 type
 * 
 * @param os C++ output stream
 * @param m GF2 element to print
 * 
 * @return std::ostream& output stream with GF2 element appended
 */
std::ostream &operator<<(std::ostream &os, gf2 const &m) 
{ 
    return os << m.e;
}

/**
 * @brief evaluates an Eigen type containing Z3 types per a Z3 model
 * 
 * @param model z3 model to evaluate with
 * @param m Eigen type containing Z3 types
 * 
 * @return Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > evaluated matrix
 */
Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > z3_eval_Eigen_type(const z3::model &model, const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m)
{
    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > ret(m.rows(), m.cols());
    for(int r = 0; r < m.rows(); r++)
        for(int c = 0; c < m.cols(); c++)
            ret(r, c) = model.eval(m(r, c).get_z3_expr());
    return ret;
}

/**
 * @brief converts a GF2 type to a plain integer type
 * 
 * @param m GF2 type
 * 
 * @return Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > integer type
 */
Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > gf2_to_int(const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m)
{
    Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > ret(m.rows(), m.cols());
    for(int r = 0; r < m.rows(); r++)
        for(int c = 0; c < m.cols(); c++)
            ret(r, c) = m(r, c).get_z3_expr().is_true();
    return ret; 
}

/**
 * @brief converts a GF2 type to a plain bool type
 * 
 * @param m GF2 type
 * @return Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > bool type
 */
Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > gf2_to_bool(const Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > &m)
{
    Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > ret(m.rows(), m.cols());
    for(int r = 0; r < m.rows(); r++)
        for(int c = 0; c < m.cols(); c++)
            ret(r, c) = m(r, c).get_z3_expr().is_true();
    return ret; 
}

/**
 * @brief Swaps two rows of a matrix
 * 
 * @param M matrix
 * @param a first row index
 * @param b second row index
 */
void gf2_swap_rows(Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &M, int a, int b)
{
    M.row(a).swap(M.row(b));
}

/**
 * @brief Adds a matrix row to another row
 * 
 * @param M matrix
 * @param r0 source row
 * @param r1 destination row
 */
void gf2_add_r0_to_r1(Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &M, int r0, int r1)
{
    uint64_t C = M.cols();
    for(uint64_t c = 0; c < C; c++)
        M(r1, c) = M(r0, c) != M(r1, c); // GF(2) addition is XNOR
}

/**
 * @brief Computes the row-reduced echelon form of a given matrix
 * 
 * @param M matrix
 * @param left_rref whether to compute the left-aligned RREF (or right-aligned RREF)
 */
void gf2_rref(Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &M, bool left_rref)
{
    int R = M.rows();
    int C = M.cols();

    if(left_rref)
    {
        // Left-aligned RREF
        int lead = 0;
        for(int r = 0; r < R; r++)
        {
            if(lead >= C)
                return;
            
            // find the row with the leader
            int i = r;
            while(M(i, lead) == 0)
            {
                i = i + 1;
                if(i == R)
                {
                    i = r;
                    lead = lead + 1;
                    if(lead == C)
                        return;
                }
            }

            // swap rows i and r (i -> r so now r is the leader)
            gf2_swap_rows(M, i, r);
            assert(M(r, lead));

            // sum out all other rows with a leader 
            for(i = 0; i < R; i++)
                if(i != r && M(i, lead))
                    gf2_add_r0_to_r1(M, r, i);
            lead = lead + 1;
        }
    }
    else
    {
        // Right-aligned RREF
        int lead = C - R;
        for(int r = 0; r < R; r++)
        {
            // find the row with the leader
            if(lead >= C)
                return;
            
            int i = r;
            while(M(i, lead) == 0)
            {
                i = i + 1;
                if(i == R)
                {
                    i = r;
                    lead = lead + 1;
                    if(lead == C)
                        return;
                }
            }

            // swap rows i and r (i -> r so now r is the leader)
            assert(M(i, lead));
            gf2_swap_rows(M, i, r);
            assert(M(r, lead));
            
            // sum out all other rows with a leader 
            for(i = 0; i < R; i++)
                if(i != r && M(i, lead))
                    gf2_add_r0_to_r1(M, r, i);
            lead = lead + 1;
        }       
    }
}

/**
 * @brief Initializes the unknown hamming code from a JSON configuration file
 * 
 * @param json_ecc_code_cfg_file JSON configuration file
 */
void unk_hamming_code::extract_ecc_code_from_cfg_file(const std::string &json_ecc_code_cfg_file)
{
    rapidjson::Document d = read_json_cfg_file(json_ecc_code_cfg_file);
   
    this->K = d["k"].GetUint64();
    this->NmK = compute_hamming_code_n_parity_bits(K);
    this->N = NmK + K;

    GT.resize(N, K);
    const rapidjson::Value &G_json_obj = d["G"];
    for(rapidjson::SizeType r = 0; r < G_json_obj.Size(); r++)
        for(rapidjson::SizeType c = 0; c < G_json_obj[r].Size(); c++)
            GT(r, c) = (bool)G_json_obj[r][c].GetUint64(); 
    
    R.resize(K, N);
    const rapidjson::Value &R_json_obj = d["R"];
    for(rapidjson::SizeType r = 0; r < R_json_obj.Size(); r++)
        for(rapidjson::SizeType c = 0; c < R_json_obj[r].Size(); c++)
            R(r, c) = (bool)R_json_obj[r][c].GetUint64(); 
    
    H.resize(NmK, N);
    const rapidjson::Value &H_json_obj = d["H"];
    for(rapidjson::SizeType r = 0; r < H_json_obj.Size(); r++)
        for(rapidjson::SizeType c = 0; c < H_json_obj[r].Size(); c++)
            H(r, c) = (bool)H_json_obj[r][c].GetUint64(); 
    
    // set the transpose members - NO NEED to add constraints to make the objects equal - z3 expressions are pointers into the AST
    this->G = GT.transpose();
    this->RT = R.transpose();
    this->HT = H.transpose();
    
    std::cout << "[INFO] Created Hamming parity matrices:" << std::endl;
    std::cout << "[INFO]     H: (" << H.rows() << ", "<< H.cols() << ")" << std::endl;
    std::cout << "[INFO]     G: (" << G.rows() << ", "<< G.cols() << ")" << std::endl;
    std::cout << "[INFO]     R: (" << R.rows() << ", "<< R.cols() << ")" << std::endl;
}

/**
 * @brief Initialize an unknown Hamming code in systematic + standard form
 * 
 * Assumes that the code parameters (e.g., length) are already set
 */
void unk_hamming_code::build_ecc_code_systematic_and_standard_form(void)
{
    std::cout << "[INFO] building hamming_code of (n: " << N << " k: " << K << " n-k: " << NmK << " d: 3) in standard form (no total order on syndromes)" << std::endl;

    // General code (not systematic, not standard form)
    // self._G = [[z3.Bool("G_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._K)]
    // self._R = [[z3.Bool("R_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._K)]
    // self._H = [[z3.Bool("H_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._NmK)]

    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > P(K, NmK);
    for(uint64_t r = 0; r < K; r++)
        for(uint64_t c = 0; c < NmK; c++)
            P(r, c) = ctx.bool_const((std::string("P_") + std::to_string(r) + "_" + std::to_string(c)).c_str());

    G.resize(K, N);
    G << Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), P;
    R.resize(K, N);
    R << Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(K, NmK);
    H.resize(NmK, N);
    H << P.transpose(), Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(NmK, NmK);
    
    // set the transpose members - NO NEED to add constraints to make the objects equal - z3 expressions are pointers into the AST
    this->GT = G.transpose();
    this->RT = R.transpose();
    this->HT = H.transpose();

    // properties of a real Hamming code - unique columns
    for(uint64_t c_base = 0; c_base < N - 1; c_base++)
    {
        std::vector< gf2 > base_column_constraints;
        base_column_constraints.reserve(N);
        for(uint64_t c = c_base + 1; c < N; c++)
        {
            std::vector< gf2 > base_column_constraints_arr;
            base_column_constraints_arr.reserve(NmK);
            for(uint64_t r = 0; r < NmK; r++)
                base_column_constraints_arr.push_back(H(r, c) == H(r, c_base));
            base_column_constraints.push_back(!z3_and_reduce(base_column_constraints_arr));
        }
        sol.add(z3_and_reduce(base_column_constraints).get_z3_expr());
    }

    // properties of a real Hamming code - nonzero columns (i.e., nonzero rows in P)
    for(uint64_t r = 0; r < K; r++)
    {
        std::vector< gf2 > row_elements;
        row_elements.reserve(NmK);
        for(uint64_t c = 0; c < NmK; c++)
            row_elements.push_back(P(r, c));
        sol.add(z3_or_reduce(row_elements).get_z3_expr());
    }

    // GH = 0
    auto GHT = G * HT;
    for(uint64_t r = 0; r < K; r++)
        for(uint64_t c = 0; c < NmK; c++)
            sol.add(!GHT(r, c).get_z3_expr());

    std::cout << "[INFO] Created Hamming parity matrices:" << std::endl;
    std::cout << "[INFO]     H: (" << H.rows() << ", "<< H.cols() << ")" << std::endl;
    std::cout << "[INFO]     G: (" << G.rows() << ", "<< G.cols() << ")" << std::endl;
    std::cout << "[INFO]     R: (" << R.rows() << ", "<< R.cols() << ")" << std::endl;
    // std::cout << "[INFO]     H:" << std::endl << H << std::endl;
    // std::cout << "[INFO]     G:" << std::endl << G << std::endl;
    // std::cout << "[INFO]     R:" << std::endl << R << std::endl;

    // check if the constraints can be satisfied
    #if 0
    std::cout << "[INFO] Solution:" << std::endl;
    // std::cout << s.to_smt2() << "\n";
    switch (sol.check()) 
    {
        case z3::unsat:   std::cout << "UNSAT\n"; break;
        case z3::unknown: std::cout << "UNKNOWN\n"; break;
        case z3::sat:     std::cout << "SAT\n"; 
        {
            z3::model m = sol.get_model();

            // we can evaluate expressions in the model.
            std::cout << "[INFO]     H:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, H)) << std::endl;
            std::cout << "[INFO]     G:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, G)) << std::endl;
            std::cout << "[INFO]     R:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, R)) << std::endl;
            break;
        }
    }
    #endif
}

/**
 * @brief Initialize an unknown Hamming code in standard form with a total order on syndromes
 * 
 * Assumes that the code parameters (e.g., length) are already set
 */
void unk_hamming_code::build_ecc_code_standard_form_with_syndromes_ordered(void)
{
    std::cout << "[INFO] building hamming_code of (n: " << N << " k: " << K << " n-k: " << NmK << " d: 3) in standard from with a total order on all syndromes" << std::endl;

    // General code (not systematic, not standard form)
    // self._G = [[z3.Bool("G_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._K)]
    // self._R = [[z3.Bool("R_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._K)]
    // self._H = [[z3.Bool("H_{0}_{1}".format(j, i)) for i in range(self._N)] for j in range(self._NmK)]

    Eigen::Matrix< gf2, Eigen::Dynamic, Eigen::Dynamic > P(K, NmK);
    for(uint64_t r = 0; r < K; r++)
        for(uint64_t c = 0; c < NmK; c++)
            P(r, c) = ctx.bool_const((std::string("P_") + std::to_string(r) + "_" + std::to_string(c)).c_str());

    G.resize(K, N);
    G << Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), P;
    R.resize(K, N);
    R << Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(K, K), Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Zero(K, NmK);
    H.resize(NmK, N);
    H << P.transpose(), Eigen::Matrix< gf2, Eigen::Dynamic /* rows */, Eigen::Dynamic /* cols */>::Identity(NmK, NmK);
    
    // set the transpose members - NO NEED to add constraints to make the objects equal - z3 expressions are pointers into the AST
    this->GT = G.transpose();
    this->RT = R.transpose();
    this->HT = H.transpose();

    // WLOG force total order of syndromes to simplify SAT problem -> forces a unique solution for a given set of syndromes
    // first K rows of H must be the data rows per standard form
    // // z3::expr_vector syndromes_as_ints(ctx);
    // // for(uint64_t r = 0; r < K; r++)
    // // {
    // //     z3::expr_vector syn_as_int(ctx);
    // //     for(uint64_t c = 0; c < NmK; c++)
    // //     {
    // //         z3::expr intval = z3::ite(HT(r, c).get_z3_expr(), ctx.int_val(1 << c), ctx.int_val(0));
    // //         syn_as_int.push_back(intval);
    // //     }
    // //     syndromes_as_ints.push_back(z3::sum(syn_as_int));
    // // }
    // // for(uint64_t ra = 0; ra < K - 1; ra++)
    // // {
    // //     uint64_t rb = ra + 1;
    // //     sol.add(syndromes_as_ints[ra] < syndromes_as_ints[rb], (std::string("total order on syndromes ") + std::to_string(ra)).c_str());
    // // }

    // This is a fancy bitwise way of saying that pairwise syndromes are totally ordered -> a < b < c < d...
    // 
    // this is SO MUCH FASTER for large codes (e.g., k=120 -> UNSAT time goes from 5s to 100ms)
    // 
    // a' || b || (a' && b) in higher bit1 || (a' && b) in higher bit2...
    // &&
    // (a == b) || (01 in higher bit1) || (01 in higher bit2)...
    // && ...
    for(uint64_t a = 0; a < K - 1; a++)
    {
        uint64_t b = a + 1;

        z3::expr_vector a_lt_b(ctx);
        z3::expr_vector a_ne_b(ctx);
        for(uint64_t bit = 0; bit < NmK; bit++)
        {
            z3::expr cond = (!P(a, bit).get_z3_expr()) || P(b, bit).get_z3_expr();
            for(uint64_t higher_bit = bit + 1; higher_bit < NmK; higher_bit++)
                cond = cond || ((!P(a, higher_bit).get_z3_expr()) && P(b, higher_bit).get_z3_expr());
            a_lt_b.push_back(cond);
            a_ne_b.push_back(P(a, bit).get_z3_expr() != P(b, bit).get_z3_expr());
        }
        sol.add(z3::mk_and(a_lt_b));
        sol.add(z3::mk_or(a_ne_b));
    }

    

    // properties of a real Hamming code - nonzero columns (i.e., nonzero rows in P)
    // no P column can be an identity column - we encode this by saying that the hamming weight of all cols is at least 1
    std::vector< int > syndrome_coeffs(NmK, 1);
    for(uint64_t r = 0; r < K; r++)
    {
        z3::expr_vector this_syndrome(ctx);
        for(uint64_t c = 0; c < NmK; c++)
            this_syndrome.push_back(P(r, c).get_z3_expr());
        sol.add(z3::pbge(this_syndrome, syndrome_coeffs.data(), 2));
    }

    // GH = 0
    auto GHT = G * HT;
    z3::expr_vector ght(ctx);
    for(uint64_t r = 0; r < K; r++)
        for(uint64_t c = 0; c < NmK; c++)
            ght.push_back(GHT(r, c).get_z3_expr());
    sol.add(!z3::mk_or(ght));
    
    std::cout << "[INFO] Created Hamming parity matrices:" << std::endl;
    std::cout << "[INFO]     H: (" << H.rows() << ", "<< H.cols() << ")" << std::endl;
    std::cout << "[INFO]     G: (" << G.rows() << ", "<< G.cols() << ")" << std::endl;
    std::cout << "[INFO]     R: (" << R.rows() << ", "<< R.cols() << ")" << std::endl;
    // std::cout << "[INFO]     H:" << std::endl << H << std::endl;
    // std::cout << "[INFO]     G:" << std::endl << G << std::endl;
    // std::cout << "[INFO]     R:" << std::endl << R << std::endl;

    // check if the constraints can be satisfied
    #if 0
    std::cout << "[INFO] Solution:" << std::endl;
    // std::cout << s.to_smt2() << "\n";
    switch (sol.check()) 
    {
        case z3::unsat:   std::cout << "UNSAT\n"; break;
        case z3::unknown: std::cout << "UNKNOWN\n"; break;
        case z3::sat:     std::cout << "SAT\n"; 
        {
            z3::model m = sol.get_model();

            // we can evaluate expressions in the model.
            std::cout << "[INFO]     H:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, H)) << std::endl;
            std::cout << "[INFO]     G:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, G)) << std::endl;
            std::cout << "[INFO]     R:" << std::endl << gf2_to_int(z3_eval_Eigen_type(m, R)) << std::endl;
            break;
        }
    }
    #endif
}

/**
 * @brief 'Injects' errors according to a particular error model
 *
 * Careful- this is the interaction between CHARGED/DISCHARGED states and the
 * actual values stored in the cells. The "error operator" is the realization of
 * what needs to happen according to the different possiblities
 *
 * @param codeword_p post-error codeword
 * @param codeword pre-correction codeword
 * @param error_mask mask for which bits can experience errors (1) or cannot (0)
 * @param eo the error model to use when injecting errors
 */
void unk_hamming_code::apply_error_operator_to_codeword(    
      Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword_p
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , const enum error_operator eo)
{
    codeword_p.resize(1, N);
    for(uint64_t i = 0; i < N; i++)
    {
        if(eo == EO_ALL_T)
            codeword_p(i) = codeword(i) && error_mask(i);
        else if(eo == EO_ALL_A)
            codeword_p(i) = codeword(i) || error_mask(i);
        else if(eo == EO_ALT_T
             || eo == EO_DATA_ALT_T_PARITY_ALT_A 
             || eo == EO_DATA_ALT_T_PARITY_T
             || eo == EO_DATA_ALT_T_PARITY_A
             || eo == EO_DATA_ALT_T_PARITY_0T
             || eo == EO_DATA_ALT_T_PARITY_1T
             || eo == EO_DATA_ALT_T_PARITY_2T
             || eo == EO_DATA_ALT_T_PARITY_3T
             || eo == EO_DATA_ALT_T_PARITY_4T
             || eo == EO_DATA_ALT_T_PARITY_5T
             || eo == EO_DATA_ALT_T_PARITY_6T
             || eo == EO_DATA_ALT_T_PARITY_7T)
        {
            if((i < K) || (eo == EO_ALT_T))
            {
                if((i % 2) == 0)
                    codeword_p(i) = codeword(i) && error_mask(i);
                else
                    codeword_p(i) = codeword(i) || error_mask(i);
            }
            else if(eo == EO_DATA_ALT_T_PARITY_ALT_A)
            {
                if((i % 2) == 0)
                    codeword_p(i) = codeword(i) || error_mask(i);
                else
                    codeword_p(i) = codeword(i) && error_mask(i);
            }
            else if(eo == EO_DATA_ALT_T_PARITY_T)
                    codeword_p(i) = codeword(i) && error_mask(i);
            else if(eo == EO_DATA_ALT_T_PARITY_A)
                    codeword_p(i) = codeword(i) || error_mask(i);
            else
            { 
                switch(eo)
                {
                    case EO_DATA_ALT_T_PARITY_0T: codeword_p(i) = (i - K <= 0) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_1T: codeword_p(i) = (i - K <= 1) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_2T: codeword_p(i) = (i - K <= 2) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_3T: codeword_p(i) = (i - K <= 3) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_4T: codeword_p(i) = (i - K <= 4) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_5T: codeword_p(i) = (i - K <= 5) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_6T: codeword_p(i) = (i - K <= 6) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    case EO_DATA_ALT_T_PARITY_7T: codeword_p(i) = (i - K <= 7) ? (codeword(i) && error_mask(i)) : (codeword(i) || error_mask(i)); break;
                    default:
                        assert(0 && "impossible");
                }                 
            }
        }
        else if(eo == EO_ALT_A
             || eo == EO_DATA_ALT_A_PARITY_ALT_T 
             || eo == EO_DATA_ALT_A_PARITY_T
             || eo == EO_DATA_ALT_A_PARITY_A
             || eo == EO_DATA_ALT_A_PARITY_0T
             || eo == EO_DATA_ALT_A_PARITY_1T
             || eo == EO_DATA_ALT_A_PARITY_2T
             || eo == EO_DATA_ALT_A_PARITY_3T
             || eo == EO_DATA_ALT_A_PARITY_4T
             || eo == EO_DATA_ALT_A_PARITY_5T
             || eo == EO_DATA_ALT_A_PARITY_6T
             || eo == EO_DATA_ALT_A_PARITY_7T)
        {
            if((i < K) || (eo == EO_ALT_A))
            {
                if((i % 2) == 0)
                    codeword_p(i) = codeword(i) || error_mask(i);
                else
                    codeword_p(i) = codeword(i) && error_mask(i);
            }
            else if(eo == EO_DATA_ALT_A_PARITY_ALT_T)
            {
                if((i % 2) == 0)
                    codeword_p(i) = codeword(i) && error_mask(i);
                else
                    codeword_p(i) = codeword(i) || error_mask(i);
            }
            else if(eo == EO_DATA_ALT_A_PARITY_T)
                    codeword_p(i) = codeword(i) && error_mask(i);
            else if(eo == EO_DATA_ALT_A_PARITY_A)
                    codeword_p(i) = codeword(i) || error_mask(i);
            else
            { 
                switch(eo)
                {
                    case EO_DATA_ALT_A_PARITY_0T: codeword_p(i) = (i - K <= 0) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_1T: codeword_p(i) = (i - K <= 1) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_2T: codeword_p(i) = (i - K <= 2) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_3T: codeword_p(i) = (i - K <= 3) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_4T: codeword_p(i) = (i - K <= 4) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_5T: codeword_p(i) = (i - K <= 5) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_6T: codeword_p(i) = (i - K <= 6) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    case EO_DATA_ALT_A_PARITY_7T: codeword_p(i) = (i - K <= 7) ? (codeword(i) || error_mask(i)) : (codeword(i) && error_mask(i)); break;
                    default:
                        assert(0 && "impossible");
                }                 
            }
        }
        else
        {
            printf("[ERROR] unsupported error operator: %s\n", error_operator_to_string.at(eo).c_str());
            assert(0 && "unsupported eo");
        }
    }    
}

// Eigen::Matrix< gf2, 1, Eigen::Dynamic > unk_hamming_code::calculate_syndrome_from_error_pattern(
//         const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_pattern)
// {

// }


/**
 * @brief Calculates the error syndrome given a single codeword
 * 
 * @param codeword ECC codeword
 * @param error_mask mask for which bits can experience errors (1) or cannot (0)
 * @param eo the error model to use when injecting errors
 * 
 * @return Eigen::Matrix< gf2, 1, Eigen::Dynamic > error syndrome
 */
Eigen::Matrix< gf2, 1, Eigen::Dynamic > unk_hamming_code::calculate_syndrome_from_codeword(
      const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , const enum error_operator eo
)
{
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword_p(1, N);
    apply_error_operator_to_codeword(codeword_p, codeword, error_mask, eo);
    return codeword_p * HT;
}

/**
 * @brief Calculates the error syndrome given an single dataword
 * 
 * @param dataword ECC dataword
 * @param error_mask mask for which bits can experience errors (1) or cannot (0)
 * @param eo the error model to use when injecting errors
 * 
 * @return Eigen::Matrix< gf2, 1, Eigen::Dynamic > error syndrome
 */
Eigen::Matrix< gf2, 1, Eigen::Dynamic > unk_hamming_code::calculate_syndrome_from_dataword(
      const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , const enum error_operator eo
)
{
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword = dataword * G;
    return calculate_syndrome_from_codeword(codeword, error_mask, eo);
}

/**
 * @brief caluculates the full corrected dataword_p given an implied ZERO
 * codeword and an error mask
 *
 * @param error_mask mask for which bits can experience errors (1) or cannot (0)
 * @param corrected_codeword ECC codeword (post-correction)
 * @param dataword_p ECC dataword (post-correction)
 * @param syndrome genreated error syndrome
 */
void unk_hamming_code::calculate_corrected_dataword_from_error_mask_given_zero_codeword(
      const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &corrected_codeword
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &dataword_p
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &syndrome
)
{
    syndrome = error_mask * HT;

    corrected_codeword.resize(1, N);
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        z3::expr_vector is_syn_pos_eq_arr(ctx);
        for(uint64_t ridx = 0; ridx < NmK; ridx++)
            is_syn_pos_eq_arr.push_back((H(ridx, bidx) == syndrome(ridx)).get_z3_expr());
        corrected_codeword(bidx) = error_mask(bidx).get_z3_expr() != z3::mk_and(is_syn_pos_eq_arr);
    }
    dataword_p = corrected_codeword * RT;
}

/**
 * @brief caluculates the full corrected dataword_p given a fixed (i.e., bool) dataword
 * 
 * @param dataword ECC dataword
 * @param error_mask mask for which bits can experience errors (1) or cannot (0)
 * @param eo the error model to use when injecting errors
 * @param codeword ECC codeword (pre-correction)
 * @param codeword_p ECC codeword (post-error)
 * @param corrected_codeword ECC codeword (post-correction)
 * @param dataword_p ECC dataword (post-correction)
 * @param syndrome genreated error syndrome
 */
void unk_hamming_code::calculate_everything(
      const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , const enum error_operator eo
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &codeword_p
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &corrected_codeword
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &dataword_p
    , Eigen::Matrix< gf2, 1, Eigen::Dynamic > &syndrome
)
{
    codeword = dataword * G;
    codeword_p.resize(1, N);
    apply_error_operator_to_codeword(codeword_p, codeword, error_mask, eo);
    syndrome = codeword_p * HT;

    corrected_codeword.resize(1, N);
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        z3::expr_vector is_syn_pos_eq_arr(ctx);
        for(uint64_t ridx = 0; ridx < NmK; ridx++)
            is_syn_pos_eq_arr.push_back((H(ridx, bidx) == syndrome(ridx)).get_z3_expr());
        corrected_codeword(bidx) = codeword_p(bidx) != z3::mk_and(is_syn_pos_eq_arr);
    }
    dataword_p = corrected_codeword * RT;
}

/**
 * @brief construct a constraint that represents the fact that a particular
 * error mask can cause a given error syndrome
 *
 * @param dataword initial ECC dataword
 * @param error_mask error mask to apply to the codeword
 * @param syndrome_idx index of the error syndrome within the codeword
 * @param eo the error model to use when injecting errors
 * 
 * @return z3::expr the constraint
 */
z3::expr unk_hamming_code::enforce_error_mask_causes_syndrome_idx(
      const Eigen::Matrix< bool, 1, Eigen::Dynamic > &dataword
    , const Eigen::Matrix< gf2, 1, Eigen::Dynamic > &error_mask
    , uint64_t syndrome_idx
    , const enum error_operator eo
)
{
    const Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome_exp = H.col(syndrome_idx);
    const Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome = calculate_syndrome_from_dataword(dataword, error_mask, eo);
    assert((uint64_t)syndrome_exp.size() == NmK);
    assert((uint64_t)syndrome.size() == NmK);

    z3::expr_vector syndrome_entries(ctx);
    for(uint64_t i = 0; i < NmK; i++)
        syndrome_entries.push_back((syndrome(i) == syndrome_exp(i)).get_z3_expr());
    return z3::mk_and(syndrome_entries);
}

/**
 * @brief computes the miscorrection profile given a fully specified ECC code
 * 
 * @param mc_profile_test_patterns List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param mc_profile_miscorrections_binarized Miscorrections observed at each bit position for each test pattern
 * @param eo the error model to use when injecting errors
 */
void unk_hamming_code::compute_miscorrection_profile(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
    , const enum error_operator eo)
{
    mc_profile_miscorrections_binarized.resize(mc_profile_test_patterns.rows(), mc_profile_test_patterns.cols());
    for(int tp_idx = 0; tp_idx < mc_profile_test_patterns.rows(); tp_idx++)
    {
        const Eigen::Matrix< bool, 1, Eigen::Dynamic, Eigen::RowMajor > &test_pattern = mc_profile_test_patterns.row(tp_idx);
        for(uint64_t data_bit_idx = 0; data_bit_idx < K; data_bit_idx++)
        {
            solver_push();

            // set up the error mask
            Eigen::Matrix< gf2, 1, Eigen::Dynamic > error_mask(1, N);
            z3::expr_vector error_mask_vars(ctx);
            for(uint64_t i = 0; i < N; i++)
            {
                error_mask(i) = ctx.bool_const(("err_mask_" + std::to_string(tp_idx) + "_" + std::to_string(data_bit_idx) + "_" + std::to_string(i)).c_str());
                error_mask_vars.push_back(error_mask(i).get_z3_expr());
            }
            
            // define an existence constraint for this particular syndrome
            z3::expr em_exists(ctx);
            if(true)
            {
                // THIS ASSUMES ONLY SYSTEMATIC ENCODING -> unknown syndrome idx
                // note that we DO NOT know the syndrome since we don't know which codeword bit maps to this data bit idx
                // the constraint is still possible but relies on the R matrix
                Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword;
                Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword_p;
                Eigen::Matrix< gf2, 1, Eigen::Dynamic > corrected_codeword;
                Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword_p;
                Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome;
                calculate_everything(test_pattern, error_mask, eo, codeword, codeword_p, corrected_codeword, dataword_p, syndrome);
                em_exists = z3::exists(error_mask_vars, (dataword_p(data_bit_idx) != test_pattern(data_bit_idx)).get_z3_expr());
                
                // constrain the error mask to avoid unnecessary state-space exploration
                for(uint64_t i = 0; i < N; i++)
                    if(!codeword(i).get_z3_expr().is_true())
                        sol.add(!error_mask(i).get_z3_expr());
            }
            else
            {
                // THIS ASSUMES SYSTEMATIC + STANDARD FORM -> data_bit_idx == syndrome_idx -> 
                // we do NOT need to know the syndrome, just its index -> this is significantly faster for the SAT solver
                // this does NOT require he R matrix
                em_exists = z3::exists(error_mask_vars, enforce_error_mask_causes_syndrome_idx(test_pattern, error_mask, data_bit_idx, eo));
            }

            // assume this is possible
            sol.add(em_exists);
            
            // solve the SAT problem
            // std::cout << "[INFO]     checking satisfiability of test pattern: " << test_pattern << " bit: " << data_bit_idx << std::endl;
            // auto start_time = std::chrono::high_resolution_clock::now();
             auto sat_result = sol.check();
            // auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
            // std::cout << "[TIME]     " << time_elapsed / std::chrono::milliseconds(1) << " ms to determine mc profile bit using SAT solver" << std::endl;

            // if possible, set miscorrections to true
            switch(sat_result) 
            {
                case z3::sat:
                    mc_profile_miscorrections_binarized(tp_idx, data_bit_idx) = true;
                    // std::cout << "MODEL " << (bool)sol.get_model().eval(em_exists) << std::endl;
                    break;
                
                case z3::unsat:
                    mc_profile_miscorrections_binarized(tp_idx, data_bit_idx) = false;
                    break;

                case z3::unknown:
                default:
                    assert(0 && "wat");
                    break;
            }
    
            // pop solver state
            solver_pop();
        }
    }
}

/**
 * @brief Get the DISCHARGED value of a particular bit position in the dataword
 * 
 * Maps a particular 'error operator' to DISCHARGED values
 * 
 * @param eo the error model to use when injecting errors
 * @param bit_idx bit index in the dataword
 *
 * @return true DISCHARGED value is 1
 * @return false DISCHARGED value is 0
 */
bool get_discharged_value_dataword(const enum error_operator eo, uint64_t bit_idx)
{
    if(eo == EO_ALL_T)  // only 1->0 errors allowed
        return false;
    else if(eo == EO_ALL_A)
        return true;
    else if(eo == EO_ALT_T
        || eo == EO_DATA_ALT_T_PARITY_ALT_A
        || eo == EO_DATA_ALT_T_PARITY_T
        || eo == EO_DATA_ALT_T_PARITY_A
        || eo == EO_DATA_ALT_T_PARITY_0T
        || eo == EO_DATA_ALT_T_PARITY_1T
        || eo == EO_DATA_ALT_T_PARITY_2T
        || eo == EO_DATA_ALT_T_PARITY_3T
        || eo == EO_DATA_ALT_T_PARITY_4T
        || eo == EO_DATA_ALT_T_PARITY_5T
        || eo == EO_DATA_ALT_T_PARITY_6T
        || eo == EO_DATA_ALT_T_PARITY_7T)
        return (bit_idx % 2 == 0) ? false : true;
    else if(eo == EO_ALT_A
        || eo == EO_DATA_ALT_A_PARITY_ALT_T
        || eo == EO_DATA_ALT_A_PARITY_T
        || eo == EO_DATA_ALT_A_PARITY_A
        || eo == EO_DATA_ALT_A_PARITY_0T
        || eo == EO_DATA_ALT_A_PARITY_1T
        || eo == EO_DATA_ALT_A_PARITY_2T
        || eo == EO_DATA_ALT_A_PARITY_3T
        || eo == EO_DATA_ALT_A_PARITY_4T
        || eo == EO_DATA_ALT_A_PARITY_5T
        || eo == EO_DATA_ALT_A_PARITY_6T
        || eo == EO_DATA_ALT_A_PARITY_7T)
        return (bit_idx % 2 == 0) ? true : false;
    else
        assert(0 && "unsupported error operator!");
}

/**
 * @brief Get the CHARGED value of a particular bit position in the dataword
 * 
 * Maps a particular 'error operator' to CHARGED values
 * 
 * @param eo the error model to use when injecting errors
 * @param bit_idx bit index in the dataword
 *
 * @return true CHARGED value is 1
 * @return false CHARGED value is 0
 */
bool get_charged_value_dataword(const enum error_operator eo, uint64_t bit_idx)
{
    return !get_discharged_value_dataword(eo, bit_idx);
}

/**
 * @brief Worker for determining the ECC function from a binarized miscorrection profile
 * 
 * @param mc_profile_test_patterns List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param mc_profile_miscorrections_binarized Miscorrections observed at each bit position for each test pattern
 * @param constrain_mc_exists whether to apply the constraint that a miscorrection is possible
 * @param constrain_mc_not_exists whether to apply the constraint that a miscorrection is impossible
 * @param eo the error model to use when injecting errors
 */
void unk_hamming_code::determine_ecc_function_from_binarized_miscorrection_profile_impl(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
    , const bool constrain_mc_exists
    , const bool constrain_mc_not_exists
    , const enum error_operator eo
)
{
    // force systematic encoding in standard form so that we reduce the search space
    std::cout << "[INFO] establishing systematic encoding constraints using eo: " << eo << std::endl;
    for(uint64_t row_idx = 0; row_idx < K; row_idx++)
        for(uint64_t data_bit_idx = 0; data_bit_idx < K; data_bit_idx++)
        {
            if(row_idx == data_bit_idx)
            {
                sol.add(GT(row_idx, data_bit_idx).get_z3_expr());
                sol.add(R(row_idx, data_bit_idx).get_z3_expr());
            }
            else
            {
                sol.add(!GT(row_idx, data_bit_idx).get_z3_expr());
                sol.add(!R(row_idx, data_bit_idx).get_z3_expr());
            }
        }
    for(uint64_t row_idx = 0; row_idx < K; row_idx++)
        for(uint64_t data_bit_idx = K; data_bit_idx < N; data_bit_idx++)
            sol.add(!R(row_idx, data_bit_idx).get_z3_expr());
    for(uint64_t col_idx = 0; col_idx < NmK; col_idx++)
        for(uint64_t parity_bit_idx = 0; parity_bit_idx < NmK; parity_bit_idx++)
        {
            if(col_idx == parity_bit_idx)
                sol.add(H(parity_bit_idx, col_idx + K).get_z3_expr());
            else
                sol.add(!H(parity_bit_idx, col_idx + K).get_z3_expr());
        }

    // establish the miscorrection profile constraints
    std::cout << "[INFO] establishing miscorrection profile constraints" << std::endl;
    assert(mc_profile_test_patterns.rows() == mc_profile_miscorrections_binarized.rows());
    assert(mc_profile_test_patterns.cols() == mc_profile_miscorrections_binarized.cols());
    for(int tp_idx = 0; tp_idx < mc_profile_test_patterns.rows(); tp_idx++)
    {
        const Eigen::Matrix< bool, 1, Eigen::Dynamic > &n_charged_test_pattern = mc_profile_test_patterns.row(tp_idx);
        const Eigen::Matrix< long, 1, Eigen::Dynamic > &miscorrections = mc_profile_miscorrections_binarized.row(tp_idx);
        
        // construct the dataword based on the one-charged test pattern in the MCP and the error model
        Eigen::Matrix< bool, 1, Eigen::Dynamic > dataword(1, K);
        for(uint64_t data_bit_idx = 0; data_bit_idx < K; data_bit_idx++)
            if(n_charged_test_pattern(data_bit_idx))
                dataword(data_bit_idx) = get_charged_value_dataword(eo, data_bit_idx);
            else
                dataword(data_bit_idx) = get_discharged_value_dataword(eo, data_bit_idx);

        // iterate over the test pattern
        for(uint64_t data_bit_idx = 0; data_bit_idx < K; data_bit_idx++)
        {
            // do not constrain any CHARGED dataword bits - they cannot present miscorrections (at least, we cannot know
            // if they are true retention errors or miscrorrections)
            if(n_charged_test_pattern(data_bit_idx))
                continue;

            // non-binary is the signal for uncertainty, so provide no constraint
            if(miscorrections(data_bit_idx) != !!miscorrections(data_bit_idx))
                continue; 

            // set up the error mask
            Eigen::Matrix< gf2, 1, Eigen::Dynamic > error_mask(1, N);
            z3::expr_vector error_mask_vars(ctx);
            for(uint64_t i = 0; i < N; i++)
            {
                if(i >= K or n_charged_test_pattern(i)) // the corresponding data bit will be charged and can fail
                {
                    error_mask(i) = ctx.bool_const(("err_mask_" + std::to_string(tp_idx) + "_" + std::to_string(data_bit_idx) + "_" + std::to_string(i)).c_str());
                    error_mask_vars.push_back(error_mask(i).get_z3_expr());
                }
                else
                    error_mask(i) = dataword(i); // this bit is discharged so canNOT fail - program it to the quiescent state and don't solve for it
            }

            // define an existence constraint for this particular syndrome
            z3::expr em_exists(ctx);
            if(false)
            {
                // THIS ASSUMES ONLY SYSTEMATIC ENCODING -> unknown syndrome idx
                // note that we DO NOT know the syndrome since we don't know which codeword bit maps to this data bit idx
                // the constraint is still possible but relies on the R matrix
                assert(0 && "UNIMPLEMENTED!");
                // dataword_extended, codeword, codeword_p, corrected_codeword, dataword_p, syndrome = self.calculate_everything(dataword, error_mask, eo)
                // em_exists = z3.Exists(error_mask_vars, dataword[data_bit_idx] != dataword_p[data_bit_idx])
            }
            else
            {
                // THIS ASSUMES SYSTEMATIC + STANDARD FORM -> data_bit_idx == syndrome_idx -> 
                // we do NOT need to know the syndrome, just its index -> this is significantly faster for the SAT solver
                // this does NOT require he R matrix
                em_exists = z3::exists(error_mask_vars, enforce_error_mask_causes_syndrome_idx(dataword, error_mask, data_bit_idx, eo));
            }

            // add possible/impossible constraint to Z3
            if(miscorrections(data_bit_idx) == 1)
            {
                // CONSTRAINT: there must exist an error mask that sets the syndrome to the identified bit
                if(constrain_mc_exists)
                    sol.add(em_exists); // ONES
            }
            else if(miscorrections(data_bit_idx) == 0)
            {
                // EXCLUSION: no error mask must be able to set the syndrome to the identified bit
                if(constrain_mc_not_exists)
                    sol.add(!em_exists); // ZEROES
            }
            else
                assert(0 && "impossible");
        }
    }
}

/**
 * @brief Determines the ECC function for a given binarized miscorrection profile
 * 
 * @param mc_profile_test_patterns List of N-charged test patterns showing which bits are charged (1) and discharged (0)
 * @param mc_profile_miscorrections_binarized Miscorrections observed at each bit position for each test pattern
 * @param G_ret solved-for generator matrix
 * @param H_ret solved-for parity-check matrix
 * @param R_ret solved-for degenerator matrix
 * @param Hs list of parityc-check matrices to *exclude* from this solve
 * @param is_rerun whether this is an incremental solve building upon a previous invocation
 * @param eo the error model to use when injecting errors
 * @param constrain_mc_exists whether to apply the constraint that a miscorrection is possible
 * @param constrain_mc_not_exists whether to apply the constraint that a miscorrection is impossible
 * 
 * @return true 
 * @return false 
 */
bool unk_hamming_code::determine_ecc_function_from_binarized_miscorrection_profile(
      const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_test_patterns
    , const Eigen::Matrix< long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > &mc_profile_miscorrections_binarized
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G_ret
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_ret
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R_ret
    , std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Hs
    , const bool is_rerun
    , const bool constrain_mc_exists
    , const bool constrain_mc_not_exists
    , const enum error_operator eo
    )
{
    if(is_rerun)
    {
        if(g_verbosity >= 1)
            std::cout << "[INFO] re-using previous constraints" << std::endl;
        solver_pop();
    }
    else
    {
        std::cout << "[INFO] no previous constraints to re-use: establishing constraints" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        determine_ecc_function_from_binarized_miscorrection_profile_impl(mc_profile_test_patterns, mc_profile_miscorrections_binarized, constrain_mc_exists, constrain_mc_not_exists, eo);
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to establish constraints" << std::endl;
    }

    // avoid finding the same matrices as previously found
    // the SET of data columns must NOT be equal to the same SET
    if(g_verbosity >= 1)
        std::cout << "[INFO] excluding already-discovered schemes" << std::endl;
    // assert(Hs.size() <= 1 && "incremental solver should only be excluding the last-found solution!");
    for(const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_exclude : Hs)
    {
        std::vector< gf2 > equalities;
        for(uint64_t row = 0; row < NmK; row++)
            for(uint64_t col = 0; col < N; col++)
                equalities.push_back(H(row, col) == H_exclude(row, col));
        sol.add(!z3_and_reduce(equalities).get_z3_expr());
    }
    
    // ensure that the SET of H columns is not the same as what has already been found
    // for H in exclude_Hs:
    //     or_expressions = []
    //     for c_base in range(self._K):
    //         coleq_expressions = []
    //         for c in range(self._K):
    //             coleq_expressions.append(z3.And([self._H[r][c] == z3.BoolVal(H[r][c_base]) for r in range(self._NmK)]))
    //         or_expressions.append(z3.Or(coleq_expressions))
    //     self._s.add(z3.Not(z3.And(or_expressions)))

    // save miscorrection profile assertions for re-use when trying to find additional solutions
    if(g_verbosity >= 1)
        std::cout << "[INFO] saving solver state for next run" << std::endl;
    solver_push();

    // check SAT
    std::cout << "[INFO] checking satisfiability of determine_ecc_function_from_binarized_miscorrection_profile" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto sat_result = sol.check();
    auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to check SAT" << std::endl;
    
    // print/return answers
    switch(sat_result) 
    {
        case z3::unsat:
        case z3::unknown:
        {
            std::cout << "[INFO] SAT failure!" << std::endl;
            return false;    
        }
        case z3::sat:
        {
            std::cout << "[INFO] SAT success!" << std::endl;
            z3::model m = sol.get_model();

            G_ret = gf2_to_bool(z3_eval_Eigen_type(m, G));
            H_ret = gf2_to_bool(z3_eval_Eigen_type(m, H));
            R_ret = gf2_to_bool(z3_eval_Eigen_type(m, R));
            return true;    
        }
        default:
            std::cout << "[ERROR] impossible SAT result: " << sat_result << std::endl;
            assert(0 && "impossible result!");
    }
}

/**
 * @brief establishes SAT solver constraints for solving for an ECC function
 * with the constraint that, given N precorrection errors, at most M
 * postcorrection errors will occur
 *
 * @param n_errs_pre_correction number of pre-correciton errors to inject
 * @param max_n_errors_post_correction maximum number of postcorrection errors to allow
 */
void unk_hamming_code::find_ecc_function_impl(int n_errs_pre_correction, int max_n_errors_post_correction)
{
    printf("[INFO] searching for ECC function that has <= %d errors post-correction given %d pre-correction\n"
        , max_n_errors_post_correction, n_errs_pre_correction);

    // generate error masks
    std::vector< int > coeffs_n(N, 1);
    std::vector< int > coeffs_em(1u << n_errs_pre_correction, 1);
    z3::expr_vector mask_vars(ctx);
    z3::expr_vector all_mask_var_constraints(ctx);
    std::vector< Eigen::Matrix< gf2, 1, Eigen::Dynamic > > error_masks_pre_permutation;
    for(uint64_t i = 1; i < (1u << n_errs_pre_correction); i++)
    {
        if((i & (i - 1)) == 0) // skip po2 because they have only 1 bit set -> correctable errors! no need to waste effort
            continue;
            
        int bitweight = hamming_weight(i);
        z3::expr_vector this_mask(ctx);
        error_masks_pre_permutation.push_back(Eigen::Matrix< gf2, 1, Eigen::Dynamic >::Zero(N));
        for(uint64_t j = 0; j < N; j++)
        {
            z3::expr bit = ctx.bool_const((std::string("mask_") + std::to_string(i) + "_" + std::to_string(j)).c_str());
            error_masks_pre_permutation.back()(j) = bit;
            mask_vars.push_back(bit);
            this_mask.push_back(bit);
        }
        all_mask_var_constraints.push_back(z3::pbeq(this_mask, coeffs_n.data(), bitweight));
    }
    for(uint64_t i = 0; i < N; i++)
    {
        z3::expr_vector this_col(ctx);
        for(const auto &m : error_masks_pre_permutation)
            this_col.push_back(m(i).get_z3_expr());

        z3::expr all_zero = !z3::mk_or(this_col);
        z3::expr n_not_zero = z3::pbeq(this_col, coeffs_em.data(), (1 << (n_errs_pre_correction - 1)) - 1);
        all_mask_var_constraints.push_back(all_zero || n_not_zero);
    }
     
    // compute the full transformation for every error mask and identify which bits are susceptible to error
    std::vector< z3::expr > per_bit_syndrome_producible(K, ctx.bool_val(0));
    for(const auto &mask : error_masks_pre_permutation)
    {
        // create a new error mask with this mask
        Eigen::Matrix< gf2, 1, Eigen::Dynamic > error_mask = mask;

        // calculate resultant codeword
        Eigen::Matrix< gf2, 1, Eigen::Dynamic > corrected_codeword;
        Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword_p;
        Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome;
        calculate_corrected_dataword_from_error_mask_given_zero_codeword(error_mask, corrected_codeword, dataword_p, syndrome);

        // log all errors
        for(uint64_t i = 0; i < K; i++)
            per_bit_syndrome_producible[i] = per_bit_syndrome_producible[i] || dataword_p(i).get_z3_expr();
    }
    z3::expr_vector per_bit_error_possibility(ctx);
    for(uint64_t i = 0; i < K; i++)
        per_bit_error_possibility.push_back(per_bit_syndrome_producible[i]);

    // create a single expression to constrain the mask vars - only need to enforce it separately if we want to find intermediate vars
    z3::expr mask_var_constraints = z3::mk_and(all_mask_var_constraints);
    // sol.add(mask_var_constraints); // only necessary if we want to actually find such a pattern

    // ENFORCE that there exists no permutation matrix SUCH THAT there are at least max_n post-correction errors
    // sol.add(z3::pbeq(per_bit_error_possibility, coeffs_n.data(), 4)); // THIS WORKS CORRECTLY
    // sol.add(z3::exists(mask_vars, mask_var_constraints && z3::pbeq(per_bit_error_possibility, coeffs_n.data(), 4)), "exists quantifier"); // THIS WORKS CORRECTLY
    sol.add(!z3::exists(mask_vars, mask_var_constraints && z3::pbge(per_bit_error_possibility, coeffs_n.data(), max_n_errors_post_correction + 1)), "exists quantifier"); // GE requires that we add one to be greater than the max
}

/**
 * @brief Searches for an ECC function with the given constraints
 * 
 * @param G_out discovered G matrix 
 * @param H_out discovered H matrix
 * @param R_out discovered R matrix
 * @param Hs_to_exclude H matrices to avoid finding again
 * @param is_rerun whether this is an incremental solve building upon a previous invocation
 * @return success or failure
 */
bool unk_hamming_code::find_ecc_function(
      Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &G_out
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_out
    , Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &R_out
    , int n_errs_pre_correction
    , int max_n_errors_post_correction
    , const std::vector< Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > > Hs_to_exclude
    , const bool is_rerun)
{
    if(is_rerun)
    {
        if(g_verbosity >= 1)
            std::cout << "[INFO] re-using previous constraints" << std::endl;
        solver_pop();
    }
    else
    {
        std::cout << "[INFO] no previous constraints to re-use: establishing constraints" << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        find_ecc_function_impl(n_errs_pre_correction, max_n_errors_post_correction);
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to establish constraints" << std::endl;
    }

    // avoid finding the same matrices as previously found -> total order on
    // matrix ensures that each combination of syndromes is unique
    if(g_verbosity >= 1)
        std::cout << "[INFO] excluding already-discovered schemes" << std::endl;
    for(size_t i = 0; i < Hs_to_exclude.size(); i++)
    {
        const Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > &H_exclude = Hs_to_exclude.at(i);
        std::vector< gf2 > equalities;
        for(uint64_t row = 0; row < NmK; row++)
            for(uint64_t col = 0; col < N; col++)
                equalities.push_back(H(row, col) == H_exclude(row, col));
        sol.add(!z3_and_reduce(equalities).get_z3_expr());
    }

    // save assertions for re-use when trying to find additional solutions
    if(g_verbosity >= 1)
        std::cout << "[INFO] saving solver state for next run" << std::endl;
    solver_push();

    // check SAT
    std::cout << "[INFO] checking satisfiability of find_ecc_function" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto sat_result = sol.check();
    auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to check SAT" << std::endl;
    
    // print/return answers
    switch(sat_result) 
    {
        case z3::unsat:
        {
            std::cout << "[INFO] UNSAT!" << std::endl;
            
            if(g_verbosity > 1)
            {
                z3::expr_vector core = sol.unsat_core();
                std::cout << "[DEBUG2] unsat core: " << core << std::endl;
                std::cout << "[DEBUG2] size: " << core.size() << std::endl;
                for (unsigned i = 0; i < core.size(); i++) 
                {
                    std::cout  << "[DEBUG2] " << core[i] << std::endl;
                }
            }
            break;
        }
        case z3::unknown:
        {
            std::cout << "[INFO] SAT failure (unknown)!" << std::endl;
            
            z3::expr_vector core = sol.unsat_core();
            std::cout << "[ERROR] unsat core: " << core << std::endl;
            std::cout << "[ERROR] size: " << core.size() << std::endl;
            for (unsigned i = 0; i < core.size(); i++) 
            {
            std::cout  << "[ERROR] " << core[i] << std::endl;
            }
            assert(0 && "bug in SAT constraints - should not occur");
        }
        case z3::sat:
        {
            std::cout << "[INFO] SAT success!" << std::endl;
            z3::model m = sol.get_model();

            G_out = gf2_to_bool(z3_eval_Eigen_type(m, G));
            H_out = gf2_to_bool(z3_eval_Eigen_type(m, H));
            R_out = gf2_to_bool(z3_eval_Eigen_type(m, R));

            return true;    
        }
        default:
            std::cout << "[ERROR] impossible SAT result: " << sat_result << std::endl;
            assert(0 && "impossible result!");
    }
    return false;
}

// skips postcorrection errors that are already known
bool unk_hamming_code::find_possible_postcorrection_errors(
      const std::unordered_set< uint64_t > &precorrection_errors
    , std::unordered_set< uint64_t > &postcorrection_errors
) const
{
    // exists a data pattern and error combination such that this bit has an error.
    z3::expr_vector unknowns(ctx);
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword(1, K);
    for(uint64_t k = 0; k < K; k++)
    {
        z3::expr bit = ctx.bool_const((std::string("dw_") + std::to_string(k)).c_str());
        dataword(k) = bit;
        unknowns.push_back(bit);
    }

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword;
    codeword = dataword * G;

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword_p(1, N);
    for(uint64_t n = 0; n < N; n++)
    {
        if(precorrection_errors.count(n) > 0)
        {
            z3::expr err = ctx.bool_const((std::string("em_") + std::to_string(n)).c_str()); 
            codeword_p(n) = codeword(n) && err; // 1 -> 0 <ONLY!> <- assumes ALL_TRUE distribution
            unknowns.push_back(err);
        }
        else
            codeword_p(n) = codeword(n); // 1 means no error possible
    }

    // Eigen::Matrix< gf2, 1, Eigen::Dynamic > error_mask(1, N);
    // for(uint64_t n = 0; n < N; n++)
    // {
    //     if(precorrection_errors.count(n) > 0)
    //     {
    //         z3::expr err = ctx.bool_const((std::string("em_") + std::to_string(n)).c_str()); 
    //         error_mask(n) = err;
    //         unknowns.push_back(err);
    //     }
    //     else
    //         error_mask(n) = 1; // 1 means no error possible
    // }


    // calculate everything
    // enum error_operator eo = EO_ALL_T;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > corrected_codeword;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword_p;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome;

    // codeword_p.resize(1, N);
    // apply_error_operator_to_codeword(codeword_p, codeword, error_mask, eo);
    // for(uint64_t i = 0; i < N; i++)
    //     codeword_p(i) = codeword(i) && error_mask(i);
    syndrome = codeword_p * HT;

    corrected_codeword.resize(1, N);
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        z3::expr_vector is_syn_pos_eq_arr(ctx);
        for(uint64_t ridx = 0; ridx < NmK; ridx++)
            is_syn_pos_eq_arr.push_back((H(ridx, bidx) == syndrome(ridx)).get_z3_expr());
        corrected_codeword(bidx) = codeword_p(bidx) != z3::mk_and(is_syn_pos_eq_arr);
    }
    dataword_p = corrected_codeword * RT;

    // run solver for each bit
    for(uint64_t dw_bit_idx = 0; dw_bit_idx < K; dw_bit_idx++)
    {
        if(postcorrection_errors.count(dw_bit_idx) > 0)
            continue;

        // save constraints
        solver_push();

        // add our constraint
        // z3::expr error_possible = z3::exists(unknowns, (dataword_p(dw_bit_idx) != dataword(dw_bit_idx)).get_z3_expr());
        z3::expr error_possible = (dataword_p(dw_bit_idx) != dataword(dw_bit_idx)).get_z3_expr();
        sol.add(error_possible);

        // check SAT
        // std::cout << "[INFO] checking satisfiability of find_possible_postcorrection_errors for dataword_bit_idx: " << dw_bit_idx << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        auto sat_result = sol.check();
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        // std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to check SAT" << std::endl;
    
        // print/return answers
        switch(sat_result) 
        {
            case z3::unsat:
                break; // no way to have an uncorrectable error in this bit
            case z3::sat:
                postcorrection_errors.insert(dw_bit_idx);
                break; // can form this bit
            case z3::unknown:
            {
                std::cout << "[INFO] SAT failure (unknown)!" << std::endl;
                
                z3::expr_vector core = sol.unsat_core();
                std::cout << "[ERROR] unsat core: " << core << std::endl;
                std::cout << "[ERROR] size: " << core.size() << std::endl;
                for (unsigned i = 0; i < core.size(); i++) 
                {
                std::cout  << "[ERROR] " << core[i] << std::endl;
                }

                return true; // ASSUME TIMEOUT
                // assert(0 && "bug in SAT constraints - should not occur");
            }
            default:
                std::cout << "[ERROR] impossible SAT result: " << sat_result << std::endl;
                assert(0 && "impossible result!");
        }

        // restore constraints
        solver_pop();
    }
    return false;
}

// return: true == timeout
bool unk_hamming_code::get_histogram_of_n_simultaneous_errors(
      const std::unordered_set< uint64_t > &precorrection_errors
    , const std::unordered_set< uint64_t > &postcorrection_errors
    , const std::unordered_set< uint64_t > &known_errors
    // , std::unordered_set< uint64_t > &histogram
    , uint64_t &max_n_unknown_posc_errors
) const
{
    // just to be safe, filter the known errors to remove non-postcorrection errors
    std::unordered_set< uint64_t > known_postcorrection_errors;
    for(const uint64_t bit_idx: known_errors)
        if(postcorrection_errors.count(bit_idx) > 0)
            known_postcorrection_errors.insert(bit_idx);

    // clearly, if all possible errors are known, we will never have errors
    if(known_postcorrection_errors.size() == postcorrection_errors.size())
    {
        max_n_unknown_posc_errors = 0;
        return false;
    }

    // exists a data pattern and error combination such that there is an N-bit error
    z3::expr_vector unknowns(ctx);
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword(1, K);
    for(uint64_t k = 0; k < K; k++)
    {
        z3::expr bit = ctx.bool_const((std::string("dw_") + std::to_string(k)).c_str());
        dataword(k) = bit;
        unknowns.push_back(bit);    
    }

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword;
    codeword = dataword * G;

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword_p(1, N);
    for(uint64_t n = 0; n < N; n++)
    {
        if(precorrection_errors.count(n) > 0)
        {
            z3::expr err = ctx.bool_const((std::string("em_") + std::to_string(n)).c_str()); 
            codeword_p(n) = codeword(n) && err; // 1 -> 0 <ONLY!> <- assumes ALL_TRUE distribution
            unknowns.push_back(err);
        }
        else
            codeword_p(n) = codeword(n); // 1 means no error possible
    }

    // calculate everything
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > corrected_codeword;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword_p;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome;

    syndrome = codeword_p * HT;

    corrected_codeword.resize(1, N);
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        z3::expr_vector is_syn_pos_eq_arr(ctx);
        for(uint64_t ridx = 0; ridx < NmK; ridx++)
            is_syn_pos_eq_arr.push_back((H(ridx, bidx) == syndrome(ridx)).get_z3_expr());
        corrected_codeword(bidx) = codeword_p(bidx) != z3::mk_and(is_syn_pos_eq_arr);
    }
    dataword_p = corrected_codeword * RT;

    // add constraints to speed up solving - errors are only possible in known bits
    z3::expr_vector unknown_postcorrection_errors(ctx);
    std::vector< int > coeffs(K, 1);
    for(uint64_t dw_bit_idx = 0; dw_bit_idx < K; dw_bit_idx++)
    {
        if(postcorrection_errors.count(dw_bit_idx) > 0 && known_postcorrection_errors.count(dw_bit_idx) == 0)
        {
            z3::expr err_this_bit = (dataword_p(dw_bit_idx) != dataword(dw_bit_idx)).get_z3_expr();
            unknown_postcorrection_errors.push_back(err_this_bit);
        }
        else
            sol.add((dataword_p(dw_bit_idx) == dataword(dw_bit_idx)).get_z3_expr()); // does not fail - hopefully speeds up the solve
    }


    // run solver for each simultaneous combination of errors
    // for(size_t n_errors = 1; n_errors <= postcorrection_errors.size() - known_postcorrection_errors.size(); n_errors++)
    for(size_t n_errors = postcorrection_errors.size() - known_postcorrection_errors.size(); n_errors > 0; n_errors--)
    {
        // std::cout << std::endl;
        // std::cout << "[DEBUG] precorrection_errors:  " << set_to_string(precorrection_errors) << std::endl;
        // std::cout << "[DEBUG] postcorrection_errors: " << set_to_string(postcorrection_errors) << std::endl;
        // std::cout << "[DEBUG] known_errors:          " << set_to_string(known_errors) << std::endl;
        // std::cout << "[DEBUG] n_unknowns:            " << std::to_string(unknown_postcorrection_errors.size()) << std::endl;
        // std::cout << "[DEBUG] now testing n_errors:  " << n_errors << std::endl;
        
        // save constraints
        solver_push();

        // add our constraint
        // z3::expr error_possible = z3::exists(unknowns, (dataword_p(dw_bit_idx) != dataword(dw_bit_idx)).get_z3_expr());
        // sol.add(z3::pbeq(unknown_postcorrection_errors, coeffs.data(), n_errors));
        sol.add(z3::exists(unknowns, z3::pbeq(unknown_postcorrection_errors, coeffs.data(), n_errors)));
        // sol.add(!z3::exists(unknowns, z3::pbeq(unknown_postcorrection_errors, coeffs.data(), n_errors)));

        // check SAT
        // std::cout << "[INFO] checking satisfiability of get_histogram_of_n_simultaneous_errors for n_errors: " << n_errors << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        auto sat_result = sol.check();
        auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
        // std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to check SAT" << std::endl;
    
        // print/return answers
        switch(sat_result) 
        {
            case z3::unsat:
                break; // no way to have an uncorrectable error in this bit
            case z3::sat:
                max_n_unknown_posc_errors = n_errors;                
                return false;
            case z3::unknown:
            {
                std::cout << "[INFO] SAT failure (unknown)!" << std::endl;
                
                z3::expr_vector core = sol.unsat_core();
                std::cout << "[ERROR] unsat core: " << core << std::endl;
                std::cout << "[ERROR] size: " << core.size() << std::endl;
                for (unsigned i = 0; i < core.size(); i++) 
                {
                std::cout  << "[ERROR] " << core[i] << std::endl;
                }

                return true; // ASSUME TIMEOUT
                // assert(0 && "bug in SAT constraints - should not occur");
            }
            default:
                std::cout << "[ERROR] impossible SAT result: " << sat_result << std::endl;
                assert(0 && "impossible result!");
        }

        // restore constraints
        solver_pop();
    }

    // no error is possible, but we had to figure that out the hard way >.>
    max_n_unknown_posc_errors = 0;
    return false;
}

// find a test pattern such that a miscorrection is possible
bool unk_hamming_code::find_beep_pattern(
          const std::unordered_set< uint64_t > &known_prec_errors
        , Eigen::Matrix< gf2_noz3, 1, Eigen::Dynamic > &dw
        , uint64_t bit_under_test
        , bool use_worst_case_constraint) const
{
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword(1, K);
    for(uint64_t k = 0; k < K; k++)
        dataword(k) = ctx.bool_const((std::string("dw_") + std::to_string(k)).c_str());

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword;
    codeword = dataword * G;
    sol.add(codeword(bit_under_test).get_z3_expr()); // the bit under test MUST be set (and MUST fail later)
    if(use_worst_case_constraint)
    {
        if(bit_under_test > 0)
            sol.add(!codeword(bit_under_test - 1).get_z3_expr()); // the bit under test MUST be set (and MUST fail later)
        if(bit_under_test < N - 1)
            sol.add(!codeword(bit_under_test + 1).get_z3_expr()); // the bit under test MUST be set (and MUST fail later)
    }

    Eigen::Matrix< gf2, 1, Eigen::Dynamic > codeword_p(1, N);
    for(uint64_t n = 0; n < N; n++)
    {
        if(known_prec_errors.count(n) > 0)
        {
            z3::expr err = ctx.bool_const((std::string("em_") + std::to_string(n)).c_str()); 
            codeword_p(n) = codeword(n) && err; // 1 -> 0 <ONLY!> <- assumes ALL_TRUE distribution
        }
        else if(n == bit_under_test)
            codeword_p(n) = 0; // the bit under test MUST fail
        else
            codeword_p(n) = codeword(n); // no errors expected anywhere else
    }

    // calculate everything
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > corrected_codeword;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > dataword_p;
    Eigen::Matrix< gf2, 1, Eigen::Dynamic > syndrome;

    syndrome = codeword_p * HT;

    corrected_codeword.resize(1, N);
    for(uint64_t bidx = 0; bidx < N; bidx++)
    {
        z3::expr_vector is_syn_pos_eq_arr(ctx);
        for(uint64_t ridx = 0; ridx < NmK; ridx++)
            is_syn_pos_eq_arr.push_back((H(ridx, bidx) == syndrome(ridx)).get_z3_expr());
        corrected_codeword(bidx) = codeword_p(bidx) != z3::mk_and(is_syn_pos_eq_arr);
    }
    dataword_p = corrected_codeword * RT;

    // add the constraint that there is a miscorrection present
    // z3::expr error_possible = z3::exists(unknowns, (dataword_p(dw_bit_idx) != dataword(dw_bit_idx)).get_z3_expr());
    z3::expr_vector miscorrections_per_bit(ctx);
    for(uint64_t bidx = 0; bidx < K; bidx++)
        miscorrections_per_bit.push_back((dataword_p(bidx) && !dataword(bidx)).get_z3_expr()); // 0->1 is definitely a miscorrection in the data-retention error model
    sol.add(z3::mk_or(miscorrections_per_bit));

    // check SAT
    // std::cout << "[INFO] checking satisfiability of find_possible_postcorrection_errors for dataword_bit_idx: " << dw_bit_idx << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    auto sat_result = sol.check();
    auto time_elapsed = std::chrono::high_resolution_clock::now() - start_time;
    // std::cout << "[TIME] " << time_elapsed / std::chrono::milliseconds(1) << " ms to check SAT" << std::endl;

    // print/return answers
    switch(sat_result) 
    {
        case z3::unknown:
        {
            // std::cout << "[INFO] SAT failure (unknown)!" << std::endl;
            
            // z3::expr_vector core = sol.unsat_core();
            // std::cout << "[ERROR] unsat core: " << core << std::endl;
            // std::cout << "[ERROR] size: " << core.size() << std::endl;
            // for (unsigned i = 0; i < core.size(); i++) 
            // {
            // std::cout  << "[ERROR] " << core[i] << std::endl;
            // }

            return false; // SILENT TIMEOUT
            // assert(0 && "bug in SAT constraints - should not occur");
        }
        case z3::unsat:
            return false;

        case z3::sat:
        {
            z3::model m = sol.get_model();
            for(uint64_t k = 0; k < K; k++)
                dw(k) = m.eval(dataword(k).get_z3_expr()).is_true();

            // we can evaluate expressions in the model.
            if(g_verbosity >= 2)
            {
                std::cout << "[DEBUG2] SAT results for BEEP test pattern:" << std::endl;
                std::cout << "[DEBUG2]     dataword:  " << gf2_to_int(z3_eval_Eigen_type(m, dataword)) << std::endl;
                std::cout << "[DEBUG2]     codeword:  " << gf2_to_int(z3_eval_Eigen_type(m, codeword)) << std::endl;
                std::cout << "[DEBUG2]     codeword_p:" << gf2_to_int(z3_eval_Eigen_type(m, codeword_p)) << std::endl;
                std::cout << "[DEBUG2]     codeword_p:" << gf2_to_int(z3_eval_Eigen_type(m, codeword_p)) << std::endl;
            }

            return true; // found a valid test pattern
        }
        
        default:
            std::cout << "[ERROR] impossible SAT result: " << sat_result << std::endl;
            assert(0 && "impossible result!");
    }
}