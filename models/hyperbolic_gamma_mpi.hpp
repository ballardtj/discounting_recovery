// Code generated by Stan version 2.18.1

#include <stan/model/model_header.hpp>

namespace hyperbolic_gamma_mpi_model_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "/QRISdata/Q0992/models/hyperbolic_gamma_mpi.stan");
    reader.add_event(102, 100, "end", "/QRISdata/Q0992/models/hyperbolic_gamma_mpi.stan");
    return reader;
}

template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic,1>
likelihood(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& phi,
               const Eigen::Matrix<T1__, Eigen::Dynamic,1>& theta,
               const std::vector<T2__>& real_data,
               const std::vector<int>& int_data, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 10;
        int Nvalid(0);
        (void) Nvalid;  // dummy to suppress unused var warning

        stan::math::fill(Nvalid, std::numeric_limits<int>::min());
        stan::math::assign(Nvalid,get_base1(int_data,1,"int_data",1));
        current_statement_begin__ = 11;
        int Nplaces(0);
        (void) Nplaces;  // dummy to suppress unused var warning

        stan::math::fill(Nplaces, std::numeric_limits<int>::min());
        stan::math::assign(Nplaces,get_base1(int_data,2,"int_data",1));
        current_statement_begin__ = 14;
        local_scalar_t__ k;
        (void) k;  // dummy to suppress unused var warning

        stan::math::initialize(k, DUMMY_VAR__);
        stan::math::fill(k,DUMMY_VAR__);
        stan::math::assign(k,get_base1(theta,1,"theta",1));
        current_statement_begin__ = 15;
        local_scalar_t__ sigma;
        (void) sigma;  // dummy to suppress unused var warning

        stan::math::initialize(sigma, DUMMY_VAR__);
        stan::math::fill(sigma,DUMMY_VAR__);
        stan::math::assign(sigma,(get_base1(phi,3,"phi",1) + (get_base1(phi,4,"phi",1) * get_base1(theta,2,"theta",1))));
        current_statement_begin__ = 19;
        validate_non_negative_index("y", "Nvalid", Nvalid);
        vector<int> y(Nvalid, 0);
        stan::math::fill(y, std::numeric_limits<int>::min());
        current_statement_begin__ = 20;
        validate_non_negative_index("p_a_logit", "Nvalid", Nvalid);
        Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  p_a_logit(static_cast<Eigen::VectorXd::Index>(Nvalid));
        (void) p_a_logit;  // dummy to suppress unused var warning

        stan::math::initialize(p_a_logit, DUMMY_VAR__);
        stan::math::fill(p_a_logit,DUMMY_VAR__);
        current_statement_begin__ = 21;
        local_scalar_t__ m_a;
        (void) m_a;  // dummy to suppress unused var warning

        stan::math::initialize(m_a, DUMMY_VAR__);
        stan::math::fill(m_a,DUMMY_VAR__);
        current_statement_begin__ = 22;
        local_scalar_t__ m_b;
        (void) m_b;  // dummy to suppress unused var warning

        stan::math::initialize(m_b, DUMMY_VAR__);
        stan::math::fill(m_b,DUMMY_VAR__);
        current_statement_begin__ = 23;
        local_scalar_t__ d_a;
        (void) d_a;  // dummy to suppress unused var warning

        stan::math::initialize(d_a, DUMMY_VAR__);
        stan::math::fill(d_a,DUMMY_VAR__);
        current_statement_begin__ = 24;
        local_scalar_t__ d_b;
        (void) d_b;  // dummy to suppress unused var warning

        stan::math::initialize(d_b, DUMMY_VAR__);
        stan::math::fill(d_b,DUMMY_VAR__);
        current_statement_begin__ = 25;
        local_scalar_t__ u_a;
        (void) u_a;  // dummy to suppress unused var warning

        stan::math::initialize(u_a, DUMMY_VAR__);
        stan::math::fill(u_a,DUMMY_VAR__);
        current_statement_begin__ = 26;
        local_scalar_t__ u_b;
        (void) u_b;  // dummy to suppress unused var warning

        stan::math::initialize(u_b, DUMMY_VAR__);
        stan::math::fill(u_b,DUMMY_VAR__);
        current_statement_begin__ = 27;
        local_scalar_t__ lp;
        (void) lp;  // dummy to suppress unused var warning

        stan::math::initialize(lp, DUMMY_VAR__);
        stan::math::fill(lp,DUMMY_VAR__);


        current_statement_begin__ = 29;
        for (int i = 1; i <= Nvalid; ++i) {

            current_statement_begin__ = 30;
            stan::math::assign(m_a, get_base1(real_data,i,"real_data",1));
            current_statement_begin__ = 31;
            stan::math::assign(m_b, get_base1(real_data,(Nplaces + i),"real_data",1));
            current_statement_begin__ = 32;
            stan::math::assign(d_a, get_base1(real_data,((2 * Nplaces) + i),"real_data",1));
            current_statement_begin__ = 33;
            stan::math::assign(d_b, get_base1(real_data,((3 * Nplaces) + i),"real_data",1));
            current_statement_begin__ = 35;
            stan::math::assign(u_a, (m_a / (1 + (k * d_a))));
            current_statement_begin__ = 36;
            stan::math::assign(u_b, (m_b / (1 + (k * d_b))));
            current_statement_begin__ = 37;
            stan::model::assign(p_a_logit, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        ((u_a - u_b) * sigma), 
                        "assigning variable p_a_logit");
            current_statement_begin__ = 39;
            stan::model::assign(y, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        get_base1(int_data,(2 + i),"int_data",1), 
                        "assigning variable y");
        }
        current_statement_begin__ = 42;
        stan::math::assign(lp, bernoulli_logit_log(y,p_a_logit));
        current_statement_begin__ = 44;
        return stan::math::promote_scalar<fun_return_scalar_t__>(transpose(stan::math::to_row_vector(stan::math::array_builder<local_scalar_t__ >().add(lp).array())));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct likelihood_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic,1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& phi,
               const Eigen::Matrix<T1__, Eigen::Dynamic,1>& theta,
               const std::vector<T2__>& real_data,
               const std::vector<int>& int_data, std::ostream* pstream__) const {
        return likelihood(phi, theta, real_data, int_data, pstream__);
    }
};

class hyperbolic_gamma_mpi_model : public prob_grad {
private:
    int Nsubj;
    int Max_obs_per_subj;
    vector<vector<double> > real_data;
    vector<vector<int> > int_data;
public:
    hyperbolic_gamma_mpi_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    hyperbolic_gamma_mpi_model(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "hyperbolic_gamma_mpi_model_namespace::hyperbolic_gamma_mpi_model";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 49;
            context__.validate_dims("data initialization", "Nsubj", "int", context__.to_vec());
            Nsubj = int(0);
            vals_i__ = context__.vals_i("Nsubj");
            pos__ = 0;
            Nsubj = vals_i__[pos__++];
            current_statement_begin__ = 50;
            context__.validate_dims("data initialization", "Max_obs_per_subj", "int", context__.to_vec());
            Max_obs_per_subj = int(0);
            vals_i__ = context__.vals_i("Max_obs_per_subj");
            pos__ = 0;
            Max_obs_per_subj = vals_i__[pos__++];
            current_statement_begin__ = 51;
            validate_non_negative_index("real_data", "Nsubj", Nsubj);
            validate_non_negative_index("real_data", "(Max_obs_per_subj * 4)", (Max_obs_per_subj * 4));
            context__.validate_dims("data initialization", "real_data", "double", context__.to_vec(Nsubj,(Max_obs_per_subj * 4)));
            validate_non_negative_index("real_data", "Nsubj", Nsubj);
            validate_non_negative_index("real_data", "(Max_obs_per_subj * 4)", (Max_obs_per_subj * 4));
            real_data = std::vector<std::vector<double> >(Nsubj,std::vector<double>((Max_obs_per_subj * 4),double(0)));
            vals_r__ = context__.vals_r("real_data");
            pos__ = 0;
            size_t real_data_limit_1__ = (Max_obs_per_subj * 4);
            for (size_t i_1__ = 0; i_1__ < real_data_limit_1__; ++i_1__) {
                size_t real_data_limit_0__ = Nsubj;
                for (size_t i_0__ = 0; i_0__ < real_data_limit_0__; ++i_0__) {
                    real_data[i_0__][i_1__] = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 52;
            validate_non_negative_index("int_data", "Nsubj", Nsubj);
            validate_non_negative_index("int_data", "(Max_obs_per_subj + 2)", (Max_obs_per_subj + 2));
            context__.validate_dims("data initialization", "int_data", "int", context__.to_vec(Nsubj,(Max_obs_per_subj + 2)));
            validate_non_negative_index("int_data", "Nsubj", Nsubj);
            validate_non_negative_index("int_data", "(Max_obs_per_subj + 2)", (Max_obs_per_subj + 2));
            int_data = std::vector<std::vector<int> >(Nsubj,std::vector<int>((Max_obs_per_subj + 2),int(0)));
            vals_i__ = context__.vals_i("int_data");
            pos__ = 0;
            size_t int_data_limit_1__ = (Max_obs_per_subj + 2);
            for (size_t i_1__ = 0; i_1__ < int_data_limit_1__; ++i_1__) {
                size_t int_data_limit_0__ = Nsubj;
                for (size_t i_0__ = 0; i_0__ < int_data_limit_0__; ++i_0__) {
                    int_data[i_0__][i_1__] = vals_i__[pos__++];
                }
            }

            // validate, data variables
            current_statement_begin__ = 49;
            current_statement_begin__ = 50;
            current_statement_begin__ = 51;
            current_statement_begin__ = 52;
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 57;
            ++num_params_r__;
            current_statement_begin__ = 58;
            ++num_params_r__;
            current_statement_begin__ = 59;
            validate_non_negative_index("k", "Nsubj", Nsubj);
            num_params_r__ += Nsubj;
            current_statement_begin__ = 61;
            ++num_params_r__;
            current_statement_begin__ = 62;
            ++num_params_r__;
            current_statement_begin__ = 63;
            validate_non_negative_index("sigma_raw", "Nsubj", Nsubj);
            num_params_r__ += Nsubj;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~hyperbolic_gamma_mpi_model() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("k_shape")))
            throw std::runtime_error("variable k_shape missing");
        vals_r__ = context__.vals_r("k_shape");
        pos__ = 0U;
        context__.validate_dims("initialization", "k_shape", "double", context__.to_vec());
        double k_shape(0);
        k_shape = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,k_shape);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable k_shape: ") + e.what());
        }

        if (!(context__.contains_r("k_scale")))
            throw std::runtime_error("variable k_scale missing");
        vals_r__ = context__.vals_r("k_scale");
        pos__ = 0U;
        context__.validate_dims("initialization", "k_scale", "double", context__.to_vec());
        double k_scale(0);
        k_scale = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,k_scale);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable k_scale: ") + e.what());
        }

        if (!(context__.contains_r("k")))
            throw std::runtime_error("variable k missing");
        vals_r__ = context__.vals_r("k");
        pos__ = 0U;
        validate_non_negative_index("k", "Nsubj", Nsubj);
        context__.validate_dims("initialization", "k", "vector_d", context__.to_vec(Nsubj));
        vector_d k(static_cast<Eigen::VectorXd::Index>(Nsubj));
        for (int j1__ = 0U; j1__ < Nsubj; ++j1__)
            k(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,k);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable k: ") + e.what());
        }

        if (!(context__.contains_r("sigma_mean")))
            throw std::runtime_error("variable sigma_mean missing");
        vals_r__ = context__.vals_r("sigma_mean");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigma_mean", "double", context__.to_vec());
        double sigma_mean(0);
        sigma_mean = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigma_mean);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_mean: ") + e.what());
        }

        if (!(context__.contains_r("sigma_sd")))
            throw std::runtime_error("variable sigma_sd missing");
        vals_r__ = context__.vals_r("sigma_sd");
        pos__ = 0U;
        context__.validate_dims("initialization", "sigma_sd", "double", context__.to_vec());
        double sigma_sd(0);
        sigma_sd = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0,sigma_sd);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_sd: ") + e.what());
        }

        if (!(context__.contains_r("sigma_raw")))
            throw std::runtime_error("variable sigma_raw missing");
        vals_r__ = context__.vals_r("sigma_raw");
        pos__ = 0U;
        validate_non_negative_index("sigma_raw", "Nsubj", Nsubj);
        context__.validate_dims("initialization", "sigma_raw", "vector_d", context__.to_vec(Nsubj));
        vector_d sigma_raw(static_cast<Eigen::VectorXd::Index>(Nsubj));
        for (int j1__ = 0U; j1__ < Nsubj; ++j1__)
            sigma_raw(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_lb_unconstrain(0,sigma_raw);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable sigma_raw: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            local_scalar_t__ k_shape;
            (void) k_shape;  // dummy to suppress unused var warning
            if (jacobian__)
                k_shape = in__.scalar_lb_constrain(0,lp__);
            else
                k_shape = in__.scalar_lb_constrain(0);

            local_scalar_t__ k_scale;
            (void) k_scale;  // dummy to suppress unused var warning
            if (jacobian__)
                k_scale = in__.scalar_lb_constrain(0,lp__);
            else
                k_scale = in__.scalar_lb_constrain(0);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  k;
            (void) k;  // dummy to suppress unused var warning
            if (jacobian__)
                k = in__.vector_lb_constrain(0,Nsubj,lp__);
            else
                k = in__.vector_lb_constrain(0,Nsubj);

            local_scalar_t__ sigma_mean;
            (void) sigma_mean;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_mean = in__.scalar_lb_constrain(0,lp__);
            else
                sigma_mean = in__.scalar_lb_constrain(0);

            local_scalar_t__ sigma_sd;
            (void) sigma_sd;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_sd = in__.scalar_lb_constrain(0,lp__);
            else
                sigma_sd = in__.scalar_lb_constrain(0);

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  sigma_raw;
            (void) sigma_raw;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma_raw = in__.vector_lb_constrain(0,Nsubj,lp__);
            else
                sigma_raw = in__.vector_lb_constrain(0,Nsubj);


            // transformed parameters
            current_statement_begin__ = 69;
            validate_non_negative_index("phi", "4", 4);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  phi(static_cast<Eigen::VectorXd::Index>(4));
            (void) phi;  // dummy to suppress unused var warning

            stan::math::initialize(phi, DUMMY_VAR__);
            stan::math::fill(phi,DUMMY_VAR__);
            current_statement_begin__ = 70;
            validate_non_negative_index("theta", "2", 2);
            validate_non_negative_index("theta", "Nsubj", Nsubj);
            vector<Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> > theta(Nsubj, (Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> (static_cast<Eigen::VectorXd::Index>(2))));
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta,DUMMY_VAR__);


            current_statement_begin__ = 73;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        k_shape, 
                        "assigning variable phi");
            current_statement_begin__ = 74;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        k_scale, 
                        "assigning variable phi");
            current_statement_begin__ = 75;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        sigma_mean, 
                        "assigning variable phi");
            current_statement_begin__ = 76;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(4), stan::model::nil_index_list()), 
                        sigma_sd, 
                        "assigning variable phi");
            current_statement_begin__ = 79;
            for (int subj = 1; subj <= Nsubj; ++subj) {

                current_statement_begin__ = 80;
                stan::model::assign(theta, 
                            stan::model::cons_list(stan::model::index_uni(subj), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                            get_base1(k,subj,"k",1), 
                            "assigning variable theta");
                current_statement_begin__ = 81;
                stan::model::assign(theta, 
                            stan::model::cons_list(stan::model::index_uni(subj), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                            get_base1(sigma_raw,subj,"sigma_raw",1), 
                            "assigning variable theta");
            }

            // validate transformed parameters
            for (int i0__ = 0; i0__ < 4; ++i0__) {
                if (stan::math::is_uninitialized(phi(i0__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: phi" << '[' << i0__ << ']';
                    throw std::runtime_error(msg__.str());
                }
            }
            for (int i0__ = 0; i0__ < Nsubj; ++i0__) {
                for (int i1__ = 0; i1__ < 2; ++i1__) {
                    if (stan::math::is_uninitialized(theta[i0__](i1__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: theta" << '[' << i0__ << ']' << '[' << i1__ << ']';
                        throw std::runtime_error(msg__.str());
                    }
                }
            }

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 69;
            current_statement_begin__ = 70;

            // model body

            current_statement_begin__ = 89;
            lp_accum__.add(normal_log<propto__>(k_shape, 0, 1));
            current_statement_begin__ = 90;
            lp_accum__.add(normal_log<propto__>(k_scale, 0, 1));
            current_statement_begin__ = 92;
            lp_accum__.add(normal_log<propto__>(sigma_mean, 0, 1));
            current_statement_begin__ = 93;
            lp_accum__.add(normal_log<propto__>(sigma_sd, 0, 1));
            current_statement_begin__ = 95;
            lp_accum__.add(gamma_log<propto__>(k, k_shape, inv(k_scale)));
            current_statement_begin__ = 96;
            lp_accum__.add(normal_log<propto__>(sigma_raw, 0, 1));
            current_statement_begin__ = 99;
            lp_accum__.add(sum(map_rect<1, likelihood_functor__>(phi, theta, real_data, int_data, pstream__)));

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("k_shape");
        names__.push_back("k_scale");
        names__.push_back("k");
        names__.push_back("sigma_mean");
        names__.push_back("sigma_sd");
        names__.push_back("sigma_raw");
        names__.push_back("phi");
        names__.push_back("theta");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(Nsubj);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(Nsubj);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(4);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(Nsubj);
        dims__.push_back(2);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "hyperbolic_gamma_mpi_model_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double k_shape = in__.scalar_lb_constrain(0);
        double k_scale = in__.scalar_lb_constrain(0);
        vector_d k = in__.vector_lb_constrain(0,Nsubj);
        double sigma_mean = in__.scalar_lb_constrain(0);
        double sigma_sd = in__.scalar_lb_constrain(0);
        vector_d sigma_raw = in__.vector_lb_constrain(0,Nsubj);
        vars__.push_back(k_shape);
        vars__.push_back(k_scale);
            for (int k_0__ = 0; k_0__ < Nsubj; ++k_0__) {
            vars__.push_back(k[k_0__]);
            }
        vars__.push_back(sigma_mean);
        vars__.push_back(sigma_sd);
            for (int k_0__ = 0; k_0__ < Nsubj; ++k_0__) {
            vars__.push_back(sigma_raw[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {
            current_statement_begin__ = 69;
            validate_non_negative_index("phi", "4", 4);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  phi(static_cast<Eigen::VectorXd::Index>(4));
            (void) phi;  // dummy to suppress unused var warning

            stan::math::initialize(phi, DUMMY_VAR__);
            stan::math::fill(phi,DUMMY_VAR__);
            current_statement_begin__ = 70;
            validate_non_negative_index("theta", "2", 2);
            validate_non_negative_index("theta", "Nsubj", Nsubj);
            vector<Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> > theta(Nsubj, (Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> (static_cast<Eigen::VectorXd::Index>(2))));
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta,DUMMY_VAR__);


            current_statement_begin__ = 73;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                        k_shape, 
                        "assigning variable phi");
            current_statement_begin__ = 74;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list()), 
                        k_scale, 
                        "assigning variable phi");
            current_statement_begin__ = 75;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(3), stan::model::nil_index_list()), 
                        sigma_mean, 
                        "assigning variable phi");
            current_statement_begin__ = 76;
            stan::model::assign(phi, 
                        stan::model::cons_list(stan::model::index_uni(4), stan::model::nil_index_list()), 
                        sigma_sd, 
                        "assigning variable phi");
            current_statement_begin__ = 79;
            for (int subj = 1; subj <= Nsubj; ++subj) {

                current_statement_begin__ = 80;
                stan::model::assign(theta, 
                            stan::model::cons_list(stan::model::index_uni(subj), stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list())), 
                            get_base1(k,subj,"k",1), 
                            "assigning variable theta");
                current_statement_begin__ = 81;
                stan::model::assign(theta, 
                            stan::model::cons_list(stan::model::index_uni(subj), stan::model::cons_list(stan::model::index_uni(2), stan::model::nil_index_list())), 
                            get_base1(sigma_raw,subj,"sigma_raw",1), 
                            "assigning variable theta");
            }

            // validate transformed parameters
            current_statement_begin__ = 69;
            current_statement_begin__ = 70;

            // write transformed parameters
            if (include_tparams__) {
            for (int k_0__ = 0; k_0__ < 4; ++k_0__) {
            vars__.push_back(phi[k_0__]);
            }
            for (int k_1__ = 0; k_1__ < 2; ++k_1__) {
                for (int k_0__ = 0; k_0__ < Nsubj; ++k_0__) {
                vars__.push_back(theta[k_0__][k_1__]);
                }
            }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "hyperbolic_gamma_mpi_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "k_shape";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "k_scale";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "k" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_mean";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_sd";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_raw" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= 4; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "phi" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_1__ = 1; k_1__ <= 2; ++k_1__) {
                for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "theta" << '.' << k_0__ << '.' << k_1__;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }


        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "k_shape";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "k_scale";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "k" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_mean";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "sigma_sd";
        param_names__.push_back(param_name_stream__.str());
        for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma_raw" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
            for (int k_0__ = 1; k_0__ <= 4; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "phi" << '.' << k_0__;
                param_names__.push_back(param_name_stream__.str());
            }
            for (int k_1__ = 1; k_1__ <= 2; ++k_1__) {
                for (int k_0__ = 1; k_0__ <= Nsubj; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "theta" << '.' << k_0__ << '.' << k_1__;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef hyperbolic_gamma_mpi_model_namespace::hyperbolic_gamma_mpi_model stan_model;

STAN_REGISTER_MAP_RECT(1, hyperbolic_gamma_mpi_model_namespace::likelihood_functor__)
