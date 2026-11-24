// Standalone NSW6 snow-moment microbenchmark on synthetic positive snow inputs.
// clang++ -std=c++17 -O3 -ffp-contract=on snow_moments.cpp -o /tmp/snow_moments
// /tmp/snow_moments > snow_moments_cpu.csv
// Benchmark scope: seven distinct empirical fits in the Roh branch of mp_nsw6.
// The directly assigned second moment X and the later fall-speed ratio are not fits.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

using Cubic = std::array<double, 4>;
using Outputs = std::array<double, 7>;
struct State { double temperature_c, x; };
constexpr Outputs orders{0, 1, 2, 3, 4, 2.25, 1.625};
constexpr std::array<double, 10> a{5.065339,-0.062659,-3.032362,0.029469,-0.000285,
                                   0.31255,0.000204,0.003199,0,-0.015952};
constexpr std::array<double, 10> b{0.476221,-0.015896,0.165977,0.007468,-0.000141,
                                   0.060366,0.000079,0.000594,0,-0.003577};
constexpr double ln10 = 2.302585092994045684;
constexpr double ln2 = 0.6931471805599453094;
constexpr double log2e = 1.442695040888963407;
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559);

Cubic coefficients(const std::array<double, 10>& c, double t) {
    return {c[0] + t * (c[1] + t * (c[4] + t * c[8])),
            c[2] + t * (c[3] + t * c[6]), c[5] + t * c[7], c[9]};
}
double horner(const Cubic& c, double n) {
    return c[0] + n * (c[1] + n * (c[2] + n * c[3]));
}

// Interpolate exp on the roots of T_(D+1) after scaling onto [-log(2)/2,log(2)/2].
// Build coefficients outside the timed region. Evaluate the monomial form with FMA.
template<int D> std::array<double, D + 1> interpolate_exp() {
    std::array<long double, D + 1> node{}, dd{}, mono{};
    const long double pi = std::acos(-1.L);
    for (int j = 0; j <= D; ++j) {
        node[j] = std::log(2.L) / 2 * std::cos(pi * (2*j+1) / (2*(D+1)));
        dd[j] = std::exp(node[j]);
    }
    for (int order = 1; order <= D; ++order)
        for (int j = D; j >= order; --j)
            dd[j] = (dd[j]-dd[j-1]) / (node[j]-node[j-order]);
    mono[0] = dd[D];
    for (int j = D-1; j >= 0; --j) {
        for (int k = D-j; k >= 1; --k) mono[k] = mono[k-1] - node[j]*mono[k];
        mono[0] = dd[j] - node[j]*mono[0];
    }
    std::array<double, D + 1> result{};
    std::copy(mono.begin(), mono.end(), result.begin());
    return result;
}

template<int D> const auto polynomial = interpolate_exp<D>();

template<int D> double reduced_exp(double y) {
    if (!(std::abs(y) < 700)) return std::exp(y);
    const double k = std::nearbyint(y*log2e);
    const double r = std::fma(-k, ln2, y);
    double p = polynomial<D>[D];
    for (int j = D-1; j >= 0; --j) p = std::fma(p, r, polynomial<D>[j]);
    // k is in [-1010,1010]; construct the exact normal binary64 power of two.
    const auto bits = std::uint64_t(int(k) + 1023) << 52;
    double scale;
    std::memcpy(&scale, &bits, sizeof(scale));
    return scale*p;
}

enum class Method { Direct, SharedLog, SharedCheb7, CubicLog, Cheb5, Cheb7, Cheb9 };

template<Method method> Outputs evaluate(const State& state) {
    const auto t = std::min(-0.1, state.temperature_c);
    const auto ac = coefficients(a, t), bc = coefficients(b, t);
    Outputs out;
    if (method == Method::Direct || !(state.x > 0)) {
        for (size_t j = 0; j < orders.size(); ++j)
            out[j] = std::exp(ln10*horner(ac, orders[j])) * std::pow(state.x, horner(bc, orders[j]));
        return out;
    }
    const double logx = std::log(state.x);
    if constexpr (method == Method::SharedLog || method == Method::SharedCheb7) {
        for (size_t j = 0; j < orders.size(); ++j) {
            const double y = ln10*horner(ac, orders[j]) + logx*horner(bc, orders[j]);
            if constexpr (method == Method::SharedLog) out[j] = std::exp(y);
            else out[j] = reduced_exp<7>(y);
        }
    } else {
        Cubic c;
        for (int p = 0; p < 4; ++p) c[p] = ln10*ac[p] + logx*bc[p];
        for (size_t j = 0; j < orders.size(); ++j) {
            const double y = horner(c, orders[j]);
            if constexpr (method == Method::CubicLog) out[j] = std::exp(y);
            if constexpr (method == Method::Cheb5) out[j] = reduced_exp<5>(y);
            if constexpr (method == Method::Cheb7) out[j] = reduced_exp<7>(y);
            if constexpr (method == Method::Cheb9) out[j] = reduced_exp<9>(y);
        }
    }
    return out;
}

template<Method m> __attribute__((noinline)) void run(
    const std::vector<State>& input, std::vector<Outputs>& output) {
    for (size_t i = 0; i < input.size(); ++i) output[i] = evaluate<m>(input[i]);
}
struct Variant {
    const char* name;
    void (*run)(const std::vector<State>&, std::vector<Outputs>&);
};
constexpr std::array<Variant,7> variants{{
    {"direct",run<Method::Direct>}, {"shared_log",run<Method::SharedLog>},
    {"cubic_log",run<Method::CubicLog>}, {"cheb5",run<Method::Cheb5>},
    {"cheb7",run<Method::Cheb7>}, {"cheb9",run<Method::Cheb9>},
    {"shared_log_cheb7",run<Method::SharedCheb7>}
}};

int main() {
    constexpr size_t count = 131072;
    std::vector<State> input(count);
    std::mt19937_64 random(20261124);
    std::uniform_real_distribution<double> temp(-75.,10.), logx(-12.,-1.);
    for (auto& s : input) s = {temp(random), std::pow(10.,logx(random))};
    std::vector<Outputs> reference(count), output(count);
    variants[0].run(input,reference);
    std::array<double,variants.size()> error{};
    for (size_t m=0; m<variants.size(); ++m) {
        variants[m].run(input,output);
        for (size_t i=0; i<count; ++i) for (size_t j=0; j<orders.size(); ++j) {
            if (!(std::isfinite(output[i][j]) && reference[i][j] > 0)) return 2;
            error[m] = std::max(error[m], std::abs(output[i][j]/reference[i][j]-1));
        }
        for (double t : {-75.,-30.,-0.1,10.}) {
            const std::vector<State> dry{{t,0.}};
            std::vector<Outputs> dry_ref(1), dry_out(1);
            variants[0].run(dry,dry_ref);
            variants[m].run(dry,dry_out);
            if (dry_ref != dry_out) return 3;
        }
    }
    // Two balanced cycles. Each method occurs twice in each timing position.
    using Clock = std::chrono::steady_clock;
    volatile double sink = 0;
    std::puts("block,position,method,ns_per_state,max_relative_error");
    for (size_t block=0; block<2*variants.size(); ++block) {
        for (size_t pos=0; pos<variants.size(); ++pos) {
            size_t m = (block < variants.size() ? block+pos : block+variants.size()-pos) % variants.size();
            const auto start = Clock::now();
            variants[m].run(input,output);
            const auto stop = Clock::now();
            const double ns = std::chrono::duration<double,std::nano>(stop-start).count()/count;
            sink += output[(block*1019+pos*37)%count][0];
            std::printf("%zu,%zu,%s,%.6f,%.10e\n",block+1,pos+1,variants[m].name,ns,error[m]);
        }
    }
    return std::isfinite(double(sink)) ? 0 : 4;
}
