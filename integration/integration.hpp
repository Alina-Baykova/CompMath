#include <vector>
#include <type_traits>
#include <algorithm>
#include <map>

const std::vector<double> LEJENDRE_ZEROS_2    = {0.577350269189626};
const std::vector<double> LEJENDRE_WEIGHTS_2  = {1.};

const std::vector<double> LEJENDRE_ZEROS_3    = {0., 0.774596669241483};
const std::vector<double> LEJENDRE_WEIGHTS_3  = {0.888888888888889, 0.55555555555555};

const std::vector<double> LEJENDRE_ZEROS_4    = {0.339981043584856, 0.86113631159405};
const std::vector<double> LEJENDRE_WEIGHTS_4  = {0.652145154862546, 0.34785484513745};


const auto Gauss2 = std::make_pair(LEJENDRE_ZEROS_2, LEJENDRE_WEIGHTS_2);
const auto Gauss3 = std::make_pair(LEJENDRE_ZEROS_3, LEJENDRE_WEIGHTS_3);
const auto Gauss4 = std::make_pair(LEJENDRE_ZEROS_4, LEJENDRE_WEIGHTS_4);


const std::map<std::size_t, std::pair<std::vector<double>, std::vector<double>>> availableQuadratures {
    {2, Gauss2},
    {3, Gauss3},
    {4, Gauss4}
};


template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Функция производит интегрирование на одном отрезке */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end){
    typename ArgumentGetter<Callable>::Argument ans = 0;
    const auto quadrature = availableQuadratures.at(N);
    const auto delta = end - start;
    if(N % 2 != 0) {
        ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[0] * (delta) / 2) * quadrature.second[0];
        for(std::size_t i = 1; i < N / 2 + 1; i++) {
            ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            ans += (delta) / 2 * func((end + start) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
        }
    }
    else {
        for(std::size_t i = 0; i < N / 2; i++) {
            ans += (delta) / 2 * func((end + start) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            ans += (delta) / 2 * func((end + start) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
        }    
    }
    return ans;
}


/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, std::size_t N>
decltype(auto) integrate(const Callable& func, const typename ArgumentGetter<Callable>::Argument& start, const typename ArgumentGetter<Callable>::Argument& end, const Dif<typename ArgumentGetter<Callable>::Argument>& dx) {
    typename ArgumentGetter<Callable>::Argument ans = 0;
    const auto quadrature = availableQuadratures.at(N);
    const auto Delta = end - start;
    const unsigned int numberOfIntervals = std::ceil(Delta / dx);
    const Dif<typename ArgumentGetter<Callable>::Argument> exactDx = Delta / numberOfIntervals;
    Dif<typename ArgumentGetter<Callable>::Argument> currStart = start;
    Dif<typename ArgumentGetter<Callable>::Argument> currEnd = start + exactDx;
    for(unsigned int j = 0; j < numberOfIntervals; j++) {
        if(N % 2 != 0) {
            auto delta = currEnd - currStart;
            ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[0] * (delta) / 2) * quadrature.second[0];
            for(std::size_t i = 1; i < N / 2 + 1; i++) {
                ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
                ans += (delta) / 2 * func((currEnd + currStart) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            }
            currStart += exactDx;
            currEnd += exactDx;
        }
        else {
            auto delta = currEnd - currStart;
            for(std::size_t i = 0; i < N / 2; i++) {
                ans += (delta) / 2 * func((currEnd + currStart) / 2 + quadrature.first[i] * (delta) / 2) * quadrature.second[i];
                ans += (delta) / 2 * func((currEnd + currStart) / 2 - quadrature.first[i] * (delta) / 2) * quadrature.second[i];
            }
            currStart += exactDx;
            currEnd += exactDx;  
        }
    }
    return ans;
}
