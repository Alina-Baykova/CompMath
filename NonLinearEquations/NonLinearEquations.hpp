#include <stdexcept>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_4
#define M_PI_4 (M_PI / 4.0)
#endif

// Функция возвращает пару: bool (было ли найдено решение) и double (E_i)
decltype(auto) keplerSolver(double ecc, double meanAnomaly, unsigned int maxIter, double tol) {
    std::pair<bool, double> ans = {false, M_PI_4 - (M_PI_4 - ecc * std::sin(M_PI_4) - meanAnomaly) / (1 - ecc * std::cos(meanAnomaly))};
    double tmp = 0.;
    for(unsigned int i = 1; i < maxIter; i++){
        tmp = ans.second - (ans.second - ecc * std::sin(ans.second) - meanAnomaly) / (1 - ecc * std::cos(ans.second));
        if(std::abs(tmp - ans.second) <= tol) {
            ans.first = true;
            return ans;
        }
        ans.second = tmp;
    }
    return ans;
}

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(   
    const Callable& func,                                             // функция F
    const RealType& tau,                                              // шаг тау
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
    const unsigned int nIteration,                                    // количество итераций
    const RealType& tolerance                                         // точность
                    ) {
    std::pair<bool, RealType> ans = {false, initialGuess + tau * func(initialGuess)};
    RealType tmp = 0.;
    for(unsigned int i = 1; i < nIteration; i++) {
        tmp = ans.second + tau * func(ans.second);
        if(std::abs(tmp - ans.second) <= tolerance) {
            ans.first = true;
            return ans;
        }
        ans.second = tmp;
    }
    return ans;
}