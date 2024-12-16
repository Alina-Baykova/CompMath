#include <Eigen/Dense>
#include <array>
#include <vector>

/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка. Я ее не заполнил */
struct RK4Table{
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = {{{0., 0., 0., 0.}, 
                                                                        {0.5, 0., 0., 0.}, 
                                                                        {0., 0.5, 0., 0.}, 
                                                                        {0., 0., 1., 0.}}};
    static constexpr std::array<double, stages> cColumn = {0, 0.5, 0.5, 1.};
    static constexpr std::array<double, stages> bString = {1./6, 1./3, 1./3, 1./6};
};

class tCube{
public:

    static constexpr unsigned int dim = 1;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg {
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{std::pow(stateAndArg.arg, 3)};
    }
};

class Oscillator {

public:

    static constexpr unsigned int dim = 2;  // размерность задачи

    using Argument = double;  // тип аргумента, тип t

    using State = Eigen::Vector<double, dim>;  // состояние

    struct StateAndArg{
        State state;
        Argument arg;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    Eigen::Vector<double, dim> calc(const StateAndArg& stateAndArg) const {
        return Eigen::Vector<double, dim>{stateAndArg.state(1), -stateAndArg.state(0)};
    } 
};

template<typename aType, typename bType, std::size_t N>
decltype(auto) dot(const std::array<aType, N> &a, const std::array<bType, N> &b) {
    bType ans = a[0] * b[0];
    for (std::size_t i = 1; i < a.size(); i++) {
        ans += a[i] * b[i];
    }
    return ans;
}


template<typename Table, typename RHS>
std::vector<typename RHS::StateAndArg> integrate(
        const typename RHS::StateAndArg &initialState,
        const typename RHS::Argument &endTime,
        double step,
        const RHS &rhs
) {
    using uType = typename RHS::StateAndArg;
    using fType = Eigen::Vector<double, RHS::dim>;
    using rType = double;
    uType uCurrent = initialState;
    std::vector<uType> uVector{uCurrent};
    std::size_t numStep = endTime / step;
    uVector.reserve(numStep);

    for (std::size_t i = 1; i <= numStep; i++) {
        std::array<fType, Table::stages> k{};
        k[0] = rhs.calc(uCurrent);
        for(unsigned int i = 1; i < Table::stages; i++) {
            k[i] = rhs.calc({uCurrent.state + step * dot<rType, fType, Table::stages>(Table::table[i], k), uCurrent.arg + step * Table::cColumn[i]});
        }
        uCurrent = {uCurrent.state + step * dot<rType, fType, Table::stages>(Table::bString, k), uCurrent.arg + step};
        uVector.push_back(uCurrent);
    }
    
    if(endTime - step * double(numStep) > 0) {
        std::array<fType, Table::stages> k{};
        k[0] = rhs.calc(uCurrent);
        for(unsigned int i = 1; i < Table::stages; i++) {
            k[i] = rhs.calc({uCurrent.state + step * dot<rType, fType, Table::stages>(Table::table[i], k), uCurrent.arg + step * Table::cColumn[i]});
        }
        uCurrent = {uCurrent.state + step * dot<rType, fType, Table::stages>(Table::bString, k), uCurrent.arg + step};
        uVector.push_back(uCurrent);
    }
    return uVector;
}