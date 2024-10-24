#include <Eigen/Dense>
#include <cmath>
#include <numeric>

const unsigned int factorial(unsigned int L) {
    if(L == 0) {
        return 1;
    }
    else {
        unsigned int l = 1;
        for(unsigned int i = 1; i < L + 1; i++) {
            l *= i;
        }
        return l;
    }
}

template<typename RealType, unsigned int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};


template<typename RealType, std::size_t N, unsigned int L>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) {
    Eigen::Vector<RealType, N> b;
    b.setZero();
    b(L - 1) = RealType(factorial(L));

    Eigen::Matrix<RealType, N, N> A;

    for(unsigned int i = 1; i < N + 1; i++) {
        for(unsigned int j = 0; j < N; j++) {
            A(i - 1, j) = std::pow(points[j], i);
        }
    }

    typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrix;
    typedef Eigen::Vector<RealType, N> AnswArray;

    Eigen::FullPivLU<RealMatrix> decomposition(A);
    AnswArray otherC = decomposition.solve(b);

    std::array<RealType, N> answer;
    
    for(unsigned int i = 0; i < N; i++) {
        answer[i] = otherC[i];
    }

    return {-std::accumulate(answer.begin(), answer.end(), RealType(0)), answer};
}