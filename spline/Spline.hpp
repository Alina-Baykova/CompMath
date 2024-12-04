#include <vector>
#include <iostream>
#include <type_traits>
#include <cassert>
#include <algorithm>


/** класс для работы с трехдиагональной матрицей **/
template<typename xType>
class ThreeDiagonalMatrix {
    unsigned int size;
    std::vector<xType> mainDiag;
    std::vector<xType> lowerDiag;
    std::vector<xType> upperDiag;

    public:
    ThreeDiagonalMatrix(const std::vector<xType>& mDiag, const std::vector<xType>& lDiag, const std::vector<xType>& uDiag) noexcept {
        assert((mDiag.size() == lDiag.size() + 1) && (mDiag.size() == uDiag.size() + 1));
        size = mDiag.size();
        mainDiag.resize(size);
        lowerDiag.resize(size - 1);
        upperDiag.resize(size - 1);
        mainDiag = mDiag;
        lowerDiag = lDiag;
        upperDiag = uDiag;
    };

    const std::vector<xType>& getMainDiag() const {
        return mainDiag;
    };

    const std::vector<xType>& getLowerDiag() const {
        return lowerDiag;
    };

    const std::vector<xType>& getUpperDiag() const {
        return upperDiag;
    };

    const unsigned int& getSize() const {
        return size;
    }
};


template<typename xType>
std::ostream& operator << (std::ostream& os, const ThreeDiagonalMatrix<xType>& matrix)
{
    unsigned int N = matrix.getSize();
    auto md = matrix.getMainDiag();
    auto ld = matrix.getLowerDiag();
    auto ud = matrix.getUpperDiag();
    for(unsigned int i = 0; i < N; i++) {
        for(unsigned int j = 0; j < N; j++) {
            if(i == j) {
                os << md[i] << "\t";
                continue;
            }
            else if(j == i + 1) {
                os << ud[i] << "\t";
                continue;
            }
            else if(j == i - 1) {
                os << ld[j] << "\t";
                continue;
            }
            else {
                os << xType(0) << "\t";
                continue;
            }
        }
        os << "\n";
    }
    return os;
}


template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

/** Функция для решения методом  прогонки **/

template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const ThreeDiagonalMatrix<mType>& matrix, const std::vector<cType>& column){
    unsigned int size = column.size();
    assert(size > 1);
    auto c = matrix.getUpperDiag();
    auto b = matrix.getMainDiag();
    auto a = matrix.getLowerDiag();
    auto d = column;
    unsigned int M = size - 1;

    std::vector<DivisType<cType, mType>> answer;
    std::vector<DivisType<cType, mType>> p;
    std::vector<DivisType<cType, mType>> q;
    
    answer.resize(size);
    p.resize(M);
    q.resize(M);

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];

    for(unsigned int i = 1; i < M; i++) {
        p[i] = (-c[i]) / (a[i-1] * p[i - 1] + b[i]);
        q[i] = (d[i] - a[i - 1] * q[i - 1]) / (a[i - 1] * p[i - 1] + b[i]);
    }

    answer.back() = (d[M] - a[M - 1] * q[M - 1]) / (a[M - 1] * p[M - 1] + b[M]);

    for(int i = M - 1; i > -1; i--) {
        answer[i] = answer[i + 1] * p[i] + q[i];
    }

    return answer;
}


template<typename xType, typename yType>
class CubicSpline {
    using DeltaXType = DiffType<xType>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;
    using Deriv3Type = DivisType<DiffType<Deriv2Type>, DeltaXType>;

    unsigned int N;
    
    std::vector<xType> x;
    std::vector<yType> u;
    std::vector<DeltaXType> h;
    std::vector<DivisType<yType, DeltaXType>> splitDiffs1;
    std::vector<DivisType<DivisType<yType, DeltaXType>, DeltaXType>> splitDiffs2MultBySix;
    std::vector<DivisType<DeltaXType, DeltaXType>> mainDiag;
    std::vector<DivisType<DeltaXType, DeltaXType>> lowerDiag;
    std::vector<DivisType<DeltaXType, DeltaXType>> upperDiag;
    std::vector<DerivType> bCoeffs;
    std::vector<Deriv2Type> cCoeffs;
    std::vector<Deriv3Type> dCoeffs;
    
    Deriv2Type leftSecondDeriv;
    Deriv2Type rightSecondDeriv;

    public:
    CubicSpline(const std::vector<xType>& points, const std::vector<yType>& values, const Deriv2Type& first, const Deriv2Type& second){
        N = points.size();

        x.resize(N);
        u.resize(N);
        h.resize(N - 1);
        mainDiag.resize(N - 2);
        lowerDiag.resize(N - 3);
        upperDiag.resize(N - 3);
        splitDiffs1.resize(N - 1);
        splitDiffs2MultBySix.resize(N - 2);
        bCoeffs.resize(N - 1);
        cCoeffs.resize(N - 2);
        dCoeffs.resize(N - 1);

        x = points;
        u = values;

        for(unsigned int i = 0; i < N - 1; i++){
            h[i] = x[i + 1] - x[i];
            splitDiffs1[i] = (u[i + 1] - u[i]) / h[i];
        }

        for(unsigned int i = 0; i < N - 2; i++){
            splitDiffs2MultBySix[i] = 6 * (splitDiffs1[i + 1] - splitDiffs1[i]) / (h[i + 1] + h[i]);
        }

        for(unsigned int i = 0; i < N - 3; i++){
            mainDiag[i] = xType(2);
            lowerDiag[i] = h[i + 1] / (h[i + 2] + h[i + 1]);
            upperDiag[i] = h[i + 1] / (h[i + 1] + h[i]);
        }

        mainDiag[N - 3] = xType(2);

        ThreeDiagonalMatrix<Deriv2Type> splineThreeDiagMatrix(mainDiag, lowerDiag, upperDiag);
        cCoeffs = solve(splineThreeDiagMatrix, splitDiffs2MultBySix);
        cCoeffs.resize(N - 1);
        cCoeffs[N - 2] = 0;
        dCoeffs[0] = cCoeffs[0] / h[0];
        bCoeffs[0] = cCoeffs[0] * h[0] / 3 + splitDiffs1[0];

        for(unsigned int i = 1; i < N - 1; i++) {
            dCoeffs[i] = (cCoeffs[i] - cCoeffs[i - 1]) / h[i];
            bCoeffs[i] = cCoeffs[i] * h[i] / 3 + cCoeffs[i - 1] * h[i] / 6 + splitDiffs1[i];
        }
    }

    const auto& getPoints() const{
        return x;
    }

    const auto& getValues() const{
        return u;
    }

    const auto& getSteps() const{
        return h;
    }

    const auto& getSplitDiffs1() const{
        return splitDiffs1;
    }

    const auto& getSplitDiffs2MultBySix() const{
        return splitDiffs2MultBySix;
    }

    const auto& getBCoeffs() const{
        return bCoeffs;
    }

    const auto& getCCoeffs() const{
        return cCoeffs;
    }

    const auto& getDCoeffs() const{
        return dCoeffs;
    }

    yType interpolate(const xType& X) const noexcept {
        if(X <= x[0]) {
            return u[1] + bCoeffs[0] * (X - x[1]) + cCoeffs[0] * std::pow((X - x[1]), 2) / 2 + dCoeffs[0] * std::pow((X - x[1]), 3) / 6;
        }
        else if(X >= x.back()) {
            return u.back() + bCoeffs.back() * (X - x.back()) + cCoeffs.back() * std::pow((X - x.back()), 2) / 2 + dCoeffs.back() * std::pow((X - x.back()), 3) / 6;
        }
        else {
            auto iter = std::lower_bound(std::begin(x), std::end(x), X, std::less<xType>());
            auto indexOfInterval = std::distance(x.begin(), iter);
            return u[indexOfInterval] + bCoeffs[indexOfInterval - 1] * (X - x[indexOfInterval]) + cCoeffs[indexOfInterval - 1] * std::pow((X - x[indexOfInterval]), 2) / 2 + dCoeffs[indexOfInterval - 1] * std::pow((X - x[indexOfInterval]), 3) / 6;
        }
    }
};
