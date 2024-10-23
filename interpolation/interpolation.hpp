#include <array>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
    std::array<yType, N> splitDifferencies;
    std::array<xType, N> points;

    public:
    NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept: splitDifferencies(values), points(points){
        for(unsigned int i = 0; i < N - 1 ; i++){
            for(int j = N - 2 - i; j >= 0; j--){
                yType splitDifferenciesTemp = (splitDifferencies[j + i + 1] - splitDifferencies[j + i])/(points[j + i + 1] - points[j]);
                splitDifferencies[j + i + 1] = splitDifferenciesTemp;
            }
        }
    }

    const std::array<xType, N>& getPoints() const{
        return points;
    }

    const std::array<yType, N>& getsplitDifferencies() const{
        return splitDifferencies;
    }

    yType interpolate(const xType& x) const noexcept{
        yType f = splitDifferencies.back();
        for(int i = N - 1; i > 0;  i--){
            f = splitDifferencies[i - 1] + (x - points[i - 1]) * f;
        }
        return f;
    }
};
