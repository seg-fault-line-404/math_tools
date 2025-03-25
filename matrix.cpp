#include <iomanip>
#include <stdexcept>

// Undefine if you only need the discrete-value matrix parts 
#define USE_POLYNOMIAL
#include <differential>
#include "matrix.hpp"


namespace Tensor {
    Matrix::Matrix(const uint32_t& W, const uint32_t& H): width(W), height(H) {
        DATA.resize(height);
        for (auto& row : DATA) {
            row.resize(width);
        }
    }

    Matrix::Matrix(const std::vector<std::vector<double>>& initList){
        height = initList.size();
        width = initList[0].size();
        for (const auto& row : initList) {
            if (row.size() != width) {
                throw std::domain_error("Rows are not the same length.");
                *this = Matrix(0,0);
                return;
            }
        }
        DATA = initList;
    }

    Matrix::Matrix(const uint32_t& S){
        Matrix base(S, S);
        for (int i = 0; i < S; ++i) {
            base.DATA[i][i] = 1.;
        }
        *this = base;
        
    }

    Matrix::~Matrix() {}

    uint32_t Matrix::rows() const { return height; }

    uint32_t Matrix::columns() const { return width; }

    std::vector<std::vector<double>> Matrix::data() const { return DATA; }

    double Matrix::determinant() {

        if (width != height) {
            throw std::domain_error("Must be a square matrix.");
            return 0.0;
        }

        if (width == 2) {
            return (DATA[0][0] * DATA[1][1] - DATA[0][1] * DATA[1][0]);
        }

        double sum = 0.0;

        for (int i = 0; i < width; ++i) {
            std::vector<std::vector<double>> sublist;
            sublist.reserve(height - 1);

            for (int j = 1; j < height; ++j) {

                std::vector<double> temp;
                temp.reserve(width - 1);
                for (int k = 0; k < width; ++k) {
                    if (k == i) continue;
                    temp.push_back(DATA[j][k]);
                }
                sublist.push_back(temp);
            }

            Matrix Submatrix(sublist);

            sum += (pow(-1, i) * DATA[0][i] * Submatrix.determinant());
        }

        return sum;
    }

    void Matrix::set(uint32_t row, uint32_t column, double data) {
        DATA[row][column] = data;
    }

    Matrix Matrix::operator+(const Matrix& other){
        if (width != other.width || height != other.height) {
            throw std::domain_error("Matrices must be the same size.\n");
            return *this;
        }

        Matrix base(*this);

        for (int i = 0; i < base.width; ++i) {              //  for i in range(0, base.width):
            for (int j = 0; j < base.height; ++j) {         //      for j in range(0, base.height):
                base.DATA[i][j] += other.DATA[i][j];        //          base.DATA[i][j] += other.DATA[i][j]
            }
        }

        return base;                                        //  return base
    }

    Matrix Matrix::operator-(const Matrix& other){
        if (width != other.width || height != other.height) {
            throw std::domain_error("Matrices must be the same size.\n");
            return *this;
        }

        Matrix base(*this);

        for (int i = 0; i < base.width; ++i) {
            for (int j = 0; j < base.height; ++j) {
                base.DATA[i][j] -= other.DATA[i][j];
            }
        }

        return base;
    }

    std::vector<double> Matrix::operator[](uint32_t ROW) {
        return DATA[ROW];
    }

    #ifdef INCLUDE_SEGFAULT_POLYNOMIAL
    std::vector<double> Matrix::eigenvalues(const double& start, const double& stop, const int& terms) {
        if (width != height) {
            throw std::domain_error("Must be a square matrix.");
            return {0.0};
        }

        polynomial ONE({1.0, 0.0});
        PolyMatrix I(width);

        for (int i = 0; i < width; ++i) I.set(i, i, ONE);

        auto mat = *this - I;
        auto eqn = mat.determinant();

        auto interval = linspace(start, stop, terms);
        auto outputspace = interval;

        for (auto& i : outputspace) {
            i = eqn.eval(i);
        }

        std::vector<std::array<double, 2>> POI;
        for (int i = 1; i < interval.size(); ++i) {
            if (outputspace[i] * outputspace[i - 1] <= 0) {
                POI.push_back({interval[i - 1], interval[i]});
            }
        }

        std::vector<double> roots;
        roots.reserve(POI.size());
        for (const auto& pair : POI) {
            roots.push_back(roots::bisectRoot(eqn, pair, 50));
        }

        return roots;
    }
    #endif


    std::vector<double> GaussSolution(Matrix mat, std::vector<double> vec) {
        if (mat.height != vec.size()) {
            throw std::domain_error("Equation results to non-unique solutions.\
                                    Matrix and vector must be the same height.");
            return {0.0};
        }

        if (mat.width != mat.height) {
            throw std::domain_error("Must be a square matrix.");
            return {0.0};
        }

        // == Make the upper triangle using magic =====
        for (uint32_t i = 0; i < vec.size(); ++i) {
            for (uint32_t j = i + 1; j < vec.size(); ++j) {
                if (!equalFloat(mat.DATA[j][i], 0.0)) {
                    const double multi = mat.DATA[j][i]/mat.DATA[i][i];
                    mat.DATA[j] = (mat.DATA[j] - (multi * mat.DATA[i]));
                    vec[j] -= (multi * vec[i]);
                }
            }
        }

        uint32_t max = vec.size()-1;
        vec[max] /= mat.DATA[max][max];

        for (int i = max - 1; i >=0; --i) {
            for (int j = i + 1; j <= max; ++j) {
                vec[i] -= (mat[i][j] * vec[j]);
            }
            vec[i] /= mat.DATA[i][i];
        }

        return vec;
    }

    // == Helper Functions ==============================

    bool equalFloat(const double& a, const double& b, const double& dx) {
        return abs(a-b) < dx;
    }

    template<typename _Ty>
    std::vector<_Ty> operator*(std::vector<_Ty> base, const double& multiplier) {
        for (_Ty& data : base) data *= multiplier;
        return base;
    }

    template<typename _Ty>
    std::vector<_Ty> operator*(const double& multiplier, std::vector<_Ty> base) {
        for (_Ty& data : base) data *= multiplier;
        return base;
    }

    
    // == Subtracts two vectors. Data must be the same size. ==================
    template<typename _Ty>
    std::vector<_Ty> operator-(std::vector<_Ty> base, const std::vector<_Ty>& other) {
        if (base.size() != other.size())  {
            throw std::domain_error("Vectors must be the same size.");
            return {0.0};
        }

        for (int i = 0; i < base.size(); ++i) base[i] -= other[i];
        
        return base;
    }
}


std::ostream& operator<<(std::ostream& out, const Tensor::Matrix& mat) {
    out << '\n';
    for (const auto& row : mat.DATA) {
        for (const auto& data : row) {
            out << data << '\t';
        }
        out << '\n';
    }
    return out;
}

#ifdef INCLUDE_SEGFAULT_POLYNOMIAL

namespace Tensor {
    PolyMatrix::PolyMatrix(const uint32_t& W, const uint32_t& H): width(W), height(H) {
        DATA.resize(height);
        for (auto& row : DATA) {
            row.resize(width);
        }
    }

    PolyMatrix::PolyMatrix(const std::vector<std::vector<polynomial>>& initList){
        height = initList.size();
        width = initList[0].size();
        for (const auto& row : initList) {
            if (row.size() != width) {
                throw std::domain_error("Rows are not the same length.");
                *this = PolyMatrix(0,0);
                return;
            }
        }
        DATA = initList;
    }

    PolyMatrix::PolyMatrix(const uint32_t& S): width(S), height(S) {
        DATA.resize(S);

        for (auto& row : DATA) row.resize(S);

        polynomial a({1.});
        
        for (int i = 0; i < S; ++i) {
            DATA[i][i] = a;
        }
    }

    PolyMatrix::PolyMatrix(const Matrix& initializer) {
        width = initializer.columns();
        height = initializer.rows();
        DATA.reserve(height);
        
        for (auto& row : initializer.data()) {
            std::vector<polynomial> temp;
            temp.reserve(row.size());
            for (const auto& data : row) {
                polynomial D({data});
                temp.push_back(D);
            }
            DATA.push_back(temp);
        }
    }

    PolyMatrix::~PolyMatrix() {}

    uint32_t PolyMatrix::rows() const { return height; }
    uint32_t PolyMatrix::columns() const { return width; }
    std::vector<std::vector<polynomial>> PolyMatrix::data() const { return DATA; }

    polynomial PolyMatrix::determinant() {

        if (width != height) {
            throw std::domain_error("Must be a square matrix.");
            polynomial fail({0.});
            return fail;
        }

        if (width == 2) {
            return (DATA[0][0] * DATA[1][1] - DATA[0][1] * DATA[1][0]);
        }

        polynomial sum({0.});

        for (int i = 0; i < width; ++i) {
            std::vector<std::vector<polynomial>> sublist;
            sublist.reserve(height - 1);

            for (int j = 1; j < height; ++j) {

                std::vector<polynomial> temp;
                temp.reserve(width - 1);
                for (int k = 0; k < width; ++k) {
                    if (k == i)  continue;
                    temp.push_back(DATA[j][k]);
                }
                sublist.push_back(temp);
            }

            PolyMatrix Submatrix(sublist);

            // std::cout << signof(pow(-1, i)) << DATA[0][i] << Submatrix << std::endl;

            sum += (pow(-1, i) * DATA[0][i] * Submatrix.determinant());
        }

        return sum;
    }

    void PolyMatrix::set(uint32_t row, uint32_t column, polynomial data) {
        DATA[row][column] = data;
    }

    PolyMatrix PolyMatrix::operator+(const PolyMatrix& other) {
        if (width != other.width || height != other.height) {
            PolyMatrix fail(1);
            throw std::domain_error("Matrices must be the same size.");
            return fail;
        }

        PolyMatrix sum(*this);

        for (uint32_t i = 0; i < height; ++i) {
            for (uint32_t j = 0; j < width; ++j) {
                sum.DATA[i][j] += other.DATA[i][j];
            }
        }

        return sum;
    }

    PolyMatrix PolyMatrix::operator-(const PolyMatrix& other) {
        if (width != other.width || height != other.height) {
            PolyMatrix fail(1);
            throw std::domain_error("Matrices must be the same size.");
            return fail;
        }

        PolyMatrix sum(*this);

        for (uint32_t i = 0; i < height; ++i) {
            for (uint32_t j = 0; j < width; ++j) {
                sum.DATA[i][j] -= other.DATA[i][j];
            }
        }

        return sum;
    }
    
    std::vector<polynomial> PolyMatrix::operator[] (uint32_t ROW) {
        return DATA[ROW];
    }

    PolyMatrix PolyMatrix::operator*(const double& other) {
        PolyMatrix base(*this);
        for (auto& row : base.DATA) {
            for (auto& data : row) {
                data *= other;
            }
        }

        return base;
    }

    PolyMatrix PolyMatrix::operator/(const double& other) {
        PolyMatrix base(*this);
        for (auto& row : base.DATA) {
            for (auto& data : row) {
                data /= other;
            }
        }

        return base;
    }

    PolyMatrix operator*(const polynomial& x, PolyMatrix mat) {
        std::vector<std::vector<polynomial>> base = mat.data();
        for (auto& row : base) {
            for (auto& data : row) {
                data *= x;
            }
        }

        return PolyMatrix(base);
    }

    PolyMatrix operator*(const PolyMatrix& mat, const polynomial& x) {
        return x * mat;
    }

    PolyMatrix operator*(const polynomial& x, const Matrix& mat) {
        PolyMatrix newMat(mat);
        return x * newMat;
    }

    PolyMatrix operator*(const Matrix& mat, const polynomial& x) {
        PolyMatrix newMat(mat);
        return x * newMat;
    }

    PolyMatrix operator+(const Matrix& base, PolyMatrix& other) {
        if (base.rows() != other.rows() || base.columns() != other.columns()) {
            throw std::domain_error("Matrices must be the same size.");
            return PolyMatrix(1);
        }

        std::vector<std::vector<polynomial>> result;
        result.resize(base.rows());

        for (uint32_t i = 0; i < base.rows(); ++i) {
            result[i].resize(base.columns());

            for (uint32_t j = 0; j < base.columns(); ++j) {
                result[i][j] = base.data()[i][j] + other.data()[i][j];
            }
        }

        PolyMatrix sum(result);
        return sum;
    }

    PolyMatrix operator+(const PolyMatrix& base, const Matrix& other) { 
        if (base.rows() != other.rows() || base.columns() != other.columns()) {
            throw std::domain_error("Matrices must be the same size.");
            return PolyMatrix(1);
        }

        std::vector<std::vector<polynomial>> result;
        result.resize(base.rows());

        for (uint32_t i = 0; i < base.rows(); ++i) {
            result[i].resize(base.columns());

            for (uint32_t j = 0; j < base.columns(); ++j) {
                result[i][j] = base.data()[i][j] + other.data()[i][j];
            }
        }

        PolyMatrix sum(result);
        return sum;
    }

    PolyMatrix operator-(const Matrix base, const PolyMatrix other) {
    if (base.rows() != other.rows() || base.columns() != other.columns()) {
        throw std::domain_error("Matrices must be the same size.");
        return PolyMatrix(1);
    }

    std::vector<std::vector<polynomial>> result;
    result.resize(base.rows());

    for (uint32_t i = 0; i < base.rows(); ++i) {
        result[i].resize(base.columns());

        for (uint32_t j = 0; j < base.columns(); ++j) {
            result[i][j] = base.data()[i][j] - other.data()[i][j];
        }
    }

    PolyMatrix sum(result);
    return sum;
}


    PolyMatrix operator-(const PolyMatrix& base, const Matrix& other) {
        if (base.rows() != other.rows() || base.columns() != other.columns()) {
            throw std::domain_error("Matrices must be the same size.");
            return PolyMatrix(1);
        }

        std::vector<std::vector<polynomial>> result;
        result.resize(base.rows());

        for (uint32_t i = 0; i < base.rows(); ++i) {
            result[i].resize(base.columns());

            for (uint32_t j = 0; j < base.columns(); ++j) {
                result[i][j] = base.data()[i][j] - other.data()[i][j];
            }
        }

        PolyMatrix sum(result);
        return sum;
    }
}

std::ostream& operator<<(std::ostream& out, const Tensor::PolyMatrix& mat) {
    out << '\n';
    for (const auto& row : mat.DATA) {
        for (const auto& data : row) {
            out << data << '\t';
        }
        out << '\n';
    }
    return out;
}
#endif