/*  TODO: 
 *  
 *  Define Matrix-Vector Multiplication
 *  Define PolyMatrix-Vector Multiplication
 *  Define Matrix-Matrix (and PolyMatrix Variants) Multiplication
 */


#pragma once

#ifndef INCLUDE_SEGFAULT_MATRIX
#define INCLUDE_SEGFAULT_MATRIX

#ifdef USE_POLYNOMIAL
    #include "polynomial"
#else
    #include <cstdint>
    #include <iostream>
    #include <vector>
    #include <minmax.h>
#endif

#include <math.h>


// == No idea why I had to do this. I just did and it worked. I mean... who cares, right?

namespace Tensor {
    class Matrix;

    #ifdef INCLUDE_SEGFAULT_POLYNOMIAL
    class PolyMatrix;
    PolyMatrix operator*(const polynomial& x, const Matrix& mat);
    PolyMatrix operator*(const polynomial& x, PolyMatrix mat);

    PolyMatrix operator-(const Matrix base, const PolyMatrix other);
    #endif
}

std::ostream& operator<<(std::ostream& out, const Tensor::Matrix& mat);

#ifdef INCLUDE_SEGFAULT_POLYNOMIAL
std::ostream& operator<<(std::ostream& out, const Tensor::PolyMatrix& mat);
#endif


// ======================================================================================

namespace Tensor {
    class Matrix {
        private:
        uint32_t width, height;
        std::vector<std::vector<double>> DATA;

        public:
        Matrix(const uint32_t& W, const uint32_t& H);
        Matrix(const std::vector<std::vector<double>>& initList);
        Matrix(const uint32_t& S);
        ~Matrix();
        
        uint32_t rows() const;
        uint32_t columns() const;
        std::vector<std::vector<double>> data() const;
        double determinant();

        void set(uint32_t row, uint32_t column, double data);
        
        std::vector<double> operator[](uint32_t ROW);

        Matrix operator+(const Matrix& other);
        Matrix operator-(const Matrix& other);

        Matrix operator*(const double& other);
        Matrix operator/(const double& other);

        std::vector<double> operator*(const std::vector<double>& other);

        #ifdef INCLUDE_SEGFAULT_VECTOR
        PolyVector operator*(const Tensor::PolyVector& vec);
        #endif

        friend std::ostream& ::operator<<(std::ostream& out, const Tensor::Matrix& mat);
        friend std::vector<double> GaussSolution(Matrix matrix, std::vector<double> vector);

        #ifdef INCLUDE_SEGFAULT_POLYNOMIAL
        std::vector<double> eigenvalues(const double& start, const double& stop, const int& terms);
        #endif
    };

    // ====================================================================================

    std::vector<double> GaussSolution(Matrix mat, std::vector<double> vec);

    // == Helper Functions ================================================================

    bool equalFloat(const double& a, const double& b, const double& dx=pow(10,-10));

    template<typename _Ty>
    std::vector<_Ty> operator*(std::vector<_Ty> base, const double& multiplier);

    template<typename _Ty>
    std::vector<_Ty> operator*(const double& multiplier, std::vector<_Ty> base);

    
    // == Subtracts two vectors. Data must be the same size. 
    template<typename _Ty>
    std::vector<_Ty> operator-(std::vector<_Ty> base, const std::vector<_Ty>& other);

};

#ifdef INCLUDE_SEGFAULT_POLYNOMIAL

namespace Tensor {
    class PolyMatrix {
        private:
        uint32_t width, height;
        std::vector<std::vector<polynomial>> DATA;

        public:
        PolyMatrix(const uint32_t& W, const uint32_t& H);
        PolyMatrix(const std::vector<std::vector<polynomial>>& initList);
        PolyMatrix(const uint32_t& S);
        PolyMatrix(const Matrix& initializer);
        ~PolyMatrix();

        uint32_t rows() const;
        uint32_t columns() const;
        std::vector<std::vector<polynomial>> data() const;
        polynomial determinant();

        void set(uint32_t row, uint32_t column, polynomial data);

        std::vector<polynomial> operator[](uint32_t ROW);

        PolyMatrix operator+(const PolyMatrix& other);
        PolyMatrix operator-(const PolyMatrix& other);

        PolyMatrix operator*(const double& other);
        PolyMatrix operator/(const double& other);

        #ifdef INCLUDE_SEGFAULT_VECTOR
        PolyVector operator*(const std::vector<double>& other);
        PolyVector operator*(const Tensor::PolyVector& vec);
        #else
        template<typename _Ty> auto operator*(const _Ty& vec) {
            try { return (*this * vec); }
            catch (auto e) {
            throw std::exception("Tried to multiply PolyMatrix to another object when\
            PolyVector is not included. Include polyvector.hpp before including matrix.hpp.")
            return _Ty; }
        }
        #endif

        friend std::ostream& ::operator<<(std::ostream& out, const Tensor::PolyMatrix& mat);
        friend PolyMatrix operator*(const polynomial& x, PolyMatrix mat);
    };

    PolyMatrix operator*(const polynomial& x, PolyMatrix mat);

    PolyMatrix operator*(const PolyMatrix& mat, const polynomial& x);

    PolyMatrix operator*(const polynomial& x, const Matrix& mat);

    PolyMatrix operator*(const Matrix& mat, const polynomial& x);
}
#endif

#endif