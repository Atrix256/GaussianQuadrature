#include <array>
#include <stdio.h>

static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_pi = 3.14159265359f;

// Define a vector as an array of floats
template<size_t N>
using TVector = std::array<float, N>;

// Define a matrix as an array of vectors
template<size_t M, size_t N>
using TMatrix = std::array<TVector<N>, M>;

#include "Legendre.h"

template <size_t NUM_SAMPLES>
void DoTests(const TVector<NUM_SAMPLES>& roots)
{
    // Integrating a function f(x) by summing weighted samples.
    // y = w_1 * f(x_1) + ... + w_n * f(x_2)
    // Gaussian Quadrature finds w's and x's that make it so order 0, order 1, ... , order m polynomials are exactly calculatable.
    // If the function is not a polynomial, or is higher order, this becomes an approximation.
    // The above comments about polynomials assume that basis functions such as power or Legendre were used.
    // Other basis functions have different statements about what they can calculate or approximate.

    // https://en.wikipedia.org/wiki/Vandermonde_matrix
    TMatrix<NUM_SAMPLES, NUM_SAMPLES> vandermonde;

    // make the constant values
    for (size_t x = 0; x < NUM_SAMPLES; ++x)
        vandermonde[0][x] = 1.0f;

    // make the values of increasing powers
    for (size_t y = 1; y < NUM_SAMPLES; ++y)
        for (size_t x = 0; x < NUM_SAMPLES; ++x)
            vandermonde[y][x] = vandermonde[y - 1][x] * roots[x];

    int ijkl = 0;
}

template <typename LAMBDA>
void TestFunction(const char* label, const LAMBDA& f, float actualValue, float a = -1.0f, float b = 1.0f)
{
    printf("%s\n", label);

    // change of interval: https://en.wikipedia.org/wiki/Gaussian_quadrature
    auto ChangeOfInterval = [a, b] (float x) -> float
    {
        return x * (b - a) / 2.0f + (a + b) / 2.0f;
    };
    float changeOfIntervalMultiplier = (b - a) / 2.0f;

    // L2
    {
        float value = 0.0f;
        for (size_t index = 0; index < _countof(L2); ++index)
            value += f(ChangeOfInterval(L2[index].x)) * L2[index].weight;
        value *= changeOfIntervalMultiplier;
        printf("  L2: %f (%f)\n", value, abs(value - actualValue));
    }

    // L3
    {
        float value = 0.0f;
        for (size_t index = 0; index < _countof(L3); ++index)
            value += f(ChangeOfInterval(L3[index].x)) * L3[index].weight;
        value *= changeOfIntervalMultiplier;
        printf("  L3: %f (%f)\n", value, abs(value - actualValue));
    }

    // L4
    {
        float value = 0.0f;
        for (size_t index = 0; index < _countof(L4); ++index)
            value += f(ChangeOfInterval(L4[index].x)) * L4[index].weight;
        value *= changeOfIntervalMultiplier;
        printf("  L4: %f (%f)\n", value, abs(value - actualValue));
    }

    // L5
    {
        float value = 0.0f;
        for (size_t index = 0; index < _countof(L5); ++index)
            value += f(ChangeOfInterval(L5[index].x)) * L5[index].weight;
        value *= changeOfIntervalMultiplier;
        printf("  L5: %f (%f)\n", value, abs(value - actualValue));
    }

    printf("\n");
}

int main(int argc, char** argv)
{

    {
        auto f = [](float x) -> float
        {
            return 1.0f;
        };
        TestFunction("y=1", f, 2.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return x;
        };
        TestFunction("y=x", f, 0.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return x * x;
        };
        TestFunction("y=x*x", f, 2.0f / 3.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return x * x * x;
        };
        TestFunction("y=x*x*x", f, 0.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return x * x * x * x;
        };
        TestFunction("y=x*x*x*x", f, 0.4f);
    }

    {
        auto f = [](float x) -> float
        {
            return x * x * x * x * x;
        };
        TestFunction("y=x*x*x*x*x", f, 0.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return x * x * x * x * x * x;
        };
        TestFunction("y=x*x*x*x*x*x", f, 2.0f / 7.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return 5.0f * x *x + 3.0f * x + 2.0f;
        };
        TestFunction("y=5x^2+3x+2", f, 22.0f / 3.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return 4.0f * x * x * x * x - 2.0f * x * x * x + 5.0f * x * x + 3.0f * x + 2.0f;
        };
        TestFunction("y=4x^4-2x^3+5x^2+3x+2", f, 134.0f / 15.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return (float)sin(x);
        };
        TestFunction("y=sin(x)", f, 0.0f);
    }

    {
        auto f = [](float x) -> float
        {
            return (float)(sin(x)*sin(x));
        };
        TestFunction("y=sin(x)*sin(x)", f, (float)(1.0f - sin(2.0f) / 2.0f));
    }

    {
        auto f = [](float x) -> float
        {
            return (float)sin(x);
        };
        TestFunction("y=sin(x) from 0 to pi", f, 2.0f, 0.0f, c_pi);
    }

    // golden Ratio LDS
    {
        TVector<2> roots;
        {
            float root = 0.0f;
            for (float &f : roots)
            {
                root = fmodf(root + c_goldenRatioConjugate, 1.0f);
                f = root;
            }
        }

        //DoTests(roots);
    }

    // Legendre polynomials
    {

    }

    return 0;
}

/*
TODO:
! try to gain better knowledge of why it works, not just how it works
* power basis polynomials
* legendre polynomials
* golden ratio polynomials?
* compare vs monte carlo. how many samples does it take to reach the error level??
* do the thing that lets you change it from being from -1 to +1, to a to b.  change of intervals. 
* could try to calculate the weights of legendre polynomials yourself too

NOTES:
* the weights multiplying the samples in, it's kind of like 1 / (N * pdf(x)). of course, the x's aren't randomly selected...
* the x's and weights are re-usable.

Links:
* all of the related things from here https://www.youtube.com/watch?v=nQZYBWB6q_k

*/
