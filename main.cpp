#include <array>
#include <stdio.h>
#include <vector>

static const float c_goldenRatioConjugate = 0.61803398875f;
static const float c_pi = 3.14159265359f;

// Define a vector as an array of floats
template<size_t N>
using TVector = std::array<float, N>;

// Define a matrix as an array of vectors
template<size_t M, size_t N>
using TMatrix = std::array<TVector<N>, M>;

struct WeightedSample
{
    float x;
    float weight;
};

struct WeightedSamples
{
    const char* label = nullptr;
    std::vector<WeightedSample> weightedSamples;
};

#include "Legendre.h"

std::vector<WeightedSamples> g_weightedSamples;

template <size_t NUM_SAMPLES>
void CalculateWeightedSamples(const TVector<NUM_SAMPLES>& roots, std::vector<WeightedSample>& weightedSamples)
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

    for (const WeightedSamples& weightedSamples : g_weightedSamples)
    {
        {
            float value = 0.0f;
            for (const WeightedSample& weightedSample : weightedSamples.weightedSamples)
                value += f(ChangeOfInterval(weightedSample.x)) * weightedSample.weight;
            value *= changeOfIntervalMultiplier;
            printf("  %s: %f (%f)\n", weightedSamples.label, value, abs(value - actualValue));
        }
    }

    printf("\n");
}

int main(int argc, char** argv)
{

    g_weightedSamples.push_back(L2);
    g_weightedSamples.push_back(L3);
    g_weightedSamples.push_back(L4);
    g_weightedSamples.push_back(L5);

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
* probably need to make an interface where you give a list of sample locations and weights and it will do the test. this works for Legendre and other things too. also makes it clearer what is being generated by eg golden ratio LDS 
! try to gain better knowledge of why it works, not just how it works
* power basis polynomials
* legendre polynomials
* golden ratio polynomials? might want symmetry on +1 and -1 though? unsure.
* compare vs monte carlo. how many samples does it take to reach the error level??
* do the thing that lets you change it from being from -1 to +1, to a to b.  change of intervals. 
* could try to calculate the weights of legendre polynomials yourself too

NOTES:
* the weights multiplying the samples in, it's kind of like 1 / (N * pdf(x)). of course, the x's aren't randomly selected...
* the x's and weights are re-usable.

Links:
* all of the related things from here https://www.youtube.com/watch?v=nQZYBWB6q_k

*/
