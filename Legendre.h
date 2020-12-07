#pragma once

// Zeroes and weights of Legendre polynomials.
// From https://www.ams.org/journals/bull/1942-48-10/S0002-9904-1942-07771-8/S0002-9904-1942-07771-8.pdf
// That document goes up to n=16

struct LegendreSample
{
    float x;
    float weight;
};

// 2 sample points - up to order 3 (cubic)
LegendreSample L2[] =
{
    {-0.577350269189626f, 1.000000000000000f},
    {0.577350269189626f, 1.000000000000000f}
};

// 3 sample points - up to order 5
LegendreSample L3[] =
{
    {-0.774596669241483f, 0.555555555555556f},
    {0.000000000000000f, 0.888888888888889f},
    {0.774596669241483f, 0.555555555555556f}
};

// 4 sample points - up to order 7
LegendreSample L4[] =
{
    {-0.861136311594053f, 0.347854845137454f},
    {-0.339981043584856f, 0.652145154862546f},
    {0.339981043584856f, 0.652145154862546f},
    {0.861136311594053f, 0.347854845137454f}
};

// 5 sample points - up to order 9
LegendreSample L5[] =
{
    {-0.906179845938664f, 0.236926885056189f},
    {-0.538469310105683f, 0.478628670499366f},
    {0.000000000000000f, 0.568888888888889f},
    {0.538469310105683f, 0.478628670499366f},
    {0.906179845938664f, 0.236926885056189f}
};
