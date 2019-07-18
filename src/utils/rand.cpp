/*
 * rand.cpp
 * Copyright (C) 2017 joseph <joseph@JMC-WORKSTATION>
 *
 * Distributed under terms of the MIT license.
 */

#include "rand.hpp"
#include <cmath>


Rand::Rand(uint32_t seed)
    : max(4294967295),min(0)
{
    init(seed);
}

Rand::~Rand()
{
}


// Re-init with a given seed
void Rand::init(const uint32_t  seed)
{
    uint32_t  i;

    mt[0] = seed;

    for ( i = 1; i < N; i++ )
    {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
    }

    index = N;
}

void Rand::Twist()
{
    uint32_t  i, x, xA;

    for ( i = 0; i < N; i++ )
    {
        x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);

        xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
    }

    index = 0;
}

// Obtain a 32-bit random number
uint32_t Rand::nextU32()
{
    uint32_t  y;
    int       i = index;

    if ( index >= N )
    {
        Twist();
        i = index;
    }

    y = mt[i];
    index = i + 1;

    y ^= (y >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}


double Rand::uniform()
{
    return (double)nextU32()/(double)max;
}


double Rand::normal()
{
    // use the Box Muller transform
    double u1 = uniform();
    double u2 = uniform();
    double norm = sqrt(-2.0*log(u1))*cos(2.0*3.14159265359*u2);
    return norm;
}
