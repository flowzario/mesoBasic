/*
 * rand.hpp
 * Copyright (C) 2017 joseph <joseph@JMC-WORKSTATION>
 *
 * Distributed under terms of the MIT license.
 *
 *
 * A Mersenne Twister Psuedo Random Number Generator Class
 *
 * Code comes from the wikipedia page:
 *
 * https://en.wikipedia.org/wiki/Mersenne_Twister
 *
 * Class methods are provided for transforming the Mersenne
 * Twister output to uniform and normal distributions. This
 * RNG is being used for its very high periodicity (2^19937)
 *
 */

#ifndef RAND_H
#define RAND_H

#include <stdint.h>

class Rand
{
    private:
        // Define MT19937 constants (32-bit RNG)
        enum
        {
            // Assumes W = 32 (omitting this)
            N = 624,
            M = 397,
            R = 31,
            A = 0x9908B0DF,

            F = 1812433253,

            U = 11,
            // Assumes D = 0xFFFFFFFF (omitting this)

            S = 7,
            B = 0x9D2C5680,

            T = 15,
            C = 0xEFC60000,

            L = 18,

            MASK_LOWER = (1ull << R) - 1,
            MASK_UPPER = (1ull << R)
        };
        uint32_t mt[N];
        uint16_t index;
        void Twist();

    public:
        Rand(uint32_t seed = 0);
        ~Rand();

        /// max value that can be generated by MT generator
        const uint32_t max; 
        
        /// min value that can be generated by MT generator
        const uint32_t min; 

        /// seeds the MT generator
        void init(const uint32_t seed);

        /// gets the next unsigned 32 bit integer from the MT generator
        uint32_t nextU32();

        /// gets a double in the range [0,1] from a uniform distribution
        double uniform();

        /// gets a double from a normal distribution with mean 0 and
        /// unit deviation
        double normal();
};

#endif /* !RAND_H */