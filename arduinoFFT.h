#ifndef arduinoFFT_h
#define arduinoFFT_h

#include <math.h>

#define FFT_FORWARD 0
#define FFT_REVERSE 1
#define FFT_WIN_TYP_HAMMING 1

class ArduinoFFT {
public:
    void Windowing(float *vData, uint16_t samples, uint8_t winType, uint8_t dir) {
        for (uint16_t i = 0; i < samples; i++) {
            vData[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (samples - 1)); // Hamming
        }
    }

    void Compute(float *vReal, float *vImag, uint16_t samples, uint8_t dir) {
        const float pi = M_PI;
        uint16_t n = samples;
        uint16_t half = n / 2;
        uint16_t j = 0;

        // Bit-reversal permutation
        for (uint16_t i = 1; i < n - 1; i++) {
            uint16_t bit = half;
            while (j >= bit) {
                j -= bit;
                bit >>= 1;
            }
            j += bit;
            if (i < j) {
                float temp = vReal[i];
                vReal[i] = vReal[j];
                vReal[j] = temp;

                temp = vImag[i];
                vImag[i] = vImag[j];
                vImag[j] = temp;
            }
        }

        // FFT
        for (uint16_t len = 2; len <= n; len <<= 1) {
            float angle = (dir == FFT_FORWARD ? -1 : 1) * 2 * pi / len;
            float wlenReal = cos(angle);
            float wlenImag = sin(angle);

            for (uint16_t i = 0; i < n; i += len) {
                float wReal = 1;
                float wImag = 0;
                for (uint16_t j = 0; j < len / 2; j++) {
                    uint16_t u = i + j;
                    uint16_t v = i + j + len / 2;
                    float tReal = wReal * vReal[v] - wImag * vImag[v];
                    float tImag = wReal * vImag[v] + wImag * vReal[v];

                    vReal[v] = vReal[u] - tReal;
                    vImag[v] = vImag[u] - tImag;
                    vReal[u] += tReal;
                    vImag[u] += tImag;

                    float tempReal = wReal * wlenReal - wImag * wlenImag;
                    wImag = wReal * wlenImag + wImag * wlenReal;
                    wReal = tempReal;
                }
            }
        }
    }

    void ComplexToMagnitude(float *vReal, float *vImag, uint16_t samples) {
        for (uint16_t i = 0; i < samples; i++) {
            vReal[i] = sqrt(vReal[i] * vReal[i] + vImag[i] * vImag[i]);
        }
    }
};

#endif
