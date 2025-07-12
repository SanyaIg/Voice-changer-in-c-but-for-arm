#include <iostream>
#include <vector>
#include <cmath>
#include <portaudio.h>
#include <string>
#include <sndfile.h>
#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 1024
#define PITCH_FACTOR 1.2f
#define GAIN 1.0f  // Lowered to avoid clipping

void pitchShift(const float* input, float* output, unsigned long frames, float factor) {
    for (unsigned long i = 0; i < frames; ++i) {
        float idx = i * factor;
        unsigned long idx0 = (unsigned long)idx;
        unsigned long idx1 = std::min(idx0 + 1, frames - 1);
        float frac = idx - idx0;
        output[i] = (1 - frac) * input[idx0] + frac * input[idx1];
    }
}

void applyWindow(float* data, unsigned long frames) {
    for (unsigned long i = 0; i < frames; ++i)
        data[i] *= 0.5f * (1 - cosf(2 * M_PI * i / (frames - 1)));
}

int findDeviceByName(const std::string& name, bool output) {
    int numDevices = Pa_GetDeviceCount();
    for (int i = 0; i < numDevices; ++i) {
        const PaDeviceInfo* info = Pa_GetDeviceInfo(i);
        std::string deviceName(info->name);
        if (deviceName.find(name) != std::string::npos) {
            if ((output && info->maxOutputChannels > 0) ||
                (!output && info->maxInputChannels > 0)) {
                return i;
            }
        }
    }
    return -1;
}

// Globals for file writing
SNDFILE* sndfile = nullptr;
SF_INFO sfinfo;

static int paCallback(const void* inputBuffer, void* outputBuffer,
                      unsigned long framesPerBuffer,
                      const PaStreamCallbackTimeInfo*,
                      PaStreamCallbackFlags, void*) {
    const float* in = (const float*)inputBuffer;
    float* out = (float*)outputBuffer;
    static float temp[FRAMES_PER_BUFFER];
    // Pitch shift and window
    pitchShift(in, temp, framesPerBuffer, PITCH_FACTOR);
    applyWindow(temp, framesPerBuffer);
    // Apply gain, clip, and output stereo (duplicate mono to both channels)
    for (unsigned long i = 0; i < framesPerBuffer; ++i) {
        float sample = temp[i] * GAIN;
        if (sample > 1.0f) sample = 1.0f;
        if (sample < -1.0f) sample = -1.0f;
        out[2*i] = sample;     // left
        out[2*i+1] = sample;   // right
    }
    // Write to WAV file (stereo)
    if (sndfile) {
        std::vector<float> stereoBuf(framesPerBuffer * 2);
        for (unsigned long i = 0; i < framesPerBuffer; ++i) {
            stereoBuf[2*i] = out[2*i];
            stereoBuf[2*i+1] = out[2*i+1];
        }
        sf_writef_float(sndfile, stereoBuf.data(), framesPerBuffer);
    }
    return paContinue;
}

int main() {
    Pa_Initialize();

    // Find BlackHole 2ch output device
    std::string virtualDeviceName = "BlackHole 2ch";
    int outputDevice = findDeviceByName(virtualDeviceName, true);
    int inputDevice = Pa_GetDefaultInputDevice();

    if (outputDevice < 0) {
        std::cerr << "Could not find BlackHole 2ch output device!" << std::endl;
        Pa_Terminate();
        return 1;
    }

    // Setup WAV file for stereo output
    sfinfo.samplerate = SAMPLE_RATE;
    sfinfo.channels = 2;
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
   sndfile = sf_open("/tmp/output.wav", SFM_WRITE, &sfinfo);
    if (!sndfile) {
        std::cerr << "Could not open output.wav for writing!" << std::endl;
        Pa_Terminate();
        return 1;
    }

    PaStreamParameters inputParams, outputParams;
    inputParams.device = inputDevice;
    inputParams.channelCount = 1;
    inputParams.sampleFormat = paFloat32;
    inputParams.suggestedLatency = Pa_GetDeviceInfo(inputParams.device)->defaultLowInputLatency;
    inputParams.hostApiSpecificStreamInfo = nullptr;

    outputParams.device = outputDevice;
    outputParams.channelCount = 2; // stereo for BlackHole 2ch
    outputParams.sampleFormat = paFloat32;
    outputParams.suggestedLatency = Pa_GetDeviceInfo(outputParams.device)->defaultLowOutputLatency;
    outputParams.hostApiSpecificStreamInfo = nullptr;

    PaStream* stream;
    Pa_OpenStream(&stream,
                  &inputParams,
                  &outputParams,
                  SAMPLE_RATE,
                  FRAMES_PER_BUFFER,
                  paNoFlag,
                  paCallback,
                  nullptr);

    Pa_StartStream(stream);

    std::cout << "Pitch shifter outputting to BlackHole 2ch and saving to output.wav. Press Enter to stop..." << std::endl;
    std::cin.get();

    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();

    if (sndfile) sf_close(sndfile);
    if (sf_error(sndfile) != SF_ERR_NO_ERROR) {
    std::cerr << "libsndfile error: " << sf_strerror(sndfile) << std::endl;
}
    std::cout << "Saved to output.wav" << std::endl;
    return 0;
}