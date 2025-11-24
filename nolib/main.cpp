// tone_asio_rtaudio.cpp
// Simple RtAudio example using ASIO on Windows to generate a sine tone.

#include <iostream>
#include <cmath>
#include "RtAudio.h"
#include <fstream>
#include <map>
#include <conio.h>
#include <chrono>
#include "rtmidi.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float empty[512];
unsigned long handleCount = 0;
struct StereoWavBuffer
{
    std::vector<float> samples;   // interleaved: L,R,L,R,...
    uint32_t sampleRate = 0;
    uint16_t channels = 0;        // original file channels (1 or 2)
    uint32_t frames = 0;          // number of frames in 'samples'
};

StereoWavBuffer snareBuffer;
StereoWavBuffer hh;
StereoWavBuffer kick;

struct Hit
{
    StereoWavBuffer& hit;
    unsigned long currentPos = 0;

    bool remove = false;
    float volume;
};
std::map<unsigned long, Hit> currentHits;


static uint32_t readUint32LE(std::ifstream& f)
{
    uint8_t b[4];
    f.read(reinterpret_cast<char*>(b), 4);
    return (uint32_t(b[0])) |
        (uint32_t(b[1]) << 8) |
        (uint32_t(b[2]) << 16) |
        (uint32_t(b[3]) << 24);
}

static uint16_t readUint16LE(std::ifstream& f)
{
    uint8_t b[2];
    f.read(reinterpret_cast<char*>(b), 2);
    return uint16_t(b[0]) | (uint16_t(b[1]) << 8);
}

// Loads a PCM 16-bit or 24-bit WAV and converts to stereo float buffer [-1,1]
StereoWavBuffer loadWavToStereoFloat(const std::string& filename)
{
    StereoWavBuffer out;

    std::ifstream f(filename, std::ios::binary);
    if (!f)
        throw std::runtime_error("Cannot open WAV file");

    // --- RIFF header ---
    char riff[4];
    f.read(riff, 4);
    if (std::strncmp(riff, "RIFF", 4) != 0)
        throw std::runtime_error("Not a RIFF file");

    (void)readUint32LE(f); // file size, ignore

    char wave[4];
    f.read(wave, 4);
    if (std::strncmp(wave, "WAVE", 4) != 0)
        throw std::runtime_error("Not a WAVE file");

    // --- read chunks until "fmt " and "data" ---
    uint16_t audioFormat = 0;
    uint16_t numChannels = 0;
    uint32_t sampleRate = 0;
    uint16_t bitsPerSample = 0;

    std::vector<uint8_t> rawData; // raw PCM bytes

    bool haveFmt = false;
    bool haveData = false;
    uint32_t dataChunkSize = 0;

    while (f && (!haveFmt || !haveData))
    {
        char chunkId[4];
        if (!f.read(chunkId, 4))
            break;

        uint32_t chunkSize = readUint32LE(f);

        if (std::strncmp(chunkId, "fmt ", 4) == 0)
        {
            audioFormat = readUint16LE(f);
            numChannels = readUint16LE(f);
            sampleRate = readUint32LE(f);
            auto rate =readUint32LE(f); // byteRate
            auto align = readUint16LE(f); // blockAlign
            bitsPerSample = readUint16LE(f);

            if (chunkSize > 16)
                f.seekg(chunkSize - 16, std::ios::cur);

            haveFmt = true;
        }
        else if (std::strncmp(chunkId, "data", 4) == 0)
        {
            dataChunkSize = chunkSize;
            rawData.resize(dataChunkSize);
            f.read(reinterpret_cast<char*>(rawData.data()), dataChunkSize);
            haveData = true;
        }
        else
        {
            f.seekg(chunkSize, std::ios::cur); // skip other chunks
        }
    }

    if (!haveFmt || !haveData)
        throw std::runtime_error("WAV missing fmt or data chunk");

    if (audioFormat != 1) // 1 = PCM
        throw std::runtime_error("Only PCM WAV supported");

    if (numChannels != 1 && numChannels != 2)
        throw std::runtime_error("Only mono or stereo WAV supported");

    if (bitsPerSample != 16 && bitsPerSample != 24)
        throw std::runtime_error("Only 16-bit or 24-bit WAV supported");

    out.sampleRate = sampleRate;
    out.channels = numChannels;

    const uint32_t bytesPerSample = bitsPerSample / 8;
    const uint32_t frameSizeBytes = bytesPerSample * numChannels;

    if (frameSizeBytes == 0)
        throw std::runtime_error("Invalid WAV frame size");

    uint32_t totalFrames = dataChunkSize / frameSizeBytes;
    out.frames = totalFrames;

    out.samples.resize(totalFrames * 2); // stereo output
    const float* outPtr = nullptr;       // just to remind ourselves type
    (void)outPtr;

    const uint8_t* p = rawData.data();
    const float scale16 = 1.0f / 32768.0f;
    const float scale24 = 1.0f / 8388608.0f; // 2^23

    for (uint32_t frame = 0; frame < totalFrames; ++frame)
    {
        float leftFloat = 0.0f;
        float rightFloat = 0.0f;

        // --- read samples for this frame, per channel ---

        if (bitsPerSample == 16)
        {
            // Little-endian signed 16-bit PCM
            for (uint16_t ch = 0; ch < numChannels; ++ch)
            {
                int16_t s = int16_t(p[0] | (p[1] << 8));
                p += 2;

                float fs = s * scale16;
                if (ch == 0) leftFloat = fs;
                else          rightFloat = fs;
            }
        }
        else // 24-bit
        {
            // Little-endian signed 24-bit PCM (3 bytes, sign bit at bit 23)
            for (uint16_t ch = 0; ch < numChannels; ++ch)
            {
                int32_t v = int32_t(p[0])
                    | (int32_t(p[1]) << 8)
                    | (int32_t(p[2]) << 16);
                p += 3;

                // Sign-extend from 24 bits to 32 bits
                if (v & 0x00800000)
                    v |= 0xFF000000;

                float fs = float(v * scale24);
                if (ch == 0) leftFloat = fs;
                else          rightFloat = fs;
            }
        }

        // Mono → copy to both channels
        if (numChannels == 1)
        {
            rightFloat = leftFloat;
        }

        out.samples[2 * frame + 0] = leftFloat;
        out.samples[2 * frame + 1] = rightFloat;
    }

    return out;
}

// Data passed to the audio callback
struct SineData
{
    double t;       // current phase [0, 1)
    double frequency;   // Hz
    double sampleRate;  // Hz
};

// RtAudio callback
int audioCallback(
    void* outputBuffer,
    void* /*inputBuffer*/,
    unsigned int nFrames,
    double t/*streamTime*/,
    RtAudioStreamStatus status,
    void* userData)
{
    if (status)
        std::cerr << "Stream underflow/overflow detected!" << std::endl;

    if (!outputBuffer || !userData)
        return 0;

    SineData* data = static_cast<SineData*>(userData);
    float* out = static_cast<float*>(outputBuffer);

      data->t;
    const double freq = data->frequency;
    const double sampleRate = data->sampleRate;
    const double phaseInc = freq / sampleRate;  // cycles per sample

    memset(out, 0, 2*nFrames*sizeof(float));
    std::vector<unsigned long> removes;
    for ( auto& hitPair : currentHits)
    {
        Hit& Hit = hitPair.second;
        auto& buffer = Hit.hit;
        auto missing = buffer.frames - Hit.currentPos;
        if (missing == 0)
        {
            removes.push_back(hitPair.first);
            continue;
        }
        int frames = nFrames < missing ? nFrames : missing;
        for (int i = 0; i < frames; ++i)
        {
            out[2 * i + 0] += Hit.volume*buffer.samples[2*Hit.currentPos]; // left
            out[2 * i + 1] += Hit.volume * buffer.samples[2 * Hit.currentPos+1]; // left; // right
            ++Hit.currentPos;
        }

    }
    for (auto& remove : removes)
        currentHits.erase(remove);

   
 
    return 0; // continue
}
uint64_t getTimeMs()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()
    ).count();
}

void midiCallback(double deltatime,
    std::vector<unsigned char>* message,
    void* userData)
{
    (void)userData;

    if (!message || message->empty())
        return;

    unsigned char status = (*message)[0];

    std::cout << "MIDI (" << deltatime << " s): ";

    // Basic decode: note on/off, CC
    unsigned char type = status & 0xF0;
    unsigned char channel = (status & 0x0F) + 1;

    if (type == 0x90 && message->size() >= 3)
    {
        unsigned char note = (*message)[1];
        unsigned char vel = (*message)[2];
        if (vel > 0)
        {
            std::cout << "Note ON  ch " << (int)channel
                << " note " << (int)note
                << " vel " << (int)vel;
            if (note==0x24)
                currentHits.insert({ handleCount ++, {snareBuffer,0,false,vel/100.0F} });
            if (note == 0x26)
                currentHits.insert({ handleCount++, {hh,0,false,vel / 100.0F} });
            if (note == 0x28)
                currentHits.insert({ handleCount++, {kick,0,false,vel / 100.0F} });
        }
        else
            std::cout << "Note OFF ch " << (int)channel
            << " note " << (int)note;
    }
    else if (type == 0x80 && message->size() >= 3)
    {
        unsigned char note = (*message)[1];
        unsigned char vel = (*message)[2];
        std::cout << "Note OFF ch " << (int)channel
            << " note " << (int)note
            << " vel " << (int)vel;
    }
    else if (type == 0xB0 && message->size() >= 3)
    {
        unsigned char cc = (*message)[1];
        unsigned char val = (*message)[2];
        std::cout << "CC      ch " << (int)channel
            << " cc " << (int)cc
            << " val " << (int)val;
    }

    std::cout << "\n";
}


int main()
{
  RtMidiIn midiIn;
    try
    {
      

        unsigned int nPorts = midiIn.getPortCount();
        std::cout << "MIDI input ports:\n";
        for (unsigned int i = 0; i < nPorts; ++i)
        {
            std::string name = midiIn.getPortName(i);
            std::cout << "  " << i << ": " << name << "\n";
        }

        if (nPorts == 0)
        {
            std::cout << "No MIDI input ports available.\n";
        
        }

        // Open first port (or let user choose)
        unsigned int portToOpen = 0;
        std::cout << "Opening port " << portToOpen << "...\n";
        midiIn.openPort(portToOpen);

        // Don't ignore anything (by default sysex/clock may be ignored)
        midiIn.ignoreTypes(false, false, false);

        // Setup callback
        midiIn.setCallback(&midiCallback);



    }
    catch (RtMidiError& e)
    {
        std::cerr << "RtMidiError: " << e.getMessage() << "\n";

    }


    memset(empty, 0, 256 * sizeof(float));
    snareBuffer = loadWavToStereoFloat("./samples/snare.wav" );
    hh = loadWavToStereoFloat("./samples/hh.wav");
 
    kick = loadWavToStereoFloat("./samples/kick/RD_K_4.wav");


    try
    {
        // Force ASIO API on Windows
        RtAudio dac(RtAudio::WINDOWS_ASIO);

        if (dac.getDeviceCount() == 0)
        {
            std::cerr << "No audio devices found (ASIO)." << std::endl;
            return 1;
        }

        // Use default output device for this API
        RtAudio::StreamParameters oParams;
        oParams.deviceId = dac.getDefaultOutputDevice();
        oParams.nChannels = 2;      // stereo
        oParams.firstChannel = 0;

        unsigned int sampleRate = 44100; // or 44100 if you prefer
        unsigned int bufferFrames = 256; // frames per buffer

        SineData data;
        data.t = 0.0;
        data.frequency = 440.0;  // A4
        data.sampleRate = static_cast<double>(sampleRate);

        RtAudio::StreamOptions options;
        // options.flags = RTAUDIO_HOG_DEVICE; // optional if you want exclusive

        dac.openStream(
            &oParams,
            nullptr,               // no input
            RTAUDIO_FLOAT32,       // 32-bit float samples
            sampleRate,
            &bufferFrames,
            &audioCallback,
            &data,
            &options);

        dac.startStream();

        std::cout << "ASIO sine tone running (440 Hz)." << std::endl;
        std::cout << "Press ENTER to stop..." << std::endl;
        auto t = getTimeMs();
        auto sh = false;
        while (true)
        {
             int key = _getch();
            if (key == 27) break;
            
            if (key == 's')
                currentHits.insert({ handleCount ++, {snareBuffer,0,false,5} });
            if (key == 'h')
                currentHits.insert({ handleCount ++, {hh,0,false,.5F} });
            if (key == 'k')
                currentHits.insert({ handleCount ++, {kick,0,false,.5F} });

          
           
        }
        if (dac.isStreamRunning())
            dac.stopStream();

        if (dac.isStreamOpen())
            dac.closeStream();

        midiIn.closePort();
    }

    catch (std::exception& e)
    {
        std::cerr << "Std exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}