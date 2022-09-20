#include "Screenshots.h"


namespace megamol {
namespace frontend_resources {

ScreenshotImageData::ScreenshotImageData(ScreenshotImageData_YUV const& yuv) {
    resize(yuv.width, yuv.height);
    std::transform(yuv.image.begin(), yuv.image.end(), image.begin(), [](ScreenshotImageData_YUV::Pixel const& yuv) {
        // https://www.fourcc.org/fccyvrgb.php
        ScreenshotImageData::Pixel rgb;
        rgb.r = static_cast<uint8_t>(
            1.164f * (static_cast<float>(yuv.y) - 16.f) + 1.596f * (static_cast<float>(yuv.v) - 128.f));
        rgb.g = static_cast<uint8_t>(1.164f * (static_cast<float>(yuv.y) - 16.f) -
                                     0.813f * (static_cast<float>(yuv.v) - 128.f) -
                                     0.391f * (static_cast<float>(yuv.u) - 128.f));
        rgb.b = static_cast<uint8_t>(
            1.164f * (static_cast<float>(yuv.y) - 16.f) + 2.018f * (static_cast<float>(yuv.u) - 128.f));
        return rgb;
    });
}


ScreenshotImageData_YUV::ScreenshotImageData_YUV(ScreenshotImageData const& rgb) {
    resize(rgb.width, rgb.height);
    std::transform(rgb.image.begin(), rgb.image.end(), image.begin(), [](ScreenshotImageData::Pixel const& rgb) {
        // https://www.fourcc.org/fccyvrgb.php
        ScreenshotImageData_YUV::Pixel yuv;
        yuv.y = static_cast<uint8_t>((0.257f * static_cast<float>(rgb.r)) + (0.504f * static_cast<float>(rgb.g)) +
                                     (0.098f * static_cast<float>(rgb.b)) + 16.f);
        yuv.u = static_cast<uint8_t>(-(0.148f * static_cast<float>(rgb.r)) - (0.291f * static_cast<float>(rgb.g)) +
                                     (0.439f * static_cast<float>(rgb.b)) + 128.f);
        yuv.v = static_cast<uint8_t>((0.439f * static_cast<float>(rgb.r)) - (0.368f * static_cast<float>(rgb.g)) -
                                     (0.071f * static_cast<float>(rgb.b)) + 128.f);
        return yuv;
    });
}


ScreenshotImageData_I420::ScreenshotImageData_I420(ScreenshotImageData_YUV const& yuv) {
    resize(yuv.width, yuv.height);
    for (int64_t h = yuv.height - 1; h >= 0; --h) {
        for (size_t w = 0; w < yuv.width; ++w) {
            auto const& pix_yuv = yuv.image[(height - 1 - h) * width + w];
            Y[h * width + w] = pix_yuv.y;
            if (h % 2 == 0 && w % 2 == 0) {
                U[0.5f * (h * width / 2 + w)] = pix_yuv.u;
                V[0.5f * (h * width / 2 + w)] = pix_yuv.v;
            }
        }
    }
}

} // namespace frontend_resources
} // namespace megamol
