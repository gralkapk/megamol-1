#include "OpenGL_GLFW_Utils.h"

#ifdef _WIN32
#include <glad/wgl.h>
#else
#include <glad/glx.h>
#endif

namespace megamol::frontend {
void request_opengl(
    int const versionMajor, int const versionMinor, bool const enableKHRDebug, bool const useCoreProfile) {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, versionMajor);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, versionMinor);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, enableKHRDebug ? GLFW_TRUE : GLFW_FALSE);
    // context profiles available since 3.2
    bool has_profiles = versionMajor > 3 || (versionMajor == 3 && versionMinor >= 2);
    if (has_profiles)
        glfwWindowHint(GLFW_OPENGL_PROFILE, useCoreProfile ? GLFW_OPENGL_CORE_PROFILE : GLFW_OPENGL_COMPAT_PROFILE);

    std::string profile_name =
        has_profiles ? ((useCoreProfile ? "Core" : "Compatibility") + std::string(" Profile")) : "";
    core::utility::log::Log::DefaultLog.WriteInfo(
        ("Requesting OpenGL " + std::to_string(versionMajor) + "." + std::to_string(versionMinor) + " " + profile_name +
            (enableKHRDebug ? ", Debug Context" : ""))
            .c_str());
}

bool load_opengl_functions(int& major, int& minor) {
    auto version = gladLoaderLoadGL();
    major = GLAD_VERSION_MAJOR(version);
    minor = GLAD_VERSION_MINOR(version);
    if (version == 0) {
        core::utility::log::Log::DefaultLog.WriteError("Failed to load OpenGL functions via glad");
        return false;
    }
#ifdef _WIN32
    if (gladLoaderLoadWGL(wglGetCurrentDC()) == 0) {
        core::utility::log::Log::DefaultLog.WriteError("Failed to load OpenGL WGL functions via glad");
        return false;
    }
#else
    Display* display = XOpenDisplay(nullptr);
    gladLoaderLoadGLX(display, DefaultScreen(display));
    XCloseDisplay(display);
#endif
}

std::vector<std::string> get_extensions(int const major) {
    std::vector<std::string> ret_ext;
    if (major < 3) {
        auto ext = std::string((char const*) glGetString(GL_EXTENSIONS));
        std::istringstream iss(ext);
        std::copy(
            std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), std::back_inserter(ret_ext));
    } else {
        GLint num_ext = 0;
        glGetIntegerv(GL_NUM_EXTENSIONS, &num_ext);

        ret_ext.resize(num_ext);
        for (GLint i = 0; i < num_ext; ++i) {
            ret_ext[i] = std::string((char const*) glGetStringi(GL_EXTENSIONS, i));
        }
    }

    core::utility::log::Log::DefaultLog.WriteInfo((
        std::string("OpenGL Context Info") + "\n\tVersion:  " + reinterpret_cast<const char*>(glGetString(GL_VERSION)) +
        "\n\tVendor:   " + reinterpret_cast<const char*>(glGetString(GL_VENDOR)) +
        "\n\tRenderer: " + reinterpret_cast<const char*>(glGetString(GL_RENDERER)) +
        "\n\tGLSL:     " + reinterpret_cast<const char*>(glGetString(GL_SHADING_LANGUAGE_VERSION)))
                                                      .c_str());

    return ret_ext;
}

void setup_khr_debug() {
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    GLuint ignorethis;
    ignorethis = 131185;
    glDebugMessageControl(GL_DEBUG_SOURCE_API, GL_DEBUG_TYPE_OTHER, GL_DONT_CARE, 1, &ignorethis, GL_FALSE);
    ignorethis = 131184;
    glDebugMessageControl(GL_DEBUG_SOURCE_API, GL_DEBUG_TYPE_OTHER, GL_DONT_CARE, 1, &ignorethis, GL_FALSE);
    ignorethis = 131204;
    glDebugMessageControl(GL_DEBUG_SOURCE_API, GL_DEBUG_TYPE_OTHER, GL_DONT_CARE, 1, &ignorethis, GL_FALSE);
    glDebugMessageControl(GL_DEBUG_SOURCE_APPLICATION, GL_DEBUG_TYPE_PUSH_GROUP, GL_DONT_CARE, 0, nullptr, GL_FALSE);
    glDebugMessageControl(GL_DEBUG_SOURCE_APPLICATION, GL_DEBUG_TYPE_POP_GROUP, GL_DONT_CARE, 0, nullptr, GL_FALSE);
    glDebugMessageCallback(opengl_debug_message_callback, nullptr);
    core::utility::log::Log::DefaultLog.WriteInfo("Enabled OpenGL debug context. Will print debug messages.");
}
} // namespace megamol::frontend
