/*
 * mmstd.volume.cpp
 *
 * Copyright (C) 2009 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "mmstd.volume.h"
#include "api/MegaMolCore.std.h"
#include "ModuleAutoDescription.h"
#include "vislib/vislibversion.h"
#include "vislib/Log.h"
#include "vislib/ThreadSafeStackTrace.h"


/*
 * mmplgPluginAPIVersion
 */
MMSTD_VOLUME_API int mmplgPluginAPIVersion(void) {
    return 100;
}


/*
 * mmplgPluginName
 */
MMSTD_VOLUME_API const char * mmplgPluginName(void) {
    return "mmstd.volume";
}


/*
 * mmplgPluginDescription
 */
MMSTD_VOLUME_API const char * mmplgPluginDescription(void) {
    return "Template for MegaMol Plugins (TODO: CHANGE this description)";
}


/*
 * mmplgCoreCompatibilityValue
 */
MMSTD_VOLUME_API const void * mmplgCoreCompatibilityValue(void) {
    static const mmplgCompatibilityValues compRev = {
        sizeof(mmplgCompatibilityValues),
        MEGAMOL_CORE_COMP_REV,
        VISLIB_VERSION_REVISION
    };
    return &compRev;
}


/*
 * mmplgModuleCount
 */
MMSTD_VOLUME_API int mmplgModuleCount(void) {
    return 0; // TODO: Implement
}


/*
 * mmplgModuleDescription
 */
MMSTD_VOLUME_API void* mmplgModuleDescription(int idx) {
    return NULL; // TODO: Implement
}


/*
 * mmplgCallCount
 */
MMSTD_VOLUME_API int mmplgCallCount(void) {
    return 0; // TODO: Implement
}


/*
 * mmplgCallDescription
 */
MMSTD_VOLUME_API void* mmplgCallDescription(int idx) {
    return NULL; // TODO: Implement
}


/*
 * mmplgConnectStatics
 */
MMSTD_VOLUME_API bool mmplgConnectStatics(int which, void* value) {
    switch (which) {

        case 1: // vislib::log
            vislib::sys::Log::DefaultLog.SetLogFileName(static_cast<const char*>(NULL), false);
            vislib::sys::Log::DefaultLog.SetLevel(vislib::sys::Log::LEVEL_NONE);
            vislib::sys::Log::DefaultLog.SetEchoTarget(new vislib::sys::Log::RedirectTarget(static_cast<vislib::sys::Log*>(value)));
            vislib::sys::Log::DefaultLog.SetEchoLevel(vislib::sys::Log::LEVEL_ALL);
            vislib::sys::Log::DefaultLog.EchoOfflineMessages(true);
            return true;

        case 2: // vislib::stacktrace
            return vislib::sys::ThreadSafeStackTrace::Initialise(
                *static_cast<const vislib::SmartPtr<vislib::StackTrace>*>(value), true);

    }
    return false;
}
